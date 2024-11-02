"""
1. connect a thorlabs K10CR1 rotator to your machine with USB

2. open a terminal in the same directory as this script, 
enter your artiq virtual environment, then run the 
NDSP server for this test with

    python launcher_multi_rotor.py

3. In another terminal in this directory in the artiq environment,
run this script with artiq:

    artiq_run example_multi_k10cr1.py --device-db=.\device_db.py
"""

# here's how you do ndsp stuff in your artiq scripts
from artiq.experiment import *
from time import sleep


class K10CR1MultiExample(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("led0")

        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")

        self.rotors = ["780_QWP","780_HWP"]
        self.degrees_list = [0]*len(self.rotors)
    
    @kernel
    def run(self):
        self.degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print(self.rotors)
        print(self.degrees_list)

        # move all of the rotors by some amount, with only one rpc call:
        self.k10cr1_ndsp.move_by(self.rotors, [20,30])
        self.led0.pulse(1*s)

        self.degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print(self.rotors)
        print(self.degrees_list)

        # you don't have to move all of the rotors we connected;
        # you can move only one by passing in only one name:
        self.k10cr1_ndsp.move_by(['780_HWP'], [-20])
        self.led0.pulse(1*s)

        self.k10cr1_ndsp.move_by(['780_QWP'], [-30])
        self.led0.pulse(1*s)

        self.degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print(self.rotors)
        print(self.degrees_list)