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

        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")

        self.rotors = ["780_QWP","780_HWP"]
    
    def run(self):
        degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print([f"{rotor} at {deg} deg." for rotor,deg in zip(self.rotors,degrees_list)])

        # move all of the rotors by some amount, with only one rpc call:
        self.k10cr1_ndsp.move_by(self.rotors, [90,80])

        # wait until all rotors are done
        i = 0
        while self.k10cr1_ndsp.is_moving(self.rotors):
            sleep(0.1)  
            i += 1
        print(i)

        degrees = self.k10cr1_ndsp.get_position(self.rotors)
        print(f"Rotator angles = {degrees} deg.")

        # you don't have to move all of the rotors we connected;
        # you can move only one by passing in only one name:
        self.k10cr1_ndsp.move_by(['780_HWP'], [-90])

        # we can wait until the previous rotor is done moving
        # (this isn't necessary, I'm just showing it as an example)
        while self.k10cr1_ndsp.is_moving(['780_HWP']):
            sleep(0.1)

        self.k10cr1_ndsp.move_by(['780_QWP'], [-80])

        while self.k10cr1_ndsp.is_moving(self.rotors):
            sleep(0.1)

        degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print([f"{rotor} at {deg} deg." for rotor,deg in zip(self.rotors,degrees_list)])
