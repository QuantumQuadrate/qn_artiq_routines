"""
1. connect a thorlabs K10CR1 rotator to your machine with USB

2. open a terminal in the same directory as this script, 
enter your artiq virtual environment, then run the 
NDSP server for this test with

    python launcher_single_rotor.py

3. In another terminal in this directory in the artiq environment,
run this script with artiq:

    artiq_run example_multi_k10cr1.py --device-db=.\device_db.py
"""

# here's how you do ndsp stuff in your artiq scripts
from artiq.experiment import *
from time import sleep


class K10CR1Example(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")
    
    def run(self):
        degrees = self.k10cr1_ndsp.get_position()
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(20)

        while self.k10cr1_ndsp.is_moving():
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position()
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(-20)

        while self.k10cr1_ndsp.is_moving():
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position()
        print(f"Rotator angle = {degrees:.2f} deg.")