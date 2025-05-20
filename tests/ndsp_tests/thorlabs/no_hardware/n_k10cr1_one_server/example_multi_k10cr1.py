"""
This example and driver used for it were intended to allow controlling multiple
Thorlabs rotators with only a single RCP call to talk to all rotators instead
of one RPC call each. However, I was never able to get this work, which seems
like a quirk of the way NDSP works in ARTIQ. - PH 2024.11.12

1. connect thorlabs K10CR1 rotators to your machine with USB, and make sure
the serial numbers in the device_db match the rotators you connected.

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


class K10CR1Example(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")
    
    def run(self):

        # Move both waveplates with only one NDSP call.
        # More efficient than the two_k10cr1_one_server example!
        degrees = self.k10cr1_ndsp.get_position(names=['780_QWP', '780_HWP'])
        print(f"Rotator angles: {degrees} deg.")

        # Move both waveplates with only one NDSP call.
        # More efficient than the two_k10cr1_one_server example!
        self.k10cr1_ndsp.move_by(degrees=[20, 30], names=['780_QWP', '780_HWP'])

        # while self.k10cr1_ndsp.is_moving('780_QWP'):
        #     sleep(0.1)

        # degrees = self.k10cr1_ndsp.get_position('780_QWP')
        # print(f"Rotator angle = {degrees:.2f} deg.")

        # self.k10cr1_ndsp.move_by(-20, '780_QWP')

        # while self.k10cr1_ndsp.is_moving('780_QWP'):
        #     sleep(0.1)

        # degrees = self.k10cr1_ndsp.get_position('780_QWP')
        # print(f"Rotator angle = {degrees:.2f} deg.")

        # # move the 780 HWP
        # degrees = self.k10cr1_ndsp.get_position('780_HWP')
        # print(f"Rotator angle = {degrees:.2f} deg.")

        # self.k10cr1_ndsp.move_by(20, '780_HWP')

        # while self.k10cr1_ndsp.is_moving('780_HWP'):
        #     sleep(0.1)

        # degrees = self.k10cr1_ndsp.get_position('780_HWP')
        # print(f"Rotator angle = {degrees:.2f} deg.")

        # self.k10cr1_ndsp.move_by(-20, '780_HWP')

        # while self.k10cr1_ndsp.is_moving('780_HWP'):
        #     sleep(0.1)

        # degrees = self.k10cr1_ndsp.get_position('780_HWP')
        # print(f"Rotator angle = {degrees:.2f} deg.")