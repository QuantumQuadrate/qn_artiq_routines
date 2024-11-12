"""
1. connect a thorlabs K10CR1 rotator to your machine with USB

2. open a terminal in the same directory as this script, 
enter your artiq virtual environment, then run the 
NDSP server for this test with

    python launcher_single_rotor.py

3. In another terminal in this directory in the artiq environment,
run this script with artiq:

    artiq_run example_multi_k10cr1.py --device-db=.\device_db.py

If you wanted to connect to multiple K10CR1 units, at least if we stick to the
way this example is structured using the K10CR1_NDSP_Driver class, you would need 
to launch a separate server for every Thorlabs rotator, which is probably not a 
good solution. A better solution would have only one server for all the K10CR1 units,
and this is what the K10CR1_Multi_NDSP_Driver tries to accomplish. See the multi
examples."""

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

        # move the 780 QWP
        degrees = self.k10cr1_ndsp.get_position('780_QWP')
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(20, '780_QWP')

        while self.k10cr1_ndsp.is_moving('780_QWP'):
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position('780_QWP')
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(-20, '780_QWP')

        while self.k10cr1_ndsp.is_moving('780_QWP'):
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position('780_QWP')
        print(f"Rotator angle = {degrees:.2f} deg.")

        # move the 780 HWP
        degrees = self.k10cr1_ndsp.get_position('780_HWP')
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(20, '780_HWP')

        while self.k10cr1_ndsp.is_moving('780_HWP'):
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position('780_HWP')
        print(f"Rotator angle = {degrees:.2f} deg.")

        self.k10cr1_ndsp.move_by(-20, '780_HWP')

        while self.k10cr1_ndsp.is_moving('780_HWP'):
            sleep(0.1)

        degrees = self.k10cr1_ndsp.get_position('780_HWP')
        print(f"Rotator angle = {degrees:.2f} deg.")