"""
Example 1 demonstrates how to control multiple Thorlabs rotators as NDSP devices
when you DO NOT need to get values from them (e.g. position) from an experiment
on the kernel. This could be relevant for a tomography experiment where the command
to check the current positions of the rotators (e.g. to save them to a database) and
then move them is done off the kernel before doing self.experiment_function() to run
your experiment on the kernel. 

For an example of how to query the device from a kernel function, see Example 2.

Instructions:

1a. connect two thorlabs K10CR1 rotators to your machine with USB
1b. in .\device_db.py, in the k10cr1_ndsp entry, 
    make sure serial numbers match those of your K10CR1 devices.

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
    def experiment_function(self):
        """some experiment that runs on the kernel"""

        self.core.reset()

        for i in range(10):
            self.led0.pulse(0.5*s)
            delay(0.1*s)
        
    def run(self):
        """
        Query the rotators to get their positions, move them by some amount,
        then do an experiment on the kernel, and finally reset the rotator positions.
        """

        degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print([f"{rotor} at {deg} deg." for rotor,deg in zip(self.rotors,degrees_list)])

        # move all of the rotors by some amount, with only one rpc call:
        self.k10cr1_ndsp.move_by(self.rotors, [90,80])

        # wait for the rotators to move. probably a more sophisticated way of doing this
        sleep(2*s)
        degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print([f"{rotor} at {deg} deg." for rotor,deg in zip(self.rotors,degrees_list)])

        # run our experiment on the kernel
        self.experiment_function()
        print("kernel experiment done")

        # move all of the rotors by some amount, with only one rpc call:
        self.k10cr1_ndsp.move_by(self.rotors, [-90,-80])

        sleep(2*s)

        degrees_list = self.k10cr1_ndsp.get_position(self.rotors)
        print([f"{rotor} at {deg} deg." for rotor,deg in zip(self.rotors,degrees_list)])