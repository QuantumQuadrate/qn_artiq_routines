"""
Example 2 demonstrates how to control multiple Thorlabs rotators as NDSP devices
when you need to get values from them (e.g. position) from an experiment on the kernel, 
where you want to interact with the device between doing things with the Sinara hardware. 
This could be relevant for a some kind of optimization or feedback experiment where you
 query to the device to get the position and then do something, e.g.e measure a photodiode
voltage with the Sampler, and then move the rotators based on the result.

For an example of how to query the device from outside a kernel function, see Example 1.

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
        self.positions = [0.0]*len(self.rotors)
    
    def get_rotator_positions(self, names: TList(TStr)) -> TList(TFloat ):
        """wrapper function"""
        positions = self.k10cr1_ndsp.get_position(names)
        return positions
    
    def rotators_are_moving(self, names: TList(TStr)) -> TBool:
        """wrapper function"""
        still_moving = self.k10cr1_ndsp.is_moving(names)
        return still_moving
         
    @kernel 
    def run(self):
        """
        Do some stuff on the kernel, query the rotators to get their positions, 
        do more stuff on the kernel, then move the rotators by some amount.
        """

        self.core.reset()

        # do some stuff with your hardware
        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2*s)

        # get the rotator positions
        self.positions = self.get_rotator_positions(self.rotors)
        print("rotators:",self.rotors)
        print("positions (deg.):",self.positions)
        delay(10*ms)

        # do something, e.g. measure a photodetector after the rotators
        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2*s)

        with sequential:
            # move all of the rotors by some amount, with only one rpc call:
            self.k10cr1_ndsp.move_by(self.rotors, [90,80])

            # wait for the rotators to move
            now = now_mu()
            t = 0
            while self.rotators_are_moving(self.rotors):
                delay(0.1*ms)
                t += now_mu() - now
            at_mu(now + t + 100)

         # do something, e.g. measure a photodetector after the rotators
        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2*s)

        # get the position of only one rotator, then move it
        self.positions = self.get_rotator_positions(["780_QWP"])
        print("rotators:","780_QWP")
        print("positions (deg.):",self.positions)
        delay(10*ms)
        self.k10cr1_ndsp.move_by(["780_QWP"], [30])

        # finally, the rest of your experiment
        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2*s)