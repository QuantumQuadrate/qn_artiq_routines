"""
1. in a terminal, run the ndsp launcher in this directory with
    python launcher.py

2. in another terminal opened in this directory, run your artiq environment then run

    artiq_run ndsp_example.py --device-db=.\device_db.py
"""

# here's how you do ndsp stuff in your artiq scripts
from artiq.experiment import *

class NDSPExample(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("led0")

        try:
            self.setattr_device("example_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")
    
    def prepare(self):
        self.some_value = 0.0

    def ndsp_return_a_number(self) -> TFloat:
        """wrapper function"""
        value = self.example_ndsp.return_a_number()
        return value

    @kernel
    def run(self):
        self.core.reset()

        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2)

        # this errors out with "error: cannot unify float with NoneType", 
        # likely due to some quirk of how NDSPs work in artiq. 
        # self.some_value = self.example_ndsp.return_a_number() 

        # same functionality but with a wrapper function:
        self.some_value = self.ndsp_return_a_number()

        print(self.some_value)
        delay(10*ms)

        for i in range(10):
            self.led0.pulse(0.2*s)
            delay(0.2)

        self.example_ndsp.close()

