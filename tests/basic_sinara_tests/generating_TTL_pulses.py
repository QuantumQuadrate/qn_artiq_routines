"""
Generating a series of short TTL pulses in Node2 to send to node1 for testing the click counter.

"""
from artiq.experiment import *
import numpy as np

class generating_TTL(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl4")

    @kernel
    def run(self):
        self.core.reset()
        self.ttl4.output()

        delay(10 * us)

        for _ in range(5000000):
            self.ttl4.pulse(20 * ns)
            delay(10 * us)

        delay (10 * ms)
        print("*****************   ALL DONE!   *****************")