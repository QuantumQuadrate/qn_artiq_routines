"""
Testing if the ttl logic gates add any jitter to the ttl signals.
I generate a train of ttl pulses at ttl12 and ttl15 simultaneously. ttl12 pulses go through the logic gate and then oscilloscope.
I use ttl_GRIN2_switch temporary as the Y-channel on the logic.

ttl15 pulses go to the oscilloscope directly and are used to trigger the scope.

I observe on the scope that the logic gate adds a 20ns delay and shortens the pulse from 100ns to about 80ns. But the jitter time is only about 1ns.

Akbar 2025-11-24

"""

from artiq.experiment import *
import logging
import numpy as np

import sys, os
# get the current working directory
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from utilities.conversions import dB_to_V_kernel as dB_to_V


class Testing_TTL_logic_gates(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

    def prepare(self):
        self.base.prepare()

    @kernel
    def run(self):
        self.core.reset()

        self.base.initialize_hardware()

        delay(1*ms)
        self.ttl_SPCM0_logic.off()
        self.ttl_SPCM1_logic.off()
        self.ttl_GRIN2_switch.on()
        delay(1 * ms)

        for i in range(100):
            with parallel:
                self.ttl_SPCM0_logic.pulse(100*ns)
                self.ttl_SPCM1_logic.pulse(100*ns)
            delay(100*ms)


        print("code done!")
