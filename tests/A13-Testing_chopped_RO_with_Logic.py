"""
Testing the chopped RO using the TTL logic gate devices to block the signals externally when the FORT is on. The SPCM gates are open
the entire time (say 20ms), but their output to the ttl input is blocked by the logic gates when the FORT is ON.

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


class Testing_chopped_RO_with_Logic(EnvExperiment):

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
        delay(1 * ms)

        for i in range(100):
            with parallel:
                self.ttl_SPCM0_logic.pulse(100*ns)
                self.ttl_SPCM1_logic.pulse(100*ns)
            delay(100*ms)


        print("code done!")
