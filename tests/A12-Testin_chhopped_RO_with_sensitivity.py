"""
Testing the chopped RO using self.ttl0._set_sensitivity(0) and (1) in DMA.
Here, I am chopping the FORT and change the ttl sensitivity both on DMA.

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


class Testing_chopped_RO_with_sensitivity(EnvExperiment):

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
        self.ttl_SPCM0_logic.on() ### turn on the ttl to let the signal pass the logic gates
        self.ttl_SPCM1_logic.on()

        delay(1 * ms)




        print("code done!")
