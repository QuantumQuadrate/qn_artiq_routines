"""
a simple experiment to print out the Sampler card values
"""
from artiq.experiment import *

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
import numpy as np


class SamplerPrint(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()
        self.sampler_list = [self.sampler0, self.sampler1, self.sampler2]
        self.buffer = np.array([0.0]*8)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        print("Print out measurement from each sampler:")
        delay(10*ms)

        for sampler in self.sampler_list:
            sampler.sample(self.buffer)
            delay(10*ms)
            print(self.buffer)
            delay(10*ms)

        print("End of Sampler Print")
