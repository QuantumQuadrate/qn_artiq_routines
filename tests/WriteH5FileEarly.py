from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class WriteH5FileEarly(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

    def prepare(self):
        self.base.prepare()

    @kernel
    def run(self):
        self.set_dataset("test",[0])

        for i in range(100):

            self.append_to_dataset("test", np.cos(2*np.pi*i/100))
            delay(1*ms)
            self.write_results() # can write each iteration

        # self.write_results() # or we can write everything at the end

        self.print_async("write h5 file test finished")