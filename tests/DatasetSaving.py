"""
For testing how saving datasets with different types and effects size of h5 file.
"""

from artiq.experiment import *
import numpy as np

import sys
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class DatasetSaving(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.setattr_argument("n_datapoints", NumberValue(100, ndecimals=0, step=1))
        # self.setattr_argument("python_data_type",StringValue('float'))
        self.setattr_argument("round_decimals", NumberValue(10, ndecimals=0, step=1))

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        np.random.seed(0)
        # python_data_type = eval(self.python_data_type)
        self.data = np.array([np.round(100*x,self.round_decimals) for x in np.random.rand(self.n_datapoints)])

    def run(self):
        self.set_dataset("random_data",self.data)

        print("Test finished")


