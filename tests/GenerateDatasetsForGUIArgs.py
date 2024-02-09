from artiq.experiment import *
import math

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer

class GenerateDatasetsFromGuiArgs(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("my_argument",NumberValue(5))

        self.base.set_datasets_from_gui_args()

    def run(self):
        pass