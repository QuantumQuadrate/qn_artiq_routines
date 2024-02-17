from artiq.experiment import *
import math

import sys, os
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')


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