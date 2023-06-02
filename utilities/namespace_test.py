from artiq.experiment import *
from utilities.BaseExperiment import *

class NamespaceTest(EnvExperiment):

    def build(self):

        assert t_MOT_loading in globals(), "oops, couldn't find the global variable t_MOT_loading"
