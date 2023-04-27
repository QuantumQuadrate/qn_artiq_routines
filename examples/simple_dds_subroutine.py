"""
A simple experiment which uses a user-defined class to perform a subroutine.

The user-defined subroutine Subsequence is given the reference of the parent
experiment when it is instantiated. In this way, it can reference hardware objects
which are initialized in the parent.
"""

from artiq.experiment import *

class Subsequence1:

    t_pulse = 200*ms

    def __init__(self, experiment):
        """
        An experiment subsequence which must be instantiated inside an experiment.

        :param experiment: a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        """
        self.exp = experiment

    @kernel
    def run(self):
        self.exp.urukul2_ch0.sw.on()
        delay(self.t_pulse)
        self.exp.urukul2_ch0.sw.off()

class ExperimentWithHierarchy(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        # initialize the hardware that we know the subroutine will use
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")

    def prepare(self):
        self.sequence1 = Subsequence1(self)

    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.urukul2_ch0.cpld.init()
        self.urukul2_ch0.init()
        self.urukul2_ch0.set_att(float(0))
        self.urukul2_ch0.set(70.0 * MHz, amplitude=0.1)

        delay(1000*ms)

        self.sequence1.run()

        print("experiment finished")



