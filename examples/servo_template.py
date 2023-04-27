"""
A template for making a subroutine servo loop using a DDS channel and a Sampler

This could be used to control the power of an AOM. The only critical part missing
from this code is in line 45 where the change in DDS amplitude is set. The code runs.
"""
from artiq.experiment import *

# todo: an ideal workflow would be to define the groups of dds channels, sampler channels,
#  and transfer functions in a cfg file.

class AOMPowerStabilizer:

    n_iterations = 5

    def __init__(self, experiment):
        """
        An experiment subsequence for reading a Sampler and adjusting Urukul output power.

        The Urukul and Sampler channels to be used are hardcode here, for now. Todo: set
        which channels are used in a cfg file.

        :param experiment: a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        """
        self.exp = experiment
        self.sample_buffer = [0.0]*8


    @kernel
    def run(self):
        """
        Run the feedback loop
        :return:
        """

        # in each iteration, measure the sampler and set the dds amplitude accordingly
        for i in range(self.n_iterations):
            with sequential:
                self.exp.sampler0.sample(self.sample_buffer)
                delta_ampl = 0 # todo adjust this based on the sampler measurement
                delay(0.5*ms)
                self.exp.urukul2_ch0.set(self.exp.dds2_0_freq,
                                         amplitude=self.exp.dds2_0_ampl+delta_ampl)


class AOMPowerStabilizerTest(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        # initialize the hardware that we know the subroutine will use
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("sampler0")

    def prepare(self):
        self.AOMservo = AOMPowerStabilizer(self)
        self.dds2_0_ampl = 0.1
        self.dds2_0_freq = 70.0 * MHz

    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.urukul2_cpld.init()
        self.urukul2_ch0.init()
        self.urukul2_ch0.set_att(float(0))
        self.urukul2_ch0.set(self.dds2_0_freq, amplitude=self.dds2_0_ampl)
        self.sampler0.init()


        start_mu = now_mu()

        at_mu(start_mu + self.core.seconds_to_mu(5*ms))
        # todo: put this inside an experiment loop
        self.AOMservo.run()


        print("experiment finished")