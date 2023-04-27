"""
A template for making a subroutine servo loop using a DDS channel and a Sampler

This could be used to control the power of an AOM for stabilizing laser intensity. The
only critical part missing from this code is in line 45 where the change in DDS amplitude
is set. The code runs.

Intended usage:
1. Instantiate the servo in the prepare method of the parent artiq experiment
2. Call get_dds_settings() in the run method of the parent, after the relevant AOMs have been set
to their nominal desired powers. For now, it is assumed that the frequency of the AOMs
being controlled with this servo do not change throughout the experiment.
3. Call run() during the parent experiment loop whenever it is desired to tune up the
AOM powers.

Preston's notes
* could be nice to have a DDS wrapper class, maybe just a NamedTuple, so that the amplitudes,
frequencies, and identifiers aren't stored in independent lists.
* should add the option to specify hardware from a config file
"""

from artiq.experiment import *

# todo: an ideal workflow would be to define the groups of dds channels, sampler channels,
#  and transfer functions in a cfg file.

class AOMPowerStabilizer:

    # Todo: add option to set which channels are used with a cfg file.

    def __init__(self, experiment, dds_names, sampler_name, sampler_channels):
        """
        An experiment subsequence for reading a Sampler and adjusting Urukul output power.

        :param experiment: a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        """
        self.exp = experiment
        self.n_iterations = 5 # number of times to adjust dds power per run() call
        self.dds_names = dds_names
        self.sampler_name = sampler_name # only one sampler allowed for now
        self.n_channels = len(self.dds_names)
        self.amplitudes = [0.0]*len(self.dds_names)
        self.frequencies = [0.0]*len(self.dds_names)
        self.dds_list = [] # will hold list of dds objects
        self.sampler_channels = sampler_channels
        self.sample_buffer = [0.0]*8

        # get hardware references by name from the parent experiment
        for i in range(self.n_channels):
            self.dds_list.append(getattr(self.exp, self.dds_names[i]))
        self.sampler = getattr(self.exp, self.sampler_name)

    @kernel
    def get_dds_settings(self):
        """
        Populates the amplitude and frequency lists
        :return:
        """
        for i in range(self.n_channels):
            self.frequencies[i], _, self.amplitudes[i] = self.dds_list[i].get()

    @kernel
    def run(self):
        """
        Run the feedback loop
        :return:
        """

        # in each iteration, measure the sampler and set the dds amplitude accordingly
        for i in range(self.n_iterations):
            with sequential:
                self.sampler.sample(self.sample_buffer)

                for ch in range(self.n_channels):

                    delta_ampl = 0 # todo adjust this based on the sampler measurement
                    delay(0.5*ms)
                    self.dds_list[ch].set(self.frequencies[ch], amplitude=self.amplitudes[ch])


class AOMPowerStabilizerTest(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        # initialize the hardware that we know the subroutine will use
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("sampler0")

    def prepare(self):
        self.AOMservo = AOMPowerStabilizer(self, ["urukul2_ch0"], "sampler0", [0])
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

        self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set
        start_mu = now_mu()

        at_mu(start_mu + self.core.seconds_to_mu(5*ms))
        # todo: put this inside an experiment loop
        self.AOMservo.run()


        print("experiment finished")