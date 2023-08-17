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

*A better design is to allow specifying which dds channels we want to stabilize,
at the experiment class level instead of hard-coding them here. We want to retain
the ability to easily run experiments that do not try to stabilize channels that
aren't used.

*All of the argument lists here can be obtained from a dictionary defined in the stabilizer
file itself, so the user doesn't need to worry about passing that stuff.

*We could separate the stabilizer into two classes: for the chip-based beams and one for the
rest since the chip based beams have to be stabilized one at a time. However, we could also
just hard-code in a block that runs the fiber AOM stabilization separately. this seems
fine.
"""

from artiq.experiment import *
import math

# todo: an ideal workflow would be to define the groups of dds channels, sampler channels,
#  and transfer functions in a cfg file.

laser_stabilizer_dict = {
    'dds_cooling_PD':
        {
            'sampler': 'sampler',
            'sampler_channel': 7, # the channel connected to the appropriate PD
            'transfer_function': lambda x : x, # converts volts to optical mW
            'setpoint': 6, # value in mW,
            'p': 0.07, # the proportionality constant
        }
}

class AOMPowerStabilizer:

    # Todo: add option to set which channels are used with a cfg file.

    def __init__(self, experiment,
                 #dds_names, sampler_name, sampler_channels, transfer_functions,
                 #setpoints, proportionals,
                 iters=10, t_meas_delay=10*ms):
        """
        An experiment subsequence for reading a Sampler and adjusting Urukul output power.

        #todo add parameter descriptions in this docstring

        :param experiment: a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        """
        self.exp = experiment
        self.n_iterations = iters # number of times to adjust dds power per run() call
        self.dds_names = dds_names # the dds channels for the AOMs to stabilize
        self.sampler_name = sampler_name # only one sampler allowed for now
        self.n_channels = len(self.dds_names)
        self.amplitudes = [0.0]*len(self.dds_names)
        self.frequencies = [0.0]*len(self.dds_names)
        self.dds_list = [] # will hold list of dds objects
        self.sampler_channels = sampler_channels
        self.xfer_funcs = transfer_functions # translates a measured voltage to the setpoint units
        self.p = proportionals # the proportionality coeffs for feedback
        self.setpoints = setpoints # one for each channel
        self.t_meas_delay = t_meas_delay
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

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def run(self):
        """
        Run the feedback loop
        :return:
        """

        # self.exp.core.break_realtime()

        # start_mu = now_mu()

        # in each iteration, measure the sampler and set the dds amplitude accordingly
        for i in range(self.n_iterations):
            # at_mu(i * self.exp.core.seconds_to_mu(10 * ms) + start_mu)
            with sequential:
                self.sampler.sample(self.sample_buffer)
                delay(10*ms)
                with parallel:  # update all the dds channels simultaneously
                    for ch in range(self.n_channels):

                        measured_power = self.xfer_funcs[ch](self.sample_buffer[self.sampler_channels[ch]])
                        # self.print(self.sample_buffer[self.sampler_channels[ch]])
                        # self.print("mW:")
                        # self.print(measured_power)
                        err = self.setpoints[ch] - measured_power
                        ampl = self.amplitudes[ch] + self.p[ch]*err

                        # todo print some warning to alert the user if we couldn't reach the setpt,
                        if ampl < 0:
                            self.amplitudes[ch] = 0.0
                        elif ampl > 0.9:
                            self.amplitudes[ch] = 0.9
                        else:  # valid amplitude for DDS
                            self.amplitudes[ch] = ampl

                        # update the dds amplitude
                        self.dds_list[ch].set(self.frequencies[ch], amplitude=self.amplitudes[ch])
                        rtio_log(self.dds_names[ch], "i", i)
            delay(self.t_meas_delay)  # the smaller this delay, the more quickly we run into an RTIO underflow error


class AOMPowerStabilizerTest(EnvExperiment):

    def build(self):
        self.setattr_device("core")

        # initialize the hardware that we know the subroutine will use, plus any others we want
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("sampler0")

        self.setattr_argument("run_with_fake_drift", BooleanValue(False))

    def prepare(self):

        # todo: eventually read conversion functions such as this from a config file
        def volts_to_optical_mW(x: TFloat) -> TFloat:
            """
            the conversion of PD voltage to cooling light power at the switchyard MOT 1 path
            """
            x += 0.011  # this accounts for a mismatch between what the Sampler reads and what
            # the multimeter that I used for the fit reads
            return -0.195395 + 17.9214 * x

        self.AOMservo = AOMPowerStabilizer(experiment=self,
                                           dds_names=["urukul0_ch1"],
                                           sampler_name="sampler0",
                                           sampler_channels=[7],
                                           transfer_functions=[volts_to_optical_mW],
                                           setpoints=[0.7],  # in mW
                                           proportionals=[0.07],
                                           iters=5, # if > x you'll underflow the rtio counter
                                           t_meas_delay=20*ms)

        # the cooling double pass AOM - this is the laser power stabilizer
        dBm = -4 # corresponds to about 70% of the max diff. eff., so we can increase power to compensate drift
        self.dds1_ampl = math.sqrt(2 * 50 * 10 ** (dBm / 10 - 3))
        self.dds1_freq = 111 * MHz

        # the cooling single pass AOM -  we'll use this to fake a power drift
        dBm = 1
        self.dds2_ampl = math.sqrt(2 * 50 * 10 ** (dBm / 10 - 3))
        self.dds2_freq = 130 * MHz

    @kernel
    def run(self):
        self.core.reset()
        self.core.break_realtime()
        self.urukul0_cpld.init()
        self.urukul0_ch2.init()  # this drives an AOM upstream (cooling single pass)
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl)
        self.urukul0_ch2.sw.on()
        self.urukul0_ch1.init()  # the dds we're going to feedback on (cooling double pass)
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch1.set(self.dds1_freq, amplitude=self.dds1_ampl)
        self.urukul0_ch1.sw.on()
        self.sampler0.init()

        self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set

        # allow AOMs to thermalize
        delay(2000*ms)

        self.core.break_realtime()

        if self.run_with_fake_drift:
            # change the upstream AOM's RF to simulate drift from the set point
            drift = 0.04  # ~ 2 dBm change
            start_mu = now_mu()
            at_mu(start_mu)
            self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl - drift)

            # allow the upstream AOM to re-thermalize, then run the feedback
            print("faking power drift with upstream AOM")

            delay(1000 * ms)
            print("running feedback")
            self.AOMservo.run()

            # reset the upstream AOM
            print("resetting the upstream AOM")
            self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl)

            delay(2000 * ms)
            # adjust the stabilizer AOM again
            print("running feedback")
            self.AOMservo.run()
        else:
            print("running feedback")
            self.AOMservo.run()
        print("experiment finished")
