"""
A subroutine feedback loop using Urukul channels and Samplers.
While ARTIQ has a Servo instrument for doing feedback, it does not allow
using the Urukul and Sampler channels for other purposes. I.e. it is for
continuous noise-eating. This code is for adjusting DDS power only when
we want to, and doesn't prevent us from using some channels on a given
Sampler for other purposes.

Intended usage:
1. Instantiate the servo in the prepare method of the parent artiq experiment
2. Call run() during the parent experiment loop whenever it is desired to tune up
the AOM powers.

See the example experiment toward the bottom of this file.
"""

from artiq.experiment import *
import math

# the dds channels for on-chip beams. these have to be fed-back to sequentially
dds_on_chip_list = ['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', 'dds_AOM_A4']#,
                   # 'dds_excitation_AL','dds_excitation_AR','dds_OP_AL','dds_OP_AR']

# group all dds channels and feedback params by sampler card
stabilizer_dict = {
    'sampler0':
        {
            'dds_cooling_PD':
                {
                    'sampler_ch': 7, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x, # converts volts to optical mW
                    'setpoint': 6, # value in mW,
                    'p': 0.07 # the proportionality constant
                } #,
        },
    'sampler1':
        {
            'dds_AOM_A1':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                },
            'dds_AOM_A2':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                },
            'dds_AOM_A3':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                },
            'dds_AOM_A4':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                },
            'dds_AOM_A5':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                },
            'dds_AOM_A6':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'setpoint': 6, # arbitrary units
                    'p': 0.07 # the proportionality constant
                }
        }
}

class FeedbackChannel:

    def __init__(self, name, dds_obj, sampler_ch, g, setpoint, p):
        """
        class which defines a DDS feedback channel

        parameters:
        'name': the name of the dds channel, e.g. dds_AOM_A1
        'dds_obj': the experiment reference to the dds object with name 'name', e.g. experiment.dds_AOM_A1
        'g': a callable which takes a float voltage and outputs a float representing mW optical power
        'setpoint': the float value in units mW that we are trying to reach
        'p': the proportional constant for feedback
        """
        self.name = name
        self.dds_obj = dds_obj
        self.sampler_ch = sampler_ch
        self.g = g # a transfer function. input V output mW
        self.setpoint = setpoint
        self.p = p
        self.frequency = 100*MHz
        self.amplitude = 0

    @kernel
    def get_dds_settings(self):
        """ get the frequency and amplitude """
        self.frequency, _, self.amplitudes = self.dds_obj.get()

    @kernel
    def feedback(self, buffer):
        """ feedback to this channel
        buffer: a list of length 8 storing a Sampler measurement
        """

        measured_power = self.g(buffer[self.sampler_ch])
        err = self.setpoint - measured_power
        ampl = self.amplitude + self.p*err

        # todo print some warning to alert the user if we couldn't reach the setpoint,
        if ampl < 0:
            self.amplitude = 0.0
        elif ampl > 0.9:
            self.amplitude = 0.9
        else:  # valid amplitude for DDS
            self.amplitude = ampl

        # update the dds amplitude
        self.dds_obj.set(frequency=self.frequency,amplitude=self.amplitude)


class AOMPowerStabilizer:

    def __init__(self, experiment, dds_names, iterations=10, t_meas_delay=10 * ms):
        """
        An experiment subsequence for reading a Sampler and adjusting Urukul output power.

        parameters:
        'experiment': a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        'dds_names': a list of the names of the dds channels to feedback to
        'iterations': integer number of feedback cycles to converge to the setpoints
        't_meas_delay': time to wait between iterations. mostly for avoiding underflow errors.
        """

        # initialized by user
        self.exp = experiment
        self.iterations = iterations # number of times to adjust dds power per run() call
        self.dds_names = dds_names # the dds channels for the AOMs to stabilize
        self.t_meas_delay = t_meas_delay
        self.sample_buffer = [0.0] * 8 # can reuse this for each sampler

        # create a dictionary of samplers and feedback objects with a structure
        # similar to stabilizer_dict
        self.channel_sampler_list = []

        for i,sampler_name in enumerate(stabilizer_dict):

            # the dds channel names associated with this sampler
            feedback_channels = stabilizer_dict[sampler_name]

            # check if any of the dds names is associated with this sampler
            if True in [dds_name in feedback_channels.keys() for dds_name in dds_names]:
                self.channel_sampler_list[i] = {'sampler': getattr(self.exp, sampler_name),
                                    'channels': []}

                # loop over the dds channels associated with this sampler
                for dds_name in feedback_channels.keys():

                    # create a feedback object for dds channels we want to feed back to
                    if dds_name in dds_names:
                        self.channel_sampler_list[i]['channels'].append(
                            FeedbackChannel(name=dds_name,
                                            dds_obj=getattr(self.exp, dds_name),
                                            sampler_ch=feedback_channels[dds_name]['sampler_ch'],
                                            g=feedback_channels[dds_name]['transfer_function'],
                                            setpoint=feedback_channels[dds_name]['setpoint'],
                                            p=feedback_channels[dds_name]['p'])
                        )

        self.on_chip_dds_list = [getattr(self.exp, dds) for dds in dds_on_chip_list]

        # channel_sampler_list looks like this:
        # [
        #   {
        #   sampler: self.sampler0,
        #   channels: [FeedbackChannel1, FeedbackChannel12, ...]
        #   },
        #   {
        #   sampler: self.sampler1,
        #   channels: [FeedbackChannel1, FeedbackChannel12, ...]
        #   }
        # ]

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def run(self):
        """
        Run the feedback loop. On exiting, this function will turn off all dds channels
         given by dds_names. If any on-chip beam dds channels are in dds_names, ALL on-chip
         dds channels will be turned off.
        :return:
        """

        # get the current frequency and amplitude for every dds channel
        # and turn each dds off
        for sampler_and_channels in self.channel_sampler_list:
            for ch in sampler_and_channels['channels']:
                ch.get_dds_settings()
                ch.dds_obj.sw.off()

        # loop over feedback iterations
        for i in range(self.iterations):

            # do feedback by sampler
            for sampler_and_channels in self.channel_sampler_list:
                sampler = sampler_and_channels['sampler'] # the sampler object
                channels = sampler_and_channels['channels'] # the Feedback channels

                with sequential:

                    free_space_channels = [ch for ch in channels if ch.name not in dds_on_chip_list]
                    on_chip_channels = [ch for ch in channels if ch.name in dds_on_chip_list]

                    [ch.dds_obj.sw.on() for ch in free_space_channels]
                    delay(10*ms)

                    # do the feedback for free-space beams
                    for i in range(n_iterations):

                        # measure the PD voltages
                        sampler.sample(self.sample_buffer)
                        delay(10 * ms)

                        with parallel:
                            for ch in free_space_channels:
                                # adjust the dds power for ch
                                ch.feedback(self.sample_buffer)

                    [ch.dds_obj.sw.off() for ch in free_space_channels]
                    delay(1*ms)

                    # if we are feeding back to any of the on-chip beams,
                    # turn off all on-chip beams
                    if len(on_chip_channels) > 1:
                        [dds.sw.off() for dds in self.on_chip_dds_list]
                        delay(1 * ms)

                    # do the feedback for on-chip beams
                    for ch in on_chip_channels:

                        with sequential: # the on-chip beams are fed back to one-at-a-time
                            for i in range(n_iterations):
                                ch.dds_obj.sw.on()
                                delay(1 * ms)

                                # measure the PD voltage
                                sampler.sample(self.sample_buffer)

                                # adjust the dds power for ch
                                ch.feedback(self.sample_buffer)

                                ch.dds_obj.sw.off()
                                delay(self.t_meas_delay)

        # feedback sequence over. turn off all the dds channels we adjusted
        for sampler_and_channels in self.channel_sampler_list:
            for ch in sampler_and_channels['channels']:
                ch.dds_obj.sw.off()

# class AOMPowerStabilizerTest(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#
#         # initialize the hardware that we know the subroutine will use, plus any others we want
#         self.setattr_device("urukul0_cpld")
#         self.setattr_device("urukul0_ch1")
#         self.setattr_device("urukul0_ch2")
#         self.setattr_device("sampler0")
#
#         self.setattr_argument("run_with_fake_drift", BooleanValue(False))
#
#     def prepare(self):
#
#         # todo: eventually read conversion functions such as this from a config file
#         def volts_to_optical_mW(x: TFloat) -> TFloat:
#             """
#             the conversion of PD voltage to cooling light power at the switchyard MOT 1 path
#             """
#             x += 0.011  # this accounts for a mismatch between what the Sampler reads and what
#             # the multimeter that I used for the fit reads
#             return -0.195395 + 17.9214 * x
#
#         self.AOMservo = AOMPowerStabilizer(experiment=self,
#                                            dds_names=["urukul0_ch1"],
#                                            sampler_name="sampler0",
#                                            sampler_channels=[7],
#                                            transfer_functions=[volts_to_optical_mW],
#                                            setpoints=[0.7],  # in mW
#                                            proportionals=[0.07],
#                                            iterations=5,  # if > x you'll underflow the rtio counter
#                                            t_meas_delay=20*ms)
#
#         # the cooling double pass AOM - this is the laser power stabilizer
#         dBm = -4 # corresponds to about 70% of the max diff. eff., so we can increase power to compensate drift
#         self.dds1_ampl = math.sqrt(2 * 50 * 10 ** (dBm / 10 - 3))
#         self.dds1_freq = 111 * MHz
#
#         # the cooling single pass AOM -  we'll use this to fake a power drift
#         dBm = 1
#         self.dds2_ampl = math.sqrt(2 * 50 * 10 ** (dBm / 10 - 3))
#         self.dds2_freq = 130 * MHz
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.core.break_realtime()
#         self.urukul0_cpld.init()
#         self.urukul0_ch2.init()  # this drives an AOM upstream (cooling single pass)
#         self.urukul0_ch2.set_att(float(0))
#         self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl)
#         self.urukul0_ch2.sw.on()
#         self.urukul0_ch1.init()  # the dds we're going to feedback on (cooling double pass)
#         self.urukul0_ch1.set_att(float(0))
#         self.urukul0_ch1.set(self.dds1_freq, amplitude=self.dds1_ampl)
#         self.urukul0_ch1.sw.on()
#         self.sampler0.init()
#
#         self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set
#
#         # allow AOMs to thermalize
#         delay(2000*ms)
#
#         self.core.break_realtime()
#
#         if self.run_with_fake_drift:
#             # change the upstream AOM's RF to simulate drift from the set point
#             drift = 0.04  # ~ 2 dBm change
#             start_mu = now_mu()
#             at_mu(start_mu)
#             self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl - drift)
#
#             # allow the upstream AOM to re-thermalize, then run the feedback
#             print("faking power drift with upstream AOM")
#
#             delay(1000 * ms)
#             print("running feedback")
#             self.AOMservo.run()
#
#             # reset the upstream AOM
#             print("resetting the upstream AOM")
#             self.urukul0_ch2.set(self.dds2_freq, amplitude=self.dds2_ampl)
#
#             delay(2000 * ms)
#             # adjust the stabilizer AOM again
#             print("running feedback")
#             self.AOMservo.run()
#         else:
#             print("running feedback")
#             self.AOMservo.run()
#         print("experiment finished")
