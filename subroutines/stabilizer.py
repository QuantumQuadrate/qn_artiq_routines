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

# for the test experiment
from utilities.BaseExperiment import BaseExperiment

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
        'g(x: float)->float': a callable that converts the measured voltage to the units of the setpoint, e.g. mW
        'setpoint': the float value we are trying to reach
        'p': the proportional constant for feedback
        """
        self.name = name
        self.dds_obj = dds_obj
        self.sampler_ch = sampler_ch
        self.g = g
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

        measured = self.g(buffer[self.sampler_ch]) # the value of the sampler ch in setpoint units
        err = self.setpoint - measured
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
        #       sampler: self.sampler0,
        #       channels: [FeedbackChannel1, FeedbackChannel2, ...]
        #   },
        #   {
        #       sampler: self.sampler1,
        #       channels: [FeedbackChannel1, FeedbackChannel2, ...]
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

class AOMPowerStabilizerTest(EnvExperiment):
    """
    this experiment assumes that the AOMs we'll feed back to have already been
    set to operate at ~ 70% efficiency by appropriate choice of default powers
    in ExperimentVariables.
    """

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_with_fake_drift", BooleanValue(False),"simulate drift")
        self.setattr_argument("drift_factor", NumberValue(0.9),"simulate drift")
        self.setattr_argument("experiment_iterations", NumberValue(20, type='int', ndecimals=0, scale=1, step=1))
        self.setattr_argument("feedback_dds_list",
                              StringValue("['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', 'dds_AOM_A4']"))
        self.setattr_argument("feedback_iterations", NumberValue(5, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")
        self.setattr_argument("t_measurement_delay", NumberValue(20*ms, unit='ms'),
                              "feedback params")

    def prepare(self):

        self.base.prepare()

        self.laser_stabilizer = AOMPowerStabilizer(experiment=self,
                                           dds_names=self.feedback_dds_list,
                                           iterations=self.feedback_iterations,
                                           t_meas_delay=self.t_measurement_delay)

        # the cooling single pass AOM -  we'll use this to fake a power drift.
        # this will suffice for feeding back to the cooling DP and the fiber AOMs,
        # but not for other beams.
        self.dds_drift = self.dds_cooling_SP
        self.dds_default_ampl = self.ampl_cooling_SP
        self.dds_default_freq = self.f_cooling_SP
        self.dds_drift_ampl = self.ampl_cooling_SP*self.drift_factor

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # allow AOMs to thermalize
        delay(2000*ms)

        self.core.break_realtime()

        if self.run_with_fake_drift:

            for i in range(self.experiment_iterations):

                # change the upstream AOM's RF to simulate drift from the set point
                start_mu = now_mu()
                at_mu(start_mu)
                self.dds_drift.set(self.dds_default_freq, amplitude=self.dds_drift_ampl)

                # allow the upstream AOM to re-thermalize, then run the feedback
                print("faking power drift with upstream AOM")

                delay(1000 * ms)
                print("running feedback")
                self.laser_stabilizer.run()

                # reset the upstream AOM
                print("resetting the upstream AOM")
                self.dds_drift.set(self.dds_drift_freq, amplitude=self.dds_drift_ampl)

                delay(2000 * ms)
                # adjust the stabilizer AOM again
                print("running feedback")
                self.laser_stabilizer.run()
        else:
            for i in range(self.experiment_iterations):
                print("running feedback")
                self.laser_stabilizer.run()
        print("experiment finished")
