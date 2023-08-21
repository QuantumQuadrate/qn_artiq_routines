"""
# todo: this should eventually replace stabilizer.py

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

Preston's notes:
todo: should be able to modify setpoints so we can optimize them as needed.
one way to do this would be to make the stabilizer dict a json file and
update the json file, e.g. after doing an optimizer routine with the AOM
powers. then whatever the voltages are after doing this routine, we update
the file. however, it might be preferable to scan over the setpoints
themselves (with proper caps on the AOM RF power), in which case, they
should be experiment variables. not obvious how to do this. the MOT can
likely be moved only using the shim fields, but the beam pairs themselves
will need to be precisely balanced in order to do good PGC. so adjusting
these setpoints will be a necessity in a matter of time.
"""

from artiq.experiment import *
import math

# group all dds channels and feedback params by sampler card
stabilizer_dict = {
    'sampler0':
        {
            'dds_cooling_PD': # signal monitored by PD0
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x, # converts volts to optical mW
                    'set_point': 6, # value in mW,
                    'p': 0.07, # the proportionality constant
                    'series': False # if series = True then these channels are fed-back to one at a time
                },
            'dds_AOM_A5': # signal monitored by PD5
                {
                    'sampler_ch': 1, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,
                    'set_point': 6, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                },
            'dds_AOM_A6': # signal monitored by PD6
                {
                    'sampler_ch': 2, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'set_point': 6, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                },
            'dds_AOM_A1': # signal monitored by Femto fW detector
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,
                    'set_point': 0.922, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                },
            'dds_AOM_A2': # signal monitored by Femto fW detector
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  # arbitrary units
                    'set_point': 0.904, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                },
            'dds_AOM_A3': # signal monitored by Femto fW detector
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,
                    'set_point': 1.436, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                },
            'dds_AOM_A4': # signal monitored by Femto fW detector
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'transfer_function': lambda x : x,  
                    'set_point': 844, # volts wrt to background
                    'p': 0.07, # the proportionality constant
                    'series': True
                }
        },
    'sampler1':
        {
        }
}

class FeedbackChannel:

    def __init__(self, name, dds_obj, buffer_index, g, set_point, p):
        """
        class which defines a DDS feedback channel

        parameters:
        'name': the name of the dds channel, e.g. dds_AOM_A1
        'dds_obj': the experiment reference to the dds object with name 'name', e.g. experiment.dds_AOM_A1
        'g(x: float)->float': a callable that converts the measured voltage to the units of the setpoint, e.g. mW
        'set_point': the float value we are trying to reach
        'p': the proportional constant for feedback
        """
        self.name = name
        self.dds_obj = dds_obj
        self.buffer_index = buffer_index
        self.g = g
        self.set_point = set_point
        self.p = p
        self.frequency = 100*MHz
        self.amplitude = 0.0

    @kernel
    def get_dds_settings(self):
        """ get the frequency and amplitude """
        self.frequency, _, self.amplitude = self.dds_obj.get()

    @kernel
    def feedback(self, buffer):
        """ feedback to this channel
        buffer: a list of length storing the measurement result from all samplers
        """

        measured = self.g(buffer[self.buffer_index]) # the value of the sampler ch in setpoint units
        err = self.set_point - measured
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


class AOMPowerStabilizer2:

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
        self.measurement = [] # a buffer to store the measurement from all sampled Samplers
        self.sampler_list = [] # stores the list of Sampler references
        self.series_channels = [] # the feedback channels that should be adjusted one-by-one
        self.parallel_channels = []

        # this block instantiates FeedbackChannel objects and groups them by series/parallel feedback flag
        for i, sampler_name in enumerate(stabilizer_dict):

            # the dds channel names associated with this sampler
            feedback_channels = stabilizer_dict[sampler_name]

            # check if any of the dds names is associated with this sampler
            if True in [dds_name in feedback_channels.keys() for dds_name in dds_names]:

                self.measurement += [0.0] * 8
                self.sampler_list.append(getattr(self.exp, sampler_name))

                # loop over the dds channels associated with this sampler
                for ch_name in feedback_channels.keys():

                    # create a feedback object for dds channels we want to feed back to
                    if ch_name in dds_names:

                        ch_params = feedback_channels[ch_name]

                        fb_channel = FeedbackChannel(name=ch_name,
                                            dds_obj=getattr(self.exp, ch_name),
                                            buffer_index=ch_params['sampler_ch'] + i*8,
                                            g=ch_params['transfer_function'],
                                            set_point=ch_params['set_point'],
                                            p=ch_params['p'],
                        )

                        if ch_params['series']:
                            self.series_channels.append(fb_channel)
                        else:
                            self.parallel_channels.append(fb_channel)

        self.background = self.measurement

        self.all_channels = self.series_channels + self.series_channels

    # create a dictionary of samplers and feedback objects with a structure
        # similar to stabilizer_dict
        # self.channel_sampler_list: list[dict] = [] #= [{}]*len(stabilizer_dict.keys())
        #
        # # todo: add a dataset to be broadcast for each dds channel so we can plot the powers
        #
        # for i,sampler_name in enumerate(stabilizer_dict):
        #
        #     # the dds channel names associated with this sampler
        #     feedback_channels = stabilizer_dict[sampler_name]
        #
        #     # check if any of the dds names is associated with this sampler
        #     if True in [dds_name in feedback_channels.keys() for dds_name in dds_names]:
        #         self.channel_sampler_list.append({'sampler': getattr(self.exp, sampler_name),
        #                             'channels': []})
        #
        #         # loop over the dds channels associated with this sampler
        #         for dds_name in feedback_channels.keys():
        #
        #             # create a feedback object for dds channels we want to feed back to
        #             if dds_name in dds_names:
        #                 self.channel_sampler_list[i]['channels'].append(
        #                     FeedbackChannel(name=dds_name,
        #                                     dds_obj=getattr(self.exp, dds_name),
        #                                     sampler_ch=feedback_channels[dds_name]['sampler_ch'],
        #                                     g=feedback_channels[dds_name]['transfer_function'],
        #                                     set_point=feedback_channels[dds_name]['set_point'],
        #                                     p=feedback_channels[dds_name]['p'])
        #                 )

        # the list of dds channels which, if included, should be adjusted in series
        # self.dds_series_list = [getattr(self.exp, dds) for dds in dds_series_list]

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def measure(self):
        """
        measure all Sampler cards used for feedback
        # todo: could add option to average for higher-res
        """
        i = 0
        for sampler in self.sampler_list:
            sampler.sample(self.measurement[i:i + 8])
            i+=1

    def measure_background(self):
        """
        measure all Sampler cards used for feedback
        # todo: could add option to average for higher-res
        """
        i = 0
        for sampler in self.sampler_list:
            sampler.sample(self.measurement[i:i + 8])
            i+=1
        self.background = self.measurement

    @kernel
    def run(self):
        """
        Run the feedback loop. On exiting, this function will turn off all dds channels
         given by dds_names. If any beams which need to be adjusted in series are in
         dds_names, all such channels will be turned off, i.e. even ones we are not
         feeding back to.
        """

        for ch in self.all_channels:
            ch.get_dds_settings()
            ch.dds_obj.sw.off()

        with sequential:

            self.measure_background() # this updates the background list
            for ch in self.parallel_channels:
                ch.dds_obj.sw.on()
            delay(1*ms)

            # do feedback on the parallel channels
            for i in range(self.iterations):
                self.measure()
                delay(1 * ms)
                with parallel: # strictly speaking, doesn't really matter if these are parallel
                    for ch in self.parallel_channels:
                        pass
                        # ch.feedback(self.measurement - self.background)
                delay(1 * ms)

            for ch in self.parallel_channels:
                ch.dds_obj.sw.off()
            delay(1*ms)

            self.measure_background()

            # do feedback on the series channels
            with sequential:
                for ch in self.series_channels:
                    ch.dds_obj.sw.on()
                    delay(1 * ms)
                    for i in range(self.iterations):
                        self.measure()
                        delay(1*ms)
                        # ch.feedback(self.measurement - self.background)
                        delay(1*ms)
                    ch.dds_obj.sw.off()