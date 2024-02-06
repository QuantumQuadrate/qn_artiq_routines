"""
A subroutine feedback loop using Urukul channels and Samplers.
While ARTIQ has a Servo instrument for doing feedback, it does not allow
using the Urukul and Sampler channels for other purposes. I.e. it is for
continuous noise-eating. This code is for adjusting DDS power only when
we want to, and doesn't prevent us from using some channels on a given
Sampler for other purposes.

Intended usage:
1. Instantiate the servo in the prepare method of the parent artiq experiment (e.g. the BaseExperiment)
2. Call run() during the parent experiment loop whenever it is desired to tune up
the AOM powers, or monitor() to, well, you get it.

See tests/AOMFeedbackTest for example usage.

Preston's notes:
todo: should set the AOMs to default frequency setting internally, rather than rely on the user to remember to do this
    in their code, as long as we also undo our setting change at the end. or, promise nothing, and do everything
    externally, i.e. in the experiment code.
todo: should be able to modify setpoints so we can optimize them as needed.
todo: add an attribute to AOMPowerStabilizers which is a list of AOMs we want to leave on after the feedback is run,
 if any. this list can then be overwritten in prepare at the experiment level, even if we give a default list in
 Base experiment
todo: make number of measurements to average specific to each feedback channel,
 because some detectors give better SNR (e.g. the I2Vs compared to the fW detector)
"""

from artiq.experiment import *
import logging
import numpy as np
import time

import sys
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.conversions import dB_to_V

"""
Group all dds channels and feedback params by sampler card

Notes about entries below. 
- All dds channel names should follow the convention "dds_description". 
THIS IS ASSUMED BY THE CODE IN AOMPowerStabilizer.
"""
stabilizer_dict = {
    'sampler0':
        {
            'dds_cooling_DP': # signal monitored by PD0
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD0_AOM_cooling_DP', # volts
                    'p': 0.00, # the proportionality constant
                    'i': 0.0, # the integral coefficient
                    'series': False, # if series = True then these channels are fed-back to one at a time
                    'dataset':'MOT_switchyard_monitor',
                    'power_dataset':'p_cooling_DP_MOT', # the dataset for the RF power in dB; see ExperimentVariables
                    # 'amplitudes': 'DummyAmplitudes',
                    't_measure_delay':1*ms, # time to wait between AOM turned on and measurement
                    'max_dB': 0
                },
            'dds_AOM_A1': # signal monitored by Thorlabs fW detector
                {
                    'sampler_ch': 7, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD1_AOM_A1', # volts
                    'p': 0.001, # the proportionality constant
                    'i': 0.000, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT1_monitor',
                    'power_dataset':'p_AOM_A1',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },
            'dds_AOM_A2': # signal monitored by Thorlabs fW detector
                {
                    'sampler_ch': 5, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD2_AOM_A2', # volts
                    'p': 0.1, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT2_monitor',
                    'power_dataset':'p_AOM_A2',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },
            'dds_AOM_A3': # signal monitored by Thorlabs fW detector
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD3_AOM_A3', # volts
                    'p': 0.2, # the proportionality constant
                    'i': 0.000, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT3_monitor',
                    'power_dataset':'p_AOM_A3',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },
            'dds_AOM_A4': # signal monitored by Thorlabs fW detector
                {
                    'sampler_ch': 4, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD4_AOM_A4', # volts
                    'p': 0.2, # the proportionality constant
                    'i': 0.0, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT4_monitor',
                    'power_dataset':'p_AOM_A4',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },
            'dds_AOM_A5': # signal monitored by PD5
                {
                    'sampler_ch': 1, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD5_AOM_A5',
                    'p': 0.08, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset':'MOT5_monitor',
                    'power_dataset':'p_AOM_A5',
                    # 'amplitudes': 'DummyAmplitudes',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },
            'dds_AOM_A6': # signal monitored by PD6
                {
                    'sampler_ch': 2, # the channel connected to the appropriate PD
                    'set_point': 'set_point_PD6_AOM_A6', # volts
                    'p': 0.1, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT6_monitor',
                    'power_dataset':'p_AOM_A6',
                    't_measure_delay':1*ms,
                    'max_dB': 0
                },

            # todo: the FORT should ideally have a dedicated feedback stage, since we will need to first
            #  feedback to motorized waveplates and then feedback to the AOM.
            'dds_FORT': # signal monitored by TTI detector connected to MM fiber
                {
                    'sampler_ch': 6, # the channel connected to the appropriate PD
                    'set_point': 'set_point_FORT_MM',
                    'p': 0.7, #
                    'i': 0.0, # the integral coefficient
                    'series': True, # setting to True because there's a bug with parallel
                    'dataset':'FORT_monitor',
                    'power_dataset':'p_FORT_loading',
                    't_measure_delay':1*ms, # time to wait between AOM turned on and measurement
                    'max_dB': 5
                },
        },
    'sampler1':
        {
        }
}


class FeedbackChannel:

    def __init__(self, name, dds_obj, buffer_index, set_point, p, i, dataset, dB_dataset, t_measure_delay,
                 error_history_length, max_dB):
        """
        class which defines a DDS feedback channel

        parameters:
        'name': the name of the dds channel, e.g. dds_AOM_A1
        'dds_obj': the experiment reference to the dds object with name 'name', e.g. experiment.dds_AOM_A1
        'g(x: float)->float': a callable that converts the measured voltage to the units of the setpoint, e.g. mW
        'set_point': the float value we are trying to reach
        'p': the proportional constant for feedback
        'error_history_length->int':, how many errors to store for the integral term
        """
        self.name = name
        self.dds_obj = dds_obj
        self.buffer_index = buffer_index
        # self.g = g
        self.set_point = set_point
        self.p = p
        self.i = i
        self.frequency = 100*MHz
        self.amplitude = 0.0
        self.value = 0.0 # the last value of the measurement
        self.value_normalized = 0.0 # self.value normalized to the set point
        self.dataset = dataset
        self.dB_dataset = dB_dataset # the name of the dataset that stores the dB RF power for the dds
        self.t_measure_delay = t_measure_delay
        self.error_history_arr = np.full(error_history_length,0.0)
        self.error_buffer = np.full(error_history_length-1,0.0)
        self.error_history_length = error_history_length
        self.cumulative_error = 0.0 # will store a sum of the
        self.max_dB = max_dB
        self.ampl_default = 0.0

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def get_dds_settings(self):
        """ get the frequency and amplitude """
        self.frequency, _, self.amplitude = self.dds_obj.get()

    @kernel
    def set_value(self, value):
        self.value = value
        self.value_normalized = value / self.set_point

    @rpc(flags={'async'})
    def update_error_history(self, error):
        """
        update the error history array.

        this is run off the kernel because the array is not updated correctly on the kernel
        :param error:
        :return:
        """
        self.error_buffer = self.error_history_arr[-(self.error_history_length - 1):]
        self.error_history_arr[:-1] = self.error_buffer
        self.error_history_arr[-1] = error
        self.cumulative_error = sum(self.error_history_arr)

    @kernel
    def feedback(self, buffer):
        """ feedback to this channel
        buffer: a list of length storing the measurement result from all samplers
        """

        # measured = self.g(buffer[self.buffer_index]) # the value of the sampler ch in setpoint units
        measured = buffer[self.buffer_index]
        err = self.set_point - measured

        # runs off the kernel, else the error array isn't correctly updated
        self.update_error_history(err)

        ampl = self.amplitude + self.p * err + self.i * self.cumulative_error

        self.set_value(measured)
        max_ampl = (2 * 50 * 10 ** (self.max_dB / 10 - 3)) ** (1 / 2)

        # todo print some warning to alert the user if we couldn't reach the setpoint,
        if ampl < 0:
            self.amplitude = 0.0
        elif ampl > max_ampl:
            self.amplitude = max_ampl
        else:  # valid amplitude for DDS
            self.amplitude = ampl

        # update the dds amplitude
        self.dds_obj.set(frequency=self.frequency,amplitude=self.amplitude)


class AOMPowerStabilizer:

    def __init__(self, experiment, dds_names, iterations=10, averages=1, leave_AOMs_on=False,
                 update_dds_settings=True, dry_run=False, open_loop_monitor_names=[]):
        """
        An experiment subsequence for reading a Sampler and adjusting Urukul output power.

        parameters:
        'experiment': a class inheriting from artiq.experiment.EnvExperiment. It is
        assumed that any hardware objects that are referenced in this class are bound to
        this instance of experiment.
        'dds_names': a list of the names of the dds channels to feedback to
        'iterations': integer number of feedback cycles to converge to the setpoints
        'averages':
        'leave_AOMs_on':
        'udpate_dds_settings':
        'dry_run':
        """

        # initialized by user
        self.exp = experiment
        self.iterations = iterations # number of times to adjust dds power per run() call
        self.dds_names = dds_names # the dds channels for the AOMs to stabilize
        self.averages = averages
        self.leave_AOMs_on = leave_AOMs_on
        self.update_dds_settings = update_dds_settings
        self.dry_run = dry_run
        self.open_loop_monitor_names = open_loop_monitor_names
        if len(open_loop_monitor_names) > 0:
            assert [x in self.dds_names for x in self.open_loop_monitor_names], \
                f"{[x for x in self.open_loop_monitor_names if x not in self.dds_names]} are not in dds_names." \
                +"Only dds channels in dds_names can be include in open_loop_monitor_channels"

        self.sampler_list = [] # stores the list of Sampler references
        self.series_channels = [] # the feedback channels that should be adjusted one-by-one
        self.parallel_channels = []
        self.open_loop_monitor_channels = []

        # this block instantiates FeedbackChannel objects and groups them by series/parallel feedback flag
        for i, sampler_name in enumerate(stabilizer_dict):

            # the dds channel names associated with this sampler
            feedback_channels = stabilizer_dict[sampler_name]

            # check if any of the dds names is associated with this sampler
            if True in [dds_name in feedback_channels.keys() for dds_name in dds_names]:

                self.sampler_list.append(getattr(self.exp, sampler_name))

                # loop over the dds channels associated with this sampler
                for ch_name in feedback_channels.keys():

                    # create a feedback object for dds channels we want to feed back to
                    if ch_name in dds_names:

                        ch_params = feedback_channels[ch_name]

                        fb_channel = FeedbackChannel(name=ch_name,
                                            dds_obj=getattr(self.exp, ch_name),
                                            buffer_index=ch_params['sampler_ch'] + i*8,
                                            # g=ch_params['transfer_function'],
                                            set_point=getattr(self.exp,ch_params['set_point']),
                                            p=ch_params['p'],
                                            i=ch_params['i'],
                                            dataset=ch_params['dataset'],
                                            dB_dataset=ch_params['power_dataset'],
                                            t_measure_delay=ch_params['t_measure_delay'],
                                            error_history_length=self.iterations,
                                            max_dB=ch_params['max_dB']
                        )

                        # make the feedback channel an attribute of the AOMPowerStabilizer instance itself
                        # so we can conveniently access the updated dds amplitudes in our experiments
                        stabilizer_ch_name = "stabilizer"+ch_name[3:]
                        setattr(self.exp, stabilizer_ch_name, fb_channel)
                        # fb_channel_ref = getattr(self.exp, stabilizer_ch_name) # use this going forward

                        if ch_params['series']:
                            # self.series_channels.append(fb_channel)
                            self.series_channels.append(getattr(self.exp, stabilizer_ch_name))
                        else:
                            # self.parallel_channels.append(fb_channel)
                            self.parallel_channels.append(getattr(self.exp, stabilizer_ch_name))

                        if ch_name in self.open_loop_monitor_names:
                            self.open_loop_monitor_channels.append(getattr(self.exp, stabilizer_ch_name))

        self.num_samplers = len(self.sampler_list)
        self.measurement_array = np.full(8 * self.num_samplers, 0.0)
        self.background_array = np.full(8 * self.num_samplers, 0.0)
        self.signal_array = np.zeros(8)
        self.bg_array = np.zeros(8)

        self.all_channels = self.parallel_channels + self.series_channels

        # for logging the measured voltages
        for ch in self.all_channels: # todo: update with the last value from the dataset
            try:
                self.exp.set_dataset(ch.dataset, [self.exp.get_dataset(ch.dataset)[-1]], broadcast=True)
            except Exception as e:
                logging.warning(e)
                self.exp.set_dataset(ch.dataset, [1.0], broadcast=True)

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def write_dds_settings(self):
        for ch in self.all_channels:
            dB = 10*(np.log10(ch.amplitude**2/(2*50)) + 3)
            self.exp.set_dataset(ch.dB_dataset, dB, broadcast=True, persist=True)

            # if ch.amplitude_setter != None:
            # ch.amplitude_setter.ampl_default = ch.amplitude
            # self.exp.dds_profiles[ch.name]["default"] = ch.amplitude
            # self.exp.__setattr__("test_amplitude",ch.amplitude)

    @kernel
    def measure(self):
        """
        measure all Sampler cards used for feedback
        # todo: could add option to average for better snr
        """
        dummy = np.full(8 * self.num_samplers, 0.0)
        measurement_array = np.full(8, 0.0)
        i = 0
        for sampler in self.sampler_list:
            for j in range(self.averages):
                sampler.sample(measurement_array)
                delay(1*ms)
                dummy[i * 8:(i + 1) * 8] += measurement_array
                delay(1 * ms)
            i += 1
        dummy /= self.averages
        self.measurement_array = dummy

    @kernel
    def measure_background(self):
        """
        measure all Sampler cards used for feedback
        # todo: could add option to average for better snr
        """

        dummy_background = np.full(8 * self.num_samplers, 0.0)
        background_array = np.full(8, 0.0)
        i = 0
        for sampler in self.sampler_list:
            for j in range(self.averages):
                sampler.sample(background_array)
                delay(1 * ms)
                # self.measurement_array[i * 8:(i + 1) * 8] += measurement_array/self.averages # outlives target error
                dummy_background[i * 8:(i + 1) * 8] += background_array  # outlives target error
                delay(1 * ms)
            i += 1
        # self.measurement_array /= self.averages
        dummy_background /= self.averages
        self.background_array = dummy_background

        self.print("self.background_array:")
        self.print(self.background_array)

    @kernel
    def open_loop_monitor(self):
        """
        make a measurement and update the datasets.

        for now, this intentionally does not respect the series flag of dds channels,
        so it should be used with only 1-series AOM at a time.
        :return:
        """
        for ch in self.open_loop_monitor_channels:
            self.measure()
            delay(1*ms)
            ch.set_value((self.measurement_array - self.background_array)[ch.buffer_index])
            self.exp.append_to_dataset(ch.dataset, ch.value_normalized)

    @kernel
    def monitor(self):
        """
        Functionally the same as run, but measurement only.
        Intended to be used to monitor beam drift in absence of feedback.
        """

        self.run(monitor_only=True)

    @kernel
    def run(self, record_all_measurements=False, monitor_only=False):
        """
        Run the feedback loop. On exiting, this function will turn off all dds channels
         given by dds_names. If any beams which need to be adjusted in series are in
         dds_names, all such channels will be turned off, i.e. even ones we are not
         feeding back to.

        :param record_all_measurements: optional, False by default. if True, every measurement point will
        be posted to the corresponding dataset, else, only post the final measurement, i.e. after the feedback
        loop has been run.
        """

        self.exp.ttl7.pulse(1*ms) # scope trigger

        # todo: up-date the set points for all channels in case they have been changed since
        #  the stabilizer was instantiated. not sure how to do this since getattr can not be
        #  used in the kernel, and I don't want to have to pass in a variable explicitly.

        # todo: modify this to only affect the DDSs we're feeding back to?
        # make sure that all of the DDSs are set the default frequency and power levels
        self.exp.named_devices.set_dds_default_settings()

        # turn off the repumps which are mixed into the cooling light
        self.exp.ttl_repump_switch.on() # block RF to the RP AOM
        # todo include pumping repump

        for ch in self.all_channels:
            ch.get_dds_settings()
            delay(1*ms)
            ch.dds_obj.sw.off()
            delay(1*ms)

        with sequential:

            # don't blind the fW detector
            self.exp.dds_FORT.sw.off()

            for ch in self.series_channels:
                ch.dds_obj.sw.off()
            delay(1 * ms)

            delay(1*ms)
            for ch in self.parallel_channels:
                ch.dds_obj.sw.on()
            delay(10*ms)

            # do feedback on the "parallel" channels
            for i in range(self.iterations):
                self.measure()
                delay(1 * ms)

                # strictly speaking, doesn't really matter if these are parallel
                for ch in self.parallel_channels:
                    delay(1*ms)
                    if not (self.dry_run or monitor_only):
                        ch.feedback(self.measurement_array - self.background_array)
                    else:
                        ch.set_value((self.measurement_array - self.background_array)[ch.buffer_index])

                    if record_all_measurements:
                        self.exp.append_to_dataset(ch.dataset, ch.value_normalized)

                delay(1 * ms)

            for ch in self.parallel_channels:
                ch.dds_obj.sw.off()
            delay(10 * ms)  # the Femto fW detector is slow

            # need to have this on
            self.exp.dds_cooling_DP.sw.on()
            delay(50 * ms)

            # self.measure_background()
            delay(1*ms)

            # do feedback on the series channels
            with sequential:
                for ch in self.series_channels:
                    ch.dds_obj.sw.on()

                    delay(ch.t_measure_delay) # allows for detector rise time
                    for i in range(self.iterations):
                        self.measure()
                        delay(1 * ms)

                        if not (self.dry_run or monitor_only):
                            ch.feedback(self.measurement_array - self.background_array)

                        else:
                            ch.set_value((self.measurement_array - self.background_array)[ch.buffer_index])
                        if record_all_measurements:
                            self.exp.append_to_dataset(ch.dataset, ch.value_normalized)

                    ## trigger for Andor Luca camera for independent verification of the measured signals
                    if self.exp.Luca_trigger_for_feedback_verification:
                        self.exp.ttl6.pulse(5 * ms)
                        delay(60*ms)
                    delay(1*ms)
                    ch.dds_obj.sw.off()

            # update the datasets with the last values if we have not already done so
            if not record_all_measurements:
                for ch in self.all_channels:
                    self.exp.append_to_dataset(ch.dataset, ch.value_normalized)

        delay(1*ms)
        # turn cooling DP back on
        self.exp.dds_cooling_DP.sw.on()

        delay(1*ms)

        if self.leave_AOMs_on:
            for ch in self.all_channels:
                ch.dds_obj.sw.on()

        delay(1 * ms)
        if self.update_dds_settings:
            self.write_dds_settings()

        # turn on the repumps which are mixed into the cooling light
        self.exp.ttl_repump_switch.on()  # block RF to the RP AOM
        # todo include pumping repump
