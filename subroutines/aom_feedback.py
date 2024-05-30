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
the AOM powers, or monitor() to measure each channel and update the monitor datasets but without actually
feeding back to the measurement.
3. See the example code block below for a typical use example and some important things to note.

class MyStableExperiment(EnvExperiment):
    ... BaseExperiment setup happens in self.build, including instantiation of an AOMPowerStabilizer

    @kernel
    def run(self):
        base.initialize_hardware()

        for n in range(self.n_measurements):

            # feedback to adjust the AOMs. this adjusts the default powers, given in DeviceAliases.DDS_DEFAULTS.
            # e.g., the default power for dds_FORT is p_FORT_loading. the set point we try to meet refers to
            # the voltage we try to reach on the detector for this particular RF setting.
            self.laser_stabilizer.run()
            self.dds_FORT.sw.on()

            # set the amplitude of the FORT to the default level. laser_stabilizer updates the p_FORT_loading
            # dataset, but we can't get_dataset on the kernel, so use the amplitude attribute that is
            # associated with a stabilizer channel as follows. The naming convention for the stabilizer
            # channel is "stabilizer_name" if your dds is named "dds_name".
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

            # now, when if we want to do PGC, we lower the FORT by some fraction of the amplitude that
            # we set with the stabilizer. here, maybe self.p_FORT_PGC = 0.4, for example.
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_PGC*self.stabilizer_FORT.amplitude)

Preston's notes:
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

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.conversions import dB_to_V
from utilities.DeviceAliases import DDS_DEFAULTS
from utilities.helper_functions import print_async


"""
Group all dds channels and feedback params by sampler card

Notes about entries below. 
- All dds channel names should follow the convention "dds_description", and the "dds_description" must be one of the
keys in DeviceAliases.DDS_DEFAULTS. THIS IS ASSUMED BY THE CODE IN AOMPowerStabilizer.
"""
stabilizer_dict = { # todo: replace series/paralle with an index specifying the order in which to run.
                    #   this will allow for giving two channels the same priority thus allowing them to be done in
                    #   parallel or series in any combination, e.g. 1 1 1 2 3 3 4, etc. This would allow us to feedback
                    #   orthogonal MOT beams in parallel in groups, thus significantly reducing the time needed for
                    #   feedback.
    'sampler0':
        {
            'dds_cooling_DP':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD0_AOM_cooling_DP'], # volts
                    'p': 0.00, # the proportionality constant
                    'i': 0.0, # the integral coefficient
                    'series': False, # if series = True then these channels are fed-back to one at a time
                    'dataset':'MOT_switchyard_monitor',
                    'power_dataset':'p_cooling_DP_MOT', # the dataset for the RF power in dB; see ExperimentVariables
                    't_measure_delay':1*ms, # time to wait between AOM turned on and measurement
                    'max_dB': 0
                },
            'dds_AOM_A1':
                {
                    'sampler_ch': 0, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD1_AOM_A1'], # volts
                    'p': 0.1, # the proportionality constant
                    'i': 0.000, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT1_monitor',
                    'power_dataset':'p_AOM_A1',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },
            'dds_AOM_A2':
                {
                    'sampler_ch': 5, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD2_AOM_A2'], # volts
                    'p': 0.1, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT2_monitor',
                    'power_dataset':'p_AOM_A2',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },
            'dds_AOM_A3':
                {
                    'sampler_ch': 3, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD3_AOM_A3'], # volts
                    'p': 0.2, # the proportionality constant
                    'i': 0.000, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT3_monitor',
                    'power_dataset':'p_AOM_A3',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },
            'dds_AOM_A4':
                {
                    'sampler_ch': 4, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD4_AOM_A4'], # volts
                    'p': 0.2, # the proportionality constant
                    'i': 0.0, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT4_monitor',
                    'power_dataset':'p_AOM_A4',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },
            'dds_AOM_A5':
                {
                    'sampler_ch': 1, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD5_AOM_A5'],
                    'p': 0.07, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset':'MOT5_monitor',
                    'power_dataset':'p_AOM_A5',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },
            'dds_AOM_A6':
                {
                    'sampler_ch': 2, # the channel connected to the appropriate PD
                    'set_points': ['set_point_PD6_AOM_A6'], # volts
                    'p': 0.08, # the proportionality constant
                    'i': 0.00, # the integral coefficient
                    'series': True,
                    'dataset': 'MOT6_monitor',
                    'power_dataset':'p_AOM_A6',
                    't_measure_delay':0.5*ms,
                    'max_dB': 0
                },

            'dds_FORT': # signal monitored by TTI detector connected to MM fiber
                {
                    'sampler_ch': 6, # the channel connected to the appropriate PD
                    'set_points': ['set_point_FORT_MM_loading', 'set_point_FORT_MM_science'],
                    'p': 0.4, # set to 0.7 if using VCA
                    'i': 0.0, # the integral coefficient
                    'series': True, # setting to True because there's a bug with parallel
                    'dataset':'FORT_monitor',
                    'power_dataset':'p_FORT_loading',
                    't_measure_delay':10*us, # time to wait between AOM turned on and measurement
                    'max_dB': 5
                },
        },
    'sampler1':
        {
        }
}


class FeedbackChannel:

    def __init__(self, stabilizer, name, dds_obj, buffer_index, set_points, p, i, frequency, amplitude, dataset,
                 dB_dataset, t_measure_delay, error_history_length, max_dB):
        """
        class which defines a DDS feedback channel

        parameters:
        'stabilizer': an instance of AOMPowerStabilizer
        'name': the name of the dds channel, e.g. dds_AOM_A1
        'dds_obj': the experiment reference to the dds object with name 'name', e.g. experiment.dds_AOM_A1
        'buffer_index': int, the Sampler channel, e.g., 3, which we measure from.
        'set_points': a list of the float value(s) we are trying to reach. an index specifies which to use, with the
            first entry in the list assumed to be the default value.
        'p': the proportional constant for feedback
        'i': the integral constant. typically zero, as integral doesn't make much sense for our feedback scheme
        'frequency': the default frequency for the dds. see DeviceAliases.py
        'dataset': the dataset which stores the process value normalized to the setpoint
        'dB_dataset': the dataset which stores the power for the dds channel in dB
        't_measure_delay': time to wait after turning on the dds before making a measurement
        'error_history_length->int':, how many errors to store for the integral term
        'max_dB': the maximum dB we can attempt to set for the dds channel
        """
        self.stabilizer = stabilizer
        self.name = name
        self.dds_obj = dds_obj
        self.buffer_index = buffer_index
        # self.g = g
        self.set_points = set_points
        self.p = p
        self.i = i
        self.frequency = frequency
        self.amplitude = amplitude
        self.amplitudes = np.zeros(len(self.set_points)) # we would rather output too no RF than too much
        self.amplitudes[0] = self.amplitude
        self.feedback_sign = 1.0
        self.value = 0.0 # the last value of the measurement
        self.value_normalized = 0.0 # self.value normalized to the set point
        self.dataset = dataset
        self.dB_dataset = dB_dataset # the name of the dataset that stores the dB RF power for the dds
        self.dB_history_dataset = dB_dataset + str("_history")
        self.t_measure_delay = t_measure_delay
        self.error_history_arr = np.full(error_history_length,0.0)
        self.error_buffer = np.full(error_history_length-1,0.0)
        self.error_history_length = error_history_length
        self.cumulative_error = 0.0 # will store a sum of the
        self.max_dB = max_dB
        self.max_ampl = (2 * 50 * 10 ** (self.max_dB / 10 - 3)) ** (1 / 2)
        self.ampl_default = 0.0

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def set_value(self, value, setpoint_index=0):
        self.value = value
        self.value_normalized = value / self.set_points[setpoint_index]

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
    def feedback(self, buffer, setpoint_index=0):
        """ feedback to this channel
        buffer: a list of length storing the measurement result from all samplers
        """

        measured = buffer[self.buffer_index]
        err = self.set_points[setpoint_index] - measured

        # runs off the kernel, else the error array isn't correctly updated
        self.update_error_history(err)

        ampl = self.amplitudes[setpoint_index] + self.feedback_sign*self.p * err + self.i * self.cumulative_error

        self.set_value(measured, setpoint_index)

        if ampl < 0:
            self.amplitudes[setpoint_index] = 0.0
        elif ampl > self.max_ampl:
            # we either overcorrected  or there is not enough laser power right now to reach the setpoint
            # and in either case we are now stuck because the slope changed sign, so we throw ourselves
            # back to the bottom of the correct side of the laser power vs RF amplitude curve
            self.amplitudes[setpoint_index] = (2 * 50 * 10 ** (-20 / 10 - 3)) ** (1 / 2)

        else:  # valid amplitude for DDS
            self.amplitudes[setpoint_index] = ampl

        if setpoint_index == 0:
            self.amplitude = self.amplitudes[0] # for seamless backwards compatibility with older code

        # update the dds amplitude
        self.dds_obj.set(frequency=self.frequency,amplitude=self.amplitudes[setpoint_index])

    @kernel
    def set_dds_to_defaults(self, setpoint_index=0):
        """
        sets the dds channel to the default frequency and to the most recent amplitude.

        Except for the first time the feedback runs in your experiment, the amplitude stored in
        self.amplitude will be the one that was set at the end of the feedback routine.
        :param setpoint_index=0: the index of the setpoint, which specifies which amplitude we should use
        """
        self.dds_obj.set(frequency=self.frequency,amplitude=self.amplitudes[setpoint_index])

    @kernel
    def run(self, monitor_only=False, record_all_measurements=False, setpoint_index=0):
        """
        Run feedback for this instance of FeedbackChannel only, and leave the DDS on

        This ignores the any other FeedbackChannels which belong to self.stabilizer, so it is the responsibility of the
        user to turn channels on/off prior to calling this.

        :param monitor_only=False:. If True, does not apply feedback but still updates datasets based on measurements
        :param record_all_measurements=False: If True, dataset will be updated with each iteration of feedback, not just
            the last value measured
        :param setpoint_index=0: the index of the setpoint to use, where 0 is the default. note that dB dataset is not
            updated for setpoint_index != 0.
        :return:
        """

        if not self.stabilizer.exp.enable_laser_feedback:
            monitor_only = True

        self.dds_obj.sw.on()
        self.set_dds_to_defaults(setpoint_index)
        delay(0.1 * ms)

        # do feedback on the "parallel" channels
        for i in range(self.stabilizer.iterations):
            self.stabilizer.measure()
            delay(0.1 * ms)

            if not (self.stabilizer.dry_run or monitor_only):
                self.feedback(self.stabilizer.measurement_array - self.stabilizer.background_array, setpoint_index)
            else:
                self.set_value((self.stabilizer.measurement_array - self.stabilizer.background_array)[
                                   self.buffer_index], setpoint_index)

            if record_all_measurements:
                self.stabilizer.exp.append_to_dataset(self.dataset, self.value_normalized)

            delay(0.1 * ms)

        # self.stabilizer.measure_background()
        # delay(1*ms)

        self.dds_obj.sw.on()

        # update the datasets with the last values if we have not already done so
        if not record_all_measurements:
            self.stabilizer.exp.append_to_dataset(self.dataset, self.value_normalized)

        delay(0.1 * ms)

        if self.stabilizer.update_dds_settings and setpoint_index == 0:
            self.stabilizer.update_dB_dataset()

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

                        fb_channel = FeedbackChannel(
                                            stabilizer=self,
                                            name=ch_name,
                                            dds_obj=getattr(self.exp, ch_name),
                                            buffer_index=ch_params['sampler_ch'] + i*8,
                                            set_points=[getattr(self.exp,sp) for sp in ch_params['set_points']],
                                            p=ch_params['p'],
                                            i=ch_params['i'],
                                            frequency=getattr(self.exp,DDS_DEFAULTS[ch_name]["frequency"]),
                                            amplitude=dB_to_V(getattr(self.exp,DDS_DEFAULTS[ch_name]["power"])),
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

        for ch in self.all_channels:
            try:
                self.exp.set_dataset(ch.dataset, [self.exp.get_dataset(ch.dataset)[-1]], broadcast=True)
                self.exp.append_to_dataset(ch.dB_history_dataset,
                                           [self.exp.get_dataset(ch.dB_dataset, archive=False)])
            except Exception as e:
                logging.warning(e)
                self.exp.set_dataset(ch.dataset, [1.0], broadcast=True)
                self.exp.set_dataset(ch.dB_history_dataset,
                                           [float(self.exp.get_dataset(ch.dB_dataset, archive=False))], broadcast=True)

    @rpc(flags={"async"})
    def print(self, x):
        print(x)

    @kernel
    def update_dB_dataset(self):
        for ch in self.all_channels:
            dB = 10*(np.log10(ch.amplitude**2/(2*50)) + 3)
            self.exp.set_dataset(ch.dB_dataset, dB, broadcast=True, persist=True)
            self.exp.append_to_dataset(ch.dB_history_dataset, dB)

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
                delay(0.1*ms)
                dummy[i * 8:(i + 1) * 8] += measurement_array
                delay(0.1 * ms)
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

        print_async("self.background_array:")
        print_async(self.background_array)

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
    def run(self, record_all_measurements=False, monitor_only=False, defaults_at_start=True):
        """
        Run the feedback loop. On exiting, this function will turn off all dds channels
         given by dds_names. If any beams which need to be adjusted in series are in
         dds_names, all such channels will be turned off, i.e. even ones we are not
         feeding back to.

        record_all_measurements: optional, False by default. if True, every measurement point will
        be posted to the corresponding dataset, else, only post the final measurement, i.e. after the feedback
        loop has been run.
        monitor_only=False: if True, measurements of the detectors are made and the monitor datasets are updated,
            but feedback is not applied. this is overridden to be True if the experiment isntance for AOMPowerStabilizer
            has enable_laser_feedback = False.
        defaults_at_start=True: if True, the dds for each FeedbackChannel is set to
            the frequency and amplitude associated with the FeedbackChannel object. This should nearly always be true,
            unless for example, in cases where one wants specifically wants to read the detector values after the dds
            amplitudes have been changed in an experiment. This happens at the end of SamplerMOTCoilAndBeamBalance.
        """
        self.exp.core.reset()

        if not self.exp.enable_laser_feedback:
            monitor_only = True

        # self.exp.ttl7.pulse(1*ms) # scope trigger
        delay(100*ms)

        # turn off the repumps which are mixed into the cooling light
        self.exp.ttl_repump_switch.on() # block RF to the RP AOM
        self.exp.dds_pumping_repump.sw.off()

        if defaults_at_start:
            for ch in self.all_channels:
                ch.set_dds_to_defaults()
                delay(0.1 * ms)

        with sequential:

            # don't blind the fW detector
            self.exp.dds_FORT.sw.off()

            for ch in self.series_channels:
                ch.dds_obj.sw.off()

            delay(0.1*ms)
            for ch in self.parallel_channels:
                ch.dds_obj.sw.on()
            delay(0.1*ms)

            # do feedback on the "parallel" channels
            for i in range(self.iterations):
                self.measure()
                delay(0.1 * ms)

                # strictly speaking, doesn't really matter if these are parallel
                for ch in self.parallel_channels:
                    delay(0.1*ms)
                    if not (self.dry_run or monitor_only):
                        ch.feedback(self.measurement_array - self.background_array)
                    else:
                        ch.set_value((self.measurement_array - self.background_array)[ch.buffer_index])

                    if record_all_measurements:
                        self.exp.append_to_dataset(ch.dataset, ch.value_normalized)

                delay(0.1 * ms)

            for ch in self.parallel_channels:
                ch.dds_obj.sw.off()
            # delay(10 * ms)  # the Femto fW detector is slow

            # need to have this on for MOT feedback
            self.exp.dds_cooling_DP.sw.on() # todo: only turn this on if the one of the FeedbackChannels depends on it
            delay(1 * ms)

            # self.measure_background()
            # delay(1*ms)

            # do feedback on the series channels
            with sequential:
                for ch in self.series_channels:
                    ch.dds_obj.sw.on()

                    delay(ch.t_measure_delay) # allows for detector rise time
                    for i in range(self.iterations):
                        self.measure()
                        delay(0.1 * ms)

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
        self.exp.dds_cooling_DP.sw.on() # todo: only turn this on if the one of the FeedbackChannels depends on it
        self.exp.ttl_repump_switch.off()  # enable RF to the RP AOM

        delay(0.1*ms)

        if self.leave_AOMs_on:
            for ch in self.all_channels:
                ch.dds_obj.sw.on()

        delay(0.1* ms)
        if self.update_dds_settings:
            self.update_dB_dataset()

