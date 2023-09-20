from artiq.experiment import *
import math

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer


class AOMFeedbackTest(EnvExperiment):
    # this experiment assumes that the AOMs we'll feed back to have already been
    # set to operate at ~ 70% efficiency by appropriate choice of default powers
    # in ExperimentVariables.

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("feedback_dds_list",
                              StringValue(
                                  "['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', 'dds_AOM_A4','dds_AOM_A4','dds_AOM_A5',"
                                  "'dds_AOM_A6','dds_cooling_DP']"))
        self.setattr_argument("reinitialize_DDSs_each_iteration", BooleanValue(False))

        self.setattr_argument("experiment_iterations", NumberValue(1, type='int', ndecimals=0, scale=1, step=1))

        # note this is ignored if running with fake drift
        self.setattr_argument("t_iteration_delay", NumberValue(1 * s, unit='s'))

        # for when, e.g., we want to monitor what a channel is doing when the feedback isn't on
        # note that this will report a near-zero value if the channel requested is not on when
        # the feedback ends. the value should be in the feedback_dds_list below.
        open_loop_group = "open loop monitor"
        self.setattr_argument("single_ch_open_loop_monitor", StringValue('dds_AOM_A1'),open_loop_group)
        self.setattr_argument("enable_open_loop_monitor", BooleanValue(False),open_loop_group)
        self.setattr_argument("t_monitor_period", NumberValue(10 * ms, unit='ms'),open_loop_group)

        self.setattr_argument("exclude_from_fake_drift_dds_list",
                              StringValue(
                                  "['dds_cooling_DP']"),"simulate drift")
        self.setattr_argument("no_feedback", BooleanValue(False),"feedback params")
        self.setattr_argument("feedback_iterations", NumberValue(10, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")
        self.setattr_argument("averages", NumberValue(10, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")

        # whether to plot all measurements or only after each feedback loop has been applied
        self.setattr_argument("record_all_measurements", BooleanValue(False), "feedback params")

        self.setattr_argument("run_with_fake_drift", BooleanValue(False), "simulate drift")
        self.setattr_argument("drift_factor", NumberValue(0.9), "simulate drift")

        self.base.set_datasets_from_gui_args()

    def prepare(self):

        self.base.prepare()

        self.timesteps = int(self.t_iteration_delay/self.t_monitor_period+0.5)

        # list of names of dds channels
        self.dds_list = eval(self.feedback_dds_list)
        self.exclude_list = eval(self.exclude_from_fake_drift_dds_list)

        if self.enable_open_loop_monitor:
            assert self.single_ch_open_loop_monitor in self.dds_list

        if self.single_ch_open_loop_monitor != '':
            self.open_loop_monitors = [self.single_ch_open_loop_monitor]
        else:
            self.open_loop_monitors = []

        # overwrites the instance created in base.prepare
        self.laser_stabilizer = AOMPowerStabilizer(experiment=self,
                                           dds_names=self.dds_list,
                                           iterations=self.feedback_iterations,
                                           averages=self.averages,
                                           leave_AOMs_on=True,
                                           update_dds_settings=False,
                                           dry_run=self.no_feedback,
                                           open_loop_monitor_names = self.open_loop_monitors
                                )

        # for running with fake drift, we need the dds references themselves
        self.dds_refs = []
        for dds_name in self.dds_list:
            if dds_name not in self.exclude_list:
                self.dds_refs.append(getattr(self, dds_name))

        self.ampl_list = [0.0] * (len(self.dds_refs))
        self.freq_list = [0.0] * (len(self.dds_refs))


        # the cooling double pass AOM -  we'll use this to fake a power drift.
        # this will suffice for feeding back to the fiber AOMs, but not to the
        # other beams.
        self.dds_drift = self.dds_cooling_DP
        self.dds_default_ampl = self.ampl_cooling_DP_MOT
        self.dds_default_freq = self.f_cooling_DP_MOT
        self.dds_drift_ampl = self.dds_default_ampl*self.drift_factor
        self.dds_drift_freq = self.dds_default_freq

    @kernel
    def run(self):
        self.base.initialize_hardware()

        self.core.break_realtime()

        if self.run_with_fake_drift:

            # first get the default settings
            i = 0
            for dds in self.dds_refs:
                freq, _, ampl = dds.get()
                delay(10*ms)
                self.ampl_list[i] = ampl
                self.freq_list[i] = freq
                i += 1

            for i in range(self.experiment_iterations):

                if self.drift_factor != 1:
                    print("faking drift")
                    j = 0
                    for dds in self.dds_refs:
                        dds.set(frequency=self.freq_list[j],amplitude=self.ampl_list[j]*self.drift_factor)
                        delay(1*ms)
                        j += 1

                self.laser_stabilizer.run(record_all_measurements=self.record_all_measurements)
                delay(self.t_iteration_delay)

        else:

            print("running feedback")
            delay(10 * ms)

            for i in range(self.experiment_iterations):

                self.laser_stabilizer.run(record_all_measurements=self.record_all_measurements,
                                          monitor_only=self.no_feedback)
                delay(10*ms)

                for t in range(self.timesteps):
                    if self.enable_open_loop_monitor:
                        self.laser_stabilizer.open_loop_monitor()
                    delay(self.t_monitor_period)

                # reset the dds settings. since we initialized the AOMPowerStabilizer
                # with update_dds_settings=False, the values used will be those used at the beginning
                if self.reinitialize_DDSs_each_iteration:
                    self.named_devices.initialize()



        print("experiment finished")