from artiq.experiment import *
import math

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer


class AOMPowerStabilizerTest(EnvExperiment):
    # this experiment assumes that the AOMs we'll feed back to have already been
    # set to operate at ~ 70% efficiency by appropriate choice of default powers
    # in ExperimentVariables.

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.setattr_argument("experiment_iterations", NumberValue(20, type='int', ndecimals=0, scale=1, step=1))

        # note this is ignored if running with fake drift
        self.setattr_argument("t_iteration_delay", NumberValue(1 * s, unit='s'))

        self.setattr_argument("feedback_dds_list",
                              StringValue(
                                  "['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', 'dds_AOM_A4','dds_AOM_A4','dds_AOM_A5',"
                                  "'dds_AOM_A6','dds_cooling_DP']"))
        self.setattr_argument("exclude_from_fake_drift_dds_list",
                              StringValue(
                                  "['dds_cooling_DP']"))
        self.setattr_argument("run_with_fake_drift", BooleanValue(False),"simulate drift")
        self.setattr_argument("drift_factor", NumberValue(0.9),"simulate drift")
        self.setattr_argument("no_feedback", BooleanValue(False),"feedback params")
        self.setattr_argument("feedback_iterations", NumberValue(5, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")
        self.setattr_argument("averages", NumberValue(5, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")

        self.base.set_datasets_from_gui_args()

    def prepare(self):

        self.base.prepare()

        # list of names of dds channels
        self.dds_list = eval(self.feedback_dds_list)
        self.exclude_list = eval(self.exclude_from_fake_drift_dds_list)

        self.laser_stabilizer = AOMPowerStabilizer(experiment=self,
                                           dds_names=self.dds_list,
                                           iterations=self.feedback_iterations,
                                           averages=self.averages,
                                           dry_run=self.no_feedback,
                                           update_dds_settings=False)

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

                # print("running feedback")
                self.laser_stabilizer.run_tuning_mode()
                delay(self.t_iteration_delay)

        else:

            print("running feedback")
            delay(10 * ms)

            for i in range(self.experiment_iterations):
                self.laser_stabilizer.run()
                delay(1*ms)

                # reset the dds settings to the experiment variables. since we initialized the AOMPowerStabiizer
                # with update_dds_settings=False, the values used will be those used at the beginning
                self.named_devices.initialize()

                delay(self.t_iteration_delay)
        print("experiment finished")