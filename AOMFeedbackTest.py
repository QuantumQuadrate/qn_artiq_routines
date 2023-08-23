from artiq.experiment import *
import math

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer2


class AOMPowerStabilizerTest2(EnvExperiment):
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
        self.setattr_argument("run_with_fake_drift", BooleanValue(False),"simulate drift")
        self.setattr_argument("drift_factor", NumberValue(0.9),"simulate drift")
        self.setattr_argument("feedback_iterations", NumberValue(5, type='int', ndecimals=0, scale=1, step=1),
                              "feedback params")
        # self.setattr_argument("t_measurement_delay", NumberValue(20*ms, unit='ms'),
        #                       "feedback params")

    def prepare(self):

        self.base.prepare()

        feedback_dds_list = eval(self.feedback_dds_list)

        self.laser_stabilizer = AOMPowerStabilizer2(experiment=self,
                                           dds_names=feedback_dds_list,
                                           iterations=self.feedback_iterations,
                                           # t_meas_delay=self.t_measurement_delay)
                                           t_meas_delay=10*ms) # not used

        # the cooling single pass AOM -  we'll use this to fake a power drift.
        # this will suffice for feeding back to the cooling DP and the fiber AOMs,
        # but not for other beams.
        self.dds_drift = self.dds_cooling_SP
        self.dds_default_ampl = self.ampl_cooling_SP
        self.dds_default_freq = self.f_cooling_SP
        self.dds_drift_ampl = self.ampl_cooling_SP*self.drift_factor
        self.dds_drift_freq = self.dds_default_freq

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
                delay(self.t_iteration_delay)
        print("experiment finished")