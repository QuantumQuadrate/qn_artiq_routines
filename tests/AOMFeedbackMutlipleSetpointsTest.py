from artiq.experiment import *
import math

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer


class AOMFeedbackMultipleSetpointsTest(EnvExperiment):
    # this experiment assumes that the AOMs we'll feed back to have already been
    # set to operate at ~ 70% efficiency by appropriate choice of default powers
    # in ExperimentVariables.

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # self.base.set_datasets_from_gui_args()

    def prepare(self):

        self.base.prepare()

        # overwrite the instance created in BaseExperiment
        self.laser_stabilizer = AOMPowerStabilizer(experiment=self,
                                                   dds_names=['dds_FORT'],
                                                   iterations=self.aom_feedback_iterations,
                                                   averages=self.aom_feedback_averages,
                                                   leave_AOMs_on=False)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        self.core.break_realtime()

        for i in range(50):
            delay(100*ms)
            if i % 2:
                self.laser_stabilizer.run()  # run the feedback loop with default set points
            else:
                self.stabilizer_FORT.run(setpoint_index=1)  # FeedbackChannel.run always leaves the AOM on

        print("experiment finished")
