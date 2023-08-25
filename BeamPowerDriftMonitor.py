"""
Runs AOMPowerStabilizer2.monitor to broadcast the dataset for each beam power.
NOTE: this is for measuring beam power drift in the absence of feedback,
i.e. this is not with feedback on.

If only measuring the MOT beams, plot with the plot_MOT_powers applet
"""
from artiq.experiment import *
import math

from utilities.BaseExperiment import BaseExperiment
from subroutines.aom_feedback import AOMPowerStabilizer2


class BeamPowerDriftMonitor(EnvExperiment):
    # this experiment assumes that the AOMs we'll feed back to have already been
    # set to operate at ~ 70% efficiency by appropriate choice of default powers
    # in ExperimentVariables.

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.setattr_argument("experiment_iterations", NumberValue(20, type='int', ndecimals=0, scale=1, step=1))
        self.setattr_argument("dds_list",
                              StringValue(
                                  "['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', 'dds_AOM_A4','dds_AOM_A4','dds_AOM_A5',"
                                  "'dds_AOM_A6','dds_cooling_DP']"))
        self.setattr_argument("t_measurement_delay", NumberValue(500*ms, unit='ms'))

    def prepare(self):

        self.base.prepare()

        dds_list = eval(self.dds_list)

        self.laser_stabilizer = AOMPowerStabilizer2(experiment=self,
                                           dds_names=dds_list,
                                           iterations=1) # not used

    @kernel
    def run(self):
        self.base.initialize_hardware()

        for i in range(self.experiment_iterations):

            self.laser_stabilizer.monitor()
            delay(self.t_measurement_delay)

    print("experiment finished")