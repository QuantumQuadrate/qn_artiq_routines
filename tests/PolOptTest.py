"""
new version based on GVS
"""


from artiq.experiment import *
import logging

import numpy as np
from numpy import array  # necessary for some override_ExperimentVariable entries

import sys, os

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions
from subroutines.aom_feedback import AOMPowerStabilizer

class Polarization_Optimization_Test(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        # it will not run correctly without this statement
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_device("core")
        self.setattr_device("sampler2")
        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")

        self.setattr_argument("max_iterations", NumberValue(4, ndecimals=1, step=1))
        self.setattr_argument("full_range", NumberValue(20, ndecimals=1, step=1))
        self.setattr_argument("sample_pts", NumberValue(9, ndecimals=1, step=1))
        self.setattr_argument("start_from_home", BooleanValue(default=False), "initialization")
        # self.setattr_argument("start_from_optimized_point", BooleanValue(default=False), "initialization")

        # n_measurements not used
        # self.setattr_argument("n_measurements", NumberValue(0, ndecimals=1, step=1))
        # self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        # self.base.prepare()

        #todo: put this in BaseExperiment.py
        self.setpoint_datasets = ["best_HWP_to_H","best_QWP_to_H"]
        self.default_setpoints = [getattr(self,dataset) for dataset in self.setpoint_datasets]


    @kernel
    def initialize_hardware(self):
        # self.base.initialize_hardware()
        self.core.reset()
        self.sampler2.init()  # for reading laser feedback
        # delay(1*s)
        if self.start_from_home:
            go_to_home(self,'780_HWP')
            go_to_home(self,'780_QWP')
        else:
            move_to_target_deg(self, name="780_HWP", target_deg=self.best_HWP_to_H)
            move_to_target_deg(self, name="780_QWP", target_deg=self.best_QWP_to_H)

    def initialize_datasets(self):
        # self.base.initialize_datasets()
        self.set_dataset("test_PDA_monitor", [0.0], broadcast=True)
        self.set_dataset("HWP_angle", [0.0], broadcast=True)
        self.set_dataset("QWP_angle", [0.0], broadcast=True)


    @kernel
    def exp_fun(self):
        delay(10*s)
        # scan in up to 2 dimensions. for each setting of the parameters, run experiment_function n_measurement times
        # iterations = 0

        max_iterations = int(self.max_iterations)  # Limit iterations to prevent infinite loops
        full_range = float(self.full_range)  # Start with a full range
        sample_pts = int(self.sample_pts)  # Number of samples per iteration

        # only half_range used in iteration. do not use full_range
        half_range = full_range / 2  # Half the full range
        measurement_count = 0  # Track number of power measurements

        if self.start_from_home:
            best_HWP, best_QWP, best_power = 0.0, 0.0, 0.0
        else:
            best_HWP, best_QWP, best_power = self.best_HWP_to_H, self.best_QWP_to_H, 0.0
        previous_hwp, previous_qwp, previous_power = 0.0, 0.0, 0.0

        power = 0.0

        # Iterative search loop
        for iteration in range(max_iterations):
            print("iteration no.", iteration) # this statement has to be here for it to run,,,, why...
            time_ite = time_to_rotate_in_ms(self, half_range)

            delay(2 * 2 * time_ite * ms + 1*s)  # rotate two waveplates

            steps = half_range * 2 / (sample_pts - 1)
            hwp_values = np.array([best_HWP - half_range + i * steps for i in range(sample_pts)])
            qwp_values = np.array([best_QWP - half_range + i * steps for i in range(sample_pts)])


            for hwp in hwp_values:
                # time to rotate hwp
                current_hwp = hwp
                time_hwp = time_to_rotate_in_ms(self, current_hwp - previous_hwp)
                delay(time_hwp * ms)
                count_down = 0

                stop_loop = False

                for qwp in qwp_values:

                    if not stop_loop:
                        current_qwp = qwp
                        # time to rotate qwp
                        time_qwp = time_to_rotate_in_ms(self, current_qwp - previous_qwp)
                        delay(time_qwp * ms)

                        self.append_to_dataset("HWP_angle", hwp)
                        self.append_to_dataset("QWP_angle", qwp)

                        with parallel:
                            move_to_target_deg(self, name="780_HWP", target_deg=hwp)
                            move_to_target_deg(self, name="780_QWP", target_deg=qwp)

                        delay(1*s)
                        power = record_PDA_power(self)  # Measure power

                        print("hwp, qwp: ", hwp, ", ", qwp, "power: ", power)

                        delay(1*s)

                        if power > best_power:  # Update best power and angles if improved
                            best_power = power
                            best_HWP = hwp
                            best_QWP = qwp

                        previous_qwp = qwp
                        measurement_count += 1  # Increment measurement counter

                        if power < previous_power:
                            count_down += 1
                            if count_down > int(sample_pts / 2):
                                stop_loop = True

                        previous_power = power

                previous_hwp = hwp
                delay(1*s)
                # iterations += 1

            half_range = steps
            sample_pts = max(3, sample_pts // 2)  # Prevent sample_pts from becoming too small
            # sample_pts must be odd number to reproduce former best parameter
            if sample_pts % 2 == 0:
                sample_pts += 1

            delay(1 * s)
            print("best_HWP, best_QWP, power = ", best_HWP,", ", best_QWP, ",", best_power)

        self.set_dataset("best_HWP_to_H", best_HWP, broadcast=True, persist=True)
        self.set_dataset("best_QWP_to_H", best_QWP, broadcast=True, persist=True)



    def run(self):
        """
        Step through the variable values defined by the scan sequences and run the experiment function.

        Because the scan variables can be any ExperimentVariable, which includes values used to initialize
        hardware (e.g. a frequency for a dds channel), the hardware is reinitialized in each step of the
        variable scan, i.e., each iteration.
        """
        self.core.reset()

        self.initialize_hardware()
        self.initialize_datasets()

        self.exp_fun()



