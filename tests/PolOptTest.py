# """
# Polarizatino Optimization Test
# """
#
# from artiq.experiment import *
#
# import numpy as np
# import sys, os
# # get the current working directory
# cwd = os.getcwd() + "\\"
# sys.path.append(cwd)
# sys.path.append(cwd+"\\repository\\qn_artiq_routines")
#
# from utilities.BaseExperiment import BaseExperiment
# from subroutines.k10cr1_functions import *
#
#
#
# class Polarization_Optimization_Test(EnvExperiment):
#
#     def build(self):
#         self.base = BaseExperiment(experiment=self)
#         self.base.build()
#         # self.setattr_argument("max_iterations", NumberValue(4, ndecimals=1, step=1))
#
#         self.base.set_datasets_from_gui_args()
#         print("build - done")
#
#
#     def prepare(self):
#         self.base.prepare()
#
#     @kernel
#     def exp_fun(self):
#         delay(1*s)
#         # start by homing the device
#         go_to_home('780_HWP')
#         go_to_home('780_QWP')
#
# ### self.
#         max_iterations = 4  # Limit iterations to prevent infinite loops
#         full_range = 20  # Start with a full range
#         sample_pts = 9  # Number of samples per iteration
#
#         current_range = full_range / 2  # Half the full range
#         measurement_count = 0  # Track number of power measurements
#
#         rotator_ave_speed = 10 # lets assume that this is 10deg/s. Actual maximum speed is 20deg/s
#
#         best_HWP, best_QWP = 0, 0
#         best_power = 0
#
#
#         # Iterative search loop
#         for iteration in range(max_iterations):
#             delay(1*s)
#             # Generate sample_pts^2 grid of HWP and QWP angles
#             hwp_values = np.linspace(best_HWP - current_range, best_HWP + current_range, sample_pts)
#             qwp_values = np.linspace(best_QWP - current_range, best_QWP + current_range, sample_pts)
#
#             # Measure power for each combination in the grid
#             best_local_HWP, best_local_QWP = best_HWP, best_QWP  # Initialize best angles
#             previous_max_power = best_power  # Store previous iteration max
#
#             previous_hwp = 0
#             previous_qwp = 0
#
#             for hwp in hwp_values:
#                 for qwp in qwp_values:
#                     current_hwp = hwp
#                     current_qwp = qwp
#
#                     time_hwp = np.round(np.absolute(current_hwp - previous_hwp) / rotator_ave_speed, 2)
#                     time_qwp = np.round(np.absolute(current_qwp - previous_qwp) / rotator_ave_speed, 2)
#
#                     longer_time = max(time_hwp, time_qwp)
#
#                     delay(longer_time* s)
#
#                     with parallel:
#                         move_to_target_deg(self, name="780_HWP", target_deg=hwp)
#                         move_to_target_deg(self, name="780_QWP", target_deg=qwp)
#
#                     measurement_count += 1  # Increment measurement counter
#
#                     delay(1*s)
#                     power = record_PDA_power(self)  # Measure power
#
#                     # measured_powers.append(power)  # Save power measurement
#                     # iterations.append(iteration + 1)  # Store iteration index
#
#                     if power > best_power:  # Update best power and angles if improved
#                         best_power = power
#                         best_local_HWP = hwp
#                         best_local_QWP = qwp
#
#                     previous_hwp = hwp
#                     previous_qwp = qwp
#
#                     iterations += 1
#
#
#             # Update the best angles and refine the search range correctly
#             best_HWP, best_QWP = best_local_HWP, best_local_QWP
#             current_range /= sample_pts
#             sample_pts = max(3, sample_pts // 2)  # Prevent sample_pts from becoming too small
#
#     def initialize_hardware(self):
#         self.base.initialize_hardware()
#
#
#     def run(self):
#         self.base.build()
#         self.base.prepare()
#         # self.initialize_hardware()
#
#         iterations = 0
#         self.set_dataset("iteration", iterations, broadcast=True)
#         self.set_dataset("test_PDA_monitor", [0.0], broadcast=True)
#
#         self.exp_fun()
#
#
#
#
#
#
#
#
#


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
        self.base = BaseExperiment(experiment=self)
        self.base.build()
        self.setattr_argument("max_iterations", NumberValue(4, ndecimals=1, step=1))
        self.setattr_argument("full_range", NumberValue(20, ndecimals=1, step=1))
        self.setattr_argument("sample_pts", NumberValue(9, ndecimals=1, step=1))

        # n_measurements not used
        self.setattr_argument("n_measurements", NumberValue(0, ndecimals=1, step=1))
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment.
        """
        self.base.prepare()


    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()
        go_to_home(self,'780_HWP')
        go_to_home(self,'780_QWP')

    def initialize_datasets(self):
        self.base.initialize_datasets()
        self.set_dataset("test_PDA_monitor", [0.0], broadcast=True)
    @rpc
    def custom_abs(self, num):
        if num < 0:
            return -num
        else:
            return num

    @rpc
    def custom_round(value, decimals=2):
        factor = 10 ** decimals
        return int(value * factor + (0.5 if value > 0 else -0.5)) / factor

    @kernel
    def exp_fun(self):
        delay(1*s)

        print("started")

        # scan in up to 2 dimensions. for each setting of the parameters, run experiment_function n_measurement times
        iterations = 0
        self.set_dataset("iteration", iterations, broadcast=True)

        max_iterations = int(self.max_iterations)  # Limit iterations to prevent infinite loops
        full_range = float(self.full_range)  # Start with a full range
        sample_pts = int(self.sample_pts)  # Number of samples per iteration

        current_range = full_range / 2  # Half the full range
        measurement_count = 0  # Track number of power measurements

        rotator_ave_speed = 10 # lets assume that this is 10deg/s. Actual maximum speed is 20deg/s

        best_HWP, best_QWP = 0.0, 0.0
        best_power = 0.0

        # Iterative search loop
        for iteration in range(max_iterations):
            # delay(10*s)

            hwp_values = np.array(
                [best_HWP - current_range + i * (2 * current_range / (sample_pts - 1)) for i in range(sample_pts)])
            qwp_values = np.array(
                [best_QWP - current_range + i * (2 * current_range / (sample_pts - 1)) for i in range(sample_pts)])

            # Measure power for each combination in the grid
            best_local_HWP, best_local_QWP = best_HWP, best_QWP  # Initialize best angles
            # previous_max_power = best_power  # Store previous iteration max

            previous_hwp = 0.0
            previous_qwp = 0.0

            current_hwp = 0.0
            current_qwp = 0.0

            for hwp in hwp_values:
                for qwp in qwp_values:
                    current_hwp = hwp
                    current_qwp = qwp

                    # time_hwp = self.custom_round(self.custom_abs(current_hwp - previous_hwp) / rotator_ave_speed, 2)
                    # time_qwp = self.custom_round(self.custom_abs(current_qwp - previous_qwp) / rotator_ave_speed, 2)
                    #
                    # longer_time = max(time_hwp, time_qwp)
                    #
                    # delay(longer_time* s)
                    print("hwp, qwp: ", hwp, ", ", qwp)
                    delay(10 * s)

                    with parallel:
                        move_to_target_deg(self, name="780_HWP", target_deg=hwp)
                        move_to_target_deg(self, name="780_QWP", target_deg=qwp)

                    measurement_count += 1  # Increment measurement counter

                    delay(1*s)
                    power = record_PDA_power(self)  # Measure power

                    print("power: ", power)
                    delay(1*s)

                    # measured_powers.append(power)  # Save power measurement
                    # iterations.append(iteration + 1)  # Store iteration index

                    if power > best_power:  # Update best power and angles if improved
                        best_power = power
                        best_local_HWP = hwp
                        best_local_QWP = qwp

                    previous_hwp = hwp
                    previous_qwp = qwp

                    iterations += 1
                    # self.append_to_dataset("iteration", iterations)


            # Update the best angles and refine the search range correctly
            best_HWP, best_QWP = best_local_HWP, best_local_QWP
            current_range /= sample_pts
            sample_pts = max(3, sample_pts // 2)  # Prevent sample_pts from becoming too small

            print(f"best_HWP, best_QWP = {best_HWP}, {best_QWP}")



    def run(self):
        """
        Step through the variable values defined by the scan sequences and run the experiment function.

        Because the scan variables can be any ExperimentVariable, which includes values used to initialize
        hardware (e.g. a frequency for a dds channel), the hardware is reinitialized in each step of the
        variable scan, i.e., each iteration.
        """
        self.core.reset()

        self.base.build()
        # self.base.prepare()
        self.initialize_hardware()
        self.initialize_datasets()

        self.exp_fun()




