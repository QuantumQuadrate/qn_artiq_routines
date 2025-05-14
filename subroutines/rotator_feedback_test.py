"""
This is a code for FORT polarization optimization with simple iterative approach.

- A few percent of the FORT light passes the parabolic mirror through a small hole.
- Then, it passes through a polarizer and collected through the MM fiber.
- This polarizer is pre-aligned to the preferred FORT polarization.
- Thus, by finding the two waveplate angles which maximize the power out from MM fiber,
- we can optimize the FORT polarization.


As long as the step size is small enough, it will always find the correct value.


"""


from artiq.experiment import *
import logging

import numpy as np
from numpy import array  # necessary for some override_ExperimentVariable entries

import time
import json

import sys, os

cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

# this is where your experiment function should live
from subroutines.experiment_functions import *
import subroutines.experiment_functions as exp_functions
from subroutines.aom_feedback import AOMPowerStabilizer

class FORTPolarizationStabilizer:

    def __init__(self, tolerance_deg, full_range, sample_pts, initialize_to_home = False, search_start_from_scratch = False, reduce_sample_pts = False):

        # self.exp = experiment

        self.tolerance_deg = tolerance_deg
        self.full_range = full_range
        self.sample_pts = sample_pts

        self.initialize_to_home = initialize_to_home
        self.search_start_from_scratch = search_start_from_scratch
        self.reduce_sample_pts = reduce_sample_pts



    @kernel
    def initialize_hardware(self):
        #todo: do this in Base Experiment
        if self.initialize_to_home:
            go_to_home(self,'852_HWP')
            go_to_home(self,'852_QWP')
            print("Both Waveplates Initialized to Home")
        else:
            # this will do nothing if it has been optimized before
            move_to_target_deg(self, name="852_HWP", target_deg=self.best_852HWP_to_max)
            move_to_target_deg(self, name="852_QWP", target_deg=self.best_852QWP_to_max)

    def initialize_datasets(self):
        # self.base.initialize_datasets()
        # self.exp.set_dataset("FORT_MM_monitor", [0.0], broadcast=True)
        # self.exp.set_dataset("FORT_APD_monitor", [0.0], broadcast=True)
        # self.exp.set_dataset("HWP_angle", [0.0], broadcast=True)
        # self.exp.set_dataset("QWP_angle", [0.0], broadcast=True)
        self.set_dataset("FORT_MM_monitor", [], broadcast=True)
        self.set_dataset("FORT_APD_monitor", [], broadcast=True)
        self.set_dataset("HWP_angle", [], broadcast=True)
        self.set_dataset("QWP_angle", [], broadcast=True)


    @kernel
    def optimization_routine(self):
        """
        Run optimization loop.
        :return:
        """

        self.core.reset()
        delay(1*s)

        # adding stabilizer
        if self.enable_laser_feedback:
            self.laser_stabilizer.run()

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()  ### turns FORT on


        # max_iterations = int(self.max_iterations)  # Limit iterations to prevent infinite loops
        tolerance = float(self.tolerance_deg)   # rather than using fixed iteration, stop when step < tolerance
        full_range = float(self.full_range)  # Start with a full range
        sample_pts = int(self.sample_pts)  # Number of samples per iteration

        # only half_range used in iteration. do not use full_range
        half_range = full_range / 2  # Half the full range
        measurement_count = 0  # Track number of power measurements

        if self.search_start_from_scratch:
            best_HWP, best_QWP, best_power = 0.0, 0.0, 0.0
        else:
            best_HWP, best_QWP, best_power = self.best_852HWP_to_max, self.best_852QWP_to_max, 0.0
        previous_hwp, previous_qwp, previous_power = 0.0, 0.0, 0.0

        power = 0.0
        power_APD = 0.0

        tolerance_satisfied = False
        iteration = 0


        # Iterative search loop
        while not tolerance_satisfied:
        # for iteration in range(max_iterations):
            print("FORT POL OPT - iteration no.", iteration)
            time_ite = time_to_rotate_in_ms(self, half_range)

            delay(2 * 2 * time_ite * ms + 1*s)  # rotate two waveplates

            steps = half_range * 2 / (sample_pts - 1)
            hwp_values = np.array([best_HWP - half_range + i * steps for i in range(sample_pts)])
            qwp_values = np.array([best_QWP - half_range + i * steps for i in range(sample_pts)])

            if steps <= tolerance:
                tolerance_satisfied = True

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
                            move_to_target_deg(self, name="852_HWP", target_deg=hwp)
                            move_to_target_deg(self, name="852_QWP", target_deg=qwp)

                        delay(1*s)
                        # power = record_PDA_power(self)  # Measure power
                        power = record_FORT_MM_power(self)
                        power_APD = record_FORT_APD_power(self)

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

            half_range = steps

            if self.reduce_sample_pts:
                sample_pts = max(3, sample_pts // 2)  # Prevent sample_pts from becoming too small
                # sample_pts must be odd number to reproduce former best parameter
                if sample_pts % 2 == 0:
                    sample_pts += 1

            iteration += 1

            delay(1 * s)
            print("best_HWP, best_QWP, best_power = ", best_HWP,", ", best_QWP, ", ", best_power)

        # move back to the best HWP, QWP
        move_to_target_deg(self, name="852_HWP", target_deg=best_HWP)
        move_to_target_deg(self, name="852_QWP", target_deg=best_QWP)

        print("full range search: ", full_range, "tolerance_deg: ", tolerance, "sample_pts: ", sample_pts)
        print("previous best_HWP, best_QWP, best_power = ", self.best_852HWP_to_max, ", ", self.best_852QWP_to_max, ", ", self.best_852_power)

        self.dds_FORT.sw.off()  ### turns FORT on

        self.set_dataset("best_852HWP_to_max", best_HWP, broadcast=True, persist=True)
        self.set_dataset("best_852QWP_to_max", best_QWP, broadcast=True, persist=True)
        self.set_dataset("best_852_power", best_power, broadcast=True, persist=True)


    @kernel
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

        self.optimization_routine()
