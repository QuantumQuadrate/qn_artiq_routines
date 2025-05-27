"""
2025.05.15
previous path and name: qn_artiq_routines/tests/PolOptTest.py
renamed to: FORT_Polarization_Optimizer.py

This is a code for FORT polarization optimization

- A few percent of the FORT light passes the parabolic mirror through a small hole.
- Then, it passes through a polarizer and collected through the MM fiber.
- This polarizer is pre-aligned to the preferred FORT polarization.
- Thus, by finding the two waveplate angles which maximize the power out from MM fiber,
- we can optimize the FORT polarization.


As long as the step size is small enough, it will always find the correct value.

Possible Questions:
Q: Why are we performing a simple scan instead of using a more advanced optimization algorithm?
    Wouldn’t a smarter algorithm find the optimal point with fewer measurements?
A: The main reason is that the waveplate rotator's rotation speed is relatively slow.
    While each measurement only takes about 1 ms, rotating the waveplate by 1 degree
    takes approximately 50 ms at best. Additionally, unless asynchronous control (e.g., using asyncio)
    is implemented, the two waveplates cannot be rotated simultaneously.

    In typical optimization algorithms, the search begins with large step sizes to explore a wide parameter space.
    However, in our case, each movement is time-consuming. As a result, using such algorithms would actually take more time overall.

Notes:
    * The maximum speed and acceleration of the rotator are 20°/s, but this performance is only achievable
    when using USB 3.0, which can supply up to 900 mA. If your USB cable or hub does not meet this power requirement,
    the speed will be limited to 10°/s. Here, I am assuming this 10deg/s.

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

class FORT_Polarization_Optimizer(EnvExperiment):
    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        # it will not run correctly without this statement
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("initialize_to_home", BooleanValue(default=False), "Initialization")

        self.setattr_argument("tolerance_deg", NumberValue(1, ndecimals=2, step=1), "optimization parameters")
        self.setattr_argument("full_range", NumberValue(90, ndecimals=1, step=1), "optimization parameters")
        self.setattr_argument("sample_pts", NumberValue(7, ndecimals=1, step=1), "optimization parameters")

        self.setattr_argument("search_start_from_scratch", BooleanValue(default=False), "optimization settings")

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

        if self.initialize_to_home:
            go_to_home(self,'852_HWP')
            go_to_home(self,'852_QWP')
            print("Both Waveplates Initialized to Home")

        # if self.which_node == 'bob':      # if I use this, it gives me llvm error! what happens in Node1??
        else:
            #### initializing the waveplates to best_852 parameters is necessary if using MM feedback
            # this will do nothing if it has been optimized before
            move_to_target_deg(self, name="852_HWP", target_deg=self.best_852HWP_to_max)
            move_to_target_deg(self, name="852_QWP", target_deg=self.best_852QWP_to_max)

    def initialize_datasets(self):
        # self.base.initialize_datasets()
        self.set_dataset("FORT_MM_monitor", [], broadcast=True)
        self.set_dataset("FORT_APD_monitor", [], broadcast=True)
        self.set_dataset("HWP_angle", [], broadcast=True)
        self.set_dataset("QWP_angle", [], broadcast=True)


    @kernel
    def optimization_routine(self):
        self.core.reset()
        delay(1*s)

        # run stabilizer
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

        # it initializes to the best waveplate parameters. thus current hwp qwp angles are:
        previous_hwp, previous_qwp, previous_power = self.best_852HWP_to_max, self.best_852QWP_to_max, 0.0

        power = 0.0
        power_APD = 0.0

        tolerance_satisfied = False
        iteration = 0


        while not tolerance_satisfied:
        # for iteration in range(max_iterations):
            print("iteration no.", iteration) # this statement has to be here for it to run,,,, why...
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
                # todo: reduce delays
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


            iteration += 1

            delay(1 * s)
            print("iteration # ", iteration," : best_HWP, best_QWP, best_power = ", best_HWP,", ", best_QWP, ", ", best_power)

        # move back to the best HWP, QWP
        move_to_target_deg(self, name="852_HWP", target_deg=best_HWP)
        move_to_target_deg(self, name="852_QWP", target_deg=best_QWP)

        print("full range: ", full_range, "tolerance_deg: ", tolerance, "sample_pts: ", sample_pts)
        print("previous best_HWP, best_QWP, best_power = ", self.best_852HWP_to_max, ", ", self.best_852QWP_to_max, ", ", self.best_852_power)
        print("difference in best_HWP, best_QWP, best_power = ", (self.best_852HWP_to_max - best_HWP), ", ", (self.best_852QWP_to_max - best_QWP), ", ", (self.best_852_power - best_power))

        self.dds_FORT.sw.off()  ### turns FORT on

        self.set_dataset("best_852HWP_to_max", best_HWP, broadcast=True, persist=True)
        self.set_dataset("best_852QWP_to_max", best_QWP, broadcast=True, persist=True)
        self.set_dataset("best_852_power", best_power, broadcast=True, persist=True)


    @kernel
    def optimization_routine_zigzag(self):
        """
        It follows the same logic as the optimization_routine, with one key difference:
        as the routine proceeds to the next HWP value, it reverses the direction in which it scans the QWP values.
        For example, after completing the scan at hwp_values[0], the waveplates are positioned at (hwp_values[0], qwp_values[last_index]).
        For the next step at hwp_values[1], the scan starts from (hwp_values[1], qwp_values[last_index])
        and proceeds backward to (hwp_values[1], qwp_values[0]). This approach minimizes the total rotation time,
        since only the HWP needs to move to the next position, avoiding a full reset of the QWP to its starting angle.
        """

        self.core.reset()
        delay(1*s)

        # run stabilizer
        if self.enable_laser_feedback:
            self.laser_stabilizer.run()

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()  ### turns FORT on

        delay(2*s)   ## to ensure the FORT AOM power output is stabilized

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
        previous_hwp, previous_qwp, previous_power = best_HWP, best_QWP, 0.0

        power = 0.0
        power_APD = 0.0

        tolerance_satisfied = False
        iteration = 0


        # Iterative search loop
        while not tolerance_satisfied:
        # for iteration in range(max_iterations):
            print("starting iteration no.", iteration) # this statement has to be here for it to run,,,, why...
            time_ite = time_to_rotate_in_ms(self, half_range)

            delay(2 * time_ite * ms + 1*s)  # rotate two waveplates - okay for full_range = 50

            steps = half_range * 2 / (sample_pts - 1)

            # power difference with steps less than 1 is negligible.
            if steps < 1:
                tolerance_satisfied = True

                steps = 1.0
                sample_pts = int(half_range/steps) * 2 + 1 # how many sample pts to cover the range
                half_range = steps * (sample_pts - 1) / 2

            hwp_values = np.array([best_HWP - half_range + i * steps for i in range(sample_pts)])
            qwp_values = np.array([best_QWP - half_range + i * steps for i in range(sample_pts)])

            # if steps <= tolerance:
            #     tolerance_satisfied = True

            for hwp_i in range(len(hwp_values)):
                # time to rotate hwp
                current_hwp = hwp_values[hwp_i]
                time_hwp = time_to_rotate_in_ms(self, current_hwp - previous_hwp)
                delay(time_hwp * ms)
                count_down = 0

                stop_loop = False


                for qwp_i in range(len(qwp_values)):

                    if hwp_i % 2 != 0:
                        qwp_i = len(qwp_values) -1 -qwp_i
                    if not stop_loop:
                        current_qwp = qwp_values[qwp_i]
                        # time to rotate qwp
                        time_qwp = time_to_rotate_in_ms(self, current_qwp - previous_qwp)
                        delay(time_qwp * ms)

                        self.append_to_dataset("HWP_angle", current_hwp)
                        self.append_to_dataset("QWP_angle", current_qwp)

                        with parallel:
                            move_to_target_deg(self, name="852_HWP", target_deg=current_hwp)
                            move_to_target_deg(self, name="852_QWP", target_deg=current_qwp)

                        delay(1*s)

                        power = record_FORT_MM_power(self)
                        power_APD = record_FORT_APD_power(self)

                        if half_range > 5:
                            delay(1*s)          # deleting this for full_range=10 scan works.

                        if power > best_power:  # Update best power and angles if improved
                            best_power = power
                            best_HWP = current_hwp
                            best_QWP = current_qwp

                        previous_qwp = current_qwp
                        measurement_count += 1  # Increment measurement counter

                        if power < previous_power:
                            count_down += 1
                            if count_down > int(sample_pts / 2):
                                stop_loop = True

                        previous_power = power

                previous_hwp = current_hwp
                delay(1*s)

            half_range = steps


            delay(1 * s)
            print("iteration # ", iteration," : best_HWP, best_QWP, best_power = ", best_HWP,", ", best_QWP, ", ", best_power)
            iteration += 1

        # move back to the best HWP, QWP
        move_to_target_deg(self, name="852_HWP", target_deg=best_HWP)
        move_to_target_deg(self, name="852_QWP", target_deg=best_QWP)

        print("previous best_HWP, best_QWP, best_power = ", self.best_852HWP_to_max, ", ", self.best_852QWP_to_max, ", ", self.best_852_power)
        # print("difference in best_HWP, best_QWP, best_power = ", (self.best_852HWP_to_max - best_HWP), ", ", (self.best_852QWP_to_max - best_QWP), ", ", (self.best_852_power - best_power))
        delay(1*s)

        self.dds_FORT.sw.off()  ### turns FORT on

        self.set_dataset("best_852HWP_to_max", best_HWP, broadcast=True, persist=True)
        self.set_dataset("best_852QWP_to_max", best_QWP, broadcast=True, persist=True)
        self.set_dataset("best_852_power", best_power, broadcast=True, persist=True)

        # self.run_feedback_and_record_ref_power()  # recording reference power

    @kernel
    def run_feedback_and_record_ref_power(self):
        """
        Uses the `run_feedback_and_record_FORT_MM_power` function from `experiment_functions.py`
        to record a reference power value for detecting FORT polarization drift consistently.

        The function is defined in `experiment_functions.py` to ensure a unified implementation
        that can be reused across all experiments.

        Note: Since the optimization is performed with the FORT continuously ON,
        the recorded power may not be directly comparable to the power measured
        immediately after the FORT is turned ON. On Bob, it takes approximately
        2 seconds for the FORT power to stabilize, with ~5% fluctuation.

        This procedure is essential for reliable polarization optimization.
        """
        self.core.reset()

        delay(1 * s)
        power = run_feedback_and_record_FORT_MM_power(self)  # in experiment_functions

        print("After feedback - best_852_power set to ", power)
        self.set_dataset("best_852_power_ref", power, broadcast=True, persist=True)

    @kernel
    def testing(self):
        abc = get_rotator_deg(self, name="852_HWP")
        print("current HWP: ", abc)

    def run(self):
        """
        Step through the variable values defined by the scan sequences and run the experiment function.

        Because the scan variables can be any ExperimentVariable, which includes values used to initialize
        hardware (e.g. a frequency for a dds channel), the hardware is reinitialized in each step of the
        variable scan, i.e., each iteration.
        """
        self.core.reset()
        self.prepare()
        self.initialize_hardware()
        self.initialize_datasets()

        self.optimization_routine_zigzag()  # running zigzag as default

        self.run_feedback_and_record_ref_power()  # recording reference power




