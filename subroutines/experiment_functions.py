from artiq.experiment import *
from artiq.coredevice.ad9910 import (
    PHASE_MODE_ABSOLUTE,
    PHASE_MODE_CONTINUOUS,
    PHASE_MODE_TRACKING,
    RAM_DEST_ASF,
    RAM_MODE_RAMPUP
)
from artiq.coredevice.urukul import CFG_MASK_NU
from artiq.language import us, ns, MHz
import logging
import numpy as np
import math
import pyvisa as visa

import os, sys
cwd = os.getcwd() + "\\"

sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.conversions import dB_to_V_kernel as dB_to_V
from subroutines.k10cr1_functions import *

"""
Table of contents:
1. Subroutine functions
2. Experiment functions
3. Diagnostic functions
"""

###############################################################################
# 1. SUBROUTINE FUNCTIONS
# These are functions that get used within various experiments. they are not
# intended to be run in a standalone fashion, e.g. from GeneralVariableScan.
# Consequently, note that the name should not end in "experiment"
###############################################################################

@kernel
def run_feedback_and_record_FORT_MM_power(self):
    """
    Used for FORT polarization check.
    * PolOptTest also uses this to record the reference value "best_852_power_ref"

    The function is defined in `experiment_functions.py` to ensure a unified implementation
    that can be reused across all experiments.

    Note: Since the optimization is performed with the FORT continuously ON,
    the recorded power may not be directly comparable to the power measured
    immediately after the FORT is turned ON. On Bob, it takes approximately
    2 seconds for the FORT power to stabilize, with ~5% fluctuation.

    This procedure is essential for reliable polarization optimization.
    """
    self.laser_stabilizer.run()

    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
    self.dds_FORT.sw.on()  ### turns FORT on
    delay(0.1 * ms)

    # todo: in record_FORT_MM_power function, power is recorded in "FORT_MM_monitor" dataset.
    #       I think I'll keep this.
    #       Just make another dataset for optimization
    power = record_FORT_MM_power(self)
    record_FORT_APD_power(self)
    self.dds_FORT.sw.off()

    return power

@kernel
def FORT_optimization_routine_zigzag(self):
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

    self.dds_FORT.sw.off()  ### turns FORT on

    self.set_dataset("best_852HWP_to_max", best_HWP, broadcast=True, persist=True)
    self.set_dataset("best_852QWP_to_max", best_QWP, broadcast=True, persist=True)
    self.set_dataset("best_852_power", best_power, broadcast=True, persist=True)

@kernel
def FORT_run_feedback_and_record_ref_power(self):
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

    delay(0.1 * ms)
    power = run_feedback_and_record_FORT_MM_power(self)  # in experiment_functions

    print("After feedback - best_852_power set to ", power)
    self.set_dataset("best_852_power_ref", power, broadcast=True, persist=True)

def FORT_polarization_check_and_optimize(self):

    power_MM = 0.0
    difference_in_power_in_perc = 0.0

    power_MM = run_feedback_and_record_FORT_MM_power(self)
    record_FORT_APD_power(self)

    # difference_in_power: normalized over the setpoint, in %.
    #todo: if best_852_power_ref = 0
    difference_in_power_in_perc = (power_MM - self.best_852_power_ref) / self.best_852_power_ref * 100

    full_range = 0
    # todo: think about the full_range
    # Only run FORT pol stabilizer if the FORT power after polarizer is less than 5% of the setpoint.
    if difference_in_power_in_perc < -2:

        self.print_async("Running FORT Polarization optimization: difference_in_power_in_perc < -5")

        FORT_not_optimized = True

        while FORT_not_optimized:

            if difference_in_power_in_perc < -30:
                full_range = 90.0
                sample_pts = 7
                # self.set_dataset("full_range", 90.0, broadcast=True)
            elif difference_in_power_in_perc < -10:
                full_range = 45.0
                sample_pts = 7
                # self.set_dataset("full_range", 45.0, broadcast=True)
            elif difference_in_power_in_perc < -5:
                full_range = 30.0
                sample_pts = 7
                # self.set_dataset("full_range", 30.0, broadcast=True)
            else:
                full_range = 20.0
                sample_pts = 5
                # self.set_dataset("full_range", 20.0, broadcast=True)

            delay(1*ms)
            self.set_dataset("full_range", full_range, broadcast=True)
            self.set_dataset("sample_pts", sample_pts, broadcast=True)


            FORT_optimization_routine_zigzag(self)

            previous_ref = self.best_852_power_ref
            FORT_run_feedback_and_record_ref_power(self)
            current_ref = self.best_852_power_ref

            difference_in_power_in_perc = (current_ref - previous_ref) / previous_ref * 100

            if difference_in_power_in_perc > -2:
                FORT_not_optimized = False
            else:
                self.set_dataset("best_852_power_ref", previous_ref, broadcast=True, persist=True)


    else:
        self.print_async("Skipping FORT Polarization optimization. difference_in_power_in_perc > -5%")


    self.dds_FORT.sw.off()



    ###############  same upto here with def run_feedback_and_record_FORT_MM_power(self)

@kernel
def FORT_power_stabilzation_test_experiment(self):

    # run the stabilizer
    self.laser_stabilizer.run()

    # Turn on the FORT with stabilized power
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
    self.dds_FORT.sw.on()
    delay(1*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.stabilizer_FORT.run(setpoint_index=1)
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

@kernel
def rotator_test_experiment(self):
    """
    ratator_test function to record the Sampler value while rotating the waveplate

    :param self:
    :return:
    """

    self.core.reset()
    delay(2* s)


    # # GVS - set target_hwp_deg as the scan_variable.
    # move_to_target_deg(self, name="780_HWP", target_deg=self.target_hwp_deg)
    # wait_move(self, "780_HWP")
    # hwp780_pos = get_rotator_position(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_pos / self.deg_to_pos, ' deg')
    # delay(1*s)
    # record_PDA_power(self)


    move_to_target_deg(self, name="780_HWP", target_deg=10)
    hwp780_deg = get_rotator_deg(self, '780_HWP')
    delay(1*s)
    self.print_async('hwp780 at ', hwp780_deg, ' deg')

    move_to_target_deg(self, name="780_HWP", target_deg=10)
    hwp780_deg = get_rotator_deg(self, '780_HWP')
    delay(1*s)
    self.print_async('hwp780 at ', hwp780_deg, ' deg')

    # # delay(2*s)
    # # move_to_target_deg(self, name="780_HWP", target_deg=350)
    #
    # go_to_home(self, '780_HWP')
    #
    # hwp780_deg = get_rotator_deg(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_deg, ' deg')

@kernel
def shot_without_measurement(self):
    """
    non-chopped first atom readout.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    :return:
    """
    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms) ## coils relaxation time

    ### set the FORT AOM to the readout settings
    if self.which_node == 'alice':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

    self.ttl_repump_switch.off() ### turn on MOT RP
    self.dds_cooling_DP.sw.on() ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(self.t_SPCM_first_shot)
    self.dds_cooling_DP.sw.off()

@kernel
def load_MOT_and_FORT(self):
    """
    The FORT loading sequence, from loading a MOT to hopefully trapping one atom

    Loads the MOT then the turns on the FORT and waits a fixed amount of time
    to load an atom. Optionally, the FORT can be turned on right away if
    FORT_on_at_start is True.

    Turns ON the following at the beginning:
        FORT AOM
        Cooling DP
        All fiber AOMs
        MOT RP

    Leaves the following OFF at the end:
        Cooling DP
        MOT RP


    :param self: the experiment instance
    :return:
    """
    #todo: [changes] getting read of FORT_on_at_MOT_start
    self.dds_FORT.sw.on()
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    delay(1 * ms)

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    # self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

    # self.ttl7.on()

    self.dds_cooling_DP.sw.on()
    self.ttl_repump_switch.off()

    delay(1*ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    self.ttl_UV.pulse(self.t_UV_pulse)

    ### wait for the MOT to load
    delay(self.t_MOT_loading)

    # zotino_stability_test(self)


    # todo: make this work for bob too.
    if self.which_node == 'alice':
        self.stabilizer_FORT.run(setpoint_index=1) # the science setpoint
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_RO)

    self.dds_cooling_DP.sw.off()
    self.ttl_repump_switch.on()
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)
    delay(1*ms)

    # self.ttl7.off()

    if self.do_PGC_after_loading and self.t_PGC_after_loading > 0:
        self.dds_cooling_DP.sw.on()
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        self.zotino0.set_dac([self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                             channels=self.coil_channels)

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_PGC)

        delay(self.t_PGC_after_loading)
        self.dds_cooling_DP.sw.off()

@kernel
def load_MOT_and_FORT_until_atom(self):
    """
    Turning on the MOT and FORT light at the same time and monitor SPCM0. Turn off the MOT as soon as an atom is trapped.

    Turns ON the following at the beginning:
        FORT AOM
        Cooling DP
        All fiber AOMs
        MOT RP

    Leaves the following OFF at the end:
        Cooling DP
        MOT RP


    :param self: the experiment instance
    :return:
    """

    if self.monitors_for_atom_loading:
        measure_Magnetometer(self)
        delay(1 * ms)
        Sampler0_test(self)
        delay(1 * ms)
        measure_coil_driver(self)

    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    ### set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
    delay(0.1 * ms)

    self.dds_cooling_DP.sw.on()  ### turn on cooling
    self.ttl_repump_switch.off()  ### turn on MOT RP

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    delay(0.1 * ms)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    self.dds_FORT.sw.off()
    delay(1*ms)

    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
    self.dds_FORT.sw.on()

    delay(1 * ms)
    self.ttl_UV.on()

    max_tries = 100  ### Maximum number of attempts before running the feedback
    atom_check_time = self.t_atom_check_time
    atom_loaded = False
    try_n = 0
    t_before_atom = now_mu() ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
    t_after_atom = now_mu()
    time_without_atom = 0.0

    self.ttl7.on()  ### for triggering oscilloscope

    while True:
        while not atom_loaded and try_n < max_tries:
            delay(100 * us)  ### Needs a delay of about 100us or maybe less
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_counter.gate_rising(atom_check_time)

            SPCM0_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_atom_check = self.ttl_SPCM1_counter.fetch_count()

            BothSPCMs_atom_check = int((SPCM0_atom_check + SPCM1_atom_check) / 2)

            try_n += 1

            if BothSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
            self.append_to_dataset("BothSPCMs_atom_check_in_loading",BothSPCMs_atom_check)

        if atom_loaded:
            # t_before_atom = t_after_atom ### I don't know why I had this! Removed and seems working fine.
            self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
            t_after_atom = now_mu()
            break  ### Exit the outer loop if an atom is loaded

        #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
        delay(0.1 * ms)
        t_no_atom = now_mu()
        time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
        self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

        ### If max_tries reached and still no atom, run feedback
        if self.enable_laser_feedback:
            # self.ttl7.pulse(10 * ms)  ### for triggering oscilloscope
            delay(0.1 * ms) ### necessary to avoid underflow

            ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            self.laser_stabilizer.run()
            self.n_feedback_per_iteration += 1
            # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
            self.dds_microwaves.sw.on()
            self.dds_FORT.sw.on()
            delay(0.1 * ms)

            try_n = 0

            # delay(10 * ms)
            # print("**************   No atom after 2 seconds. Running feedback   ***************")
            # delay(10 * ms)

    self.ttl_UV.off()

    ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)

    self.ttl_repump_switch.on()  ### turn off MOT RP
    self.dds_cooling_DP.sw.off()  ### turn off cooling

    delay(1 * ms)
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

    ###########  PGC on the trapped atom  #############
    if self.do_PGC_after_loading:
        if self.which_node == 'alice':
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        elif self.which_node == 'bob':
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)
        ### set the cooling DP AOM to the PGC settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        self.ttl_repump_switch.off()  ### turn on MOT RP
        self.dds_cooling_DP.sw.on()  ### turn on cooling
        delay(10 * us)
        # self.dds_AOM_A5.sw.off()
        # self.dds_AOM_A6.sw.off()
        delay(self.t_PGC_after_loading)  ### this is the PGC time
    ###################################################
    self.ttl7.off()

    ### I don't know what this SPCM0_FORT_science is used for. Set to 0 for now:
    self.SPCM0_FORT_science = 0
    # t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    # self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)

    ### saving the atom loading time for each loaded atom.
    self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
    self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
    delay(1 * ms)
    self.append_to_dataset("atom_loading_wall_clock", now_mu())  ### just to plot Atom_loading_time vs actual time in analysis
    self.n_atom_loaded_per_iteration += 1
    delay(1 * ms)

@kernel
def load_atom_smooth_FORT(self):
    """
    The same as load_MOT_and_FORT, but with FORT ramped down to science setpoint smoothly to cool the atom.
    Turning on the MOT and FORT light at the same time and monitor SPCM0. Turn off the MOT as soon as an atom is trapped.

    Works fine.
    """

    if self.monitors_for_atom_loading:
        measure_Magnetometer(self)
        delay(1 * ms)
        Sampler0_test(self)
        delay(1 * ms)
        measure_coil_driver(self)

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    ### set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    self.dds_cooling_DP.sw.on()
    self.ttl_repump_switch.off()

    delay(1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    delay(10*us)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    self.ttl_UV.on()

    self.ttl7.on()

    ############################################  Turning on FORT in RAM mode
    delay(1*ms)

    ### Enabling RAM playback, not playing yet.
    self.dds_FORT.set_cfr1(ram_enable=1,
                           ram_destination=RAM_DEST_ASF)

    ### turns on the FORT to loading set point
    self.dds_FORT.set_profile_ram(
        start=self.FORT_total_points // 2,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP)
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    self.dds_FORT.sw.on()
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    ########################################  End of Turning on FORT

    ### wait for the MOT to load
    delay(self.t_MOT_loading)

    self.ttl_UV.off()
    self.ttl7.off()

    self.dds_cooling_DP.sw.off()
    self.ttl_repump_switch.on()
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate
    delay(1 * ms)

@kernel
def load_until_atom_smooth_FORT(self):
    """
    The same as load_MOT_and_FORT_until_atom, but with FORT ramped down to science setpoint smoothly to cool the atom.
    """

    if self.monitors_for_atom_loading:
        measure_Magnetometer(self)
        delay(1 * ms)
        Sampler0_test(self)
        delay(1 * ms)
        measure_coil_driver(self)

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    ### set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    self.dds_cooling_DP.sw.on()
    self.ttl_repump_switch.off()

    delay(1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    delay(10*us)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    self.ttl_UV.on()

    self.ttl7.on()

    self.dds_FORT.sw.off()  ### to lose any trapped atom and trap a new atom

    ############################################  Turning on FORT in RAM mode
    delay(1*ms)
    ### Enabling RAM playback, not playing yet.
    self.dds_FORT.set_cfr1(ram_enable=1,
                           ram_destination=RAM_DEST_ASF)

    ### turns on the FORT to loading set point
    self.dds_FORT.set_profile_ram(
        start=self.FORT_total_points - 1,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP)
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    self.dds_FORT.sw.on()
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    ########################################  End of Turning on FORT

    max_tries = 100  ### Maximum number of attempts before running the feedback
    atom_check_time = self.t_atom_check_time
    atom_loaded = False
    try_n = 0
    t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
    t_after_atom = now_mu()
    time_without_atom = 0.0

    while True:
        while not atom_loaded and try_n < max_tries:
            delay(100 * us)  ### Needs a delay of about 100us or maybe less
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_counter.gate_rising(atom_check_time)

            SPCM0_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_atom_check = self.ttl_SPCM1_counter.fetch_count()

            BothSPCMs_atom_check = int((SPCM0_atom_check + SPCM1_atom_check) / 2)

            try_n += 1

            if BothSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
            self.append_to_dataset("BothSPCMs_atom_check_in_loading", BothSPCMs_atom_check)

        if atom_loaded:
            # t_before_atom = t_after_atom ### I don't know why I had this! Removed and seems working fine.
            self.set_dataset("time_without_atom", 0.0, broadcast=True)  ### resetting time_without_atom when we load an atom
            t_after_atom = now_mu()
            break  ### Exit the outer loop if an atom is loaded

        #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
        delay(0.1 * ms)
        t_no_atom = now_mu()
        time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
        self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

        ### If max_tries reached and still no atom, run feedback
        if self.enable_laser_feedback:
            ### exit RAM
            self.dds_FORT.set_cfr1(ram_enable=0)
            self.dds_FORT.cpld.io_update.pulse_mu(8)

            delay(0.1 * ms)  ### necessary to avoid underflow
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            self.laser_stabilizer.run()
            self.n_feedback_per_iteration += 1
            delay(0.1 * ms)

            try_n = 0

            ############################################  Turning on FORT in RAM mode
            delay(1 * ms)
            ### Enabling RAM playback, not playing yet.
            self.dds_FORT.set_cfr1(ram_enable=1,
                                   ram_destination=RAM_DEST_ASF)

            ### turns on the FORT to loading set point
            self.dds_FORT.set_profile_ram(
                start=self.FORT_total_points // 2,
                end=self.FORT_total_points - 1,
                step=self.FORT_step_size,
                profile=7,
                mode=RAM_MODE_RAMPUP)
            self.dds_FORT.cpld.io_update.pulse_mu(8)
            self.dds_FORT.sw.on()
            self.dds_FORT.cpld.io_update.pulse_mu(8)
            ########################################  End of Turning on FORT


    self.ttl_UV.off()

    ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)

    self.ttl_repump_switch.on()  ### turn off MOT RP
    self.dds_cooling_DP.sw.off()  ### turn off cooling

    delay(1 * ms)
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

    # self.ttl7.on()  ### for triggering oscilloscope

    # ###########  PGC on the trapped atom  #############
    # if self.do_PGC_after_loading:
    #     if self.which_node == 'alice':
    #         self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    #     elif self.which_node == 'bob':
    #         self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)
    #     ### set the cooling DP AOM to the PGC settings
    #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
    #     self.ttl_repump_switch.off()  ### turn on MOT RP
    #     self.dds_cooling_DP.sw.on()  ### turn on cooling
    #     delay(10 * us)
    #     # self.dds_AOM_A5.sw.off()
    #     # self.dds_AOM_A6.sw.off()
    #     delay(self.t_PGC_after_loading)  ### this is the PGC time
    # ###################################################
    self.ttl7.off()

    ### I don't know what this SPCM0_FORT_science is used for. Set to 0 for now:
    self.SPCM0_FORT_science = 0
    # t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    # self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)

    ### saving the atom loading time for each loaded atom.
    self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
    self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
    delay(1 * ms)
    self.append_to_dataset("atom_loading_wall_clock", now_mu())  ### just to plot Atom_loading_time vs actual time in analysis
    self.n_atom_loaded_per_iteration += 1
    delay(1 * ms)

@kernel
def load_MOT_and_FORT_until_atom_recycle(self):
    """
    Before attempting to load MOT and FORT, it checks if there is already an atom in the FORT. If not, then it turns on the MOT and FORT
    light at the same time and monitor SPCM0. Turn off the MOT as soon as an atom is trapped.

    Turns ON the following at the beginning:
        FORT AOM
        Cooling DP
        All fiber AOMs
        MOT RP

    Leaves the following OFF at the end:
        Cooling DP
        MOT RP

    :param self: the experiment instance
    :return:
    """

    ### First check if there is already an atom in the FORT based on RO2
    delay(100 * us)
    if self.measurement > 0:
        if self.BothSPCMs_RO2/self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True

            ###########  PGC on the trapped atom  #############
            if self.do_PGC_after_loading:
                ### Set the coils to PGC setting
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)
                delay(0.4 * ms)
                ### set the cooling DP AOM to the PGC settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                delay(0.1 * ms)
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                self.dds_cooling_DP.sw.on()  ### turn on cooling
                self.ttl_repump_switch.off()  ### turn on MOT RP
                delay(self.t_PGC_after_loading)  ### this is the PGC time
                self.dds_cooling_DP.sw.off()  ### turn off cooling
                self.ttl_repump_switch.on()  ### turn off MOT RP
            ###################################################
        else:
            atom_loaded = False
    else:
        atom_loaded = False


    ### load an atom if atom_loaded = False
    if not atom_loaded:

        if self.monitors_for_atom_loading:
            measure_Magnetometer(self)
            delay(1*ms)
            Sampler0_test(self)
            delay(1*ms)
            measure_coil_driver(self)

        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(0.1 * ms)
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()

        delay(1 * ms)
        self.ttl_UV.on()

        max_tries = 100  ### Maximum number of attempts before running the feedback
        # atom_check_time = 20 * ms
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu() ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        while True:
            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)

                BothSPCMs_atom_check = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)

                try_n += 1

                if BothSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

                ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
                self.append_to_dataset("BothSPCMs_atom_check_in_loading",BothSPCMs_atom_check)
                delay(1*ms)

            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
                t_after_atom = now_mu()
                break  ### Exit the outer loop if an atom is loaded

            #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
            delay(0.1 * ms)
            t_no_atom = now_mu()
            time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
            self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms) ### necessary to avoid underflow

                # delay(0.1 * ms)
                # self.ttl7.pulse(100 * us)  ### for triggering oscilloscope
                # delay(0.1 * ms)

                ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                self.laser_stabilizer.run()
                self.n_feedback_per_iteration += 1
                # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
                self.dds_microwaves.sw.on()
                self.dds_FORT.sw.on()
                delay(0.1 * ms)

                try_n = 0

                # delay(10 * ms)
                # print("**************   No atom after 2 seconds. Running feedback   ***************")
                # delay(10 * ms)

        self.ttl_UV.off()

        ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(0.4 * ms)

        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling

        delay(1 * ms)
        delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

        ###########  PGC on the trapped atom  #############
        if self.do_PGC_after_loading:
            self.ttl_repump_switch.off()  ### turn on MOT RP
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            delay(10 * us)
            if self.which_node == 'alice':
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[0])
            elif self.which_node == 'bob':
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)
            ### set the cooling DP AOM to the PGC settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
            # self.dds_AOM_A5.sw.off()
            # self.dds_AOM_A6.sw.off()
            delay(self.t_PGC_after_loading)  ### this is the PGC time
        ###################################################

        ### I don't know what this SPCM0_FORT_science is used for. Set to 0 for now:
        self.SPCM0_FORT_science = 0
        # t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
        # self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)

        ### saving the atom loading time for each loaded atom.
        self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
        self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
        delay(1 * ms)
        self.append_to_dataset("atom_loading_wall_clock", now_mu()) ### just to plot Atom_loading_time vs actual time in analysis
        self.n_atom_loaded_per_iteration += 1
        delay(10 * ms)

@kernel
def load_MOT_and_FORT_for_Luca_scattering_measurement(self):
    """
    Modified version of load_MOT_and_FORT for imaging scattering in the chamber with the Luca
    and with an APD. Assumes the APD is connected to the Sampler card and channel specified in
    the code below.

    I typically run this with the MonitorFORTWithLuca experiment.

    The following shots are taken with the Luca:
    1. FORT scattering after feedback to loading setpoint
    2. FORT + MOT beam scattering

    :param self: the experiment instance
    :return:
    """

    self.dds_FORT.sw.on()

    if self.MOT_repump_off or self.MOT_light_off: # this is useful so that the SPCM0_RO1 dataset is all background
        self.ttl_repump_switch.on()  # turns the RP AOM off
    if self.MOT_light_off:
        self.dds_cooling_DP.sw.off()

    if not self.FORT_on_at_MOT_start:
        self.dds_FORT.sw.off()
    else:
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    # record FORT scattering with Luca and record Raman scattering from SM fiber
    self.ttl_Luca_trigger.pulse(5 * ms) # FORT loading scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_loading = self.ttl_SPCM0.count(t_gate_end)
    delay(10*us)

    self.APD_FORT_volts_loading = 0.0
    for i in range(self.APD_averages):
        self.sampler1.sample(self.APD_buffer)
        self.APD_FORT_volts_loading += self.APD_buffer[7]
        delay(1*ms)
    self.APD_FORT_volts_loading /= self.APD_averages

    delay(1 * ms)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    # self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

    if not self.MOT_light_off:
        self.dds_cooling_DP.sw.on()

        delay(1*ms)

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

    delay(1*ms) # if this delay is not here, the following line setting the dac doesn't execute

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    # wait for the MOT to load
    delay(self.t_MOT_loading/2)
    self.ttl_Luca_trigger.pulse(5 * ms) # total scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_and_MOT = self.ttl_SPCM0.count(t_gate_end)
    delay(1 * ms)
    delay(self.t_MOT_loading/2)

    # turn on the dipole trap and wait to load atoms if it is not already on
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    if not self.FORT_on_at_MOT_start:
        delay_mu(self.t_FORT_loading_mu)
    delay(1*ms)
    if not self.no_feedback:
        self.stabilizer_FORT.run(setpoint_index=1) # set to  science setpoint

    self.dds_cooling_DP.sw.off()
    delay(1*ms)
    # record FORT scattering with Luca and record Raman scattering from SM fiber
    self.ttl_Luca_trigger.pulse(5 * ms)  # FORT loading scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)
    delay(10*us)

    self.APD_FORT_volts_science = 0.0
    for i in range(self.APD_averages):
        self.sampler1.sample(self.APD_buffer)
        self.APD_FORT_volts_science += self.APD_buffer[7]
        delay(1*ms)

    self.APD_FORT_volts_science /= self.APD_averages
    delay(100*ms)  # this is 100 ms to ensure the Luca takes the next shot

    if not self.MOT_light_off:
        self.dds_cooling_DP.sw.on()

    if self.do_PGC_after_loading and self.t_PGC_after_loading > 0:

        if not self.no_feedback:
            self.dds_FORT.set(frequency=self.f_FORT,
                              amplitude=self.stabilizer_FORT.amplitudes[1])

        self.zotino0.set_dac([self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                             channels=self.coil_channels)

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_PGC)
        delay(self.t_PGC_after_loading)

    self.dds_cooling_DP.sw.off()

@kernel
def first_shot_chopped(self):
    """
    chopped first atom readout.

    # Turns on:
    #     Cooling DP
    #     MOT RP
    #     All 6 fiber AOMs
    #
    # Turns off at the end:
    #     Cooling DP
    #     MOT RP

    :return:
    """
    ### set the coils to PGC settings

    n_RO_total = 500
    self.SPCM0_test_RO = 0

    # self.zotino0.set_dac(
    #     [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
    #     channels=self.coil_channels)
    # delay(0.4 * ms) ## coils relaxation time
    #
    # ### set the FORT AOM to the readout settings
    # self.dds_FORT.set(frequency=self.f_FORT,
    #                   amplitude=self.stabilizer_FORT.amplitudes[1])
    #
    # ### set the cooling DP AOM to the readout settings
    # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
    #                         amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)
    #
    # delay(10 * us)
    # self.dds_AOM_A1.sw.on()
    # self.dds_AOM_A2.sw.on()
    # self.dds_AOM_A3.sw.on()
    # self.dds_AOM_A4.sw.on()
    # self.dds_AOM_A5.sw.on()
    # self.dds_AOM_A6.sw.on()
    #
    # delay(10 * us)
    # self.dds_cooling_DP.sw.off()  ### turn off cooling
    # self.ttl_repump_switch.on()  ### turn off MOT RP

    self.ttl7.on()
    delay(100 * us)

    for n_RO in range(n_RO_total):
        t1 = now_mu()

        self.dds_FORT.sw.off()  ### Turn off FORT
        # self.ttl_repump_switch.off()  ### turn on MOT RP
        # self.dds_cooling_DP.sw.on()  ### Turn on cooling

        # at_mu(t1 + 200)
        # self.ttl_SPCM0_counter.gate_rising(800 * ns)



        at_mu(t1 + 700)
        self.dds_FORT.sw.on()  ### Turn on FORT
        # self.dds_cooling_DP.sw.off()  ### turn off cooling
        # self.ttl_repump_switch.on()  ### turn off MOT RP

        # SPCM0_test_RO_n = self.ttl_SPCM0_counter.fetch_count()

        delay(5 * us)
        # self.SPCM0_test_RO += SPCM0_test_RO_n

    self.ttl7.off()

@kernel
def first_shot(self):
    """
    non-chopped first atom readout.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    :return:
    """
    ### set the FORT AOM to the readout settings
    if self.which_node == 'alice':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms) ## coils relaxation time

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

    self.ttl_repump_switch.off() ### turn on MOT RP
    self.dds_cooling_DP.sw.on() ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.which_node != 'alice':  # edge counters only enabled on Alice gateware so far
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(t_gate_end)
        self.BothSPCMs_RO1 = int(self.SPCM0_RO1) #todo: connect a working SPCM1
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()
    else:

        if self.use_chopped_readout:
            ro_dma_handle = self.core_dma.get_handle("first_chopped_readout")
            delay(10 * ms)

            self.ttl_repump_switch.off()  # turns the RP AOM on
            self.dds_cooling_DP.sw.off()  # the chop sequence likes to turn the FORT off

            delay(1 * ms)
            # self.ttl7.pulse(100 * us)
            self.dds_FORT.sw.on()  # the chop sequence likes to turn the FORT off

            delay(0.1 * ms)

            # we want to initiate the chop playback and read in detector clicks while the chop sequence is playing.
            # the ttl.gate_rising(duration) function is equivalent to:
            #     ttl._set_sensitivity(1)
            #     delay(duration)
            #     ttl._set_sensitivity(0)
            #     return now_mu()
            #
            # we want the dma playback to happen during the gating, so we call the _set_sensitivity functions directly

            self.ttl_SPCM0._set_sensitivity(1)
            self.core_dma.playback_handle(ro_dma_handle)
            self.ttl_SPCM0._set_sensitivity(0)
            self.SPCM0_RO1 = self.ttl_SPCM0.count(now_mu())

        else:
            delay(1 * ms)
            # self.ttl7.pulse(100 * us)
            delay(1 * ms)
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_first_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_first_shot)

            self.SPCM0_RO1 = self.ttl_SPCM0_counter.fetch_count()
            self.SPCM1_RO1 = self.ttl_SPCM1_counter.fetch_count()
            self.BothSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)
            delay(0.1 * ms)
            self.dds_cooling_DP.sw.off() ### turn off cooling
            self.ttl_repump_switch.on() ### turn off MOT RP
            delay(10 * us)

@kernel
def second_shot(self):
    """
    non-chopped second atom readout.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    warning: assumes the fiber AOMs are already on, which is usually the case
    :return:
    """
    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)  ## coils relaxation time


    ### set the FORT AOM to the readout settings
    if self.which_node == 'alice':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)



    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.which_node != 'alice':  # edge counters only enabled on Alice gateware so far
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)
        self.BothSPCMs_RO2 = int(self.SPCM0_RO2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()

    else:
        if self.use_chopped_readout:
            # rtio_log("chop_RO_counter", 0)
            # rtio_log("chop_RO_dma", 0)
            delay(10*ms)
            ro_dma_handle = self.core_dma.get_handle("second_chopped_readout")
            delay(10 * ms)

            # todo set RO coils here?

            self.ttl_repump_switch.off()  # turns the RP AOM on
            self.dds_cooling_DP.sw.off()  # the chop sequence likes to turn the FORT off

            delay(1 * ms)
            # self.ttl7.pulse(100 * us)
            self.dds_FORT.sw.on()  # the chop sequence likes to turn the FORT off

            delay(10 * ms)

            # we want to initiate the chop playback and read in detector clicks while the chop sequence is playing.
            # the edge_counter.gate_rising(duration) function is equivalent to:
            #         edge_counter.set_config(
            #             count_rising=count_rising,
            #             count_falling=count_falling,
            #             send_count_event=False,
            #             reset_to_zero=True)
            #         delay_mu(duration_mu)
            #         edge_counter.set_config(
            #             count_rising=False,
            #             count_falling=False,
            #             send_count_event=True,
            #             reset_to_zero=False)
            #
            # we want the dma playback to happen during the gating, so we call the _set_sensitivity functions directly

            now = now_mu()
            rtio_log("chop_RO_counter", 1)
            self.ttl_SPCM0_counter.set_config(
                    count_rising=True,
                    count_falling=False,
                    send_count_event=False,
                    reset_to_zero=True)
            delay(1*us)
            rtio_log("chop_RO_dma", 1)
            self.core_dma.playback_handle(ro_dma_handle) # not sure if I need extra delay here
            rtio_log("chop_RO_dma", 0)
            at_mu(now + self.core.seconds_to_mu(self.t_SPCM_second_shot+10*us))
            self.ttl_SPCM0_counter.set_config(
                    count_rising=False,
                    count_falling=False,
                    send_count_event=True,
                    reset_to_zero=False)
            self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
            rtio_log("chop_RO_counter", 0)
            delay(10*ms)

        else:
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

            self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
            self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
            self.BothSPCMs_RO2 = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)

            delay(0.1 * ms)
            self.dds_cooling_DP.sw.off() ### turn off cooling
            self.ttl_repump_switch.on()  ### turn off MOT RP
            delay(10 * us)

    ### set the FORT AOM back to loading setting
    if self.which_node == 'alice':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[0])
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

@kernel
def first_shot_smooth_FORT(self):
    """
    non-chopped first atom readout with smooth FORT

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    :return:
    """
    # self.ttl7.on()

    self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg | 1 << CFG_MASK_NU + 0)  # Mask the DDS channel which is on RAM mode
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)
    self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg & ~(1 << CFG_MASK_NU + 0))  # Unmask the DDS channel which is on RAM mode


    ### Configure the RAM to playback the first half (ramping down to science point)
    self.dds_FORT.set_profile_ram(
        start=0,
        end=len(self.FORT_amplitudes_list) // 2 - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP)
    self.dds_FORT.cpld.io_update.pulse_mu(8)

    delay(self.t_FORT_ramp)

    # self.ttl7.off()


    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms) ## coils relaxation time



    self.ttl_repump_switch.off() ### turn on MOT RP
    self.dds_cooling_DP.sw.on() ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.which_node != 'alice':  # edge counters only enabled on Alice gateware so far
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(t_gate_end)
        self.BothSPCMs_RO1 = int(self.SPCM0_RO1) #todo: connect a working SPCM1
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()
    else:
        # self.ttl7.pulse(100 * us)
        delay(1 * ms)
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_first_shot)
            self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_first_shot)

        self.SPCM0_RO1 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO1 = self.ttl_SPCM1_counter.fetch_count()
        self.BothSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off() ### turn off cooling
        self.ttl_repump_switch.on() ### turn off MOT RP
        delay(10 * us)

@kernel
def second_shot_smooth_FORT(self):
    """
    non-chopped second atom readout with smooth FORT.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    warning: assumes the fiber AOMs are already on, which is usually the case
    :return:
    """
    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)  ## coils relaxation time

    # self.ttl7.on()

    ### set the cooling DP AOM to the readout settings
    self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg | 1 << CFG_MASK_NU + 0)  # Mask the DDS channel which is on RAM mode
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)
    self.dds_FORT.cpld.cfg_write(self.dds_FORT.cpld.cfg_reg & ~(1 << CFG_MASK_NU + 0))  # Unmask the DDS channel which in on RAM mode


    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(0.1 * ms)

    # self.ttl7.off()

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.which_node != 'alice':  # edge counters only enabled on Alice gateware so far
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)
        self.BothSPCMs_RO2 = int(self.SPCM0_RO2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()

    else:
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
            self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

        self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
        self.BothSPCMs_RO2 = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)

        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off() ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        delay(10 * us)

    ### Ramping up the FORT to loading amplitude
    self.dds_FORT.set_profile_ram(
        start=self.FORT_total_points // 2,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP)
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    delay(self.t_FORT_ramp)

    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    ### exit RAM
    self.dds_FORT.set_cfr1(ram_enable=0)
    self.dds_FORT.cpld.io_update.pulse_mu(8)

@kernel
def atom_parity_shot(self):
    """
    non-chopped atom readout used in atom parity oscillation experiments. Uses the second_shot parameters.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    warning: assumes the fiber AOMs are already on, which is usually the case
    :return:
    """
    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)  ## coils relaxation time


    ### set the FORT AOM to the readout settings
    if self.which_node == 'alice':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    elif self.which_node == 'bob':
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)



    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.which_node != 'alice':  # edge counters only enabled on Alice gateware so far
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)
        self.BothSPCMs_parity_RO = int(self.SPCM0_RO2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()

    else:
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
            self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

        self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
        self.BothSPCMs_parity_RO = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)

        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        delay(10 * us)

@kernel
def record_chopped_readout(self, readout_duration: TFloat, label: TStr):
    """

    :param self:
    :return:
    """
    # todo. copy OP chopping of the FORT but also chop the SPCM gate
    # todo: still uses ttl_SPCM_gate which was removed on 2025-01-21. Use normal ttl.count without a switch.

    n_chop_cycles = int(readout_duration / self.t_RO_chop_period + 0.5)
    assert n_chop_cycles >= 1, "t_SPCM_first_shot should be > t_RO_chop_period"

    RO_pulse = self.t_RO_chop_period * self.duty_cycle_RO
    FORT_pulse = self.t_RO_chop_period * self.duty_cycle_FORT

    self.core.reset()

    start = now_mu()
    period_mu = self.core.seconds_to_mu(self.t_RO_chop_period)

    # does this need to be set within the dma context manager?
    RO_pulse_length_mu = self.core.seconds_to_mu(RO_pulse)
    FORT_pulse_length_mu = self.core.seconds_to_mu(FORT_pulse)
    FORT_on_mu = self.core.seconds_to_mu(0.0)
    RO_on_mu = self.core.seconds_to_mu(self.t_RO_chop_offset)
    gate_on_mu = self.core.seconds_to_mu(self.t_RO_gate_offset)

    with self.core_dma.record(label):

        # if not self.chopped_RO_light_off:
        for i in range(n_chop_cycles):
            at_mu(start + i * period_mu + FORT_on_mu)
            self.dds_FORT.sw.off()
            delay_mu(RO_pulse_length_mu)
            self.dds_FORT.sw.on()
            self.dds_cooling_DP.sw.on()

            at_mu(start + i * period_mu + gate_on_mu)
            self.ttl_SPCM_gate.off()  # unblocks the SPCM output ### remove: SPCM gate no longer exists
            at_mu(start + i * period_mu + RO_on_mu)
            self.dds_cooling_DP.sw.on()
            at_mu(start + i * period_mu + gate_on_mu + RO_pulse_length_mu)
            self.ttl_SPCM_gate.on()  # blocks the SPCM output ### remove: SPCM gate no longer exists
            at_mu(start + i * period_mu + RO_on_mu + RO_pulse_length_mu)
            self.dds_cooling_DP.sw.off()

            # cooling light doesn't seem synced up with SPCM gating based on photon count rate compared to RO_light_off case
            # with parallel:
            #     self.ttl_SPCM_gate.off() # unblocks the SPCM output ### remove: SPCM gate no longer exists
            #     self.dds_cooling_DP.sw.on()
            # delay_mu(RO_pulse_length_mu)
            # with parallel:
            #     self.dds_cooling_DP.sw.off()
            #     self.ttl_SPCM_gate.on()  # blocks the SPCM output ### remove: SPCM gate no longer exists
        # else:
        #     for i in range(n_chop_cycles):
        #         at_mu(start + i * period_mu + FORT_on_mu)
        #         self.dds_FORT.sw.off()
        #         delay_mu(RO_pulse_length_mu)
        #         self.dds_FORT.sw.on()
        #         at_mu(start + i * period_mu + gate_on_mu)
        #         self.ttl_SPCM_gate.off()  # unblocks the SPCM output ### remove: SPCM gate no longer exists
        #         delay_mu(RO_pulse_length_mu)
        #         self.ttl_SPCM_gate.on()  # blocks the SPCM output ### remove: SPCM gate no longer exists

@kernel
def record_chopped_blow_away(self):
    """

    :param self:
    :return:
    """

    # todo: change OP -> BA

    n_chop_cycles = int(self.t_blowaway/self.t_BA_chop_period + 0.5)
    # self.print_async("blowaway cycles:", n_chop_cycles)

    assert n_chop_cycles >= 1, "t_blowaway should be > t_BA_chop_period"

    BA_pulse = self.t_BA_chop_period * 0.35
    FORT_pulse = self.t_BA_chop_period - BA_pulse


    self.core.reset()

    with self.core_dma.record("chopped_blow_away"):

        start = now_mu()
        period_mu = self.core.seconds_to_mu(self.t_BA_chop_period)
        BA_pulse_length_mu = self.core.seconds_to_mu(BA_pulse)
        BA_on_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_pulse_length_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_on_mu = self.core.seconds_to_mu(BA_pulse)

        self.dds_FORT.sw.off()
        delay_mu(BA_pulse_length_mu)

        if not self.blowaway_light_off:
            self.dds_cooling_DP.sw.on()
        else:
            self.dds_cooling_DP.sw.off()

        for i in range(n_chop_cycles):
            at_mu(start+i*period_mu+FORT_on_mu)
            self.dds_FORT.sw.on()
            delay_mu(FORT_pulse_length_mu)
            self.dds_FORT.sw.off()
            at_mu(start+i*period_mu+BA_on_mu)
            # self.dds_cooling_DP.sw.on() # the cooling AOM seems incredibly slow so I'm just leaving it on the whole time
            delay_mu(BA_pulse_length_mu)
            # self.dds_cooling_DP.sw.off()
        self.dds_FORT.sw.on()

@kernel
def chopped_blow_away(self):

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    self.ttl_repump_switch.on()  # turns off the RP AOM

    delay(0.1 * ms)

    # set coils for blowaway
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_blowaway, self.AZ_top_volts_blowaway,
         self.AX_volts_blowaway, self.AY_volts_blowaway],
        channels=self.coil_channels)
    delay(0.3 * ms)

    with sequential:

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_blowaway, amplitude=self.ampl_cooling_DP_MOT)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        delay(0.1 * ms)

        if self.blowaway_light_off:
            self.dds_cooling_DP.sw.off()
            self.dds_AOM_A6.sw.off()
        else:
            # just turn the AOM up all the way. as long as we're 'saturating' the blowaway, it's okay if this doesn't
            # always give the same optical power
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-7.0))
            self.dds_AOM_A6.sw.on()
            self.dds_cooling_DP.sw.on()

    self.core_dma.playback_handle(ba_dma_handle)


    self.dds_cooling_DP.sw.off() ### Turns off cooling DP
    self.ttl_repump_switch.on() ### turns off MOT RP
    self.dds_AOM_A6.sw.off() ### turns off fiber AOM6

    delay(0.1 * ms)

    ### reset AOM RF powers
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
    delay(0.1*ms)

@kernel
def record_chopped_optical_pumping(self):
    """

    :param self:
    :return:
    """

    n_chop_cycles = int(self.t_pumping/self.t_OP_chop_period + 0.5)
    assert n_chop_cycles >= 1, "t_pumping should be > t_OP_chop_period"

    # todo: use duty cycle ExperimentVariables
    OP_pulse = self.t_OP_chop_period * self.duty_cycle_OP
    FORT_pulse = self.t_OP_chop_period * self.duty_cycle_FORT

    n_depump_chop_cycles = int(self.t_depumping /self.t_OP_chop_period)

    self.core.reset()

    with self.core_dma.record("chopped_optical_pumping"):

        # hardcoded for offsets for now, but this gives decent separation between the FORT and OP
        # pulses
        start = now_mu()
        period_mu = self.core.seconds_to_mu(self.t_OP_chop_period)

        OP_pulse_length_mu = self.core.seconds_to_mu(OP_pulse)
        FORT_pulse_length_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_on_mu = self.core.seconds_to_mu(0.0)
        OP_on_mu = self.core.seconds_to_mu(self.t_OP_chop_offset)

        if not (self.pumping_light_off or self.D1_off_in_OP_phase):
            for i in range(n_chop_cycles):
                at_mu(start+i*period_mu+FORT_on_mu)
                self.dds_FORT.sw.off()
                delay_mu(OP_pulse_length_mu)
                self.dds_FORT.sw.on()
                at_mu(start+i*period_mu+OP_on_mu)
                self.dds_D1_pumping_DP.sw.on()
                delay_mu(OP_pulse_length_mu)
                self.dds_D1_pumping_DP.sw.off()
        else:
            for i in range(n_chop_cycles):
                at_mu(start + i * period_mu + FORT_on_mu)
                self.dds_FORT.sw.off()
                delay_mu(OP_pulse_length_mu)
                self.dds_FORT.sw.on()

        if n_depump_chop_cycles > 0:

            # turn off the pumping repump
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            self.dds_pumping_repump.sw.off()
            start = now_mu() + self.core.seconds_to_mu(0.5 * us)

            if not (self.pumping_light_off or self.D1_off_in_depump_phase):
                for i in range(n_depump_chop_cycles):
                    at_mu(start + i * period_mu + FORT_on_mu)
                    self.dds_FORT.sw.off()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_FORT.sw.on()
                    at_mu(start + i * period_mu + OP_on_mu)
                    self.dds_D1_pumping_DP.sw.on()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_D1_pumping_DP.sw.off()
            else:
                for i in range(n_depump_chop_cycles):
                    at_mu(start + i * period_mu + FORT_on_mu)
                    self.dds_FORT.sw.off()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_FORT.sw.on()

@kernel
def chopped_optical_pumping(self):
    """
    Optical pumping with FORT chopped.

    Turns on:
        Fiber AOMs 5 and 6
        GRIN1and2 DDS
        GRIN1 AOM

    Turns off at the end:
        D1 DP
        Pumping Repumper
        Fiber AOMs 5 and 6
        GRIN1and2 DDS
        GRIN1 AOM
        cooling DP

    :param self:
    :return:
    """

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    if self.t_depumping + self.t_pumping > 3*ms:
        delay(2000 * us)  # we need extra slack

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light
    delay(1 * us)
    self.dds_AOM_A1.sw.off()
    self.dds_AOM_A2.sw.off()
    self.dds_AOM_A3.sw.off()
    self.dds_AOM_A4.sw.off()
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    delay(1 * us)

    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))
    delay(1 * us)

    ### Tunring on pumping RP:
    self.dds_pumping_repump.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(1*ms)

    ### so that D1 can pass
    self.GRIN1and2_dds.sw.on()
    self.ttl_GRIN1_switch.off()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.4 * ms)  # coil relaxation time

    self.core_dma.playback_handle(op_dma_handle)
    delay(self.t_depumping)
    self.dds_D1_pumping_DP.sw.off() ### turning off D1 DP
    self.dds_pumping_repump.sw.off() ### turning off pumping RP

    delay(2*us)
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    delay(100 * us)

    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
    delay(1*ms)
    self.GRIN1and2_dds.sw.off()
    self.ttl_GRIN1_switch.on()
    # self.dds_cooling_DP.sw.off()

@kernel
def optical_pumping(self):
    """
    optical pumping without chopping the FORT

    :param self:
    :return:
    """

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light

    ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
    # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(self.p_pumping_repump_A5))
    # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(self.p_pumping_repump_A6))
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(1*ms)

    ### so that D1 can pass
    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
    self.GRIN1and2_dds.sw.on()
    self.ttl_GRIN1_switch.off()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.4 * ms)  # coil relaxation time

    ### Optical pumping phase
    self.dds_D1_pumping_DP.set(frequency=self.f_D1_pumping_DP, amplitude=dB_to_V(self.p_D1_pumping_DP))
    delay (10 * us)
    self.dds_D1_pumping_DP.sw.on()
    self.dds_pumping_repump.sw.on()
    delay(self.t_pumping)
    self.dds_D1_pumping_DP.sw.off()
    delay(self.t_depumping)
    self.dds_pumping_repump.sw.off()


    delay(2*us)

    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

    delay(100 * us)

    # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
    delay(1*ms)
    self.GRIN1and2_dds.sw.off()
    self.ttl_GRIN1_switch.on()

@kernel
def measure_FORT_MM_fiber(self):
    # ALICE & BOB: both use Sampler1 - ch7
    # BOB: IF FORT feedback use APD, make sure to change MM smapler ch & APD sampler ch in BaseExperiment.py

    measurement_buf = np.array([0.0]*8)
    measurement = 0.0
    avgs = 50
    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        measurement += measurement_buf[self.FORT_MM_sampler_ch]
        delay(0.1*ms)
    measurement /= avgs
    self.append_to_dataset("FORT_MM_science_volts", measurement)

@kernel
def measure_GRIN1(self):

    """
    used for monitring GRIN1
    GRIN1_sampler_ch = 4 defined at "BaseExperiment.py"

    can be used for monitoring exitation power

    """
    measurement_buf = np.array([0.0]*8)
    measurement = 0.0
    avgs = 10

    # self.dds_FORT.sw.off() ### Why turnning off FORT?
    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_cooling_DP.sw.off()

    # todo: D1 feedback
    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))

    # UNBLOCKING GRIN1
    self.ttl_GRIN1_switch.off()
    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(0.0))
    self.GRIN1and2_dds.sw.on()
    # self.ttl_excitation_switch.off()

    #turning D1 ON
    self.dds_D1_pumping_DP.sw.on()

    delay(0.1*ms)

    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        measurement += measurement_buf[self.GRIN1_sampler_ch]
        delay(0.1*ms)
    measurement /= avgs
    self.append_to_dataset("GRIN1_D1_monitor", measurement)

    # turning D1 off & repump on
    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
    self.dds_D1_pumping_DP.sw.off()
    # self.ttl_repump_switch.off()  # turns the RP AOM on
    self.ttl_exc0_switch.off() # turns EXC0 AOM on

    delay(0.1*ms)
    measurement = 0.0

    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        measurement += measurement_buf[self.GRIN1_sampler_ch]
        delay(0.1*ms)
    measurement /= avgs
    self.append_to_dataset("GRIN1_EXC_monitor", measurement)


    # turning EXC OFF
    self.GRIN1and2_dds.sw.off()
    # self.ttl_excitation_switch.on()
    self.ttl_GRIN1_switch.on()
    self.ttl_exc0_switch.on()

    delay(0.1 * ms)

@kernel
def measure_REPUMP(self):
    """
    used for monitring PUMPING REPUMP power

    This is in end_measurement

    AOM1: Sampler0, 7
    AOM2: Sampler0, 5

    """
    measurement_buf = np.array([0.0]*8)
    measurement1 = 0.0 # Repump 1
    measurement2 = 0.0 # Repump 2

    avgs = 10

    # self.dds_FORT.sw.off() ### Why turnning of FORT?
    self.ttl_repump_switch.off()  # turns the RP AOM on
    self.dds_pumping_repump.sw.off()
    self.dds_cooling_DP.sw.off()

    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()

    delay(0.1 * ms)

    # Repump 1
    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[7] # Repump 1
        delay(0.1 * ms)
        measurement2 += measurement_buf[5] # Repump 2


    measurement1 /= avgs
    measurement2 /= avgs

    self.append_to_dataset("REPUMP1_monitor", measurement1)
    self.append_to_dataset("REPUMP2_monitor", measurement2)

    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_AOM_A1.sw.off()
    self.dds_AOM_A2.sw.off()

    delay(0.1 * ms)

@kernel
def measure_PUMPING_REPUMP(self):
    """
    used for monitring REPUMP power
    PUMPING_REPUMP1_monitor, PUMPING_REPUMP2_monitor defined

    This is in end_measurement

    Pumping Repumper is sent to AOM5 & 6

    AOM5: Sampler0, 1
    AOM6: Sampler0, 2

    """
    measurement_buf = np.array([0.0]*8)
    measurement1 = 0.0 # Repump 5
    measurement2 = 0.0 # Repump 6

    avgs = 10

    # self.dds_FORT.sw.off() ### Why tunning of FORT?
    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_cooling_DP.sw.off()

    self.dds_pumping_repump.sw.on() # turns of PR AOM

    # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-8.0))
    # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-8.0))

    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(0.1 * ms)

    # Repump 1
    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[1] # Repump 5
        # delay(0.1 * ms)
        measurement2 += measurement_buf[2] # Repump 6


    measurement1 /= avgs
    measurement2 /= avgs

    self.append_to_dataset("PUMPING_REPUMP1_monitor", measurement1)
    self.append_to_dataset("PUMPING_REPUMP2_monitor", measurement2)

    self.dds_pumping_repump.sw.off()  # turns the PR AOM off

    # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)


    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

    delay(0.1 * ms)

@kernel
def zotino_stability_test(self):
    '''
    Zotino Stability test function to see if Zotino voltage output drifts.

    this is called in load_MOT_and_FORT

    * for testing purposes,
    1) AZ_bottom_volts_MOT
        zotino_test_1_Zotino_channel:
        zotino_test_1_Sampler_channel:

    # Defined in BaseExperiment.py
    # zotino_test_1_Zotino_channel = 6  # Zotino 0 - ch6
    # zotino_test_2_Zotino_channel = 7  # Zotino 0 - ch7

    # Defined in BaseExperiment.py
    # self.experiment.set_dataset("zotino_test1_monitor", [0.0], broadcast=True)
    # self.experiment.set_dataset("zotino_test2_monitor", [0.0], broadcast=True)

    '''

    zotino_test1_Sampler_channel = 5 # Sampler 1 - ch5
    zotino_test2_Sampler_channel = 6  # Sampler 1 - ch6

    measurement_buf = np.array([0.0]*8)
    measurement1 = 0.0 #  1
    measurement2 = 0.0 #  2

    avgs = 50

    # Repump 1
    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[zotino_test1_Sampler_channel]
        delay(0.1 * ms)
        measurement2 += measurement_buf[zotino_test2_Sampler_channel]


    measurement1 /= avgs
    measurement2 /= avgs

    self.append_to_dataset("zotino_test1_monitor", measurement1)
    self.append_to_dataset("zotino_test2_monitor", measurement2)

    delay(0.1 * ms)

@kernel
def Sampler0_test(self):
    '''
    I have conencted 5V signal from the TENMA DC power supply to Sample0, ch0. Sampler0 is used to feedback to the AOMs.
    I am going to see if the measurement of this channel fluctuates especially when atom loading is bad, or not.

    '''

    ch = 0 # Sampler 0 - ch0
    measurement_buf = np.array([0.0]*8)
    measurement1 = 0.0
    avgs = 10

    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[ch]

    measurement1 /= avgs

    self.append_to_dataset("Sampler0_test", measurement1)

    delay(0.1 * ms)

@kernel
def measure_coil_driver(self):
    '''
    I have connected "monitor out" of the coil drivers to Sampler2 ch4, 5, 6, and 7, to monitor the output
    of the coils while loading MOT and atom to see if there is any correlation between bad atom loading
    and coil driver outputs.

    '''

    avgs = 10

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)
    delay(100 * ms)

    measurement_buf = np.array([0.0] * 8)
    coilZ_bottom = 0.0
    coilZ_top = 0.0
    coilX = 0.0
    coilY = 0.0

    for i in range(avgs):
        self.sampler2.sample(measurement_buf)
        coilZ_bottom += measurement_buf[4]
        coilZ_top += measurement_buf[7]
        coilX += measurement_buf[6]
        coilY += measurement_buf[5]

        delay(0.1 * ms)
    coilZ_bottom /= avgs
    coilZ_top /= avgs
    coilX /= avgs
    coilY /= avgs

    self.append_to_dataset("coil_driver_AZ_bottom_MOT", coilZ_bottom)
    self.append_to_dataset("coil_driver_AZ_top_MOT", coilZ_top)
    self.append_to_dataset("coil_driver_AX_MOT", coilX)
    self.append_to_dataset("coil_driver_AY_MOT", coilY)

    delay(0.1 * ms)

@kernel
def measure_MOT_end(self):
    """
    used for monitring MOT power
    MOT_end_monitor1 defined

    This is in end_measurement

    AOM1: Sampler0, 7
    AOM2: Sampler0, 5
    AOM3: Sampler0, 3
    AOM4: Sampler0, 4
    AOM5: Sampler0, 1
    AOM6: Sampler0, 2

    """
    ao_s1 = 7
    ao_s2 = 5
    ao_s3 = 3
    ao_s4 = 4
    ao_s5 = 1
    ao_s6 = 2

    avgs = 50


    self.dds_FORT.sw.off()
    self.ttl_repump_switch.on()  # turns the RP AOM on
    self.dds_cooling_DP.sw.on()

    delay(0.1 * ms)

    ### MOT1 & MOT2 & MOT5
    measurement_buf = np.array([0.0] * 8)
    measurement1 = 0.0  # MOT1
    measurement2 = 0.0  # MOT2
    measurement3 = 0.0  # MOT5

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A5.sw.on()

    delay(0.1 * ms)

    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[ao_s1]  # MOT1
        delay(0.1 * ms)
        measurement2 += measurement_buf[ao_s2]  # MOT2
        measurement3 += measurement_buf[ao_s5]  # MOT5

    measurement1 /= avgs
    measurement2 /= avgs
    measurement3 /= avgs

    self.append_to_dataset("MOT1_end_monitor", measurement1)
    self.append_to_dataset("MOT2_end_monitor", measurement2)
    self.append_to_dataset("MOT5_end_monitor", measurement3)

    delay(0.1 * ms)

    # # Repump
    #
    # self.ttl_repump_switch.off()  # turns the RP AOM on
    # self.dds_cooling_DP.sw.off()
    #
    # delay(0.1 * ms)
    #
    # measurement_buf = np.array([0.0] * 8)
    # measurement1 = 0.0  # MOT1
    # measurement2 = 0.0  # MOT2
    # measurement3 = 0.0  # MOT5
    #
    # for i in range(avgs):
    #     self.sampler0.sample(measurement_buf)
    #     delay(0.1 * ms)
    #     measurement1 += measurement_buf[ao_s1]  # MOT1
    #     delay(0.1 * ms)
    #     measurement2 += measurement_buf[ao_s2]  # MOT2
    #     measurement3 += measurement_buf[ao_s5]  # MOT5
    #
    # measurement1 /= avgs
    # measurement2 /= avgs
    # measurement3 /= avgs
    #
    #
    # self.append_to_dataset("REPUMP1_monitor", measurement1)
    # self.append_to_dataset("REPUMP2_monitor", measurement2)
    # self.append_to_dataset("REPUMP5_monitor", measurement3)
    #
    # self.dds_AOM_A1.sw.off()
    # self.dds_AOM_A2.sw.off()
    # self.dds_AOM_A5.sw.off()
    #
    # delay(0.1 * ms)


    ### MOT3 & MOT4 & MOT6

    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_cooling_DP.sw.on()

    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(0.1 * ms)

    measurement_buf = np.array([0.0] * 8)
    measurement1 = 0.0
    measurement2 = 0.0
    measurement3 = 0.0

    delay(0.1 * ms)

    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[ao_s3]
        delay(0.1 * ms)
        measurement2 += measurement_buf[ao_s4]
        measurement3 += measurement_buf[ao_s6]

    measurement1 /= avgs
    measurement2 /= avgs
    measurement3 /= avgs

    self.append_to_dataset("MOT3_end_monitor", measurement1)
    self.append_to_dataset("MOT4_end_monitor", measurement2)
    self.append_to_dataset("MOT6_end_monitor", measurement3)


    # # repump
    #
    # self.ttl_repump_switch.off()  # turns the RP AOM on
    # self.dds_cooling_DP.sw.off()
    #
    # measurement_buf = np.array([0.0] * 8)
    # measurement1 = 0.0
    # measurement2 = 0.0
    # measurement3 = 0.0
    #
    # delay(0.1 * ms)
    #
    # for i in range(avgs):
    #     self.sampler0.sample(measurement_buf)
    #     delay(0.1 * ms)
    #     measurement1 += measurement_buf[ao_s3]
    #     delay(0.1 * ms)
    #     measurement2 += measurement_buf[ao_s4]
    #     measurement3 += measurement_buf[ao_s6]
    #
    # measurement1 /= avgs
    # measurement2 /= avgs
    # measurement3 /= avgs
    #
    # self.append_to_dataset("REPUMP3_monitor", measurement1)
    # self.append_to_dataset("REPUMP4_monitor", measurement2)
    # self.append_to_dataset("REPUMP6_monitor", measurement3)
    #
    #
    # self.dds_AOM_A3.sw.off()
    # self.dds_AOM_A4.sw.off()
    # self.dds_AOM_A6.sw.off()
    #
    # self.ttl_repump_switch.on()  # turns the RP AOM off


    delay(0.1 * ms)

@kernel
def measure_Magnetometer(self):
    ### x,y, and z axes are connected to Sampler2 Ch1,2, and 3, respectively.
    avgs = 10

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)
    delay(1*ms)

    measurement_buf = np.array([0.0]*8)
    MagnetometerX = 0.0
    MagnetometerY = 0.0
    MagnetometerZ = 0.0

    for i in range(avgs):
        self.sampler2.sample(measurement_buf)
        MagnetometerX += measurement_buf[self.Magnetometer_X_ch]
        MagnetometerY += measurement_buf[self.Magnetometer_Y_ch]
        MagnetometerZ += measurement_buf[self.Magnetometer_Z_ch]
        delay(0.1*ms)
    MagnetometerX /= avgs
    MagnetometerY /= avgs
    MagnetometerZ /= avgs
    self.append_to_dataset("Magnetometer_MOT_X", MagnetometerY * 350) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_MOT_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_MOT_Z", MagnetometerZ * 350)

    #####################################  Measure in the OP phase
    ### Set the coils to OP setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(1 * ms)

    measurement_buf = np.array([0.0] * 8)
    MagnetometerX = 0.0
    MagnetometerY = 0.0
    MagnetometerZ = 0.0

    for i in range(avgs):
        self.sampler2.sample(measurement_buf)
        MagnetometerX += measurement_buf[self.Magnetometer_X_ch]
        MagnetometerY += measurement_buf[self.Magnetometer_Y_ch]
        MagnetometerZ += measurement_buf[self.Magnetometer_Z_ch]
        delay(0.1 * ms)
    MagnetometerX /= avgs
    MagnetometerY /= avgs
    MagnetometerZ /= avgs
    self.append_to_dataset("Magnetometer_OP_X", MagnetometerY * 350) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_OP_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_OP_Z", MagnetometerZ * 350)

    #####################################  Measure with Zotino set to zero V
    ### Turn off all the coils
    self.zotino0.set_dac(
        [0.0, 0.0, 0.0, 0.0],
        channels=self.coil_channels)
    delay(100 * ms)

    measurement_buf = np.array([0.0] * 8)
    MagnetometerX = 0.0
    MagnetometerY = 0.0
    MagnetometerZ = 0.0

    for i in range(avgs):
        self.sampler2.sample(measurement_buf)
        MagnetometerX += measurement_buf[self.Magnetometer_X_ch]
        MagnetometerY += measurement_buf[self.Magnetometer_Y_ch]
        MagnetometerZ += measurement_buf[self.Magnetometer_Z_ch]
        delay(0.1 * ms)

    MagnetometerX /= avgs
    MagnetometerY /= avgs
    MagnetometerZ /= avgs
    self.append_to_dataset("Magnetometer_Zero_X", MagnetometerY * 350) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_Zero_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_Zero_Z", MagnetometerZ * 350)
    delay(0.1*ms)

@kernel
def end_measurement(self):
    """
    End the measurement by setting datasets and deciding whether to increment the measuement index
    :param self:
    :return measurement: TInt32, the measurement index
    """

    ### update the datasets
    self.set_dataset(self.measurements_progress, 100 * self.measurement / self.n_measurements, broadcast=True)

    self.append_to_dataset('SPCM0_RO1_current_iteration', self.SPCM0_RO1)
    self.append_to_dataset('SPCM1_RO1_current_iteration', self.SPCM1_RO1)
    delay(1 * ms)
    self.append_to_dataset('SPCM0_RO2_current_iteration', self.SPCM0_RO2)
    self.append_to_dataset('SPCM1_RO2_current_iteration', self.SPCM1_RO2)
    delay(1 * ms)
    self.append_to_dataset('BothSPCMs_RO1_current_iteration', self.BothSPCMs_RO1)
    self.append_to_dataset('BothSPCMs_RO2_current_iteration', self.BothSPCMs_RO2)
    delay(1 * ms)
    self.SPCM0_RO1_list[self.measurement] = self.SPCM0_RO1
    self.SPCM1_RO1_list[self.measurement] = self.SPCM1_RO1
    self.SPCM0_RO2_list[self.measurement] = self.SPCM0_RO2
    self.SPCM1_RO2_list[self.measurement] = self.SPCM1_RO2

    self.BothSPCMs_RO1_list[self.measurement] = self.BothSPCMs_RO1
    self.BothSPCMs_RO2_list[self.measurement] = self.BothSPCMs_RO2
    self.atom_loading_time_list[self.measurement] = self.atom_loading_time

    self.append_to_dataset("SPCM0_FORT_science", self.SPCM0_FORT_science)

    ### now done at the beginning of the experiment for FORT POL stabilization
    # delay(1*ms)
    # measure_FORT_MM_fiber(self)

    delay(1*ms)
    measure_GRIN1(self)
    delay(1*ms)
    measure_PUMPING_REPUMP(self)
    delay(1*ms)

    # if self.which_node == "alice":
    #     measure_Magnetometer(self)
    #     delay(1*ms)
    #     Sampler0_test(self)
    #     delay(1*ms)
    #     measure_coil_driver(self)

    # measure_MOT_end(self)
    # delay(1*ms)
    measure_REPUMP(self)
    delay(1*ms)

    advance = 1
    if self.__class__.__name__ != 'ExperimentCycler':
        if self.require_atom_loading_to_advance:
            if not self.SPCM0_RO1/self.t_SPCM_first_shot > self.single_atom_threshold:
                advance = 0

        if self.require_D1_lock_to_advance:
            t1_D1_checked = now_mu()
            while True:
                self.ttl_D1_lock_monitor.sample_input()
                delay(0.1 * ms)
                laser_locked = int(1 - self.ttl_D1_lock_monitor.sample_get()) ## this is 1 when the laser is locked, it is 0 otherwise.

                if laser_locked:
                    advance = 1
                    self.set_dataset("time_without_D1", 0.0, broadcast=True)  ### resetting time_without_D1 when the laser is locked
                    delay(10*us)
                    break
                else:
                    t2_D1_checked = now_mu()
                    time_without_D1 = self.core.mu_to_seconds(t2_D1_checked - t1_D1_checked)
                    self.set_dataset("time_without_D1", time_without_D1, broadcast=True)
                    delay(2 * s)

    if advance:
        self.measurement += 1
        delay(1 * ms)
        self.append_to_dataset('SPCM0_RO1', self.SPCM0_RO1)
        self.append_to_dataset('SPCM1_RO1', self.SPCM1_RO1)
        self.append_to_dataset('SPCM0_RO2', self.SPCM0_RO2)
        self.append_to_dataset('SPCM1_RO2', self.SPCM1_RO2)
        self.append_to_dataset('BothSPCMs_RO1', self.BothSPCMs_RO1)
        self.append_to_dataset('BothSPCMs_RO2', self.BothSPCMs_RO2)
        self.append_to_dataset('SPCM0_test_RO', self.SPCM0_test_RO)
        delay(1 * ms)

@rpc(flags={"async"})
def set_RigolDG1022Z(frequency: TFloat, vpp: TFloat, vdc: TFloat):
    """
    Set the frequency, V_DC, and V_pp of a connected Rigol DG1000 series function generator.

    See https://www.batronix.com/pdf/Rigol/ProgrammingGuide/DG1000Z_ProgrammingGuide_EN.pdf

    :return TFloat: the frequency of the Rigol DG1022Z in Hz
    """
    # NAME = 'USB0::0x1AB1::0x0642::DG1ZA252402452::INSTR'

    rm = visa.ResourceManager()
    instruments = rm.list_resources()
    for instr in instruments:
        if 'DG' in instr:
            NAME = instr
            break

    funcgen = rm.open_resource(NAME, timeout=20, chunk_size=1024000)
    funcgen.timeout = 1000  # in ms I think

    # Make sure mod is off
    funcgen.write(r":SOUR1:MOD OFF")

    # Set the frequency
    freq = int(frequency)
    funcgen.write(f":SOUR1:FREQ {freq}")
    actual_freq = float(funcgen.query("SOUR1:FREQ?"))
    assert actual_freq == freq, "Oops! The device frequency is not set to the requested value!"
    print(f"f_carrier: {actual_freq} Hz")

    # Set the Vpp
    funcgen.write(f":SOUR1:VOLT {vpp}")
    actual_vpp = float(funcgen.query("SOUR1:VOLT?"))
    assert actual_vpp == vpp, "Oops! The device V_pp is not set to the requested value!"
    print(f"Vpp: {actual_vpp} V")

    # Set the dc offset
    funcgen.write(f":SOUR1:VOLT:OFFS {vdc}")
    actual_vdc = float(funcgen.query("SOUR1:VOLT:OFFS?"))
    assert actual_vdc == vdc, "Oops! The device V_DC is not set to the requested value!"
    print(f"Vdc: {actual_vdc} V")

###############################################################################
# 2. EXPERIMENT FUNCTIONS
# These are the experiments we run, and the name of each should end with
# experiment in order to have it show up in GeneralVariableScans
###############################################################################

@kernel
def MOT_loading_experiment(self):
    self.core.reset()

    for measurement in range(self.n_measurements):
        self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs
        load_MOT_and_FORT(self)

@kernel
def atom_loading_experiment(self):
    """
    The most basic two-readout single atom experiment.

    Load a MOT, load a single atom, readout, wait self.t_delay_between_shots, readout again.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    self.require_D1_lock_to_advance = False # override experiment variable

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)
        delay(0.1*ms)

        first_shot(self)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        # self.dds_cooling_DP.sw.off()
        delay(self.t_delay_between_shots)
        # self.dds_cooling_DP.sw.on()

        # rtio_log("2nd_shot_block",1)
        second_shot(self)
        # rtio_log("2nd_shot_block", 0)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)

    # self.dds_FORT.sw.off()

@kernel
def atom_loading_2_experiment(self):
    """
    Simple atom loading experiment based on load_MOT_and_FORT_until_atom:
    turn on MOT and FORT until an atom is detected. Then turn off MOT.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.ttl7.pulse(100 * us)  ### for triggering oscilloscope
        # delay(0.1 * ms)

        # load_MOT_and_FORT(self)
        load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

        delay(1*ms)
        first_shot(self)
        delay(1 * ms)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)
        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

    self.dds_FORT.sw.off()

@kernel
def atom_loading_3_experiment(self):
    """
    Same as atom_loading_2_experiment but with smooth FORT ramp to science setpoint.
    No errors. But two issues:
    1- changing the settings on cooling dds affects the FORT output! It starts the RAM on FORT!
    2- I don't get atoms, the background is lower that expected. So, something is off!

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    delay(1 * ms)

    ####################### preparing and writing the RAM
    ### Exit RAM Mode
    self.dds_FORT.set_cfr1(ram_enable=0)
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    delay(0.1 * ms)
    prepare_FORT_RAM_profile(self)
    delay(1 * ms)

    ### Configures the RAM playback engine
    self.dds_FORT.set_profile_ram(
        start=0,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP,
    )

    ### write the data onto RAM
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    self.dds_FORT.write_ram(self.FORT_amplitudes_list)
    ####################### End of preparing and writing the RAM


    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.ttl7.pulse(100 * us)  ### for triggering oscilloscope
        # delay(0.1 * ms)

        # load_atom_smooth_FORT(self)
        load_until_atom_smooth_FORT(self)

        delay(1*ms)
        first_shot_smooth_FORT(self)
        # first_shot(self)
        delay(1 * ms)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)
        second_shot_smooth_FORT(self)
        # second_shot(self)

        end_measurement(self)

    ### Exit RAM Mode
    self.dds_FORT.set_cfr1(ram_enable=0)
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    delay(1 * ms)

    self.dds_FORT.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def atom_loading_optimizer_experiment(self):
    """
    To be used in GVO to optimize atom loading time.
    Simple atom loading experiment to optimize atom loading time.
    Use this in GVO with atom_loading_time_cost. This measures the atom loading time based on RO1 rather in FORT loading.
    Does not try to load an atom forever. It stops after a few trials.
    Otherwise, we cannot use this function to optimize atom loading time; with wrong shims settings, for example,
    it would keep trying and never finishes the iteration.

    turn on MOT and FORT until an atom is detected. Then turn off MOT to dissipate, and do RO1 and RO2.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(1 * ms)
        self.dds_FORT.sw.off()

        ###################################################### atom loading attempt
        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(0.1 * ms)
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        delay(1*ms)


        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()

        delay(1 * ms)
        self.ttl_UV.on()

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        atom_loaded = False
        try_n = 0
        t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()

        while not atom_loaded and try_n < max_tries:
            delay(100 * us)  ### Needs a delay of about 100us or maybe less
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_counter.gate_rising(atom_check_time)

            SPCM0_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_atom_check = self.ttl_SPCM1_counter.fetch_count()

            BothSPCMs_atom_check = int((SPCM0_atom_check + SPCM1_atom_check) / 2)

            try_n += 1

            if BothSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
            self.append_to_dataset("BothSPCMs_atom_check_in_loading", BothSPCMs_atom_check)

        self.ttl_UV.off()

        if not atom_loaded:
            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms)  ### necessary to avoid underflow
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                self.laser_stabilizer.run()
                self.dds_FORT.sw.on()
                delay(0.1 * ms)

        ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
        delay(1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(0.4 * ms)

        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling

        delay(1 * ms)
        delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate


        # ###########  PGC on the trapped atom  #############
        # if self.do_PGC_after_loading:
        #     if self.which_node == 'alice':
        #         self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        #     elif self.which_node == 'bob':
        #         self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)
        #     ### set the cooling DP AOM to the PGC settings
        #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        #     self.ttl_repump_switch.off()  ### turn on MOT RP
        #     self.dds_cooling_DP.sw.on()  ### turn on cooling
        #     delay(10 * us)
        #     # self.dds_AOM_A5.sw.off()
        #     # self.dds_AOM_A6.sw.off()
        #     delay(self.t_PGC_after_loading)  ### this is the PGC time
        # ###################################################

        ### saving the atom loading time for each loaded atom.
        self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
        delay(1 * ms)
        ###################################################### end of atom loading attempt


        first_shot(self)
        delay(1 * ms)

        if self.BothSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(100 * us)  ### Needs a delay of about 100us or maybe less
            atom_loaded = True
            t_after_atom = now_mu()
            self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
        else:
            self.atom_loading_time = 10e9
            ### Just a large number to show no atom loading. This is compared to the typical 0.5 second atom loading.

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(5*ms)
        second_shot(self)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def atom_loading_for_optimization_experiment(self):
    """
    Simple atom loading experiment based on load_MOT_and_FORT_until_atom. We use this function with additional stages
    compared to atom_loading_2_experiment to scan and optimize different parameters like 852 waveplates.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable


    ######### to scan 852 waveplates with t_FORT_drop = 10us for example and find max retention (low T)
    # delay(1*ms)
    # move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
    # move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)
    #
    # delay(5*ms)
    # self.core.reset()
    #
    # position_852_HWP = get_rotator_deg(self, name="852_HWP")
    # position_852_QWP = get_rotator_deg(self, name="852_QWP")
    #
    # delay(5 * ms)
    # self.core.reset()
    # self.print_async("position_852_HWP: ", position_852_HWP)
    # self.print_async("position_852_QWP:", position_852_QWP)
    ###################################################################################################


    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.ttl7.pulse(100 * us)  ### for triggering oscilloscope
        # delay(0.1 * ms)

        # load_MOT_and_FORT(self)
        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)

        delay(1*ms)
        first_shot(self)
        delay(1 * ms)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()




        first_shot_chopped(self)
        delay(1*ms)





        ################## to see if RO heats atoms
        # for i in range(5):
        #     shot_without_measurement(self)
        ###########################################

        # delay(self.t_delay_between_shots)
        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

    self.dds_FORT.sw.off()

@kernel
def atom_loading_for_PGC_optimization_experiment(self):
    """
    Simple atom loading experiment based on load_MOT_and_FORT_until_atom and an extra RO or PGC between the shots
    to optimize PGC settings without affeting loading.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable


    if self.enable_laser_feedback:
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.ttl7.pulse(100 * us)  ### for triggering oscilloscope
        # delay(0.1 * ms)

        # load_MOT_and_FORT(self)
        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)

        delay(1*ms)
        first_shot(self)
        delay(1 * ms)

        ###########  PGC on the trapped atom  #############
        ### Set the coils to PGC setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(0.4 * ms)
        ### set the cooling DP AOM to the PGC settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(10 * us)
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP
        delay(20 * ms)  ### this is the PGC time
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        ###################################################

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()


        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

    self.dds_FORT.sw.off()

@kernel
def trap_frequency_experiment(self):
    """
    For spectroscopy of the trap vibrational frequencies with the Rigol D11022Z function generator.

    One way to use this is with GeneralVariableScan and scan_variable1 = f_Rigol_modulation,
    where f_Rigol_modulation must be within the boundaries Rigol_carrier_frequency +/- Rigol_FM_deviation
    specified in ExpermientVariables (which in turn must be set on the Rigol D11022Z).

    :param self: an experiment instance.
    :return:
    """

    set_RigolDG1022Z(frequency=self.Rigol_carrier_frequency,
                     vpp=self.Rigol_V_pp,
                     vdc=self.Rigol_V_DC)

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0

    self.require_D1_lock_to_advance = False # override experiment variable

    # self.set_dataset(self.SPCM0_rate_dataset,
    #                  [0.0],
    #                  broadcast=True)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        #TODO: just set the rigol frequency using pyvisa

        if self.enable_laser_feedback:
            self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        # set the FORT AOM to the science setting. this is only valid if we have run
        # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT
        self.dds_FORT.set(frequency=self.f_FORT,
                                amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        if not self.no_first_shot:
            first_shot(self)

        ##############################################################################
        # modulate the FORT
        ##############################################################################

        self.FORT_mod_switch.on()  # toggle the modulation to the VCA
        delay(10 * ms)
        self.FORT_mod_switch.off()
        delay(1 * ms)

        second_shot(self)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def microwave_Rabi_experiment(self):
    """
    Microwave and optical pumping experiment.
    This experiment is used for any experiment which involves optical pumping and a microwave rotations with a constant
    microwave amplitude over a specified (possibly zero) duration. For example, this can be used for
    - microwave Rabi oscillations to test the optical pumping fidelity
    - microwave spectroscopy (useful for zeroing the magnetic field by finding the resonances for different ground state
    transitions |F=1,m>->|F=2,m'>)
    - depumping measurements (by using a non-zero depump time)

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            self.laser_stabilizer.run()
            # bug -- microwave dds is off after AOM feedback; not clear why yet. for now, just turn it back on
            self.dds_microwaves.sw.on()

        load_MOT_and_FORT(self)

        delay(0.1 * ms)

        first_shot(self)

        # todo: do we need to pump into F=2? Remove see what happens.
        # self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        # delay(1 * ms) # leave the repump on so atoms are left in F=2
        # self.ttl_repump_switch.on()  # turns the MOT RP AOM off
        delay (1 * ms)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        ### With chopped pumping:
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)
            delay(1*ms)

        # ### with cw pumping:
        # if self.t_pumping > 0.0:
        #     delay (10 * us)
        #     optical_pumping(self)
        #     delay(1*ms)

        ############################
        # microwave phase
        ############################

        if self.t_microwave_pulse > 0.0:
            # self.ttl_repump_switch.on()  # turns off the RP AOM
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
            delay(5 * us)

            self.zotino0.set_dac(
                [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
                 self.AX_volts_microwave, self.AY_volts_microwave],
                channels=self.coil_channels)
            delay(0.3*ms)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(100 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            chopped_blow_away(self)


        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.dds_microwaves.sw.off()

@kernel
def microwave_Rabi_2_experiment(self):
    """
    Microwave and optical pumping experiment based on load_MOT_and_FORT_until_atom(self).

    This experiment is used for any experiment which involves optical pumping and a microwave rotations with a constant
    microwave amplitude over a specified (possibly zero) duration. For example, this can be used for
    - microwave Rabi oscillations to test the optical pumping fidelity
    - microwave spectroscopy (useful for zeroing the magnetic field by finding the resonances for different ground state
    transitions |F=1,m>->|F=2,m'>)
    - depumping measurements (by using a non-zero depump time)

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)

    # self.ttl7.on()
    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    # self.ttl7.off()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)

        delay(1 * ms)

        first_shot(self)

        # todo: do we need to pump into F=2? Remove see what happens.
        # self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        # delay(1 * ms) # leave the repump on so atoms are left in F=2
        # self.ttl_repump_switch.on()  # turns the MOT RP AOM off
        delay (1 * ms)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        ### With chopped pumping:
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)
            delay(1*ms)

        # ### with cw pumping:
        # if self.t_pumping > 0.0:
        #     delay (10 * us)
        #     optical_pumping(self)
        #     delay(1*ms)

        ############################
        # microwave phase
        ############################

        if self.t_microwave_pulse > 0.0:
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
            delay(2 * us)

            # ### Changing the bias field for testing.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(0.5*ms)

            # self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            chopped_blow_away(self)


        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.dds_microwaves.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def microwave_Ramsey_00_experiment(self):
    """
    Ramsey experiment with two pi/2 MW pulses with a variable time delay in the middle.
    This experiment can be used to measure the T2* of the clock qubit.

    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)


    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)

        delay(1 * ms)

        first_shot(self)

        delay (1 * ms)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        ### With chopped pumping:
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)
            delay(1*ms)

        ############################
        # microwave phase
        ############################
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        ### first pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse)
        self.ttl_microwave_switch.on()

        ### use only one of the following:

        ### without spin echo:
        delay(self.t_delay_between_shots)

        # ### with spin echo. This does not work well. Needs debugging.
        # delay(self.t_delay_between_shots/2)
        # self.ttl_microwave_switch.off()
        # delay(self.t_microwave_00_pulse)
        # self.ttl_microwave_switch.on()
        # delay(self.t_delay_between_shots/2)


        ### second pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse)
        self.ttl_microwave_switch.on()

        delay(2*us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)


        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.dds_microwaves.sw.off()

@kernel
def microwave_Ramsey_11_experiment(self):
    """
    Ramsey experiment with two pi/2 MW pulses with a variable time delay in the middle.

    This experiment can be used to measure the T2* of the mF=+1 and mF'=+1 qubit.

    1- atom loading
    2- OP
    3- MW transfer from mF=0 to mF'=1
    4- two MW pi/2 pulses between mF=+1 and mF'=+1 with a variable delay in between.
    5- Blowaway F=2 states and RO

    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)


    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)

        delay(1 * ms)

        first_shot(self)
        delay (1 * ms)

        ############################ optical pumping phase - pumps atoms into F=1,m_F=0
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)
            delay(1*ms)

        ############################ microwave 1: transfer mF=0 to mF'=1
        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8*self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        self.ttl_microwave_switch.off()
        delay(self.t_microwave_01_pulse)
        self.ttl_microwave_switch.on()
        delay(2 * us)

        ############################ microwave 2: Ramsey between mF=1 and mF'=1
        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        ### first pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse)
        self.ttl_microwave_switch.on()

        ### use only one of the following:

        ### without spin echo:
        delay(self.t_delay_between_shots)

        # ### with spin echo. This does not work well. Needs debugging.
        # delay(self.t_delay_between_shots/2)
        # self.ttl_microwave_switch.off()
        # delay(self.t_microwave_00_pulse)
        # self.ttl_microwave_switch.on()
        # delay(self.t_delay_between_shots/2)


        ### second pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse)
        self.ttl_microwave_switch.on()

        delay(2*us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.dds_microwaves.sw.off()

@kernel
def microwave_freq_scan_experiment(self):
    """
    This experiment scans the microwave frequency to find the transitions from F=1 to F=2 manifold btw different Zeeman levels.
    This is the first step for atom mapping. The code works, but the results are not very good. It is better to
    use microwave_freq_scan_with_photons_experiment to find these transitions.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- Uses MOT RP to pump atoms into F=2 manifold
    3- Uses microwave pi pulse to transfer population, say from (F=2, mF=0), to F=1
    4- Blow away F=2 atoms
    5- Measure retention

    * If microwave is off resonance, the retention is zero
    * If microwave is on resonance, the retention is ~1/5 = 20%.

    * another method for pumping added later to use pumping repumper instead of MOT RP
    * the base retention in this case would be 100% off-resonance

    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)


    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)
        delay(10 * us)

        first_shot(self)

        ### Turn on the bias field for OP
        self.zotino0.set_dac([self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP,
                              self.AX_volts_OP, self.AY_volts_OP], channels=self.coil_channels)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(0.1 * ms)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.
        delay(0.5 * ms) ### coils relaxation time



        ### ******************** Use either Method A or B for pumping:
        ############################ The pumping phase - Method A - pumps atoms into F=1 manifold
        self.dds_pumping_repump.sw.on() ### turning on pumping RP
        delay(1 * ms)  # leave the repump on
        self.dds_pumping_repump.sw.off() ### turning off pumping RP
        delay(1 * ms)

        # ############################ The pumping phase - Method B - pumps atoms into F=2 manifold
        # self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        # delay(1 * ms)  # leave the repump on so atoms are left in F=2
        # self.ttl_repump_switch.on()  # turns the MOT RP AOM off
        # delay(1 * ms)

        ### ************************************************************


        ### Change the bias field for microwave pulse
        self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
                              self.AX_volts_microwave, self.AY_volts_microwave], channels=self.coil_channels)
        delay(0.5 * ms)


        ############################ microwave phase
        if self.t_microwave_pulse > 0.0:
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
            delay(5 * us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(10 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.dds_microwaves.sw.off()

@kernel
def microwave_freq_scan_with_photons_experiment(self):
    """
    This experiment scans the microwave frequency to find the transitions from F=1 to F=2 manifold btw different Zeeman levels.
    This is the first step for atom mapping. Works well and we found the resonance usign this experiment.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- pump atom into F=1, mF=0
    3- excite atom to F=0,mF=0 from which it decays to F=1 manifold
    3- Uses microwave pi pulse to transfer population from one of the Zeeman states in F=1, say mF=0, to F=2
    4- Blow away F=2 atoms
    5- Measure retention

    * If microwave is off resonance, the retention is 100%
    * If microwave is on resonance, the retention drops to ~1/3

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    record_chopped_blow_away(self)
    delay(100 * ms)

    record_chopped_optical_pumping(self)
    delay(100 * ms)
    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)



    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1 * ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(0.1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)
        delay(10 * us)

        first_shot(self)

        delay(1*ms)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        ### With chopped pumping:
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)
            delay(1 * ms)

        # ########################################################
        # # lower level optical pumping and excitation sequence to optimize for speed
        # ########################################################
        # delay(1 * us)
        # self.dds_AOM_A1.sw.off()
        # self.dds_AOM_A2.sw.off()
        # self.dds_AOM_A3.sw.off()
        # self.dds_AOM_A4.sw.off()
        # self.dds_AOM_A5.sw.off()
        # self.dds_AOM_A6.sw.off()
        #
        # delay(1 * us)
        #
        # ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # ### use ttl_excitation to swith on/off D1 or Exc light
        # self.GRIN1and2_dds.sw.on()
        #
        # delay(1 * ms)
        #
        # ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
        # ############################
        # ### optical pumping phase - pumps atoms into F=1,m_F=0
        # ############################
        #
        # if self.t_pumping > 0.0:
        #     self.ttl_repump_switch.on()  # turns off the MOT RP AOM
        #     self.ttl_exc0_switch.on()  # turns off the excitation
        #     self.dds_cooling_DP.sw.off()  # no cooling light
        #     delay(1 * us)
        #
        #     ### set coils for pumping
        #     self.zotino0.set_dac(
        #         [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        #         channels=self.coil_channels)
        #     delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.
        #
        #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
        #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        #     delay(1 * us)
        #
        #     ### Tunring on pumping RP:
        #     self.dds_pumping_repump.sw.on()
        #     self.dds_AOM_A5.sw.on()
        #     self.dds_AOM_A6.sw.on()
        #
        #     delay(1 * ms)
        #
        #     self.ttl_GRIN1_switch.off()  ### Turn on GRIN1 AOM
        #     delay(10 * us)
        #
        #     self.core_dma.playback_handle(op_dma_handle)
        #     delay(self.t_depumping)
        #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
        #     self.dds_pumping_repump.sw.off()  ### turning off pumping RP
        #
        #     delay(2 * us)
        #     self.dds_AOM_A5.sw.off()
        #     self.dds_AOM_A6.sw.off()
        #     delay(100 * us)
        #
        #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        #     delay(1 * ms)
        #
        #     self.ttl_GRIN1_switch.on()  ### Turn off GRIN1 AOM
        #     delay(10 * us)


        if self.p_excitation > -39.0:

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            ### Change the bias field during excitation and detection
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
                channels=self.coil_channels)

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(1 * ms)

            t1 = now_mu()

            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation

            at_mu(t1 + int(self.gate_start_offset_mu))

            ######## Using the edge_counter (works well):
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
            SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()
            #############################################

            delay(15 * us)
            self.ttl_exc0_switch.on()  # block Excitation
            delay(0.1 * ms)  # ttl count consumes all the RTIO slack.


        ############################ microwave phase to transfer population from F=1 to F=2
        if self.t_microwave_pulse > 0.0:
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8*self.stabilizer_FORT.amplitudes[1])
            delay(5*us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(5 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(6 * ms) ### hopefully to avoid underflow.
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_freq_scan_2_with_photons_experiment(self):
    """
    The same as microwave_freq_scan_with_photons_experiment but with smooth microwave pulse rather than square.

    This experiment scans the microwave frequency to find the transitions from F=1 to F=2 manifold btw different Zeeman levels.
    This is the first step for atom mapping. Works well and we found the resonance usign this experiment.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- pump atom into F=1, mF=0
    3- excite atom to F=0,mF=0 from which it decays to F=1 manifold
    3- Uses microwave pi pulse to transfer population from one of the Zeeman states in F=1, say mF=0, to F=2
    4- Blow away F=2 atoms
    5- Measure retention

    * If microwave is off resonance, the retention is 100%
    * If microwave is on resonance, the retention drops to ~1/3

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    record_chopped_blow_away(self)
    delay(100 * ms)

    record_chopped_optical_pumping(self)
    delay(100 * ms)
    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    ############################################ microwave amplitude ramp setting
    self.ramp_time = 2 * us  # ramp time from first to max amplitude
    self.N = 100  # number of amplitude points
    amp_start, amp_end = 0.0, 1.0  # start and end points in scale from 0 to 1

    ### Gaussian function
    x_vals = [-3.0 + 3.0 * i / (self.N - 1) for i in range(self.N)]  ### a list from -3 to 0
    raw = [math.exp(-0.5 * x * x) for x in x_vals]
    g_min, g_max = raw[0], raw[-1]
    norm = [(r - g_min) / (g_max - g_min) for r in raw]

    ### scale into respective ramp according to amplitude. For some reason, we have to generate the list in reverse.
    ### For example, for ramping up amp_points should be decreasing.
    amp_points = [amp_end + n * (amp_start - amp_end) for n in norm]

    ### some data conversion thats needed for RAM
    amplitudes_arr = np.zeros(self.N, dtype=np.int32)
    self.dds.amplitude_to_ram(amp_points, amplitudes_arr)  ### updates arr according to the amp_points profile
    self.amplitudes_list = list(amplitudes_arr)

    ### This is calculation of steps based on above parameters
    self.step_ticks = int((self.ramp_time / self.N) / (4 * ns))  ### AD9910 can update the RAM every 4 clock cycles, i.e. every 4ns.
    ### step_ticks=25, for example, mean the dds is updated every 25*4ns = 100ns.
    ##############################################################################

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1 * ms)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)
        delay(10 * us)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        delay(1 * ms)

        ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
        ############################
        ### optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        if self.t_pumping > 0.0:
            self.ttl_repump_switch.on()  # turns off the MOT RP AOM
            self.ttl_exc0_switch.on()  # turns off the excitation
            self.dds_cooling_DP.sw.off()  # no cooling light
            delay(1 * us)

            ### set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            delay(1 * us)

            ### Tunring on pumping RP:
            self.dds_pumping_repump.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            delay(1 * ms)

            self.ttl_GRIN1_switch.off()  ### Turn on GRIN1 AOM
            delay(10 * us)

            self.core_dma.playback_handle(op_dma_handle)
            delay(self.t_depumping)
            self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            self.dds_pumping_repump.sw.off()  ### turning off pumping RP

            delay(2 * us)
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(100 * us)

            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            delay(1 * ms)

            self.ttl_GRIN1_switch.on()  ### Turn off GRIN1 AOM
            delay(10 * us)


        ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

        ### Change the bias field during excitation and detection
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            channels=self.coil_channels)

        self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
        delay(1 * ms)

        t1 = now_mu()

        self.dds_FORT.sw.off()  ### turns FORT off

        at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
        self.dds_FORT.sw.on()  ### turns FORT on

        at_mu(t1 + int(self.t_excitation_offset_mu))
        self.ttl_GRIN2_switch.off()  # turns on excitation

        at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
        self.ttl_GRIN2_switch.on()  # turns off excitation

        at_mu(t1 + int(self.gate_start_offset_mu))

        ######## Using the edge_counter (works well):
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
            self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
        SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
        SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()
        #############################################

        delay(15 * us)
        self.ttl_exc0_switch.on()  # block Excitation
        delay(0.1 * ms)  # ttl count consumes all the RTIO slack.


        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(0.1 * ms)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

        ############################ microwave phase to transfer population from F=1 to F=2
        if self.t_microwave_pulse > 0.0:
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
            delay(5 * us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(10 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(6 * ms) ### hopefully to avoid underflow.
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    delay(10*ms)
    self.dds_FORT.sw.off()
    delay(1*ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_map01_map11_experiment(self):
    """
    Using the microwave frequencies found with microwave_freq_scan_with_photons_experiment, this experiment maps from
    F=1,mF=0 to F=2,mF=1, then maps F=2,mF=1 to F=1,mF=1. We should be able to transfer >90% of
    population and find the microwave resonance frequencies more accurately than with microwave_freq_scan_with_photons_experiment.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave pi pulse to transfer population from F=1,mF=0 to F=2,mF=1
    4- Apply another microwave pulse to do Rabi oscillation between F=2,mF=1 and F=1,mF=1
    5- Blow away F=2 manifold
    6- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    record_chopped_blow_away(self)
    delay(100 * ms)

    record_chopped_optical_pumping(self)
    delay(100 * ms)
    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        load_MOT_and_FORT_until_atom_recycle(self)
        delay(10 * us)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        delay(1 * ms)

        ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
        ############################
        ### optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        if self.t_pumping > 0.0:
            self.ttl_repump_switch.on()  # turns off the MOT RP AOM
            self.ttl_exc0_switch.on()  # turns off the excitation
            self.dds_cooling_DP.sw.off()  # no cooling light
            delay(1 * us)

            ### set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            delay(1 * us)

            ### Tunring on pumping RP:
            self.dds_pumping_repump.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            delay(1 * ms)

            self.ttl_GRIN1_switch.off()  ### Turn on GRIN1 AOM
            delay(10 * us)

            self.core_dma.playback_handle(op_dma_handle)
            delay(self.t_depumping)
            self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            self.dds_pumping_repump.sw.off()  ### turning off pumping RP

            delay(2 * us)
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(100 * us)

            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            delay(1 * ms)

            self.ttl_GRIN1_switch.on()  ### Turn off GRIN1 AOM
            delay(10 * us)

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1
        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
        delay(5 * us)

        if self.t_microwave_01_pulse > 0.0:
            delay(5 * us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_01_pulse)
            self.ttl_microwave_switch.on()
            delay(5 * us)

        ############################ microwave phase to transfer population from F=2,mF=1 to F=1,mF=1
        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10 * us)

        if self.t_microwave_11_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_11_pulse)
            self.ttl_microwave_switch.on()

        delay(5 * us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(6 * ms)  ### hopefully to avoid underflow.
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    delay(10 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_map01_map11_CORPSE_experiment(self):
    """
    Using the microwave frequencies found with microwave_freq_scan_with_photons_experiment, this experiment maps from
    F=1,mF=0 to F=2,mF=1, then maps F=2,mF=1 to F=1,mF=1. It uses CORPSE pulses for both mapping. CORPSE is going to help mittigate the
    effect of resonance fluctuation.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave CORPSE pulses to transfer population from F=1,mF=0 to F=2,mF=1
    4- Uses microwave CORPSE pulses to transfer population from F=2,mF=1 to F=1,mF=1
    5- Blow away F=2 manifold
    6- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_TRACKING)
    t_microwave_01_CORPSE1 = 420./180. * self.t_microwave_01_pulse
    t_microwave_01_CORPSE2 = 300./180. * self.t_microwave_01_pulse
    t_microwave_01_CORPSE3 = 60./180. * self.t_microwave_01_pulse

    t_microwave_11_CORPSE1 = 420. / 180. * self.t_microwave_11_pulse
    t_microwave_11_CORPSE2 = 300. / 180. * self.t_microwave_11_pulse
    t_microwave_11_CORPSE3 = 60. / 180. * self.t_microwave_11_pulse

    ### overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    record_chopped_blow_away(self)
    delay(100 * ms)

    record_chopped_optical_pumping(self)
    delay(100 * ms)
    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        delay(10 * us)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        delay(1 * ms)

        ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
        ############################
        ### optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        if self.t_pumping > 0.0:
            self.ttl_repump_switch.on()  # turns off the MOT RP AOM
            self.ttl_exc0_switch.on()  # turns off the excitation
            self.dds_cooling_DP.sw.off()  # no cooling light
            delay(1 * us)

            ### set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            delay(1 * us)

            ### Tunring on pumping RP:
            self.dds_pumping_repump.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            delay(1 * ms)

            self.ttl_GRIN1_switch.off()  ### Turn on GRIN1 AOM
            delay(10 * us)

            self.core_dma.playback_handle(op_dma_handle)
            delay(self.t_depumping)
            self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            self.dds_pumping_repump.sw.off()  ### turning off pumping RP

            delay(2 * us)
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(100 * us)

            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            delay(1 * ms)

            self.ttl_GRIN1_switch.on()  ### Turn off GRIN1 AOM
            delay(10 * us)


        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1
        ### set coils for microwave
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            channels=self.coil_channels)
        delay(1 * ms)  # coil relaxation time.

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
        delay(5 * us)

        phi = 0.0
        if self.t_microwave_01_pulse > 0.0:
            self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves), phase = 0.0)
            self.ttl_microwave_switch.off()
            delay(t_microwave_01_CORPSE1)
            # self.ttl_microwave_switch.on()
            # delay(2 * us)

            self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.5 + phi)
            # self.ttl_microwave_switch.off()
            delay(t_microwave_01_CORPSE2)
            # self.ttl_microwave_switch.on()
            # delay(2 * us)

            self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0 + 0)
            # self.ttl_microwave_switch.off()
            delay(t_microwave_01_CORPSE3)
            self.ttl_microwave_switch.on()
            delay(1 * us)

        ############################ microwave phase to transfer population from F=2,mF=1 to F=1,mF=1
        if self.t_microwave_11_pulse > 0.0:
            self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)
            delay(100 * us)
            self.ttl_microwave_switch.off()
            delay(t_microwave_11_CORPSE1)
            self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.5)
            delay(t_microwave_11_CORPSE2)
            self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)
            delay(t_microwave_11_CORPSE3)
            self.ttl_microwave_switch.on()

        delay(10 * us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(6 * ms)  ### hopefully to avoid underflow.
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    delay(10 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def single_photon_experiment(self):
    """
    This experiment pumps the atom into F=1,m=0 then excites it to F=0,0.
    This sequence is repeated multiple times, but only one SPCM is monitored,
    so we can not use this result to verify single photons. We can use it to
    make sure we are only getting one or zero clicks after each excitation
    attempt, and that the one click events only occur when there is an atom
    loaded.

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    # overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0  # not used in this function
    self.SPCM0_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0
    # SPCM0_SinglePhoton_array = [0]
    # rtio_log("2nd_shot_block", 0) # todo: delete. for debugging.

    # self.set_dataset(self.SPCM0_rate_dataset,
    #                  [0.0],
    #                  broadcast=True)

    record_chopped_optical_pumping(self)
    delay(100 * ms)

    if self.verify_OP_in_photon_experiment:
        if self.t_blowaway > 0.0:
            record_chopped_blow_away(self)
            delay(100 * ms)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10 * ms)
        self.dds_microwaves.sw.on()
        delay(100 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    self.measurement = 0
    while self.measurement < self.n_measurements:

        SPCM0_SinglePhoton_array = [0] * self.max_excitation_cycles
        SPCM1_SinglePhoton_array = [0] * self.max_excitation_cycles

        if self.enable_laser_feedback:
            self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

            # bug -- microwave dds is off after AOM feedback; not clear why yet. for now, just turn it back on
            if self.verify_OP_in_photon_experiment:
                self.dds_microwaves.sw.on()

        delay(10 * ms)
        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        ### set coils to PGC settings
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        ### set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        ### set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        delay(0.2 * ms)  ### for coil relaxation

        if not self.no_first_shot:
            first_shot(self)
        delay(1 * ms)
        self.ttl_repump_switch.on()  # turns the RP AOM off

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        # self.ttl7.pulse(10 * us)  # in case we want to look at signals on an oscilloscope

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1 * us)

        # self.ttl_SPCM_gate.on()  # blocks the SPCM output ### remove: SPCM gate no longer exists
        loop_start_mu = now_mu()

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time

        # this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()


        for excitation_cycle in range(self.max_excitation_cycles):  # todo: revert later

            delay(0.5 * ms)

            # low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            # todo: make sure this is consistent with any updates in chopped_optical_pumping function
            ############################
            # optical pumping phase - pumps atoms into F=1,m_F=0
            ############################

            # todo: D1 feedback
            if self.t_pumping > 0.0:
                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(0.0))
                # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))
                # self.ttl_repump_switch.on()  # turns off the MOT RP AOM

                if not self.pumping_light_off:
                    self.dds_pumping_repump.sw.on()

                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(.1 * ms)  # maybe this can be even shorter

                # GRIN1 on
                # self.ttl_excitation_switch.off()
                self.ttl_GRIN1_switch.off()

                with sequential:

                    delay(1 * us)

                    self.core_dma.playback_handle(op_dma_handle)
                    delay(self.t_depumping)

                    self.dds_D1_pumping_DP.sw.off()
                    self.dds_pumping_repump.sw.off()  # turn the repump back on

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                # self.ttl_excitation_switch.on()
                self.ttl_GRIN1_switch.on()

                ############################
                # microwave phase - ONLY USED FOR VERIFYING OP.
                ############################

                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                ############################
                # blow-away phase - push out atoms in F=2 only
                ############################

                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################
            # excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            ############################
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.ttl_exc0_switch.off()  # turns excitation0 AOM on

            t1 = now_mu()

            self.dds_FORT.sw.off() ### turns FORT off
            at_mu(t1 + 150 + int(self.t_photon_collection_time / ns) + 150)
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + 150)
            self.ttl_GRIN1_switch.off() ### turns excitation on

            at_mu(t1 + 150 + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN1_switch.on()  ### turns excitation off

            at_mu(t1 + 150 + int(self.gate_start_offset_mu))
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)

            SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

            delay(10 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            SPCM0_SinglePhoton_array[excitation_cycle] = SPCM0_SinglePhoton
            SPCM1_SinglePhoton_array[excitation_cycle] = SPCM1_SinglePhoton
            delay(0.1 * ms)  # ttl count consumes all the RTIO slack.


            ############################
            # recooling phase
            ############################

            # # todo: use a specific detuning for this stage?
            delay(1 * ms)
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                delay(0.4 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * ms)

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

        self.GRIN1and2_dds.sw.off() ## turns off excitation dds

        delay(1 * ms)

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        delay(1 * ms)

        with sequential:

            self.ttl_repump_switch.off()  # turns the RP AOM on

            # take the second shot
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)
            delay(0.2 * ms)

            second_shot(self)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        end_measurement(self)
        for val in SPCM0_SinglePhoton_array:
            self.append_to_dataset('SPCM0_SinglePhoton', val)
        for val in SPCM1_SinglePhoton_array:
            self.append_to_dataset('SPCM1_SinglePhoton', val)

        delay(10 * ms)

    self.dds_FORT.sw.off()

@kernel
def single_photon_experiment_atom_loading_advance(self):
    """
    This experiment pumps the atom into F=1,m=0 then excites it to F=0,0.
    This sequence is repeated multiple times. We do an atom readout after each pumping/excitation attempt.

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    # overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window

    SPCM0_every_exc_RO_array = [0]

    record_chopped_optical_pumping(self)
    delay(100*ms)

    if self.verify_OP_in_photon_experiment:
        if self.t_blowaway > 0.0:
            record_chopped_blow_away(self)
            delay(100*ms)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10 * ms)
        self.dds_microwaves.sw.on()
        delay(100 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    self.measurement = 0  # advanced in end_measurement
    tries = 0

    while self.measurement < self.n_measurements:

        SPCM0_SinglePhoton_array = [0] * self.max_excitation_cycles
        SPCM1_SinglePhoton_array = [0] * self.max_excitation_cycles

        SPCM0_every_exc_RO_array = [0] * self.max_excitation_cycles

        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles)]

        self.laser_stabilizer.run()

        self.ttl_exc0_switch.on()  # turns off the excitation

        atom_loaded = False

        while not atom_loaded:
            delay(1*ms)
            load_MOT_and_FORT(self)

            delay(0.1 * ms)
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            ### set the FORT AOM to the readout settings
            self.dds_FORT.set(frequency=self.f_FORT,
                              amplitude=self.stabilizer_FORT.amplitudes[1])

            ### set the cooling DP AOM to the readout settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                    amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

            delay(0.4 * ms) ### coil relaxation time


            first_shot(self)

            # tries to load an atom 20 times before running laser feedback again.
            if self.require_atom_loading_to_advance_in_single_photon_exp:
                # if no atoms!!
                if not self.SPCM0_RO1/self.t_SPCM_first_shot > self.single_atom_threshold:
                    if tries < 20:
                        tries += 1
                    else:               # limit: 5% atom loading
                        self.print_async("ERROR: atom loading is below 5%. Stuck at measurement : ", self.measurement, "/ ", self.n_measurements, " --Running feedback--")
                        self.laser_stabilizer.run()
                        tries = 0
                        # todo: pause and resume the experiment? pause => schedule another experiment with higher priority => solve issue => override dataset;
                else:
                    # if atom loaded, initialize tries = 0
                    tries = 0
                    atom_loaded = True
            else:
                atom_loaded = True

        delay(1 * ms)


        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*us)

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time

        # this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()


        for excitation_cycle in range(self.max_excitation_cycles): # todo: revert later

            delay(1*ms)

            # low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            ############################
            # optical pumping phase - pumps atoms into F=1,m_F=0
            ############################

            # todo: D1 feedback
            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0)) ### set to 5V for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.dds_pumping_repump.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                # self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.dds_pumping_repump.sw.off()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                # self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1

                ############################
                # microwave phase - ONLY USED FOR VERIFYING OP.
                ############################

                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                ############################
                # blow-away phase - push out atoms in F=2 only
                ############################

                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################
            # excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            ############################
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(1 * ms)

            t1 = now_mu()

            self.dds_FORT.sw.off() ### turns FORT off

            at_mu(t1 + 150 + int(self.t_photon_collection_time / ns) + 150)
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1+150)
            self.ttl_GRIN1_switch.off() # turns on excitation

            at_mu(t1 + 150 + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN1_switch.on()  # turns off excitation

            at_mu(t1 + 150 + int(self.gate_start_offset_mu))

            ######### Using the edge_counter (works well):
            # with parallel:
            #     self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
            #     self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
            # # delay(15 * us)
            # SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

            ######### Using ttl.count (works well):
            # with parallel:
            #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
            #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
            # SPCM0_SinglePhoton = self.ttl_SPCM0.count(t_end_SPCM0)
            # SPCM1_SinglePhoton = self.ttl_SPCM1.count(t_end_SPCM1)

            ######### Using time stamping and counting the photons (works well):
            SPCM0_click_counter = 0
            SPCM1_click_counter = 0

            with parallel:
                t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

            ### timestamping SPCM0 events
            while SPCM0_click_counter < max_clicks:
                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                if SPCM0_click_time == -1.0:
                    break
                SPCM0_timestamps[excitation_cycle][SPCM0_click_counter] = self.core.mu_to_seconds(SPCM0_click_time - t1)
                SPCM0_click_counter += 1

            ### timestamping SPCM1 events
            while SPCM1_click_counter < max_clicks:
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                if SPCM1_click_time == -1.0:
                    break
                SPCM1_timestamps[excitation_cycle][SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time - t1)
                SPCM1_click_counter += 1

            ### counting by adding up the timestamps
            SPCM0_SinglePhoton = SPCM0_click_counter
            SPCM1_SinglePhoton = SPCM1_click_counter

            delay(15 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            SPCM0_SinglePhoton_array[excitation_cycle] = SPCM0_SinglePhoton
            SPCM1_SinglePhoton_array[excitation_cycle] = SPCM1_SinglePhoton
            delay(0.1 * ms)  # ttl count consumes all the RTIO slack.



            ############################
            # recooling phase
            ############################

            t_SPCM_recool_and_shot_mu = int(self.t_SPCM_recool_and_shot / ns)

            # # todo: use a specific detuning for this stage?
            delay(1*ms)
            if self.t_recooling > 0 or self.record_every_shot:

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                delay(0.4*ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                delay(1 * us)
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()



                if self.record_every_shot:

                    self.ttl_SPCM0_counter.gate_rising(t_SPCM_recool_and_shot_mu * ns)
                    SPCM0_every_exc_RO = self.ttl_SPCM0_counter.fetch_count()
                    SPCM0_every_exc_RO_array[excitation_cycle] = SPCM0_every_exc_RO

                    delay(1*ms)

                else:
                    delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1*ms)

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time


        self.GRIN1and2_dds.sw.off()

        delay(1*ms)

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()


        delay(1*ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)

        delay(10 * ms)

        for val in SPCM0_SinglePhoton_array:
            self.append_to_dataset('SPCM0_SinglePhoton', val)
        for val in SPCM1_SinglePhoton_array:
            self.append_to_dataset('SPCM1_SinglePhoton', val)
        for val in SPCM0_every_exc_RO_array:
            self.append_to_dataset('SPCM0_every_exc_RO', val)

        delay(1 * ms)
        for i in range(len(SPCM0_timestamps)):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        for i in range(len(SPCM1_timestamps)):
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])

        delay(10*ms)

    delay(10 * ms)
    self.dds_FORT.sw.off()

@kernel
def single_photon_experiment_2_atom_loading_advance(self):
    """
    This is based on load_MOT_and_FORT_until_atom for fast atom loading. Checks for atom after each excitation.

    This experiment pumps the atom into F=1,m=0 then excites it to F=0,0.
    This sequence is repeated multiple times. We do an atom readout after each pumping/excitation attempt.

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    ### overwritten below but initialized here so they are always initialized
    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0
    SPCM0_SinglePhoton = 0
    SPCM1_SinglePhoton = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window

    SPCM0_every_exc_RO_array = [0]

    record_chopped_optical_pumping(self)
    delay(100*ms)

    if self.verify_OP_in_photon_experiment:
        if self.t_blowaway > 0.0:
            record_chopped_blow_away(self)
            delay(100*ms)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10 * ms)
        self.dds_microwaves.sw.on()
        delay(100 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        delay(1 * ms)

        SPCM0_SinglePhoton_array = [0] * self.max_excitation_cycles
        SPCM1_SinglePhoton_array = [0] * self.max_excitation_cycles

        SPCM0_every_exc_RO_array = [0] * self.max_excitation_cycles

        ### the first element of SPCM0_timestamps will be the reference time t1.
        SPCM0_timestamps = [[-1.0] * (max_clicks + 1) for _ in range(self.max_excitation_cycles)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles)]

        self.ttl_exc0_switch.on()  # turns off the excitation

        load_MOT_and_FORT_until_atom(self)
        delay(1 * ms)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*us)

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time

        # this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):

            delay(1*ms)

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            ############################
            ### optical pumping phase - pumps atoms into F=1,m_F=0
            ############################

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0)) ### set to 5V for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.dds_pumping_repump.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                # self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.dds_pumping_repump.sw.off()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                # self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1

                ############################
                # microwave phase - ONLY USED FOR VERIFYING OP.
                ############################

                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                ############################
                # blow-away phase - push out atoms in F=2 only
                ############################

                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################
            # excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            ############################
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(1 * ms)

            t1 = now_mu()

            self.dds_FORT.sw.off() ### turns FORT off

            at_mu(t1 + 150 + int(self.t_photon_collection_time / ns) + 150)
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1+150)
            self.ttl_GRIN1_switch.off() # turns on excitation

            at_mu(t1 + 150 + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN1_switch.on()  # turns off excitation

            at_mu(t1 + 150 + int(self.gate_start_offset_mu))

            ######### Using the edge_counter (works well):
            # with parallel:
            #     self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
            #     self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
            # # delay(15 * us)
            # SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

            ######### Using ttl.count (works well):
            # with parallel:
            #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
            #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
            # SPCM0_SinglePhoton = self.ttl_SPCM0.count(t_end_SPCM0)
            # SPCM1_SinglePhoton = self.ttl_SPCM1.count(t_end_SPCM1)

            ######### Using time stamping and counting the photons (works well):
            SPCM0_click_counter = 0
            SPCM1_click_counter = 0

            with parallel:
                t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

            ### timestamping SPCM0 events
            ### We use the element of SPCM0_timestamps to keep t1 for reference
            SPCM0_timestamps[excitation_cycle][0] = self.core.mu_to_seconds(t1)
            while SPCM0_click_counter < max_clicks:
                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                if SPCM0_click_time == -1.0:
                    break
                SPCM0_timestamps[excitation_cycle][SPCM0_click_counter + 1] = self.core.mu_to_seconds(SPCM0_click_time)
                SPCM0_click_counter += 1

            ### timestamping SPCM1 events
            while SPCM1_click_counter < max_clicks:
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                if SPCM1_click_time == -1.0:
                    break
                SPCM1_timestamps[excitation_cycle][SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                SPCM1_click_counter += 1

            ### counting by adding up the timestamps
            SPCM0_SinglePhoton = SPCM0_click_counter
            SPCM1_SinglePhoton = SPCM1_click_counter

            delay(15 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            SPCM0_SinglePhoton_array[excitation_cycle] = SPCM0_SinglePhoton
            SPCM1_SinglePhoton_array[excitation_cycle] = SPCM1_SinglePhoton
            delay(0.1 * ms)  # ttl count consumes all the RTIO slack.


            ############################
            # readout to see if the atom survived
            ############################
            # delay(1*ms)

            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            delay(0.4*ms)

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.off()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            delay(1 * us)
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
            SPCM0_every_exc_RO = self.ttl_SPCM0_counter.fetch_count()
            SPCM0_every_exc_RO_array[excitation_cycle] = SPCM0_every_exc_RO

            ### stopping the excitation cycle after the atom is lost
            if SPCM0_every_exc_RO / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * ms)

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                break

            delay(1*ms)

            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(1*ms)

            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(0.4 * ms)  # coil relaxation time


        self.GRIN1and2_dds.sw.off()

        delay(1*ms)

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()


        delay(1*ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)

        delay(10 * ms)

        ### only the elements in range [0:excitation_cycle] contain non-zero values because the loop exits after
        ### the atom is lost.
        for val in SPCM0_SinglePhoton_array[0:excitation_cycle + 1]:
            self.append_to_dataset('SPCM0_SinglePhoton', val)
        for val in SPCM1_SinglePhoton_array[0:excitation_cycle + 1]:
            self.append_to_dataset('SPCM1_SinglePhoton', val)
        for val in SPCM0_every_exc_RO_array[0:excitation_cycle + 1]:
            self.append_to_dataset('SPCM0_every_exc_RO', val)

        delay(1 * ms)
        for i in range(excitation_cycle + 1):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        for i in range(excitation_cycle + 1):
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        delay(10*ms)

    delay(10 * ms)
    self.dds_FORT.sw.off()

@kernel
def single_photon_experiment_3_atom_loading_advance(self):
    """
    This is based on load_MOT_and_FORT_until_atom. Does not check for atom after each excitation attempt:

    for excitation_cycle in range(self.max_excitation_cycles):
        O.P.
        three excitation attempts, for example
        cooling (~5ms)
        R.O. every 5 cycles, for example
        if atom detected -> continue the excitation_cycle loop
        else: break the for loop, record n_excitation_cycles, and go to atom loading.

    Then we can plot n_excitation_cycles (multiples of 5) as a function of excitation attempt or cooling time, etc.
    Since there is no RO after each excitation, all data is assumed to be with_atom; there is no distinction between
    with_atom and no_atom.

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    BothSPCMs_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(100*ms)

    if self.verify_OP_in_photon_experiment:
        if self.t_blowaway > 0.0:
            record_chopped_blow_away(self)
            delay(100*ms)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(10 * ms)
        self.dds_microwaves.sw.on()
        delay(100 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        BothSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_MOT_and_FORT_until_atom(self)
        delay(1 * ms)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):

            delay(1000 * us)

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:

                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0)) ### set to 5V for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.dds_pumping_repump.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.dds_pumping_repump.sw.off()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1
                delay(10 * us)

                ############ microwave phase - ONLY USED FOR VERIFYING OP.
                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                ############ blow-away phase - push out atoms in F=2 only
                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(1 * ms)

            for excitation_attempt in range(self.n_excitation_attempts):

                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                ### timestamping SPCM0 events
                while SPCM0_click_counter < max_clicks:
                    SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                    if SPCM0_click_time == -1.0:
                        break
                    SPCM0_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)
                    SPCM0_click_counter += 1

                ### timestamping SPCM1 events
                while SPCM1_click_counter < max_clicks:
                    SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                    if SPCM1_click_time == -1.0:
                        break
                    SPCM1_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                    SPCM1_click_counter += 1

                # at_mu(t1 + 30000)
                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                delay(1 * us)
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                delay(1 * us)
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)


            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(0.4*ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                delay(1 * us)
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                BothSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = BothSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if BothSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        delay(5 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in BothSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('BothSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        delay(1*ms)

    delay(15 * ms)

@kernel
def atom_photon_partity_1_experiment(self):
    """
    A simple parity oscillation experiment. It uses a few excitation_attempts, but not excitation_cycle. So, goes
    back to atom loading after say 5 excitation attempts. However, it recycles atoms. So effectively, we are
    attempting excitation after about 100ms.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- Blow away F=2 and measure retention
    7- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    BothSPCMs_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        delay(0.1 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(0.1 * ms)

        load_MOT_and_FORT_until_atom_recycle(self)
        delay(1 * ms)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        delay(1000 * us)

        ############################### optical pumping phase - pumps atoms into F=1,m_F=0
        if self.t_pumping > 0.0:

            self.ttl_repump_switch.on()  # turns off the MOT RP AOM
            self.ttl_exc0_switch.on()  # turns off the excitation
            self.dds_cooling_DP.sw.off()  # no cooling light
            delay(1 * us)

            ### set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            delay(1 * us)

            ### Tunring on pumping RP:
            self.dds_pumping_repump.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

            delay(1 * ms)

            self.ttl_GRIN1_switch.off()
            delay(10 * us)

            self.core_dma.playback_handle(op_dma_handle)
            delay(self.t_depumping)
            self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            self.dds_pumping_repump.sw.off()  ### turning off pumping RP

            delay(2 * us)
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(100 * us)

            self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            delay(1 * ms)

            self.ttl_GRIN1_switch.on()
            delay(10 * us)

            ############ microwave phase - ONLY USED FOR VERIFYING OP.
            if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                delay(5 * us)
                self.ttl_microwave_switch.off()
                delay(self.t_microwave_pulse)
                self.ttl_microwave_switch.on()
                delay(5 * us)
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            ############ blow-away phase - push out atoms in F=2 only
            if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                chopped_blow_away(self)

        ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

        self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
        delay(5 * ms)

        for excitation_attempt in range(self.n_excitation_attempts):

            t1 = now_mu()

            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation

            at_mu(t1 + int(self.gate_start_offset_mu))

            ######### Using the edge_counter (works well):
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
            SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

            delay(20 * us)  ### 20us is not enough

            if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                if self.BothSPCMs_RO1/self.t_SPCM_first_shot > self.single_atom_threshold:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_11_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(100 * us)
                    chopped_blow_away(self)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('BothSPCMs_parity_RO', self.BothSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)

                break

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1*ms)
        self.GRIN1and2_dds.sw.off()

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)

        delay(5 * ms)

    delay(15 * ms)

@kernel
def atom_photon_partity_2_experiment(self):
    """
    A simple parity oscillation experiment. To speed up the experiment, I reuse the atom and
    do OP after say 5 excitation attempts and try excitation again. I repeate the loop as long as
    there is an atom measured every atom_check_every_n.
    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- Blow away F=2 and measure retention
    7- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    # BothSPCMs_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        # BothSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        SPCM0_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        SPCM1_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)

        delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_MOT_and_FORT_until_atom_recycle(self)
        delay(1 * ms)

        first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.BothSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(1000 * us)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:

                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.dds_pumping_repump.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.dds_pumping_repump.sw.off()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM


            for excitation_attempt in range(self.n_excitation_attempts):
                delay(200 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))

                ######### Using the edge_counter (works well):
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
                SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

                delay(20 * us)

                if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.8 * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_11_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(200 * us)
                    chopped_blow_away(self)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('BothSPCMs_parity_RO', self.BothSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(100*us)
            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                delay(1 * us)
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                delay(1 * us)
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)


            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(0.4 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                delay(1 * us)
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
                delay(0.1 * ms)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

                delay(1*ms)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                delay(1 * us)
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

                ### stopping the excitation cycle after the atom is lost
                if BothSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

                else:
                    atom_loaded = False

            excitation_cycle +=1


        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1*ms)
        self.GRIN1and2_dds.sw.off()

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        self.measurement -= 1
        end_measurement(self)


        delay(5 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python stops the loop one count earlier.
        # for val in BothSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
        #     self.append_to_dataset('BothSPCMs_RO_atom_check', val)

        # delay(1 * ms)
        # for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(15 * ms)

@kernel
def rotator_test_experiment(self):
    """
    A simple experiment to test using GVS to control the waveplates. You can scan FORT waveplates and see the effect on MM fiber on scope.

    :param self:
    :return:
    """

    self.core.reset()

    if self.enable_laser_feedback:
        self.laser_stabilizer.run()

    self.dds_FORT.sw.on() ### Tunrning on the FORT to see the effect of waveplate rotation on scope
    delay(1* ms)


    ### GVS - set target_hwp_deg as the scan_variable.
    move_to_target_deg(self, name="852_HWP", target_deg=self.target_hwp_deg)
    # wait_move(self, "852_HWP")
    hwp852_pos = get_rotator_position(self, '852_HWP')
    delay(1*s)
    self.print_async('hwp852 at ', hwp852_pos / self.deg_to_pos, ' deg')
    delay(1*s)

    # record_PDA_power(self)


    # move_to_target_deg(self, name="780_HWP", target_deg=10)
    # hwp780_deg = get_rotator_deg(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_deg, ' deg')


    # # delay(2*s)
    # # move_to_target_deg(self, name="780_HWP", target_deg=350)
    #
    # go_to_home(self, '780_HWP')
    #
    # hwp780_deg = get_rotator_deg(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_deg, ' deg')

@kernel
def prepare_FORT_RAM_profile(self):
    """
    Prepares FORT RAM waveform to be used to ramp down and ramp up the FORT smoothly.
    """

    FORT_steps_rise = 50  # number of amplitude points during the rise and fall time
    FORT_amp_high = self.stabilizer_FORT.amplitudes[0] # low and high amplitudes in scale from 0 to 1
    FORT_amp_low = self.stabilizer_FORT.amplitudes[1]

    # print(FORT_amp_high)
    # delay(10*ms)
    # print(FORT_amp_low)
    # delay(10 * ms)

    FORT_step_ticks = int(
        (self.t_FORT_ramp / FORT_steps_rise) / (
                4 * ns))  ### Step size

    FORT_dwell_time = 100*us

    FORT_steps_dwell = int(
        FORT_steps_rise / self.t_FORT_ramp * FORT_dwell_time)  # Number of dwell points

    FORT_x_vals = [i / (FORT_steps_rise - 1) for i in range(FORT_steps_rise)]
    FORT_norm = [6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3 for x in FORT_x_vals] # quintic smoothstep function

    FORT_amp_points_rise = [FORT_amp_low + n * (FORT_amp_high - FORT_amp_low) for n in FORT_norm]
    FORT_amp_points_fall = [FORT_amp_points_rise[i] for i in range(FORT_steps_rise - 1, -1, -1)]

    ### The full waveform.
    FORT_amp_points = (
            FORT_amp_points_fall +
            [FORT_amp_low] * FORT_steps_dwell +
            FORT_amp_points_rise
    )

    ### some data conversion needed for RAM
    FORT_amplitudes_arr = [0] * len(FORT_amp_points)
    self.dds_FORT.amplitude_to_ram(FORT_amp_points, FORT_amplitudes_arr)
    FORT_amplitudes_list = list(FORT_amplitudes_arr)

    ### This is calculation of steps based on above parameters
    FORT_total_points = len(FORT_amplitudes_list)

    self.FORT_step_size = FORT_step_ticks
    self.FORT_total_points = FORT_total_points
    self.FORT_amplitudes_list = FORT_amplitudes_list

    delay(1*ms)

@kernel
def Testing_FORT_ramp_experiment(self):
    """
    Simply to test the FORT ramping down and up. Prepares the RAM amplitudes on kernel. Works well.
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    delay(1*ms)

    prepare_FORT_RAM_profile(self)

    self.dds_FORT.set_frequency(self.f_FORT)  ### takes 0.7us
    self.dds_FORT.set_att(0.0)  ### takes 1.6us

    self.dds_FORT.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
    self.dds_FORT.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings

    ### Configures the RAM playback engine
    self.dds_FORT.set_profile_ram(
        start=0,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP,
    )

    ### write the data onto RAM
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    self.dds_FORT.write_ram(self.FORT_amplitudes_list)

    # ### This is not required. It is just for testing if we can ramp the dds without interupting the power.
    # ### But the moment I enable RAM below, it turns off the dds for about 2ms before playing RAM.
    # self.dds_FORT.sw.on()
    # self.dds_FORT.cpld.io_update.pulse_mu(8)
    # delay(3 * ms)


    ### Enabling RAM playback, not playing yet.
    self.dds_FORT.set_cfr1(ram_enable=1,
                           ram_destination=RAM_DEST_ASF)

    self.ttl7.on()

    # delay(1*ms)

    ### Configure the RAM to playback the second half. This effectively turns on the FORT to loading set point.
    self.dds_FORT.set_profile_ram(
        start=self.FORT_total_points // 2,
        end=self.FORT_total_points - 1,
        step=self.FORT_step_size,
        profile=7,
        mode=RAM_MODE_RAMPUP)
    self.dds_FORT.cpld.io_update.pulse_mu(8)

    ### Running the RAM
    self.dds_FORT.sw.on()
    self.dds_FORT.cpld.io_update.pulse_mu(8)
    delay(10*ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(1 * ms)

        delay(1*ms)

        ### Configure the RAM to playback the first half
        self.dds_FORT.set_profile_ram(
            start=0,
            end=len(self.FORT_amplitudes_list) // 2 - 1,
            step=self.FORT_step_size,
            profile=7,
            mode=RAM_MODE_RAMPUP)
        self.dds_FORT.cpld.io_update.pulse_mu(8)

        delay(2 * ms)  # Leave at-least enough time to cover up the first ram time, plus any extra time we want.

        ### turning the FORT off and on to simulate FORT drop, for example, in atom temperature measruement.
        self.dds_FORT.sw.off()
        delay(1*ms)
        self.dds_FORT.sw.on()
        delay(2*ms)

        ### Configure the RAM to playback the second half
        self.dds_FORT.set_profile_ram(
            start=self.FORT_total_points // 2,
            end=self.FORT_total_points - 1,
            step=self.FORT_step_size,
            profile=7,
            mode=RAM_MODE_RAMPUP)
        self.dds_FORT.cpld.io_update.pulse_mu(8)

        self.ttl7.off()

        delay(5 * ms)
        self.measurement += 1

    delay(1 * ms)
    ### Exit RAM Mode
    self.dds_FORT.set_cfr1(ram_enable=0)
    self.dds_FORT.cpld.io_update.pulse_mu(8)


    ### If we want to change frequency or amplitude after exiting RAM mode, we should do it like below
    ### with profile=0. The first time we change the amplitude the FORT gets off for about 2us.
    ### The second time, the FORT stays ON all the time, just the amplitude changes.
    # self.dds_FORT.cpld.set_profile(0)
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.2)
    delay(2*ms)
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=0.3)

    delay(10*ms)

    self.dds_FORT.sw.off()


#
# ###############################################################################
# # 3. DIAGNOSTIC FUNCTIONS
# # These are functions that are used for various tests, but are not typical
# # experiments. You might not want to run them with GeneralVariableScan, but you
# # can import these functions from other files. The advantage of having them in
# # experiment_functions is that they can benefit from the other functions defined
# # here.
# ###############################################################################

@kernel
def FORT_monitoring_with_Luca_experiment(self):
    """
    A modified version of atom_loading_experiment for monitoring the FORT scattering.
    This is most easily run with tests/MonitorFORTPowerWithLuca.py

    Load a MOT, load a single atom, readout, wait self.t_delay_between_shots, readout again.

    For this experiment, we were still using the MM fiber monitor (i.e. after the polarizer) to
    feed back to the FORT power. We used this experiment to characterize how well the signal
    from an APD monitoring the scattered 852 nm light in the chamber correlated with the signal
    seen by the Luca in an ROI restricted to a section of the parabolic mirror tube.

    The assumed connection of the Sampler card and channel used to monitor the APD is given
    in load_MOT_and_FORT_for_Luca_scattering_measurement.

    ----
    Alternatively, if you are using the APD to feedback to the FORT:
    The MM fiber is monitored with the channel specified in load_MOT_and_FORT_for_Luca_scattering_measurement.

    The purpose of this experiment is to be able to compare the normalized FORT power recorded
    by the dds_FORT feedback channel, the scattering seen by the camera, and the voltage of
    the detector monitoring the MM fiber after the polarizer. Rotating the 852 nm
    motorized waveplates, e.g. with the APT program, we can check whether the scattering in the
    chamber seen with the Luca correlates with polarization. This is a way to vet our APD
    feedback scheme.

    ----
    For analysis, see FORT feedback/monitor_FORT_scattering_and_Raman_light.ipynb

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()

    self.set_dataset("SPCM0_FORT_loading", [0.0], broadcast=True)
    self.set_dataset("SPCM0_FORT_and_MOT", [0.0], broadcast=True)
    self.set_dataset("SPCM0_FORT_science", [0.0], broadcast=True)
    self.set_dataset("APD_FORT_volts_loading", [0.0], broadcast=True)
    self.set_dataset("APD_FORT_volts_science", [0.0], broadcast=True)

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0

    self.require_D1_lock_to_advance = False # override experiment variable
    self.require_atom_loading_to_advance = False # override experiment variable

    # self.set_dataset(self.SPCM0_rate_dataset,
    #                  [0.0],
    #                  broadcast=True)

    delay(100*ms)
    now = now_mu()
    for i in range(self.warm_up_shots):
        at_mu(now+i*self.core.seconds_to_mu(500*ms))
        self.ttl_Luca_trigger.pulse(5 * ms)
        # now = now_mu()

    self.core.wait_until_mu(now_mu())

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            self.laser_stabilizer.run(monitor_only=self.no_feedback)  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT_for_Luca_scattering_measurement(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        # set the FORT AOM to the science setting. this is only valid if we have run
        # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT
        if not self.no_feedback:
            self.dds_FORT.set(frequency=self.f_FORT,
                                    amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        # take the first shot
        if not self.MOT_light_off:
            self.dds_cooling_DP.sw.on()
        with parallel:
            self.ttl_Luca_trigger.pulse(5 * ms)
            t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(t_gate_end)
        delay(1 * ms)
        self.dds_cooling_DP.sw.off()

        delay(self.t_delay_between_shots)

        # take the second shot
        if not self.MOT_light_off:
            self.dds_cooling_DP.sw.on()
        # with parallel:
            # self.ttl_Luca_trigger.pulse(5 * ms)
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)

        # update experiment-specific datasets:
        self.append_to_dataset("SPCM0_FORT_loading", self.SPCM0_FORT_loading)
        self.append_to_dataset("SPCM0_FORT_and_MOT", self.SPCM0_FORT_and_MOT)
        self.append_to_dataset("SPCM0_FORT_science", self.SPCM0_FORT_science)
        self.append_to_dataset("APD_FORT_volts_loading", self.APD_FORT_volts_loading)
        self.append_to_dataset("APD_FORT_volts_science", self.APD_FORT_volts_science)

    self.dds_FORT.sw.off()

@kernel
def atom_state_mapping(self):
    """ by Eunji """
    detuning = self.f_microwaves_detuning

    ### mapping |1,+1> to |2,1>
    self.dds_microwaves.set(frequency=self.f_microwaves_dds + 2 * detuning, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)

    self.ttl_microwave_switch.off()
    delay(self.t_pi_microwave_pulse)  # todo: change the pulse time
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

    ### mapping |2,+1> to |1, 0>
    self.dds_microwaves.set(frequency=self.f_microwaves_dds + 1 * detuning, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)

    self.ttl_microwave_switch.off()
    delay(self.t_pi_microwave_pulse)  # todo: change the pulse time
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

    ### mapping |1,-1> to |2,0>
    self.dds_microwaves.set(frequency=self.f_microwaves_dds - detuning, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)

    self.ttl_microwave_switch.off()
    delay(self.t_pi_microwave_pulse)  # todo: change the pulse time
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

@kernel
def atom_rotation(self, t_pulse):
    """ by Eunji """
    ### rotating |1,0> and  |2,0>
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)

    self.ttl_microwave_switch.off()
    delay(t_pulse)
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

@kernel
def atom_rotation_x(self):
    """
    At the atom state mapping stage:
        |1,-1> mapped to |2,0>
        |1,+1> mapped to |1,0>
    Here, we rotate the atoms using microwave pi/2 pulse
    for measurement in X basis

    phase set to 0.0

    by Eunji
    """
    ### rotating |1,0> and |2,0>
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

    self.ttl_microwave_switch.off()
    delay(self.t_pi_microwave_pulse/2)
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

@kernel
def atom_rotation_y(self):
    """
    At the atom state mapping stage:
        |1,-1> mapped to |2,0>
        |1,+1> mapped to |1,0>
    Here, we rotate the atoms using microwave pi/2 pulse
    for measurement in X basis

    phase set to np.pi/2

    by Eunji
    """
    # phase in SI unit (radians)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves), phase=np.pi/2)

    self.ttl_microwave_switch.off()
    delay(self.t_pi_microwave_pulse/2)
    self.ttl_microwave_switch.on()
    delay(0.1 * ms)

@kernel
def atom_photon_tomography_experiment(self):
    """
    by Eunji
    """

    # todo: SPCM0_RO_atom_check should be replaced by BothSPCMs_RO_atom_check. See single_photon_experiment_3_atom_loading_advance.

    self.core.reset()
    delay(1 * ms)

    # K10CR1 set the waveplates
    with parallel:      # note: this does not make two wavplates to rotate at the same time.
        # GVS variable - hwp_move_to_deg, qwp_move_to_deg

        #todo: atom state projection via photon measurement - 2D scan

        move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
        move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
        # kernel waits until this job is done.
        # might have to add delay here to avoid underflow error.
        # ex) delay(time_to_rotate_in_ms(self, current_hwp - previous_hwp) * ms)


    self.core.reset()  # to guarantee positive slack after moving the waveplates

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    SPCM0_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(100 * ms)

    record_chopped_blow_away(self)
    delay(100 * ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on()
    delay(100 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()
        # bug -- microwave dds is off after AOM feedback; not clear why yet. for now, just turn it back on
        self.dds_microwaves.sw.on()

    self.measurement = 0  # advances in end_measurement

    SPCM0_total_click_counter = 0
    SPCM1_total_click_counter = 0

    while self.measurement < self.n_measurements:

        SPCM0_RO_atom_check_array = [0] * int(self.max_excitation_cycles / self.atom_check_every_n)
        tStamps_t1 = [0.0] * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        delay(40 * ms)  ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_MOT_and_FORT_until_atom(self)
        delay(1 * ms)

        first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use ttl_excitation to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 1  ### just for initialization.

        # to keep everything consistent to single_photon_experiment, I'll just use these parameters
        # rather than getting rid of the for_loop

        self.max_excitation_cycles = 1   # one excitation cycle.
        self.n_excitation_attempts = 10  # increase the number; increasing the success probability of getting a photon out

        for excitation_cycle in range(self.max_excitation_cycles):

            delay(1000 * us)

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:

                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation,
                                        amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.dds_pumping_repump.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                # self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.dds_pumping_repump.sw.off()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                # self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1

                ############ microwave phase - readout
                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(0.1 * ms)

                ############ blow-away phase - push out atoms in F=2 only
                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(1 * ms)

            for excitation_attempt in range(self.n_excitation_attempts):

                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + 50)
                self.ttl_GRIN1_switch.off()  # turns on excitation

                at_mu(t1 + 50 + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN1_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                ############################ photon state measurement - use count or timestamp

                # Notes on timestamping:
                # .timestamp_mu(up_to_timestamp_mu): timestamp (in machine units) of the first event received; -1 on timeout.
                # up_to_timestamp_mu: timestamp up to which execution is blocked
                # up to which input events are guaranteed to be taken into account.
                #  (Events with later timestamps might still be registered if they are already available.)

                ### timestamping SPCM0 events
                while SPCM0_click_counter < max_clicks:
                    SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                    if SPCM0_click_time == -1.0:
                        break
                    SPCM0_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)
                    SPCM0_click_counter += 1
                    SPCM0_total_click_counter += 1  # total counts received during n_measurements

                ### timestamping SPCM1 events
                while SPCM1_click_counter < max_clicks:
                    SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                    if SPCM1_click_time == -1.0:
                        break
                    SPCM1_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                    SPCM1_click_counter += 1
                    SPCM1_total_click_counter += 1  # total counts received during n_measurements

                # at_mu(t1 + 30000)
                tStamps_t1[
                    excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation


            ############################ atom state measurement

            ############################
            # microwave phase - transfer |1,-1> -> |2,+1> with uw+RF
            ############################

            # coils already set to OP.

            # map the states to |1,0> and |2,0>
            atom_state_mapping(self)

            # rotate the states to measure in different basis
            if self.atom_rotation_to_x:
                atom_rotation_x(self)
            elif self.atom_rotation_to_y:
                atom_rotation_y(self)


            ############################
            # blow-away phase - push out atoms in F=2 only
            ############################

            if self.t_blowaway > 0.0:
                chopped_blow_away(self)


            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()

        delay(0.1 * ms)

        ############################
        # readout phase
        ############################


        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        delay(1 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in SPCM0_RO_atom_check_array[0:int(excitation_cycle / self.atom_check_every_n)]:
            self.append_to_dataset('SPCM0_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1) * self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)



        delay(1 * ms)


    delay(15 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)
    self.dds_microwaves.sw.off()


    # saves the total clicks in each SPCM0 and SPCM1 - applet
    self.append_to_dataset('SPCM0_total_click_counter', SPCM0_total_click_counter)
    self.append_to_dataset('SPCM1_total_click_counter', SPCM1_total_click_counter)

@kernel
def Pulse_microwave_smooth(self, MW_freq):
    """
    This uses dds RAM profiles prepared in BaseExperiment to generate smooth MW pulses. Works but blocks
    other dds channels. Needs fix.
    Use it in the experiment like:
    Pulse_microwave_smooth(self, self.f_microwaves_00_dds)

    """
    ### predefine these 3 parameters. Needed but not used.
    MW_step_size = 1
    MW_total_points = 1
    MW_amplitudes_list = [0]
    MW_pulse_length = 0.0

    if MW_freq == self.f_microwaves_00_dds:
        MW_step_size = self.MW_00_step_size
        MW_total_points = self.MW_00_total_points
        MW_amplitudes_list = self.MW_00_amplitudes_list
        MW_pulse_length = self.t_microwave_00_pulse

    if MW_freq == self.f_microwaves_01_dds:
        MW_step_size = self.MW_01_step_size
        MW_total_points = self.MW_01_total_points
        MW_amplitudes_list = self.MW_01_amplitudes_list
        MW_pulse_length = self.t_microwave_01_pulse

    if MW_freq == self.f_microwaves_11_dds:
        MW_step_size = self.MW_11_step_size
        MW_total_points = self.MW_11_total_points
        MW_amplitudes_list = self.MW_11_amplitudes_list
        MW_pulse_length = self.t_microwave_11_pulse

    self.dds_microwaves.set_frequency(MW_freq) ### takes 0.7us
    self.dds_microwaves.set_att(0.0) ### takes 1.6us


    if MW_pulse_length > 0.0:

        self.dds_microwaves.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
        self.dds_microwaves.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings
        ### the two lines above take 0.7us

        ### Configures the RAM playback engine. Takes 1.3us
        self.dds_microwaves.set_profile_ram(
            start=0,
            end=MW_total_points - 1,
            step=MW_step_size,
            profile=0,
            mode=RAM_MODE_RAMPUP,
        )

        self.dds_microwaves.cpld.set_profile(0) ### takes 0.4us
        # self.dds.cpld.io_update.pulse_mu(8)
        self.ttl_microwave_switch.off()  ### takes 0us

        self.dds_microwaves.write_ram(MW_amplitudes_list)  ### write the data onto RAM.
        ### Takes 33us with 60 MW_amp_points. Takes 48us with 90 points. Takes 96us with 180 points.

        ### Enabling RAM playback, not playing yet. takes 0.7us
        self.dds_microwaves.set_cfr1(
            internal_profile=0,
            ram_enable=1,
            ram_destination=RAM_DEST_ASF,
        )

        self.ttl7.on() ## for triggering scope
        self.dds_microwaves.sw.on()
        self.dds_microwaves.cpld.io_update.pulse_mu(8)  ### This runs the RAM

        delay(MW_pulse_length)  ### keep the delay as ramp time

        ### shutting off. Takes 0.7us
        self.dds_microwaves.set_cfr1(ram_enable=0)
        self.dds_microwaves.cpld.io_update.pulse_mu(8)
        self.dds_microwaves.sw.off()

        self.core.reset()
        self.urukul2_cpld.init()
        self.urukul2_ch0.cpld.init()
        self.urukul2_ch1.cpld.init()
        self.urukul2_ch2.cpld.init()
        self.urukul2_ch3.cpld.init()

        self.ttl_microwave_switch.on()
        self.ttl7.off()

@kernel
def track_1_microwave_transition_experiment(self):
    """
    Starting with testing the dataset read.
    1- Read the dataset and set it to f01. Use getattr to avoid "is already define message"
    2- Use that as the starting point to find the new res. f.
    3- update the dataset. We never want to run the subroutine before exiting the kernel. GVS exits the kernel between each scan.
    So, the dataset gets updated with the new res freq which can be used for the next run of the subroutine.

    """
    self.core.reset()
    delay(1*ms)



    # try:
    #     f01 = self.f01_resonance
    # except ValueError:
    #     f01 = self.get_dataset("f_microwaves_01_dds")

    # self.print_async(f01)
    # delay(10*ms)


    # self.default_setpoints = np.array([getattr(self, dataset) for dataset in self.setpoint_datasets])

    # f = getattr(self, self.f_microwaves_01_dds)

    # self.set_dataset("f_microwaves_01_dds_test", 0.0, broadcast=True)
    #
    # self.f_microwaves_01_dds = 333.0 * MHz
    # self.append_to_dataset("f_microwaves_01_dds", 0.0)


