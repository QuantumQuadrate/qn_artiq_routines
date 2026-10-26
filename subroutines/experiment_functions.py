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
from utilities.helper_functions import _rand_shift, _count_threshold_crossings_in_loading_window

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
def run_feedback_and_record_FORT_MM_power(self, record_power = True):
    """
    Function:
        1. Runs feedback to everything in list (6 MOT AOMs and 3 FORT setpoints)
        2. records FORT MM and APD powers

    * IF you want to only run feedback and disable recording, set "record_power = False"

    """
    self.stabilizer_FORT.run(setpoint_index=2)  # FORT science holding setpoint
    delay(0.1*ms)
    self.stabilizer_FORT.run(setpoint_index=1)  # FORT science setpoint
    delay(0.1*ms)
    self.laser_stabilizer.run() # 6 MOT AOMs and FORT loading setpoint
    delay(0.1*ms)

    ## if self.laser_stabilizer.run() is in the last sequence, it will leave the FORT at loading setpoint.

    ## record FORT MM and APD powers
    if record_power:
        ### The three stabilizer runs above budget their time with fixed
        ### delays and never call break_realtime, so any host work that
        ### outruns those budgets erodes slack permanently. Nothing below is
        ### timing-critical relative to the feedback, so resync the cursor
        ### before driving the FORT switch rather than inheriting whatever
        ### slack the feedback left. (Observed as RTIOUnderflow on channel 39
        ### with ~-5.9 ms slack under the DRTIO master gateware, which costs
        ### slightly more CPU per RTIO event than the standalone build.)
        self.core.break_realtime()
        self.dds_FORT.sw.on()  ### turns FORT on
        delay(0.1*ms)

        record_FORT_MM_power(self)
        record_FORT_APD_power(self)

    ### Feedback is wall-clock work: three stabilizer runs, each doing
    ### blocking sampler reads and dataset writes, budgeted only with fixed
    ### 0.1 ms delays. Whatever the host overruns is slack lost for good,
    ### and with record_power=False nothing above restores it, so the
    ### CALLER's next RTIO event is submitted in the past -- observed as
    ### RTIOUnderflow at -5.58 ms on channel 56 (dds_microwaves.sw)
    ### immediately after this returned. Restoring here covers every call
    ### site at once; no caller can rely on a deterministic cursor across a
    ### feedback run, and the branch above has restored since 2026-09-10.
    self.core.break_realtime()

@kernel
def tune_shims_for_atom_loading(self) -> TInt32:
    """
    Helper function for automatic shim tuning during atom loading.
    Scan AX_volts_MOT and AY_volts_MOT within +/-0.3 V of their current values.
    Keep MOT+FORT on continuously (assumed already on when this is called).
    Use single_atom_threshold_for_loading as threshold (counts/s).
    Update AX_volts_MOT and AY_volts_MOT datasets ONLY if strictly better than baseline.
    Returns 1 if updated, 0 otherwise.

    It checks atom loading at shim values with grid points that vary slightly in each trial.
    Otherwise, the grid points will always be limited to a few shim values and the code toggles
    between these values only.

    Akbar 2025-12-17
    """
    self.core.break_realtime()

    # --- config ---
    gate_time = 20 * ms
    n_samples = 50                  # 1s = gate_time * n_samples evaluation window
    target_crossings = 8           # early exit threshold

    coarse_step = 0.15
    fine_step   = 0.05
    fine_span   = 0.10              # search +/- 0.10 V around best in fine pass

    # --- anchor (safety box around current shims) ---
    ax0 = self.AX_volts_MOT
    ay0 = self.AY_volts_MOT
    ax_min, ax_max = ax0 - 0.3, ax0 + 0.3
    ay_min, ay_max = ay0 - 0.3, ay0 + 0.3

    # If current point is already good enough, do nothing
    self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, ax0, ay0],
                         channels=self.coil_channels)
    delay(3 * s) ### warm up
    baseline = _count_threshold_crossings_in_loading_window(
        self, self.single_atom_threshold_for_loading, gate_time, n_samples)
    if baseline >= target_crossings:
        return 0

    best_atoms = baseline
    best_ax, best_ay = ax0, ay0

    # --- per-call grid dithers (change EVERY call, even hours later) ---
    coarse_shift_x = _rand_shift(self, coarse_step, 1)
    coarse_shift_y = _rand_shift(self, coarse_step, 2)
    fine_shift_x   = _rand_shift(self, fine_step,   3)
    fine_shift_y   = _rand_shift(self, fine_step,   4)

    # --- coarse scan: dx,dy in [-0.3, +0.3] step 0.15, with dithers ---
    # points per axis: 5  (i=0..4)
    for ix in range(5):
        dx = (-0.30 + ix * coarse_step) + coarse_shift_x
        ax = ax0 + dx
        if ax < ax_min: ax = ax_min
        elif ax > ax_max: ax = ax_max

        for iy in range(5):
            dy = (-0.30 + iy * coarse_step) + coarse_shift_y
            ay = ay0 + dy
            if ay < ay_min: ay = ay_min
            elif ay > ay_max: ay = ay_max

            delay(1 * ms)

            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, ax, ay],
                                 channels=self.coil_channels)
            delay(1 * ms)

            atoms = _count_threshold_crossings_in_loading_window(
                self, self.single_atom_threshold_for_loading, gate_time, n_samples)

            # EARLY EXIT: accept immediately if "good enough"
            if atoms >= target_crossings:
                self.AX_volts_MOT = ax
                self.AY_volts_MOT = ay
                self.set_dataset("AX_volts_MOT", ax, broadcast=True, persist=True)
                self.set_dataset("AY_volts_MOT", ay, broadcast=True, persist=True)
                self.print_async("Target loading reached at ", atoms)
                delay(1 * ms)
                return 1

            if atoms > best_atoms:
                self.print_async("Best atom loading crossing threshold so far (coarse scan): ", atoms)
                best_atoms, best_ax, best_ay = atoms, ax, ay

    # --- fine scan around best: +/- fine_span with step 0.05, with dithers ---
    # points per axis: 5  (i=0..4) covering [-0.10, -0.05, 0, +0.05, +0.10] + shift
    for ix in range(5):
        dx = (-fine_span + ix * fine_step) + fine_shift_x
        ax = best_ax + dx
        if ax < ax_min: ax = ax_min
        elif ax > ax_max: ax = ax_max

        for iy in range(5):
            dy = (-fine_span + iy * fine_step) + fine_shift_y
            ay = best_ay + dy
            if ay < ay_min: ay = ay_min
            elif ay > ay_max: ay = ay_max

            delay(1 * ms)

            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, ax, ay],
                                 channels=self.coil_channels)
            delay(1 * ms)

            atoms = _count_threshold_crossings_in_loading_window(
                self, self.single_atom_threshold_for_loading, gate_time, n_samples)

            if atoms >= target_crossings:
                self.AX_volts_MOT = ax
                self.AY_volts_MOT = ay
                self.set_dataset("AX_volts_MOT", ax, broadcast=True, persist=True)
                self.set_dataset("AY_volts_MOT", ay, broadcast=True, persist=True)
                self.print_async("Target loading reached at ", atoms)
                delay(1 * ms)
                return 1

            if atoms > best_atoms:
                self.print_async("Best atom loading crossing threshold so far (fine scan): ", atoms)
                best_atoms, best_ax, best_ay = atoms, ax, ay

    # --- update only if improved ---
    if (best_atoms > baseline) and ((best_ax != ax0) or (best_ay != ay0)):
        self.AX_volts_MOT = best_ax
        self.AY_volts_MOT = best_ay
        self.set_dataset("AX_volts_MOT", best_ax, broadcast=True, persist=True)
        self.set_dataset("AY_volts_MOT", best_ay, broadcast=True, persist=True)
        delay(1 * ms)
        return 1

    return 0

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

    delay(0.01 * ms)

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
    # delay(0.1 * ms)
    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    self.dds_cooling_DP.sw.on()
    self.ttl_repump_switch.off()

    delay(1*ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    ### wait for the MOT to load
    delay(self.t_MOT_loading)

    # # todo: make this work for bob too.
    # if self.which_node == 'alice':
    #     self.stabilizer_FORT.run(setpoint_index=1) # the science setpoint
    # elif self.which_node == 'bob':
    #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_RO)

    self.dds_cooling_DP.sw.off()
    self.ttl_repump_switch.on()
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)
    delay(1*ms)

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)  ## coils relaxation time

    ### Lower the FORT to science setpoint
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


    ###########  PGC on the trapped atom  #############
    if self.do_PGC_after_loading:
        ### set the cooling DP AOM to the PGC settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        self.ttl_repump_switch.off()  ### turn on MOT RP
        self.dds_cooling_DP.sw.on()  ### turn on cooling
        delay(10 * us)
        if self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(self.t_PGC_after_loading)  ### this is the PGC time
        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling
    ###################################################

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
        Sampler_test(self)
        delay(1 * ms)

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

    if self.which_node == "alice":
        delay(1 * ms)
        self.zotino0.set_dac([3.5], self.UV_trig_channel)

    max_tries = 100  ### Maximum number of attempts before running the feedback
    atom_check_time = self.t_atom_check_time
    atom_loaded = False
    try_n = 0
    t_before_atom = now_mu() ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
    t_after_atom = now_mu()
    time_without_atom = 0.0

    AllSPCMs_atom_check_loaded = 0  ### for initilization
    AllSPCMs_atom_check_not_loaded = 0

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    while True:
        while not atom_loaded and try_n < max_tries:
            delay(100 * us)  ### Needs a delay of about 100us or maybe less
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

            # SPCM0_atom_check = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_atom_check = self.ttl_SPCM1_counter.fetch_count()

            # AllSPCMs_atom_check = int((SPCM0_atom_check + SPCM1_atom_check) / 2)
            AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                  self.ttl_SPCM1_counter.fetch_count() + \
                                  self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                  self.ttl_SPCM1_OtherNode_counter.fetch_count())

            try_n += 1

            ### To save only one photon counts of unloaded case for each loaded atom. Otherwise, the unloaded counts
            ### would overwhelm the dataset.
            if try_n == 1:
                AllSPCMs_atom_check_not_loaded = AllSPCMs_atom_check

            if AllSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True
                AllSPCMs_atom_check_loaded = AllSPCMs_atom_check


        if atom_loaded:
            # t_before_atom = t_after_atom ### I don't know why I had this! Removed and seems working fine.
            self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
            t_after_atom = now_mu()

            ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
            self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_loaded)
            self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
            delay(1 * ms)
            break  ### Exit the outer loop if an atom is loaded

        #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
        delay(0.1 * ms)
        t_no_atom = now_mu()
        time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
        self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

        ### If max_tries reached and still no atom, run feedback
        if self.enable_laser_feedback:
            # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
            # delay(0.1 * ms)
            # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
            delay(0.1 * ms) ### necessary to avoid underflow

            ### Same as load_until_atom_smooth_FORT_recycle: the slack is eaten
            ### by the set_dataset RPC above, NOT by the DRTIO crossing -- both
            ### destinations measured the same 1 us minimum lead, so Node1 is
            ### exposed here too. This branch has already given up on the
            ### current loading attempt, so nothing here is timing-critical.
            self.core.break_realtime()

            ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            run_feedback_and_record_FORT_MM_power(self, record_power = False)
            self.n_feedback_per_iteration += 1
            # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
            self.dds_microwaves.sw.on()
            self.dds_FORT.sw.on()
            delay(0.1 * ms)

            try_n = 0

    ### The loading loop above drains a few us of slack per iteration:
    ### host RPCs advance wall clock while the timeline cursor does not
    ### (set_dataset(broadcast) measured at ~4 us through this master,
    ### plus print_async and the per-try CPU work), and nothing in the
    ### loop ever restores it. After a few hundred tries the deficit
    ### reaches milliseconds and the first RTIO event after the loop --
    ### this Zotino write, channel 33 -- is submitted in the past.
    ### Seen at -0.51 ms and -3.58 ms. The loop has just ended and
    ### nothing here is timing-critical, so restore the slack.
    self.core.break_realtime()

    self.zotino0.set_dac([0.0], self.UV_trig_channel)
    delay(100*us)

    ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(0.4 * ms)

    self.ttl_repump_switch.on()  ### turn off MOT RP
    self.dds_cooling_DP.sw.off()  ### turn off cooling

    delay(1 * ms)
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

    ### Lower the FORT to science setpoint
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


    ###########  PGC on the trapped atom  #############
    if self.do_PGC_after_loading:
        ### set the cooling DP AOM to the PGC settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        delay(1*us)
        if self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        self.ttl_repump_switch.off()  ### turn on MOT RP
        self.dds_cooling_DP.sw.on()  ### turn on cooling
        delay(10 * us)

        delay(self.t_PGC_after_loading)  ### this is the PGC time
        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling
    ###################################################
    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

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
        if self.AllSPCMs_RO2/self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True

            ### Lower the FORT to science setpoint
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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
                if not self.PGC_and_RO_with_on_chip_beams:
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
            Sampler_test(self)
            delay(1*ms)

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

        if self.which_node == 'alice':
            delay(1 * ms)
            self.zotino0.set_dac([3.5], self.UV_trig_channel)

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu() ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        AllSPCMs_atom_check_loaded = 0 ### for initilization
        AllSPCMs_atom_check_not_loaded = 0

        while True:
            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

                # AllSPCMs_atom_check = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)
                AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                           self.ttl_SPCM1_counter.fetch_count() + \
                                           self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                           self.ttl_SPCM1_OtherNode_counter.fetch_count())
                try_n += 1

                ### To save only one photon counts of unloaded case for each loaded atom. Otherwise, the unloaded counts
                ### would overwhelm the dataset.
                if try_n==1:
                    AllSPCMs_atom_check_not_loaded = AllSPCMs_atom_check

                if AllSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True
                    AllSPCMs_atom_check_loaded = AllSPCMs_atom_check


            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
                t_after_atom = now_mu()

                ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_loaded)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
                delay(1 * ms)
                break  ### Exit the outer loop if an atom is loaded

            #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
            delay(0.1 * ms)
            t_no_atom = now_mu()
            time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
            self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms) ### necessary to avoid underflow

                ### Same as load_until_atom_smooth_FORT_recycle: the slack is
                ### eaten by the set_dataset RPC above, NOT by the DRTIO
                ### crossing -- both destinations measured the same 1 us
                ### minimum lead, so Node1 is exposed here too. This branch has
                ### already given up on the current loading attempt, so nothing
                ### here is timing-critical.
                self.core.break_realtime()

                ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                run_feedback_and_record_FORT_MM_power(self, record_power = False)
                self.n_feedback_per_iteration += 1
                # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
                self.dds_microwaves.sw.on()
                self.dds_FORT.sw.on()
                delay(0.1 * ms)

                try_n = 0

        ### The loading loop above drains a few us of slack per iteration:
        ### host RPCs advance wall clock while the timeline cursor does not
        ### (set_dataset(broadcast) measured at ~4 us through this master,
        ### plus print_async and the per-try CPU work), and nothing in the
        ### loop ever restores it. After a few hundred tries the deficit
        ### reaches milliseconds and the first RTIO event after the loop --
        ### this Zotino write, channel 33 -- is submitted in the past.
        ### Seen at -0.51 ms and -3.58 ms. The loop has just ended and
        ### nothing here is timing-critical, so restore the slack.
        self.core.break_realtime()

        if self.which_node == 'alice':
            self.zotino0.set_dac([0.0], self.UV_trig_channel)
        delay(100*us)

        ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(0.4 * ms)

        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling

        delay(1 * ms)
        delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

        ### Lower the FORT to science setpoint
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


        ###########  PGC on the trapped atom  #############
        if self.do_PGC_after_loading:
            ### set the cooling DP AOM to the PGC settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
            if self.PGC_and_RO_with_on_chip_beams:
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
            self.ttl_repump_switch.off()  ### turn on MOT RP
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            # delay(10 * us)

            delay(self.t_PGC_after_loading)  ### this is the PGC time
            self.ttl_repump_switch.on()  ### turn off MOT RP
            self.dds_cooling_DP.sw.off()  ### turn off cooling
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
def load_until_atom_smooth_FORT_recycle(self):
    """
    Based on load_MOT_and_FORT_until_atom_recycle but lowering FORT smoothly to Science set point instead of step function.

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
    ### This will miss the first measurement of each iteration (except the first iteration)
    ### and results in the first RO1 to be less than threshold. But it is OK.
    delay(100 * us)
    if self.AllSPCMs_RO2/self.t_SPCM_second_shot > self.single_atom_threshold:
        atom_loaded = True

        # ### Lower the FORT to science setpoint
        # if self.which_node == 'alice':
        #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        # elif self.which_node == 'bob':
        #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        ###########  PGC on the trapped atom  #############
        if self.do_PGC_after_loading:
            ### Set the coils to PGC setting
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)
            delay(1 * ms)
            ### set the cooling DP AOM to the PGC settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            delay(0.1 * ms)
            if self.PGC_and_RO_with_on_chip_beams:
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

            self.dds_cooling_DP.sw.on()  ### turn on cooling
            self.ttl_repump_switch.off()  ### turn on MOT RP
            delay(self.t_PGC_after_loading)  ### this is the PGC time
            self.dds_cooling_DP.sw.off()  ### turn off cooling
            self.ttl_repump_switch.on()  ### turn off MOT RP
        ###################################################
    else:
        atom_loaded = False

    ### load an atom if atom_loaded = False
    if not atom_loaded:

        if self.monitors_for_atom_loading:
            measure_Magnetometer(self)
            delay(1*ms)
            Sampler_test(self)
            delay(1*ms)
            measure_GRIN1(self)
            delay(1 * ms)

        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        ### set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(1 * ms)

        self.dds_cooling_DP.sw.on()  ### turn on cooling
        self.ttl_repump_switch.off()  ### turn on MOT RP
        delay(0.1 * ms)

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        delay(0.1 * ms)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
        self.dds_FORT.sw.on()

        delay(1 * ms)
        self.zotino0.set_dac([3.5], self.UV_trig_channel)

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu() ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        AllSPCMs_atom_check_loaded = 0  ### for initilization
        AllSPCMs_atom_check_not_loaded = 0

        shim_tune_runs = 0

        while True:
            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

                AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())

                try_n += 1

                ### To save only one photon counts of unloaded case for each loaded atom. Otherwise, the unloaded counts
                ### would overwhelm the dataset.
                if try_n == 1:
                    AllSPCMs_atom_check_not_loaded = AllSPCMs_atom_check

                if AllSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True
                    AllSPCMs_atom_check_loaded = AllSPCMs_atom_check


            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
                t_after_atom = now_mu()

                ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_loaded)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
                delay(1 * ms)
                break  ### Exit the outer loop if an atom is loaded

            #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
            delay(0.1 * ms)
            t_no_atom = now_mu()
            time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
            self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

            ### If max_tries reached and atom loading is too bad, run shim tuning.
            t_shim_tune_trigger = 10.0  # seconds
            max_shim_tune_runs = 3  # prevents infinite tuning loop

            if self.tune_shims_when_loading_is_bad and (time_without_atom > t_shim_tune_trigger) and (shim_tune_runs < max_shim_tune_runs):
                self.print_async("Atom loading is bad. Tuning X and Y shims.")
                tune_shims_for_atom_loading(self)
                shim_tune_runs += 1

                ### Shim tuning is host-heavy for the same reason; restore
                ### before going back round the loading loop.
                self.core.break_realtime()

                ### restart the loading attempt with the (possibly) improved shims
                try_n = 0
                t_before_atom = now_mu()
                delay(0.1 * ms)
                continue

            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms) ### necessary to avoid underflow

                ### NOT a DRTIO penalty. measure_drtio_event_cost (2026-09-10)
                ### showed master and satellite behave identically: both accept
                ### a 1 us lead and both underflow at 0.5 us. The slack is spent
                ### BEFORE this call -- the set_dataset RPC above and the loop's
                ### host-side work advance wall clock while the cursor stands
                ### still. The observed failure was 0.586 ms behind, ~500x past
                ### that floor, so no plausible delay() would have covered it.
                ### Nothing here is timing-critical (the loop has given up and
                ### feedback is wall-clock work), so restore the slack instead.
                self.core.break_realtime()

                # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
                # delay(0.1 * ms)
                # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

                ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                run_feedback_and_record_FORT_MM_power(self, record_power = False)
                self.n_feedback_per_iteration += 1
                # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
                self.dds_microwaves.sw.on()
                self.dds_FORT.sw.on()
                delay(0.1 * ms)

                try_n = 0

        ### The loading loop above drains a few us of slack per iteration:
        ### host RPCs advance wall clock while the timeline cursor does not
        ### (set_dataset(broadcast) measured at ~4 us through this master,
        ### plus print_async and the per-try CPU work), and nothing in the
        ### loop ever restores it. After a few hundred tries the deficit
        ### reaches milliseconds and the first RTIO event after the loop --
        ### this Zotino write, channel 33 -- is submitted in the past.
        ### Seen at -0.51 ms and -3.58 ms. The loop has just ended and
        ### nothing here is timing-critical, so restore the slack.
        self.core.break_realtime()

        self.zotino0.set_dac([0.0], self.UV_trig_channel)
        delay(100*us)

        ### Set the coils to PGC setting even when we don't want PGC. Effectively, this is turning off coils.
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(1 * ms)

        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling

        delay(1 * ms)
        delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

        ### Lower the FORT to science setpoint
        FORT_ramp1_smoothstep(self, direction="down")
        # if self.which_node == 'alice':
        #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        #     FORT_ramp1_smoothstep(self, direction="down")
        # elif self.which_node == 'bob':
        #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)
        #     FORT_ramp1_smoothstep(self, direction="down")

        ###########  PGC on the trapped atom  #############
        if self.do_PGC_after_loading:
            ### set the cooling DP AOM to the PGC settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
            self.ttl_repump_switch.off()  ### turn on MOT RP
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            # delay(10 * us)
            if self.PGC_and_RO_with_on_chip_beams:
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
            delay(self.t_PGC_after_loading)  ### this is the PGC time
            self.ttl_repump_switch.on()  ### turn off MOT RP
            self.dds_cooling_DP.sw.off()  ### turn off cooling
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
def first_shot(self):
    """
    first atom readout.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP

    :return:
    """

    ### set the FORT AOM to the science settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(0.1 * ms)

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    if not self.PGC_and_RO_with_on_chip_beams:
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    if self.use_chopped_readout:
        chopped_RO_handle = self.core_dma.get_handle("chopped_RO")
        delay(10 * ms)
        self.ttl_repump_switch.off()  ### turn on MOT RP
        ### DMA playback
        self.core_dma.playback_handle(chopped_RO_handle)
        delay(10 * ms)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(now_mu())
        self.SPCM1_RO1 = self.ttl_SPCM1.count(now_mu())
        self.AllSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        delay(10 * us)

    else:
        # self.ttl_repump_switch.off()  ### turn on MOT RP
        # self.dds_cooling_DP.sw.on()  ### Turn on cooling
        # delay(0.1 * ms)
        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
        # delay(1 * ms)
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_first_shot)
            self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_first_shot)
            self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_first_shot)
            self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_first_shot)

        self.SPCM0_RO1 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO1 = self.ttl_SPCM1_counter.fetch_count()
        self.SPCM0_OtherNode_RO1 = self.ttl_SPCM0_OtherNode_counter.fetch_count()
        self.SPCM1_OtherNode_RO1 = self.ttl_SPCM1_OtherNode_counter.fetch_count()
        self.AllSPCMs_RO1 = self.SPCM0_RO1 + self.SPCM1_RO1 + self.SPCM0_OtherNode_RO1 + self.SPCM1_OtherNode_RO1
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off() ### turn off cooling
        self.ttl_repump_switch.on() ### turn off MOT RP
        delay(10 * us)

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

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
    delay(1 * ms)  ## coils relaxation time

    ### set the FORT AOM to the readout settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
    # FORT_ramp2_smoothstep(self, direction="up")

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    if not self.PGC_and_RO_with_on_chip_beams:
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
    else:
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
    delay(0.1 * ms)

    if self.use_chopped_readout:
        chopped_RO_handle = self.core_dma.get_handle("chopped_RO")
        delay(10 * ms)
        ### DMA playback
        self.core_dma.playback_handle(chopped_RO_handle)
        delay(10 * ms)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(now_mu())
        self.SPCM1_RO2 = self.ttl_SPCM1.count(now_mu())
        self.AllSPCMs_RO2 = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)
        delay(0.1 * ms)
        # self.dds_cooling_DP.sw.off()  ### turn off cooling
        # self.ttl_repump_switch.on()  ### turn off MOT RP
        # delay(10 * us)

    # if self.use_chopped_readout:
    #     # rtio_log("chop_RO_counter", 0)
    #     # rtio_log("chop_RO_dma", 0)
    #     delay(10*ms)
    #     ro_dma_handle = self.core_dma.get_handle("second_chopped_readout")
    #     delay(10 * ms)
    #
    #     # todo set RO coils here?
    #
    #     self.ttl_repump_switch.off()  # turns the RP AOM on
    #     self.dds_cooling_DP.sw.off()  # the chop sequence likes to turn the FORT off
    #
    #     delay(1 * ms)
    #     # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
    #     # delay(0.1 * ms)
    #     # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
    #     self.dds_FORT.sw.on()  # the chop sequence likes to turn the FORT off
    #
    #     delay(10 * ms)
    #
    #     # we want to initiate the chop playback and read in detector clicks while the chop sequence is playing.
    #     # the edge_counter.gate_rising(duration) function is equivalent to:
    #     #         edge_counter.set_config(
    #     #             count_rising=count_rising,
    #     #             count_falling=count_falling,
    #     #             send_count_event=False,
    #     #             reset_to_zero=True)
    #     #         delay_mu(duration_mu)
    #     #         edge_counter.set_config(
    #     #             count_rising=False,
    #     #             count_falling=False,
    #     #             send_count_event=True,
    #     #             reset_to_zero=False)
    #     #
    #     # we want the dma playback to happen during the gating, so we call the _set_sensitivity functions directly
    #
    #     now = now_mu()
    #     rtio_log("chop_RO_counter", 1)
    #     self.ttl_SPCM0_counter.set_config(
    #             count_rising=True,
    #             count_falling=False,
    #             send_count_event=False,
    #             reset_to_zero=True)
    #     delay(1*us)
    #     rtio_log("chop_RO_dma", 1)
    #     self.core_dma.playback_handle(ro_dma_handle) # not sure if I need extra delay here
    #     rtio_log("chop_RO_dma", 0)
    #     at_mu(now + self.core.seconds_to_mu(self.t_SPCM_second_shot+10*us))
    #     self.ttl_SPCM0_counter.set_config(
    #             count_rising=False,
    #             count_falling=False,
    #             send_count_event=True,
    #             reset_to_zero=False)
    #     self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
    #     rtio_log("chop_RO_counter", 0)
    #     delay(10*ms)

    else:
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
            self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)
            self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)
            self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)

        self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
        self.SPCM0_OtherNode_RO2 = self.ttl_SPCM0_OtherNode_counter.fetch_count()
        self.SPCM1_OtherNode_RO2 = self.ttl_SPCM1_OtherNode_counter.fetch_count()
        self.AllSPCMs_RO2 = self.SPCM0_RO2 + self.SPCM1_RO2 + self.SPCM0_OtherNode_RO2 + self.SPCM1_OtherNode_RO2

        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off() ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        delay(10 * us)

    # ### set the FORT AOM back to loading setting
    # if self.which_node == 'alice':
    #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[0])
    # elif self.which_node == 'bob':
    #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

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
    delay(1 * ms)  ## coils relaxation time

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

    delay(10 * us)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    if not self.PGC_and_RO_with_on_chip_beams:
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
    else:
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
    delay(10 * us)

    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(10 * us)

    with parallel:
        self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)
        self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)

    self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
    self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
    self.SPCM0_OtherNode_RO2 = self.ttl_SPCM0_OtherNode_counter.fetch_count()
    self.SPCM1_OtherNode_RO2 = self.ttl_SPCM1_OtherNode_counter.fetch_count()
    self.AllSPCMs_parity_RO = self.SPCM0_RO2 + self.SPCM1_RO2 + self.SPCM0_OtherNode_RO2 + self.SPCM1_OtherNode_RO2

    delay(0.1 * ms)
    self.dds_cooling_DP.sw.off()  ### turn off cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP
    delay(10 * us)

@kernel
def recooling_after_first_shot(self):
    """
    Recooling after first shot. Improves the atom temperature.
    At the moment, this is a Node2 specific code, using only MOT 2, 4 for recooling.

    Eunji
    """
    delay(1 * ms)
    self.dds_AOM_A1.sw.off()   ## 1D PGC with beam 2,4
    # self.dds_AOM_A2.sw.off()
    self.dds_AOM_A3.sw.off()   ## 1D PGC with beam 2,4
    # self.dds_AOM_A4.sw.off()
    # self.dds_AOM_A5.sw.off()  ## already off if RO only done with on-chip beams.
    # self.dds_AOM_A6.sw.off()  ## already off if RO only done with on-chip beams.

    delay(10 * us)
    self.dds_cooling_DP.sw.on()  ### turn on cooling
    self.ttl_repump_switch.off()  ### turn on MOT RP

    delay(self.t_recooling_after_first_shot)

    self.dds_cooling_DP.sw.off()  ### turn off cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP
    delay(10 * us)

    ###todo: decide to turn off or on
    ## Turning AOM 1 and 3 back ON so they are in the same state as after first_shot,
    ## to minimize disruptions to the current code.
    # self.dds_AOM_A1.sw.on()
    # self.dds_AOM_A3.sw.on()

@kernel
def record_chopped_RO(self):
    ### Akbar 2025-11-25
    t_FORT_offset = self.core.seconds_to_mu(100 * ns)  ### delay wrt "t0" in each cycle
    t_FORT_OFF = self.core.seconds_to_mu(self.chop_RO_FORT_off)  ### the length of the FORT OFF time

    t_cooling_offset = self.core.seconds_to_mu(self.chop_RO_pulse_offset)  ### delay wrt "t0" in each cycle
    t_cooling_length = self.core.seconds_to_mu(self.chop_RO_pulse_length)  ### duration of the cooling pulse

    t_gate_offset = self.core.seconds_to_mu(self.chop_RO_gate_offset)  ### delay wrt t0 for gate to open
    t_gate_length = self.core.seconds_to_mu(self.chop_RO_gate_length)  ### gate open time

    with self.core_dma.record("chopped_RO"):
        for i in range(4000): ### number of cycles to repeat
            t0 = now_mu()

            at_mu(t0 + t_FORT_offset)
            self.dds_FORT.sw.off()
            at_mu(t0 + t_FORT_offset + t_FORT_OFF)
            self.dds_FORT.sw.on()

            at_mu(t0 + t_cooling_offset)
            # self.dds_MW_RF.sw.on()  ### a dummy dds to turn on instead of cooling
            self.dds_cooling_DP.sw.on()  ### Turn on cooling
            at_mu(t0 + t_cooling_offset + t_cooling_length)
            # self.dds_MW_RF.sw.on()  ### a dummy dds to turn off instead of cooling
            self.dds_cooling_DP.sw.off()  ### Turn off cooling

            at_mu(t0 + t_gate_offset)
            self.ttl0._set_sensitivity(1)
            self.ttl1._set_sensitivity(1)

            at_mu(t0 + t_gate_offset + t_gate_length)
            self.ttl0._set_sensitivity(0)
            self.ttl1._set_sensitivity(0)

            delay(0.5 * us)

@kernel
def chopped_RO(self):
    """
    chopped readout for testing the effect of it on retention, for example.
    Akbar 2025-11-26
    """

    ### Set the coils to PGC_optimization setting to scan the values
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization, self.AX_volts_PGC_optimization,
         self.AY_volts_PGC_optimization],
        channels=self.coil_channels)
    delay(1 * ms)

    ### set the cooling DP AOM to the PGC optimization setting to scan the values
    ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)

    ### set the FORT AOM to the science settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

    self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    delay(1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    if not self.PGC_and_RO_with_on_chip_beams:
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
    else:
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
    delay(0.1 * ms)

    chopped_RO_handle = self.core_dma.get_handle("chopped_RO")
    delay(10 * ms)
    self.ttl_repump_switch.off()  ### turn on MOT RP
    ### DMA playback
    self.core_dma.playback_handle(chopped_RO_handle)
    delay(10 * ms)
    self.SPCM0_test_RO = self.ttl_SPCM0.count(now_mu())
    delay(0.1 * ms)
    self.dds_cooling_DP.sw.off()  ### turn off cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP

    self.zotino0.set_dac([0.0], self.Osc_trig_channel)
    delay(1 * ms)

    ### Set the coils to PGC setting:
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC,
         self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(1 * ms)

    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
    delay(1 * ms)

@kernel
def record_chopped_blow_away(self):
    """Record the chopped blow-away sequence as DMA."""

    n_chop_cycles = int(self.t_blowaway / self.t_BA_chop_period + 0.5)
    assert n_chop_cycles >= 1, "t_blowaway should be >= t_BA_chop_period"

    BA_pulse = self.t_BA_chop_period * 0.35
    period_mu = self.core.seconds_to_mu(self.t_BA_chop_period)
    BA_pulse_length_mu = self.core.seconds_to_mu(BA_pulse)

    with self.core_dma.record("chopped_blow_away"):

        self.ttl_repump_switch.on()  ### Repump AOM OFF

        ### Set blow-away coil voltages using forward scheduling
        self.zotino0.write_dac(self.coil_channels[0], self.AZ_bottom_volts_blowaway)
        self.zotino0.write_dac(self.coil_channels[1], -self.AZ_bottom_volts_blowaway)
        self.zotino0.write_dac(self.coil_channels[2], self.AX_volts_blowaway)
        self.zotino0.write_dac(self.coil_channels[3], self.AY_volts_blowaway)
        delay(1.5 * us)
        self.zotino0.load()

        delay(1 * ms)  ### Physical coil-relaxation time

        ### Configure blow-away light
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_blowaway, amplitude=self.ampl_cooling_DP_MOT)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()

        if self.blowaway_light_off:
            self.dds_cooling_DP.sw.off()
            self.dds_AOM_A6.sw.off()
        else:
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-7.0))
            self.dds_AOM_A6.sw.on()
            self.dds_cooling_DP.sw.on()

        ### Chopped FORT sequence
        start = now_mu()

        for i in range(n_chop_cycles):
            at_mu(start + i * period_mu)
            self.dds_FORT.sw.off()

            at_mu(start + i * period_mu + BA_pulse_length_mu)
            self.dds_FORT.sw.on()

        at_mu(start + n_chop_cycles * period_mu)
        self.dds_FORT.sw.on()

        ### Turn off blow-away light
        self.dds_cooling_DP.sw.off()
        self.ttl_repump_switch.on()
        self.dds_AOM_A6.sw.off()

        ### Restore normal RF settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

@kernel
def record_CW_optical_pumping(self):
    """Record the CW optical-pumping DMA trace "CW_optical_pumping".

    The D1 optical pumping of every single-node experiment that pumps with D1,
    on both nodes since the 2026-10-20 Node1 rewire, after which Node1's GRIN
    AOMs are wired like Node2's: dds_GRIN1 carries the D1 pumping tone to the GRIN1 AOM (gate
    ttl_GRIN1_switch), dds_GRIN2 the excitation tone to the GRIN2 AOM, and the
    D1 light comes from an external source gated by ttl_D1_pumping_switch.
    Both nodes deliver the pumping repump through fiber AOMs A5/A6 only:
    Node1's extra A1-A4 repump (commented out below) was dropped after the
    user's 2026-10-21 test on Node1 showed A5/A6 alone pump fine.

    Named record_CW_optical_pumping_node2, recording "CW_optical_pumping_node2",
    until 2026-10-23. The old Node1 versions, the non-DMA CW_optical_pumping_*
    and the chopped and without-coil variants are in
    experiment_functions_obsolete.py.
    """
    with self.core_dma.record("CW_optical_pumping"):
        ### Turn off unused light
        with parallel:
            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            self.ttl_exc0_switch.on()
            self.ttl_D1_pumping_switch.on()

        ### Configure pumping-repump delivery AOMs
        ### These must remain sequential if they share the same Urukul SPI bus
        # if self.which_node == "alice":
        #     ### Node1 also delivers the pumping repump through A1-A4
        #     self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=dB_to_V(-5.0))
        #     self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=dB_to_V(-5.0))
        #     self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=dB_to_V(-5.0))
        #     self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

        # if self.which_node == "alice":
        #     with parallel:
        #         self.dds_AOM_A1.sw.on()
        #         self.dds_AOM_A2.sw.on()
        #         self.dds_AOM_A3.sw.on()
        #         self.dds_AOM_A4.sw.on()
        # delay(5 * us)
        with parallel:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        ### Program the coil voltages forward in time
        self.zotino0.write_dac(self.coil_channels[0], self.AZ_bottom_volts_OP)
        self.zotino0.write_dac(self.coil_channels[1], -self.AZ_bottom_volts_OP)
        self.zotino0.write_dac(self.coil_channels[2], self.AX_volts_OP)
        self.zotino0.write_dac(self.coil_channels[3], self.AY_volts_OP)
        delay(1.5 * us)
        self.zotino0.load()

        ### Physical coil-relaxation time
        # delay(0.395 * ms)
        delay(self.t_OP_coil_relax_time)

        ### Configure GRIN1 DDS
        self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        self.dds_GRIN1.sw.on()
        delay(2*us)

        ### Optical pumping
        with parallel:
            self.ttl_pumping_repump_switch.off()
            self.ttl_D1_pumping_switch.off()
            self.ttl_GRIN1_switch.off()

        delay(self.t_pumping)

        ### Turn off D1 pumping light
        with parallel:
            self.ttl_D1_pumping_switch.on()
            self.ttl_GRIN1_switch.on()

        delay(self.t_depumping)

        ### Turn off all pumping RF/light paths
        with parallel:
            self.ttl_pumping_repump_switch.on()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            self.dds_GRIN1.sw.off()
        # if self.which_node == "alice":
        #     delay(5 * us)
        #     with parallel:
        #         self.dds_AOM_A1.sw.off()
        #         self.dds_AOM_A2.sw.off()
        #         self.dds_AOM_A3.sw.off()
        #         self.dds_AOM_A4.sw.off()

        ### Restore feedback amplitudes while RF switches are off
        # if self.which_node == "alice":
            # self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=self.stabilizer_AOM_A1.amplitude)
            # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=self.stabilizer_AOM_A2.amplitude)
            # self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)
            # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=self.stabilizer_AOM_A4.amplitude)
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        # delay(1.7 * us)  ### to match the timing with node1

@kernel
def record_recooling(self):
    with self.core_dma.record("recooling"):

        ### Set PGC coil voltages using forward scheduling
        self.zotino0.write_dac(self.coil_channels[0], self.AZ_bottom_volts_PGC)
        self.zotino0.write_dac(self.coil_channels[1], -self.AZ_bottom_volts_PGC)
        self.zotino0.write_dac(self.coil_channels[2], self.AX_volts_PGC)
        self.zotino0.write_dac(self.coil_channels[3], self.AY_volts_PGC)
        delay(1.5 * us)
        self.zotino0.load()

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        delay(0.4 * ms)  ### coil relaxation time

        self.dds_cooling_DP.sw.on()
        self.ttl_repump_switch.off()

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()

        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
        else:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(self.t_recooling)

        self.dds_cooling_DP.sw.off()
        self.ttl_repump_switch.on()

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

@kernel
def record_excitation_and_photon_collection(self):
    collection_mu = self.core.seconds_to_mu(self.t_photon_collection_time)
    excitation_pulse_mu = self.core.seconds_to_mu(self.t_excitation_pulse)

    with self.core_dma.record("excitation_and_photon_collection"):
        self.ttl_exc0_switch.off()
        delay(5 * us)

        ### t1 is located here during playback
        with parallel:
            with sequential:
                self.dds_FORT.sw.off()
                delay_mu(50 + collection_mu)
                self.dds_FORT.sw.on()

            with sequential:
                delay_mu(int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()
                delay_mu(excitation_pulse_mu)
                self.ttl_GRIN2_switch.on()
                delay_mu(1000)
                self.ttl_exc0_switch.on()

            with sequential:
                delay_mu(int(self.gate_start_offset_mu))
                with parallel:
                    self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

@kernel
def record_parity_mapping(self):
    with self.core_dma.record("parity_mapping"):

        ###todo: if FORT power drifts a lot involving feedcback, FORT power setting should be outside dma
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2*us)

        t_mapping = now_mu()
        self.ttl_microwave_switch.off()

        at_mu(t_mapping + int(self.t_microwave_11_pulse / ns))
        self.ttl_microwave_switch.on()

        if self.t_MW_RF_pulse > 0:
            at_mu(t_mapping + 500 + int(self.t_microwave_11_pulse / ns))
            self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
            self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

            # at_mu(t_mapping + 1500 + int(self.t_microwave_11_pulse / ns)
            delay(2 * us)
            with parallel:
                self.ttl_microwave_switch.off()
                self.dds_MW_RF.sw.on()
            delay(self.t_MW_RF_pulse)
            # at_mu(t_mapping + 1500 + int(self.t_microwave_11_pulse / ns) + int(self.t_MW_RF_pulse / ns))
            with parallel:
                self.ttl_microwave_switch.on()
                self.dds_MW_RF.sw.off()

        delay(2*us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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
    avgs = 5

    # self.dds_FORT.sw.off() ### Why turnning off FORT?
    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_cooling_DP.sw.off()

    # todo: D1 feedback
    # self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))

    # UNBLOCKING GRIN1
    self.ttl_GRIN1_switch.off()
    # self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(0.0))
    self.dds_GRIN1.sw.on()

    #turning D1 ON. Since 2026-10-19 Node1's D1 AOM runs from an external RF
    # source switched by ttl_D1_pumping_switch (active low: off() = D1 on),
    # not by the Urukul channel dds_GRIN2.
    self.ttl_D1_pumping_switch.off()

    delay(0.1*ms)

    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        measurement += measurement_buf[self.GRIN1_sampler_ch]
        delay(0.1*ms)
    measurement /= avgs
    self.append_to_dataset("GRIN1_D1_monitor", measurement)

    # turning D1 off & repump on
    # self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
    self.ttl_D1_pumping_switch.on()  # D1 off
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
    self.dds_GRIN1.sw.off()
    self.ttl_GRIN1_switch.on()
    self.ttl_exc0_switch.on()

    delay(0.1 * ms)

@kernel
def measure_REPUMP(self):
    """
    This is for monitring MOT REPUMP power

    Used in end_measurement

    It uses feedback_channels.json to find which which sampler channel is used for AOM 5 and 6.
    This is defined in BaseExperiment.

    """
    measurement_buf = np.array([0.0]*8)
    measurement1 = 0.0 # RP on AOM5
    measurement2 = 0.0 # RP on AOM6

    avgs = 10

    self.ttl_repump_switch.off()  # turns on RP AOM
    self.ttl_pumping_repump_switch.on() ### turn off pumping RP
    self.dds_cooling_DP.sw.off() ## turn off cooling

    delay(0.1 * ms)

    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(0.1 * ms)

    for i in range(avgs):
        self.sampler0.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[self.RP_AOM5_ch] # AOM5
        delay(0.1 * ms)
        measurement2 += measurement_buf[self.RP_AOM6_ch] # AOM6

    measurement1 /= avgs
    measurement2 /= avgs

    self.append_to_dataset("REPUMP1_monitor", measurement1)
    self.append_to_dataset("REPUMP2_monitor", measurement2)

    self.ttl_repump_switch.on()  # turns the RP AOM off
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

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

    self.ttl_pumping_repump_switch.off() # turns on PR AOM

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

    self.ttl_pumping_repump_switch.on()  # turns the PR AOM off

    # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)


    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

    delay(0.1 * ms)

@kernel
def Sampler_test(self):
    '''
    I have conencted 1V signal from the linear DC power supply to different sampler channels.
    I am going to see if the measurement of these channels fluctuate.

    In the first block, I save only a single channel from each Sampler. For this, you need to change the dataset in BaseExperiment
    to self.experiment.set_dataset("Sampler0_test", [0.0], broadcast=True), and the same for Sampler1&2.

    In the 2nd block below, I am saving all the Sampler channels in the datasets to monitor all channels.
    Since most of the lasers, if not all, are off during this measurement, the Sampler channels (except those connected to DC
    power supply for monitoring) should show zero volts.
    For this block, you should use self.experiment.set_dataset("Sampler0_test", [np.zeros(8, dtype=float)], broadcast=True)
    and the same for Sampler1&2.

    '''
    #### --------------- Use only one of the two following blocks ------------------

    # ####     saving only one channel from each Sampler ##################
    # Sampler0_ch = 2 # Sampler 0 - ch2
    # Sampler1_ch = 5 # Sampler 1 - ch5
    # Sampler2_ch = 0 # Sampler 2 - ch0
    #
    # Sampler0_measurement_buf = np.array([0.0]*8)
    # Sampler1_measurement_buf = np.array([0.0]*8)
    # Sampler2_measurement_buf = np.array([0.0]*8)
    #
    # Sampler0_measurement = 0.0
    # Sampler1_measurement = 0.0
    # Sampler2_measurement = 0.0
    # avgs = 5
    #
    # delay(50 * us)
    # for i in range(avgs):
    #     self.sampler0.sample(Sampler0_measurement_buf)
    #     delay(20 * us)
    #     self.sampler1.sample(Sampler1_measurement_buf)
    #     delay(20 * us)
    #     self.sampler2.sample(Sampler2_measurement_buf)
    #
    #     delay(20 * us)
    #     Sampler0_measurement += Sampler0_measurement_buf[Sampler0_ch]
    #     Sampler1_measurement += Sampler1_measurement_buf[Sampler1_ch]
    #     Sampler2_measurement += Sampler2_measurement_buf[Sampler2_ch]
    #     delay(20 * us)
    #
    # Sampler0_measurement /= avgs
    # Sampler1_measurement /= avgs
    # Sampler2_measurement /= avgs
    #
    # self.append_to_dataset("Sampler0_test", Sampler0_measurement)
    # self.append_to_dataset("Sampler1_test", Sampler1_measurement)
    # self.append_to_dataset("Sampler2_test", Sampler2_measurement)



    ####     saving all the channels from the Samplers ##################

    ### set the cooling DP AOM to the MOT settings to monitor the MOT powers in MOT 1,2,5, and 6.
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
    delay(1 * ms)

    self.dds_cooling_DP.sw.on()  ### turn on cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP
    delay(0.1 * ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.off()
    self.dds_AOM_A4.sw.off()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    delay(0.1 * ms)

    ### Set the coils to 1V. So Sampler_test also tests coil driver output in combination with Zotino.
    # self.dds_FORT.sw.on()

    self.zotino0.set_dac(
        [1.0, 1.0, 1.0, 1.0], channels=self.coil_channels)
    delay(1.0 * ms)

    Sampler0_measurement_buf = np.array([0.0] * 8)
    Sampler1_measurement_buf = np.array([0.0] * 8)
    Sampler2_measurement_buf = np.array([0.0] * 8)


    delay(50 * us)

    self.sampler0.sample(Sampler0_measurement_buf)
    delay(20 * us)
    self.sampler1.sample(Sampler1_measurement_buf)
    delay(20 * us)
    self.sampler2.sample(Sampler2_measurement_buf)
    delay(20 * us)

    self.append_to_dataset("Sampler0_test", Sampler0_measurement_buf)
    self.append_to_dataset("Sampler1_test", Sampler1_measurement_buf)
    self.append_to_dataset("Sampler2_test", Sampler2_measurement_buf)

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
    avgs = 1

    ### Re-arm slack before the first output. end_measurement calls this
    ### right after its set_dataset(measurements_progress) RPC, which is
    ### SYNCHRONOUS: the kernel waits for the host while the timeline stands
    ### still, and when that wait outran the delays after it, the coil write
    ### below underflowed -- RID 38968, -550 us on zotino0, after six good
    ### readings. break_realtime only ever moves the cursor forward, and
    ### nothing here is timed against an earlier event.
    self.core.break_realtime()

    #####################################  Measure with Zotino set to zero V
    ### Turn off all the coils
    self.zotino0.set_dac(
        [0.0, 0.0, 0.0, 0.0],
        channels=self.coil_channels)
    delay(2 * ms)

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
    self.append_to_dataset("Magnetometer_Zero_X", MagnetometerZ * 350) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_Zero_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_Zero_Z", MagnetometerY * 350)
    delay(0.1*ms)

    #####################################  Measure in the OP phase
    ### Set the coils to OP setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(2 * ms)

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
    self.append_to_dataset("Magnetometer_OP_X", MagnetometerZ * 350)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_OP_Y", MagnetometerX * 350)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_OP_Z", MagnetometerY * 350)

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)
    delay(2 * ms)

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
    self.append_to_dataset("Magnetometer_MOT_X", MagnetometerZ * 350)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_MOT_Y", MagnetometerX * 350)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_MOT_Z", MagnetometerY * 350)
    delay(1*ms)

@kernel
def test_magnetometer_experiment(self):
    self.core.reset()
    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(100*ms)
        measure_Magnetometer_Node2(self)
        # measure_Magnetometer_Node2_dummy_step(self)
        delay(100*ms)

        self.measurement += 1
        self.set_dataset(self.measurements_progress, 100 * self.measurement / self.n_measurements, broadcast=True)

@kernel
def measure_Magnetometer_Node2(self):
    ### x,y, and z axes are connected to Sampler2 Ch1,2, and 3, respectively.
    avgs = 1
    mGauss_per_Volt = 1 #625

    #####################################  Measure with Zotino set to zero V
    ### Turn off all the coils
    self.zotino0.set_dac(
        [0.0, 0.0, 0.0, 0.0],
        channels=self.coil_channels)
    delay(200 * ms)

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
    self.append_to_dataset("Magnetometer_Zero_X", MagnetometerX * mGauss_per_Volt) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_Zero_Y", MagnetometerY * mGauss_per_Volt) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_Zero_Z", MagnetometerZ * mGauss_per_Volt)
    delay(0.1*ms)

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
    self.append_to_dataset("Magnetometer_OP_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_OP_Y", MagnetometerY * mGauss_per_Volt)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_OP_Z", MagnetometerZ * mGauss_per_Volt)

    # #### I added Mag690 magnetometer (borrowed from Josiah) temporary to test our magnetometer.
    # #### It is connected to Sampler1_ch0 to 2. To be deleted soon.
    # measurement_buf = np.array([0.0] * 8)
    # MagnetometerX = 0.0
    # MagnetometerY = 0.0
    # MagnetometerZ = 0.0
    #
    # for i in range(avgs):
    #     self.sampler1.sample(measurement_buf)
    #     MagnetometerX += measurement_buf[0]
    #     MagnetometerY += measurement_buf[1]
    #     MagnetometerZ += measurement_buf[2]
    #     delay(0.1 * ms)
    #
    # MagnetometerX /= avgs
    # MagnetometerY /= avgs
    # MagnetometerZ /= avgs
    # self.append_to_dataset("Magnetometer_Mag690_Zero_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 100 mG
    # self.append_to_dataset("Magnetometer_Mag690_Zero_Y", MagnetometerY * mGauss_per_Volt)  ###
    # self.append_to_dataset("Magnetometer_Mag690_Zero_Z", MagnetometerZ * mGauss_per_Volt)
    # delay(0.1 * ms)

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)
    delay(10 * ms)

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
    self.append_to_dataset("Magnetometer_MOT_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_MOT_Y", MagnetometerY * mGauss_per_Volt)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_MOT_Z", MagnetometerZ * mGauss_per_Volt)
    delay(1*ms)

@kernel
def measure_Magnetometer_Node2_X_Y_Z(self):
    ### x,y, and z axes are connected to Sampler2 Ch1,2, and 3, respectively.
    ### to test each coils separately
    ## Zero: 3V @ AZ_bot
    ## OP: 3V @ AZ_top
    ## MOT: 3V @ AX
    ## TEST: 3V @ AY


    avgs = 1
    mGauss_per_Volt = 1 #625

    ### Re-arm slack before the first output. end_measurement calls this
    ### right after its set_dataset(measurements_progress) RPC, which is
    ### SYNCHRONOUS: the kernel waits for the host while the timeline stands
    ### still, and when that wait outran the delays after it, the coil write
    ### below underflowed -- RID 38968, -550 us on zotino0, after six good
    ### readings. break_realtime only ever moves the cursor forward, and
    ### nothing here is timed against an earlier event.
    self.core.break_realtime()

    #####################################  Measure with Zotino set to 3V
    self.zotino0.set_dac(
        [3.0, 0.0, 0.0, 0.0],
        channels=self.coil_channels)
    delay(200 * ms)

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
    self.append_to_dataset("Magnetometer_Zero_X", MagnetometerX * mGauss_per_Volt) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_Zero_Y", MagnetometerY * mGauss_per_Volt) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_Zero_Z", MagnetometerZ * mGauss_per_Volt)
    delay(0.1*ms)

    #####################################  Measure in the OP phase
    ### Set the coils to OP setting
    self.zotino0.set_dac(
        [0.0, 3.0, 0.0, 0.0],
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
    self.append_to_dataset("Magnetometer_OP_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_OP_Y", MagnetometerY * mGauss_per_Volt)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_OP_Z", MagnetometerZ * mGauss_per_Volt)

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [0.0, 0.0, 3.0, 0.0],
        channels=self.coil_channels)
    delay(10 * ms)

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
    self.append_to_dataset("Magnetometer_MOT_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_MOT_Y", MagnetometerY * mGauss_per_Volt)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_MOT_Z", MagnetometerZ * mGauss_per_Volt)
    delay(1*ms)

    #####################################
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [0.0, 0.0, 0.0, 3.0],
        channels=self.coil_channels)
    delay(10 * ms)

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
    self.append_to_dataset("Magnetometer_TEST_X", MagnetometerX * mGauss_per_Volt)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_TEST_Y", MagnetometerY * mGauss_per_Volt)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_TEST_Z", MagnetometerZ * mGauss_per_Volt)
    delay(1*ms)

@kernel
def end_measurement(self):
    """
    End the measurement by setting datasets and deciding whether to increment the measuement index
    :param self:
    :return measurement: TInt32, the measurement index
    """
    in_health_check = self.in_health_check

    ### update the datasets
    self.set_dataset(self.measurements_progress, 100 * self.measurement / self.n_measurements, broadcast=True)

    self.append_to_dataset('SPCM0_RO1_current_iteration', self.SPCM0_RO1)
    self.append_to_dataset('SPCM1_RO1_current_iteration', self.SPCM1_RO1)
    self.append_to_dataset('SPCM0_OtherNode_RO1_current_iteration', self.SPCM0_OtherNode_RO1)
    self.append_to_dataset('SPCM1_OtherNode_RO1_current_iteration', self.SPCM1_OtherNode_RO1)
    delay(1 * ms)
    self.append_to_dataset('SPCM0_RO2_current_iteration', self.SPCM0_RO2)
    self.append_to_dataset('SPCM1_RO2_current_iteration', self.SPCM1_RO2)
    self.append_to_dataset('SPCM0_OtherNode_RO2_current_iteration', self.SPCM0_OtherNode_RO2)
    self.append_to_dataset('SPCM1_OtherNode_RO2_current_iteration', self.SPCM1_OtherNode_RO2)
    delay(1 * ms)
    self.append_to_dataset('AllSPCMs_RO1_current_iteration', self.AllSPCMs_RO1)
    self.append_to_dataset('AllSPCMs_RO2_current_iteration', self.AllSPCMs_RO2)
    delay(1 * ms)
    ### Alternating RO
    self.append_to_dataset('AllSPCMs_alternating_RO_alice_current_iteration', self.AllSPCMs_alternating_RO_alice)
    self.append_to_dataset('AllSPCMs_alternating_RO_bob_current_iteration', self.AllSPCMs_alternating_RO_bob)
    delay(1*ms)


    self.SPCM0_RO1_list[self.measurement] = self.SPCM0_RO1
    self.SPCM1_RO1_list[self.measurement] = self.SPCM1_RO1
    self.SPCM0_RO2_list[self.measurement] = self.SPCM0_RO2
    self.SPCM1_RO2_list[self.measurement] = self.SPCM1_RO2

    self.SPCM0_OtherNode_RO1_list[self.measurement] = self.SPCM0_OtherNode_RO1
    self.SPCM1_OtherNode_RO1_list[self.measurement] = self.SPCM1_OtherNode_RO1
    self.SPCM0_OtherNode_RO2_list[self.measurement] = self.SPCM0_OtherNode_RO2
    self.SPCM1_OtherNode_RO2_list[self.measurement] = self.SPCM1_OtherNode_RO2
    delay(1 * ms)
    self.AllSPCMs_RO1_list[self.measurement] = self.AllSPCMs_RO1
    self.AllSPCMs_RO2_list[self.measurement] = self.AllSPCMs_RO2
    self.atom_loading_time_list[self.measurement] = self.atom_loading_time

    ### now done at the beginning of the experiment for FORT POL stabilization
    # delay(1*ms)
    # measure_FORT_MM_fiber(self)

    # if self.which_node == "alice":
    #     # delay(1 * ms)
    #     # measure_GRIN1(self)
    #     # delay(1 * ms)
    #     # measure_PUMPING_REPUMP(self)
    #     # delay(1 * ms)
    #     # measure_Magnetometer(self)
    #     # delay(1*ms)
    #     # Sampler_test(self)
    #     # delay(1*ms)
    #     # measure_MOT_end(self)
    #     # delay(1*ms)
    #     # measure_REPUMP(self)
    # # elif self.which_node == "bob":
    # #     measure_Magnetometer_Node2(self)

    if self.monitor_magnetometer_in_end_measurement:
        if self.which_node == "alice":
            measure_Magnetometer(self)
        else:
            # measure_Magnetometer_Node2(self)
            measure_Magnetometer_Node2_X_Y_Z(self)

    delay(1*ms)

    self.append_to_dataset('AllSPCMs_alternating_RO_alice', self.AllSPCMs_alternating_RO_alice)
    self.append_to_dataset('AllSPCMs_alternating_RO_bob', self.AllSPCMs_alternating_RO_bob)
    delay(1*ms)

    self.measurement += 1
    delay(1 * ms)
    if not in_health_check:  ## advance and in_health_check are different type so can't be mixed.
        self.append_to_dataset('SPCM0_RO1', self.SPCM0_RO1)
        self.append_to_dataset('SPCM1_RO1', self.SPCM1_RO1)
        self.append_to_dataset('SPCM0_RO2', self.SPCM0_RO2)
        self.append_to_dataset('SPCM1_RO2', self.SPCM1_RO2)
        delay(1 * ms)
        self.append_to_dataset('SPCM0_OtherNode_RO1', self.SPCM0_OtherNode_RO1)
        self.append_to_dataset('SPCM1_OtherNode_RO1', self.SPCM1_OtherNode_RO1)
        self.append_to_dataset('SPCM0_OtherNode_RO2', self.SPCM0_OtherNode_RO2)
        self.append_to_dataset('SPCM1_OtherNode_RO2', self.SPCM1_OtherNode_RO2)

        self.append_to_dataset('AllSPCMs_RO1', self.AllSPCMs_RO1)
        self.append_to_dataset('AllSPCMs_RO2', self.AllSPCMs_RO2)
        delay(1 * ms)
    else:
        self.append_to_dataset('SPCM0_RO1_in_health_check', self.SPCM0_RO1)
        self.append_to_dataset('SPCM1_RO1_in_health_check', self.SPCM1_RO1)
        self.append_to_dataset('SPCM0_RO2_in_health_check', self.SPCM0_RO2)
        self.append_to_dataset('SPCM1_RO2_in_health_check', self.SPCM1_RO2)
        delay(1 * ms)
        self.append_to_dataset('SPCM0_OtherNode_RO1_in_health_check', self.SPCM0_OtherNode_RO1)
        self.append_to_dataset('SPCM1_OtherNode_RO1_in_health_check', self.SPCM1_OtherNode_RO1)
        self.append_to_dataset('SPCM0_OtherNode_RO2_in_health_check', self.SPCM0_OtherNode_RO2)
        self.append_to_dataset('SPCM1_OtherNode_RO2_in_health_check', self.SPCM1_OtherNode_RO2)

        self.append_to_dataset('AllSPCMs_RO1_in_health_check', self.AllSPCMs_RO1)
        self.append_to_dataset('AllSPCMs_RO2_in_health_check', self.AllSPCMs_RO2)

    ### Re-arm the timeline for the next measurement.
    ###
    ### This function spends real time on 25 host RPCs: 24 append_to_dataset
    ### calls, which are async, plus the measurements_progress set_dataset,
    ### which is NOT -- it returns a value, so it blocks for a full round
    ### trip, every measurement -- while advancing the timeline by only ~10 ms
    ### of delay(). When the host is busy (every broadcast dataset also goes
    ### out to the dashboard and to every subscribed applet) the real time
    ### spent exceeds the timeline advanced, and the loop wraps with the
    ### cursor BEHIND the wall clock. The next RTIO event then underflows,
    ### inside whatever the caller does first rather than here, which is why
    ### this reads as a bug in the atom loading routine:
    ###
    ###   RTIOUnderflow ... channel 36, slack -6144984 mu
    ###     load_MOT_and_FORT_until_atom_recycle, at the dds_FORT.set that
    ###     lowers the FORT to the science setpoint
    ###
    ### -6.14 ms there, after the delay(10 * ms) at the top of
    ### atom_loading_2_experiment's measurement loop, means this function
    ### returned about 16 ms behind.
    ###
    ### Placed here rather than at the 47 call sites. core.break_realtime()
    ### only ever moves the cursor FORWARD -- "if the time cursor is already
    ### after that position, this function does nothing" -- so on a healthy
    ### run, standalone included, this is a no-op, and it only acts in the
    ### case that would otherwise have raised. It also makes the dozen
    ### hand-placed `delay(5 * ms) ### hopefully to avoid underflow.` lines
    ### after other end_measurement calls redundant; they are left alone
    ### because they now serve as plain inter-measurement spacing.
    self.core.break_realtime()

@kernel
def FORT_ramp1_smoothstep(self, direction="down"):
    """
    #todo: this could be written as dma. But, needs to be re-recorded after laser feedback

    For ramping FORT from loading setpoint to science and vice versa.
    Smoothly ramp FORT power using a quintic smoothstep profile. If t_FORT_ramp is too short (<1ms), it uses less
    number of steps to avoid Underflow errors. This can handle any t_FORT_ramp, from 1us to 10ms, for example.

    direction: "down" or "up"
    """
    #### With t_FORT_ramp as a control variable to scan
    # assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"
    #
    # p_high = self.stabilizer_FORT.amplitudes[0]
    # p_low = self.stabilizer_FORT.amplitudes[1]
    # n_steps_max = 2000
    # step_delay_min = 10 * us
    #
    # ### Choose step count so delay >= step_delay_min, but not more than n_steps_max
    # n_steps = int(self.t_FORT_ramp / step_delay_min)
    # if n_steps > n_steps_max:
    #     n_steps = n_steps_max
    # elif n_steps < 1:
    #     n_steps = 1  # safety in extreme case
    #
    # step_delay = self.t_FORT_ramp / n_steps
    #
    # for i in range(n_steps):
    #     x = i / (n_steps - 1) if n_steps > 1 else 1.0
    #     smoothstep = 6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3
    #
    #     if direction == "down":
    #         p_FORT = p_high - smoothstep * (p_high - p_low)
    #     else:
    #         p_FORT = p_low + smoothstep * (p_high - p_low)
    #
    #     delay(step_delay)
    #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=p_FORT)



    #### With t_FORT_step as a control variable to scan
    assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"

    p_high = self.stabilizer_FORT.amplitudes[0]
    p_low = self.stabilizer_FORT.amplitudes[1]
    t_ramp = 10 * ms
    step_delay_min = 10 * us

    if self.t_FORT_step < step_delay_min:
        self.t_FORT_step = step_delay_min

    n_steps = int(t_ramp / self.t_FORT_step)
    if n_steps < 1:
        n_steps = 1  # safety in extreme case

    for i in range(n_steps):
        x = i / (n_steps - 1) if n_steps > 1 else 1.0
        smoothstep = 6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3

        if direction == "down":
            p_FORT = p_high - smoothstep * (p_high - p_low)
        else:
            p_FORT = p_low + smoothstep * (p_high - p_low)

        delay(self.t_FORT_step)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=p_FORT)

@kernel
def FORT_ramp2_smoothstep(self, direction="down"):
    """
    #todo: this could be written as dma. But, needs to be re-recorded after laser feedback

    For ramping FORT from science setpoint to holding (microwave) and vice versa.
    Smoothly ramp FORT power using a quintic smoothstep profile. If t_FORT_ramp is too short (<1ms), it uses less
    number of steps to avoid Underflow errors. This can handle any t_FORT_ramp, from 1us to 10ms, for example.

    direction: "down" or "up"
    """

    assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"

    p_high = self.stabilizer_FORT.amplitudes[1]
    # p_low = self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1]
    p_low = self.stabilizer_FORT.amplitudes[2]
    n_steps_max = 2000
    step_delay_min = 10 * us

    ### Choose step count so delay >= step_delay_min, but not more than n_steps_max
    n_steps = int(self.t_FORT_ramp / step_delay_min)
    if n_steps > n_steps_max:
        n_steps = n_steps_max
    elif n_steps < 1:
        n_steps = 1  # safety in extreme case

    step_delay = self.t_FORT_ramp / n_steps

    for i in range(n_steps):
        x = i / (n_steps - 1) if n_steps > 1 else 1.0
        smoothstep = 6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3

        if direction == "down":
            p_FORT = p_high - smoothstep * (p_high - p_low)
        else:
            p_FORT = p_low + smoothstep * (p_high - p_low)

        delay(step_delay)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=p_FORT)

###############################################################################
# 2. EXPERIMENT FUNCTIONS
# These are the experiments we run, and the name of each should end with
# experiment in order to have it show up in GeneralVariableScans
###############################################################################
@kernel
def test_ttl_pulse_experiment(self):
    """
    testing TTL on scope
                                             Node1              Node2
    ttl_node1_ output1/input1:            ttl 12(out)         ttl 8 (in)
    ttl_node2_ input1/output1:            ttl 11(in)          ttl 15 (out)
    ttl_node1_ output2/input2:            ttl 15(out)         ttl 9 (in)

    """
    self.core.reset()
    self.measurement = 0

    if self.which_node == 'alice':
        self.ttl11.input()
        self.ttl12.output()
        self.ttl15.output()
    else:
        self.ttl8.input()
        self.ttl15.output()
        self.ttl9.input()

    while self.measurement < self.n_measurements:
        delay(1 * s)
        # with parallel:
        #     self.ttl8.sample_input()
        #     self.ttl9.sample_input()
        #
        # with parallel:
        #     readout1 = int(self.ttl8.sample_get())
        #     readout2 = int(self.ttl9.sample_get())

        print("-----n_measurements: ", self.measurement ,"-------")
        print("pulsed ")

        # self.ttl15.on()
        # self.ttl8.sample_input()
        # readout1 = int(self.ttl8.sample_get())
        # delay(2*s)
        # self.ttl15.off()
        # self.ttl8.sample_input()
        # readout1 = int(self.ttl8.sample_get())
        # delay(5*s)
        # delay(6 * s)

        # print("Reading out ttl12: ", readout1)
        # print("Reading out ttl15: ", readout2)
        delay(5 * s)

        self.measurement += 1

        # delay(100*ns)

        # self.ttl8.sample_input()
        # readout = int(self.ttl8.sample_get())  ## this is 1 when the laser is locked, it is 0 otherwise.
        # delay(1*s)
        # print('ttl15 on, readout with ttl8 = ', readout)


        # # delay(5*s)
        # delay(80*ns)
        # # delay(100*ns)
        #
        # self.ttl8.sample_input()
        #  = int(self.ttl8.sample_get())  ## this is 1 when the laser is locked, it is 0 otherwise.
        # delay(1*s)
        # print('ttl15 off, readout with ttl8 = ', readout)

        # delay(1 * s)

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

    delay(10 * ms)



    ### Akbar 2025-11-25
    ### record and get handle
    record_chopped_RO(self)
    delay(10 * ms)
    self.core.break_realtime()



    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
            ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            delay(0.1 * ms)
            run_feedback_and_record_FORT_MM_power(self)

        load_MOT_and_FORT(self)
        delay(0.1*ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

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
    Simple atom loading experiment.

    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.tune_852_waveplates_to_target_in_experiment:
        delay(100 * ms)
        move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
        move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)
        delay(10 * ms)

        self.core.reset()

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # # delay(0.1 * ms)
        # # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
        #
        # # load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        #
        # # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # # delay(0.1 * ms)
        # # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
        #
        # delay(1*ms)
        # first_shot(self)
        # delay(1 * ms)
        #
        # if self.t_recooling_after_first_shot > 0:
        #     recooling_after_first_shot(self)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        # delay(1*ms)
        # two_node_alternating_shot(self)
        # delay(1 * ms)

        delay(self.t_delay_between_shots)
        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def blowaway_fidelity_measurement_experiment(self):
    """
    An experiment to measure the fidelity of blowaway.


    """

    self.core.reset()

    self.require_D1_lock_to_advance = False

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # ############################### pump into F=1
        # ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
        # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        # self.dds_AOM_A5.sw.on()
        # self.dds_AOM_A6.sw.on()
        #
        # self.ttl_pumping_repump_switch.off()
        # # delay(self.t_depumping)
        # delay(1000*us)
        # self.ttl_pumping_repump_switch.on()
        # self.dds_AOM_A5.sw.off()
        # self.dds_AOM_A6.sw.off()
        # delay(100 * us)
        # self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        # self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        # delay(100 * us)
        # #############################

        ############################### pump into F=2
        ### Since the atoms after loading are already in F=2 (mostly), we first pump into F=1 with 100us fixed time.
        ### Then this part shows how long it takes to pump into F=2.
        self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        # delay(self.t_depumping)  # leave the repump on so atoms are left in F=2
        delay(1000 * us)
        self.ttl_repump_switch.on()  # turns the MOT RP AOM off
        delay(0.1 * ms)
        #############################

        ######################## blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

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
        self.zotino0.set_dac([3.5], self.UV_trig_channel)

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
                self.ttl_SPCM0_OtherNode_counter.gate_rising(atom_check_time)
                self.ttl_SPCM1_OtherNode_counter.gate_rising(atom_check_time)

            # SPCM0_atom_check = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_atom_check = self.ttl_SPCM1_counter.fetch_count()
            # SPCM0_OtherNode_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
            # SPCM1_OtherNode_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

            AllSPCMs_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                      self.ttl_SPCM1_counter.fetch_count() + \
                                      self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                      self.ttl_SPCM1_OtherNode_counter.fetch_count())
            try_n += 1

            if AllSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
            self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check)

        delay(1 * ms)
        self.zotino0.set_dac([0.0], self.UV_trig_channel)
        delay(100 * us)

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
        #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        #     ### set the cooling DP AOM to the PGC settings
        #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
        #     self.ttl_repump_switch.off()  ### turn on MOT RP
        #     self.dds_cooling_DP.sw.on()  ### turn on cooling
        #     delay(10 * us)
        #     if self.PGC_and_RO_with_on_chip_beams:
        #             self.dds_AOM_A5.sw.off()
        #             self.dds_AOM_A6.sw.off()
        #     delay(self.t_PGC_after_loading)  ### this is the PGC time
        # ###################################################

        ### saving the atom loading time for each loaded atom.
        self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
        delay(1 * ms)
        ###################################################### end of atom loading attempt


        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
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
    Simple atom loading experiment using load_until_atom_smooth_FORT_recycle. We use this function with additional stages
    compared to atom_loading_2_experiment to scan and optimize different parameters like 852 waveplates.
    For PGC optimization use atom_loading_for_PGC_optimization_experiment below.

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    self.set_dataset('SPCM0_test_RO', [0], broadcast=True)
    self.set_dataset('SPCM1_test_RO', [0], broadcast=True)
    self.set_dataset('SPCM0_OtherNode_test_RO', [0], broadcast=True)
    self.set_dataset('SPCM1_OtherNode_test_RO', [0], broadcast=True)
    self.set_dataset('AllSPCMs_test_RO', [0], broadcast=True)




    ######### to scan 852 waveplates with t_FORT_drop = 10us for example and find max retention (low T)
    delay(1*ms)
    move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
    move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)

    delay(5*ms)
    self.core.reset()

    position_852_HWP = get_rotator_deg(self, name="852_HWP")
    position_852_QWP = get_rotator_deg(self, name="852_QWP")

    delay(5 * ms)
    self.core.reset()
    self.print_async("position_852_HWP: ", position_852_HWP)
    self.print_async("position_852_QWP:", position_852_QWP)
    ###################################################################################################


    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

        # load_MOT_and_FORT(self)
        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)

        delay(1*ms)
        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # ###########  PGC on the trapped atom to optimize coils and cooling_DP_PGC, for example #############
        # ### Set the coils to PGC_optimization setting:
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization, self.AX_volts_PGC_optimization, self.AY_volts_PGC_optimization],
        #     channels=self.coil_channels)
        # delay(1 * ms)
        # ### set the cooling DP AOM to the PGC settings
        # ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
        # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)
        # self.ttl_repump_switch.off()  ### turn on MOT RP
        # self.dds_cooling_DP.sw.on()  ### turn on cooling
        # delay(10 * us)
        # if self.PGC_and_RO_with_on_chip_beams:
        #     self.dds_AOM_A5.sw.off()
        #     self.dds_AOM_A6.sw.off()
        # delay(self.t_PGC_after_loading)  ### this is the PGC time
        # self.ttl_repump_switch.on()  ### turn off MOT RP
        # self.dds_cooling_DP.sw.off()  ### turn off cooling
        # delay(10*us)
        #
        # ### Set the coils to PGC setting:
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC,
        #      self.AY_volts_PGC],
        #     channels=self.coil_channels)
        # delay(0.4 * ms)
        # ###################################################


        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()


        # ################## to see if RO heats atoms
        # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
        # delay(5*us)
        # test_shot(self)
        # delay(10*us)
        # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        # ###########################################


        # ################## to optimize readout settings
        # ### Set the coils to PGC_optimization setting:
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization, self.AX_volts_PGC_optimization,
        #      self.AY_volts_PGC_optimization],
        #     channels=self.coil_channels)
        # delay(1 * ms)
        #
        # ### set the cooling DP AOM to the PGC optimization settings
        # ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
        # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)
        #
        # # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
        # FORT_ramp2_smoothstep(self, direction="down") ### ramping FORT to p_FORT_holding
        #
        # delay(5 * us)
        # test_shot(self)
        # delay(10 * us)
        #
        # # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        # FORT_ramp2_smoothstep(self, direction="up")
        # ###########################################


        # delay(self.t_delay_between_shots)
        second_shot(self)

        self.append_to_dataset('SPCM0_test_RO', self.SPCM0_test_RO)
        self.append_to_dataset('SPCM1_test_RO', self.SPCM1_test_RO)
        self.append_to_dataset('SPCM0_OtherNode_test_RO', self.SPCM0_OtherNode_test_RO)
        self.append_to_dataset('SPCM1_OtherNode_test_RO', self.SPCM1_OtherNode_test_RO)
        self.append_to_dataset('AllSPCMs_test_RO', self.AllSPCMs_test_RO)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)


    # self.dds_FORT.sw.off()

@kernel
def atom_loading_for_PGC_optimization_experiment(self):
    """
    Simple atom loading experiment using load_until_atom_smooth_FORT_recycle and an extra RO or PGC between the shots
    to optimize PGC settings without affecting loading.

    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable

    ######### to scan 852 waveplates with t_FORT_drop = 10us for example and find max retention (low T)
    delay(1 * ms)
    move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
    move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)

    delay(5 * ms)
    self.core.reset()


    if self.enable_laser_feedback:
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

        # load_MOT_and_FORT(self)
        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)

        delay(1*ms)
        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ###########  PGC on the trapped atom to optimize coils and cooling_DP_PGC, for example #############
        ### Set the coils to PGC_optimization setting:
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization, self.AX_volts_PGC_optimization,
             self.AY_volts_PGC_optimization],
            channels=self.coil_channels)
        delay(1 * ms)
        ### set the cooling DP AOM to the PGC optimization settings
        ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)
        self.ttl_repump_switch.off()  ### turn on MOT RP
        self.dds_cooling_DP.sw.on()  ### turn on cooling
        delay(10 * us)
        if self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(self.t_PGC_after_loading)  ### this is the PGC time
        self.ttl_repump_switch.on()  ### turn off MOT RP
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        delay(10 * us)

        ### Set the coils to PGC setting:
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC,
             self.AY_volts_PGC],
            channels=self.coil_channels)
        delay(1 * ms)
        ###################################################

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()


        second_shot(self)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def beam_balancing_with_atoms_experiment(self):
    """
    The idea is to use the loaded atom to balance the MOT beam pairs one by one.
    1- Load an atom
    2- Do optical pumping (not neccesary)
    3- Turn on MOT 1 and 3 for a time t such that retention drops to 50%
    4- Use MOT1 power as is (constant) and change MOT3 power (the opposite beam to MOT1) to maximize retention.

    Repeat this for MOT 2 and 4, and MOT5 and 6.
    """

    self.core.reset()

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    ### Recording costs CPU time, not timeline. Restore the slack that the
    ### delay(100 * ms) after the chopped-pumping recording (archived
    ### 2026-10-23; never played back) used to give the first event.
    self.core.break_realtime()
    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)

        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # ############################# optical pumping phase - pumps atoms into F=1,m_F=0
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        # ### with cw pumping:
        # if self.t_pumping > 0.0:
        #     delay(10 * us)
        #     CW_optical_pumping_node1(self)
        #     delay(10 * us)

        delay(1*ms)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*ms)

        FORT_ramp1_smoothstep(self, direction="down")

        ### Changing power of fiber AOM3
        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=dB_to_V(self.p_excitation))

        if self.t_blowaway > 0:
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_blowaway, amplitude=self.ampl_cooling_DP_MOT)
            delay(10*us)
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            self.ttl_repump_switch.off()  ### turn on MOT RP
            with parallel:
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A3.sw.on()
            delay(self.t_blowaway)
            with parallel:
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A3.sw.off()

        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)


        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

    delay(10 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)

@kernel
def beam_balancing_with_atoms2_experiment(self):
    """
    In beam_balancing_with_atoms_experiment the retention was not sensitive enough to beam imbalance. So, here, I am dropping the
    FORT for 10us before measuring retention to see if imbalance heats up the atoms. With this we are basically
    mimizing the temperature, which we care about. With t_recooling=50ms, I see a nice trend in atom retention vs beam imbalance.

    The idea is to use the loaded atom to balance the MOT beam pairs one by one.
    1- Load an atom, and do RO1.
    2- Turn on MOT 1 and 3 for a time t (say 20ms similar to RO time).
    3- Keep MOT1 power as is (constant) and change MOT3 power (the opposite beam to MOT1) to maximize retention after 10us FORT
    drop.

    Repeat this for MOT 2 and 4, and MOT5 and 6.
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)

        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1*ms)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*ms)


        ### Changing power of fiber AOM3
        self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=dB_to_V(self.p_AOM_A3_optimization))
        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=dB_to_V(self.p_AOM_A3_optimization))

        # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=dB_to_V(self.p_AOM_A2_optimization))
        # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=dB_to_V(self.p_AOM_A2_optimization))

        if self.t_recooling > 0:
            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_blowaway, amplitude=self.ampl_cooling_DP_MOT)
            delay(10*us)
            self.dds_cooling_DP.sw.on()  ### turn on cooling
            self.ttl_repump_switch.off()  ### turn on MOT RP
            with parallel:
                # self.dds_AOM_A1.sw.on()
                # self.dds_AOM_A3.sw.on()

                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A4.sw.on()

            delay(self.t_recooling)
            with parallel:
                # self.dds_AOM_A1.sw.off()
                # self.dds_AOM_A3.sw.off()

                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A4.sw.off()

        self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=self.stabilizer_AOM_A1.amplitude)
        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)

        # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=self.stabilizer_AOM_A2.amplitude)
        # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=self.stabilizer_AOM_A4.amplitude)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(10*us)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

    delay(10 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)

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

    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0

    # delay(1 * ms)
    # move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
    # move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)
    #
    # delay(5 * ms)
    # self.core.reset()

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    # node1_epoch, node1_duration_mu, node1_ptr = CW_OP_node1_handle
    # node2_epoch, node2_duration_mu, node2_ptr = CW_OP_node2_handle
    #
    # self.print_async("Node 1 DMA duration:", node1_duration_mu, "mu")
    # self.print_async("Node 2 DMA duration:", node2_duration_mu, "mu")
    # self.print_async("Node 2 minus Node 1:", node2_duration_mu - node1_duration_mu, "mu")

    self.core.break_realtime()

    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    # delay(1 * ms)
    # move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    # move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    # delay(10 * ms)
    self.core.reset()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        #todo: check node1 pumping repump efficiency;
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # delay(1 * ms)
        # FORT_ramp2_smoothstep(self, direction="down")
        # delay(2 * us)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     # delay(1 * ms)
        #     self.core.break_realtime()

        # ### 60Hz triggering
        # phase_delay_mu = self.core.seconds_to_mu(self.dummy_variable)
        # t_sequence_mu = wait_for_ac_trigger(self, phase_delay_mu)
        # at_mu(t_sequence_mu)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)

        ##todo: this is hardcoded delay to acount for coil drift- set it smaller so that mapping is tuned correclty
        delay(10*us)
        # delay(self.dummy_variable)
        # delay(2*ms)

        ############################
        # microwave phase
        ############################

        if self.t_microwave_pulse > 0.0:
            # if self.which_node == "alice":
            #     #todo: does this really help?
            #     ### Changing the bias field
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            #     delay(1 * ms)

            #todo: does this really help?
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2]) ### lower the FORT power
            delay(2 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            delay(2 * us)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

    # delay(10*ms)
    # self.dds_FORT.sw.off()
    delay(1 * ms)
    self.dds_microwaves.sw.off()
    self.dds_FORT.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def microwave_Ramsey_00_experiment(self):
    """
    Ramsey experiment with two pi/2 MW pulses with a variable time delay in the middle.
    This experiment can be used to measure the T2* of the clock or 01 qubit.
    """

    self.core.reset()

    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0

    # if self.t_pumping > 0.0:
    #     record_chopped_optical_pumping(self)
    #     delay(100*ms)
    #
    # record_chopped_blow_away(self)
    #
    # ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1*ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)

                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)

            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1*ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        ############################
        # microwave phase
        ############################
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave], channels=self.coil_channels)
        # delay(1 * ms)

        # FORT_ramp2_smoothstep(self, direction="down")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])

        ### first pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse/2)
        self.ttl_microwave_switch.on()

        ################## use only one of the following:
        ### without spin echo:
        delay(self.t_delay_between_shots)

        # ### with spin echo. This does not work well. Needs debugging.
        # delay(self.t_delay_between_shots/2)
        # self.ttl_microwave_switch.off()
        # delay(self.t_microwave_00_pulse)
        # self.ttl_microwave_switch.on()
        # delay(self.t_delay_between_shots/2)

        ########################################

        ### second pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_pulse/2)
        self.ttl_microwave_switch.on()

        delay(2*us)
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        # FORT_ramp2_smoothstep(self, direction="up")

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

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

    # record_chopped_blow_away(self)
    #
    # ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)

                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)

            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################ optical pumping phase - pumps atoms into F=1,m_F=0
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1*ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        # ### Changing the bias field
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)
        # delay(1 * ms)

        ############################ microwave 1: transfer mF=0 to mF'=1
        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding*self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        # FORT_ramp2_smoothstep(self, direction="down")
        delay(2 * us)

        self.ttl_microwave_switch.off()
        delay(self.t_microwave_01_pulse)
        self.ttl_microwave_switch.on()
        # delay(2 * us)

        ############################ microwave 2: Ramsey between mF=1 and mF'=1
        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        ### first pi/2 pulse
        self.ttl_microwave_switch.off()
        delay(self.t_microwave_11_pulse/2)
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
        delay(self.t_microwave_11_pulse/2)
        self.ttl_microwave_switch.on()

        delay(2*us)
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2*us)

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

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
def microwave_Rabi_2_and_EXC_experiment(self):
    """
    Microwave and optical pumping and Excitation experiment

    Sequence:
        OP => EXC (~n attempts) => MW00 => BA

    This experiment is used to characterize the excitation efficiency and optimize the power/pulse.

    """

    self.core.reset()

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_until_atom_smooth_FORT_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                first_shot(self)   # starts with FORT science setpoint;

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            first_shot(self)   # starts with FORT science setpoint;

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5 * us)
        self.dds_GRIN2.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.dds_GRIN1.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)
        self.dds_GRIN2.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10*us)

        for excitation_attempt in range(self.n_excitation_attempts):  ##### testing with multiple excitations

            slack = now_mu() - self.core.get_rtio_counter_mu()
            if slack < 1e5:
                # self.print_async("slack added in measurement:", self.measurement)
                self.core.break_realtime()

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(5 * us)

            t1 = now_mu()
            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on 50ns after photon collection ended

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation after t_excitation_offset_mu

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
            self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

        # self.core.break_realtime()
        if self.t_microwave_pulse > 0.0:
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(2 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            delay(2 * us)

        self.core.break_realtime()

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms) ### hopefully to avoid underflow.

    delay(10 * ms)
    self.dds_FORT.sw.off()
    delay(1 * ms)
    self.dds_microwaves.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

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

    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ### Turn on the bias field for OP
        self.zotino0.set_dac([self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP,
                              self.AX_volts_OP, self.AY_volts_OP], channels=self.coil_channels)

        self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(0.1 * ms)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.
        delay(0.5 * ms) ### coils relaxation time



        ### ******************** Use either Method A or B for pumping:
        ############################ The pumping phase - Method A - pumps atoms into F=1 manifold
        ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        self.ttl_pumping_repump_switch.off() ### turning on pumping RP
        delay(3 * ms)  # leave the repump on
        self.ttl_pumping_repump_switch.on() ### turning off pumping RP
        delay(1 * ms)
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        # ############################ The pumping phase - Method B - pumps atoms into F=2 manifold
        # self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        # delay(1 * ms)  # leave the repump on so atoms are left in F=2
        # self.ttl_repump_switch.on()  # turns the MOT RP AOM off
        # delay(1 * ms)

        ### ************************************************************


        ### Change the bias field for microwave pulse
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave], channels=self.coil_channels)
        # delay(0.5 * ms)


        ############################ microwave phase
        if self.t_microwave_pulse > 0.0:
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(5 * us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(10 * us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])


        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

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
    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0

    # SPCM0_SinglePhoton = 0
    # SPCM1_SinglePhoton = 0
    # SPCM0_OtherNode_SinglePhoton = 0
    # SPCM1_OtherNode_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    delay(1 * ms)

    self.core.break_realtime()

    self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs.
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        # ### 60Hz triggering
        # phase_delay_mu = self.core.seconds_to_mu(self.dummy_variable)
        # t_sequence_mu = wait_for_ac_trigger(self, phase_delay_mu)
        # at_mu(t_sequence_mu)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            # delay(5.5 * us) ##todo: this is hardcoded delay to acount for coil drift- set it smaller so that mapping is tuned correctly

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1

        ### Changing the bias field for microwave
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)

        # delay(2*ms)

        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(2 * us)

        # FORT_ramp2_smoothstep(self, direction="down")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        ##### Mapping from |1,0> to |2,1>
        if self.t_microwave_01_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_01_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

        ### The FORT stays at holding through both pulses: an extra up/down to
        ### science between them cost Node1 0.1-0.2 retention (RIDs 39332-39339).
        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        ##### Mapping from |2,1> to |1,1>
        if self.t_microwave_11_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_11_pulse)
            self.ttl_microwave_switch.on()

        delay(5 * us)
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        self.core.break_realtime()
        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        # delay(6 * ms)  ### hopefully to avoid underflow.
        # self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        # self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        # self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton)
        # self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton)
        delay(1 * ms)

    delay(10 * ms)
    # self.dds_FORT.sw.off()
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
    record_CW_optical_pumping(self)
    # record_chopped_optical_pumping(self)
    # delay(100 * ms)
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        ### optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

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
        #     self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
        #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        #     delay(1 * us)
        #
        #     ### Tunring on pumping RP:
        #     self.ttl_pumping_repump_switch.off()
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
        #     self.dds_GRIN2.sw.off()  ### turning off D1 DP
        #     self.ttl_pumping_repump_switch.on()  ### turning off pumping RP
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


        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1
        ### set coils for microwave
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
        #     channels=self.coil_channels)
        # delay(1 * ms)  # coil relaxation time.

        # FORT_ramp2_smoothstep(self, direction="down")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

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
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

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
def microwave_map00_map0m1_experiment(self):
    """
    This experiment maps from F=1,mF=0 to F=2,mF=0, then maps F=2,mF=0 to F=1,mF=-1.

    1- loads an atom using load_MOT_and_FORT_until_atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave pi pulse to transfer population from F=1,mF=0 to F=2,mF=0
    4- Apply another microwave pulse to transfer from F=2,mF=0 to F=1,mF=-1
    5- Blow away F=2 manifold
    6- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0
    # SPCM0_SinglePhoton = 0
    # SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0
    #
    # record_chopped_optical_pumping(self)
    # delay(100 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            # delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        # ### Changing the bias field
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)
        # delay(1 * ms)

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=0
        self.dds_microwaves.set(frequency=self.f_microwaves_00_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        # FORT_ramp2_smoothstep(self, direction="down")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        if self.t_microwave_00_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_00_pulse)
            self.ttl_microwave_switch.on()

        ############################ microwave phase to transfer population from F=2,mF=0 to F=1,mF=-1
        self.dds_microwaves.set(frequency=self.f_microwaves_m10_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(2 * us)

        if self.t_microwave_m10_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_m10_pulse)
            self.ttl_microwave_switch.on()

        # delay(5 * us)
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        # delay(6 * ms)  ### hopefully to avoid underflow.
        # self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        # self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    delay(10 * ms)
    # self.dds_FORT.sw.off()
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_map01_MWRFm11_experiment(self):
    """
    This experiment tests the MW+RF two-photon transition from F=2,mF=1 to F=1,mF=-1. Can be used to find the transition
    resonance or Rabi frequency

    1- loads an atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave pi pulse to transfer population from F=1,mF=0 to F=2,mF=1
    4- Apply microwave + RF pulses to do Rabi oscillation between F=1,mF=-1 and F=2,mF=1
    6- Blow away F=2 manifold
    7- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0
    # SPCM0_SinglePhoton = 0
    # SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1 * ms)
    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    delay(1 * ms)
    self.core.break_realtime()

    self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    delay(1 * ms)

    self.dds_MW_RF.sw.off()  ### turn off RF

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave], channels=self.coil_channels)
        # delay(1 * ms)
        # delay(5*us)

        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        # FORT_ramp2_smoothstep(self, direction="down")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        if self.t_microwave_01_pulse > 0.0:
            # delay(5 * us)
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_01_pulse)
            self.ttl_microwave_switch.on()
            delay(0.5 * us)

        ############################ microwave + RF phase to transfer population from F=2,mF=1 to F=1,mF=-1
        if self.t_MW_RF_pulse>0:
            self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
            self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)
            delay(2 * us)

            with parallel:
                self.ttl_microwave_switch.off() ### turn on MW
                self.dds_MW_RF.sw.on()  ### turn on RF

            delay(self.t_MW_RF_pulse)

            with parallel:
                self.ttl_microwave_switch.on() ### turn off MW
                self.dds_MW_RF.sw.off()  ### turn off RF

        # delay(5*us)
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(5*us)

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        self.core.break_realtime()

        # delay(6 * ms)  ### hopefully to avoid underflow.
        # self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        # self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        # delay(1 * ms)

    self.core.break_realtime()
    # self.dds_FORT.sw.off()
    self.dds_MW_RF.sw.off()  ### turn off RF
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_Ramsey_MWRFm11_experiment(self):
    """
    Ramsey experiment with two pi/2 MW+RF pulses with a variable time delay in the middle.
    This experiment can be used to measure the T2* of F=1,mF=-1 and F=2,mF=1 qubit.

    1- loads an atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave pi pulse to transfer population from F=1,mF=0 to F=2,mF=1
    4- Applies two MW + RF pi/2 pulse with a variable delay in between
    5- Blow away F=2 manifold
    6- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0
    # SPCM0_SinglePhoton = 0
    # SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.
    self.dds_MW_RF.sw.off()  ### turn off RF

    delay(1 * ms)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        # ### Changing the bias field
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)
        # delay(1 * ms)

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1
        self.dds_microwaves.set(frequency=self.f_microwaves_01_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        # FORT_ramp2_smoothstep(self, direction="down")
        delay(2 * us)

        if self.t_microwave_01_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_01_pulse)
            self.ttl_microwave_switch.on()

        ############################################### Ramsey phase
        ####### First MW+RF pi/2 pulse
        if self.t_MW_RF_pulse>0:
            self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
            delay(2 * us)

            with parallel:
                self.ttl_microwave_switch.off() ### turn on MW
                self.dds_MW_RF.sw.on()  ### turn on RF

            delay(self.t_MW_RF_pulse/2)

            with parallel:
                self.ttl_microwave_switch.on() ### turn off MW
                self.dds_MW_RF.sw.off()  ### turn off RF


        delay(self.t_delay_between_shots)

        ####### Second MW+RF pi/2 pulse
        if self.t_MW_RF_pulse > 0:
            with parallel:
                self.ttl_microwave_switch.off()  ### turn on MW
                self.dds_MW_RF.sw.on()  ### turn on RF

            delay(self.t_MW_RF_pulse/2)

            with parallel:
                self.ttl_microwave_switch.on()  ### turn off MW
                self.dds_MW_RF.sw.off()  ### turn off RF

        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        # delay(6 * ms)  ### hopefully to avoid underflow.
        # self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        # self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    self.core.break_realtime()
    # self.dds_FORT.sw.off()
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def microwave_map00_RF01_map00_experiment(self):
    """
    This experiment tries to find the RF transition rate and frequency around 2MHz.

    1- loads an atom
    2- Pumps the atom into F=1,mF=0
    3- Uses microwave pi pulse to transfer population from F=1,mF=0 to F=2,mF=0
    4- Uses RF on-resonance to transfer population from F=2,mF=0 to F=2,mF=1
    5- Uses microwave pi pulse to transfer population from F=2,mF=0 back to F=1,mF=0
    6- Blow away F=2 manifold
    7- Measure retention

    """

    self.core.reset()
    delay(1 * ms)

    ### overwritten below but initialized here so they are always initialized
    # self.SPCM0_RO1 = 0
    # self.SPCM0_RO2 = 0
    # self.SPCM1_RO1 = 0
    # self.SPCM1_RO2 = 0
    # SPCM0_SinglePhoton = 0
    # SPCM1_SinglePhoton = 0

    self.n_feedback_per_iteration = 2
    self.n_atom_loaded_per_iteration = 0

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    # record_chopped_optical_pumping(self)
    # delay(100 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        delay(1 * ms)

    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")

    self.core.break_realtime()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(10 * us)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(10 * us)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        # ### With chopped pumping:
        # if self.t_pumping > 0.0:
        #     chopped_optical_pumping(self)
        #     delay(1 * ms)

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay(10 * us)
            # CW_optical_pumping_node2(self)
            self.core_dma.playback_handle(CW_OP_handle)
            delay(10 * us)

        # ### Changing the bias field
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)
        # delay(1 * ms)

        # FORT_ramp2_smoothstep(self, direction="down")

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=0
        self.dds_microwaves.set(frequency=self.f_microwaves_00_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        if self.t_microwave_00_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_00_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

        ############################ RF phase to transfer population from F=2,mF=0 to F=2,mF=1
        if self.t_MW_RF_pulse>0.0:
            self.dds_MW_RF.sw.on()  ### turn on RF
            delay(self.t_MW_RF_pulse)
            self.dds_MW_RF.sw.off()  ### turn off RF
            delay(2 * us)

        ############################ microwave phase to transfer population from F=2,mF=0 to F=1,mF=0
        if self.t_microwave_00_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_00_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)
        # FORT_ramp2_smoothstep(self, direction="up")


        ############################ blow-away phase - push out atoms in F=2 only
        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        # delay(6 * ms)  ### hopefully to avoid underflow.
        # self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
        # self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
        delay(1 * ms)

    self.core.break_realtime()
    # self.dds_FORT.sw.off()
    delay(1 * ms)
    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)
    self.dds_microwaves.sw.off()

@kernel
def single_photon_experiment_9_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_8_atom_loading_advance

    - optical pumping: record_CW_optical_pumping, with the OP coils (the
      without-coil version this note once named is in
      experiment_functions_obsolete.py)
    - waiting "0.4ms" for OP coil to relax at the beginning of excitation cycles and after atom check
    - inside fast cycles, use 0.4ms relaxation time as before.
    - short_timestamps_length increased to int(self.max_excitation_cycles * 0.5)
    - Drain all additional edges, including when the timestamp buffer is already full
    """

    self.core.reset()
    delay(1 * ms)

    #### recording DMA
    record_CW_optical_pumping(self)
    # record_CW_optical_pumping_node2_without_coil(self)
    record_recooling(self)
    record_excitation_and_photon_collection(self)

    #### getting DMA handles
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")
    # CW_OP_node2_without_coil_handle = self.core_dma.get_handle("CW_optical_pumping_node2_without_coil")
    recooling_dma_handle = self.core_dma.get_handle("recooling")
    excitation_dma_handle = self.core_dma.get_handle("excitation_and_photon_collection")

    delay(1 * ms)
    self.core.break_realtime()
    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    # short_timestamps_length = int(self.max_excitation_cycles * 0.1)
    short_timestamps_length = int(self.max_excitation_cycles * 0.5)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    delay(100*ms)

    while self.measurement < self.n_measurements:
        Any_SPCM_click_counter = 0

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * short_timestamps_length
        SPCM0_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_timestamps = [-1.0]  * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0]  * short_timestamps_length

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # delay(1 * ms)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(5 * us)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5 * us)
        self.dds_GRIN2.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.dds_GRIN1.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)
        self.dds_GRIN2.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        # ### set coils for pumping - this is inside the recorded optical pumping function
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        #     channels=self.coil_channels)
        # # delay(5 * ms)
        # delay(0.4 * ms)

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                # delay(10 * us)
                self.core_dma.playback_handle(CW_OP_handle)
                # self.core_dma.playback_handle(CW_OP_node2_without_coil_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                dma_start_mu = now_mu()
                self.core_dma.playback_handle(excitation_dma_handle)

                t1 = dma_start_mu + self.core.seconds_to_mu(5 * us)
                t_end_SPCM = t1 + int(self.gate_start_offset_mu) + self.core.seconds_to_mu(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM)

                if SPCM0_click_time > 0.0 or SPCM1_click_time > 0.0 or  SPCM0_OtherNode_click_time > 0.0 or SPCM1_OtherNode_click_time > 0.0:
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0.0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)

                        if SPCM1_click_time > 0.0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)

                        if SPCM0_OtherNode_click_time > 0.0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)

                        if SPCM1_OtherNode_click_time > 0.0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)

                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

                    ### Drain all additional edges, including when the timestamp buffer is already full
                    if SPCM0_click_time > 0: self.ttl_SPCM0.count(t_end_SPCM)
                    if SPCM1_click_time > 0: self.ttl_SPCM1.count(t_end_SPCM)
                    if SPCM0_OtherNode_click_time > 0: self.ttl_SPCM0_OtherNode.count(t_end_SPCM)
                    if SPCM1_OtherNode_click_time > 0: self.ttl_SPCM1_OtherNode.count(t_end_SPCM)

                delay(5 * us)

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0  and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.core_dma.playback_handle(recooling_dma_handle)
                delay(5*us)

                # if not (excitation_cycle + 1) % self.atom_check_every_n == 0:
                ## set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4*ms)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                delay(5*us)

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)
                delay(0.4 * ms)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

                ### todo: this is inside recooling block. I am recooling every cycle so do not have set coils here

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)

                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                SPCM0_OtherNode_RO_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
                SPCM1_OtherNode_RO_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

                AllSPCMs_RO_atom_check = int(SPCM0_RO_atom_check + SPCM1_RO_atom_check + SPCM0_OtherNode_RO_atom_check + SPCM1_OtherNode_RO_atom_check)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    break

                delay(10 * us)
                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                # delay(5*ms)
                delay(0.4 * ms)

            delay(5 * us)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        self.core.break_realtime()

        self.dds_GRIN1.sw.off()
        self.dds_GRIN2.sw.off()
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

        # delay(5 * ms)
        self.core.break_realtime()

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        i = 0
        while i < len(tStamps_t1) and tStamps_t1[i] > 0.0:
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])
            i += 1

        ### Number of compact timestamp rows saved during this measurement
        self.append_to_dataset('n_photon_events', Any_SPCM_click_counter)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def atom_photon_parity_11_AllSPCM_experiment(self):
    """
    *** compatible with either node1 or node2
    *** DMA

    Same sequence as atom_photon_parity_10_AllSPCM_experiment, plus complete draining of
    leftover SPCM edges, ported from single_photon_experiment_9_atom_loading_advance_AllSPCM:

    - ttl.timestamp_mu() pops only ONE edge per channel from the RTIO input FIFO. Any extra
      edge in the same gate (double click, afterpulse, dark count) stays buffered and is
      returned by a later attempt's timestamp_mu() as a false early count, whose stale
      timestamp makes at_mu(mapping_click_time + t_start_MW_mapping_mu) underflow unless
      t_start_MW_mapping_mu is kept large.
    - parity_10 drained only after a valid single-detector event (detector >= 0). Here every
      channel that clicked is drained in ALL outcomes, including multi-detector coincidences
      (detector == -1), so no edge survives into the next excitation attempt.
    - The parity readout itself was never contaminated: atom_parity_shot counts on the
      separate ttl_SPCM*_counter channels, not on the timestamp channels drained here.
    - Goal: with false early counts gone, t_start_MW_mapping_mu can be reduced.

    The sequence:
    while self.measurement < self.n_measurements:
        O.P.
        excitation
        MW mapping
        Blow away
        Atom_parity_RO
        measurement += 1
        excitation_cycle += 1

    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    self.core.reset()
    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    #### recording DMA
    record_chopped_blow_away(self)
    record_CW_optical_pumping(self)
    # record_recooling(self)
    record_excitation_and_photon_collection(self)
    record_parity_mapping(self)

    #### getting DMA handles
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_handle = self.core_dma.get_handle("CW_optical_pumping")
    # recooling_dma_handle = self.core_dma.get_handle("recooling")
    excitation_dma_handle = self.core_dma.get_handle("excitation_and_photon_collection")
    parity_mapping_dma_handle = self.core_dma.get_handle("parity_mapping")
    delay(10 * ms)

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    SPCM0_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0
    angle_780_QWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_until_atom_smooth_FORT_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                self.ttl_microwave_switch.on()  ### close the switch
                delay(20*us)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            self.ttl_microwave_switch.on()  ### close the switch
            delay(20*us)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### GRIN1 RF (D1 tone) and GRIN2 RF (excitation tone) stay on for the entire excitation + OP loop;
        ### the GRIN1 and GRIN2 switches turn the D1 and Exc light on and off
        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5 * us)
        self.dds_GRIN2.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.dds_GRIN1.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)
        self.dds_GRIN2.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                # delay(10 * us)
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_handle)
            # t_after_OP_dma = now_mu()

            # if self.which_node == "alice":
            #     ### Changing the bias field.
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            #     delay(1 * ms)

            #
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            # delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                # slack = now_mu() - self.core.get_rtio_counter_mu()
                # if slack < 1e5:
                #     # self.print_async("slack added in measurement:", self.measurement)
                #     self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)

                dma_start_mu = now_mu()
                self.core_dma.playback_handle(excitation_dma_handle)

                t1 = dma_start_mu + self.core.seconds_to_mu(5 * us)
                t_end_SPCM = t1 + int(self.gate_start_offset_mu) + self.core.seconds_to_mu(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM)

                detector = -1
                mapping_click_time = 0

                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 0, SPCM0_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 1, SPCM1_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 2, SPCM0_OtherNode_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:
                    detector, mapping_click_time = 3, SPCM1_OtherNode_click_time

                ### Drain all additional edges, including when the timestamp buffer is already full.
                ### With no valid single-detector event (detector == -1, e.g. two SPCMs clicked in the
                ### same gate), parity_10 never drained, so the leftover edges reappeared in the next
                ### attempt as false early counts. The valid-event path drains after atom_parity_shot
                ### below, keeping the time-critical at_mu() mapping jump free of extra input reads.
                if detector < 0:
                    if SPCM0_click_time > 0: self.ttl_SPCM0.count(t_end_SPCM)
                    if SPCM1_click_time > 0: self.ttl_SPCM1.count(t_end_SPCM)
                    if SPCM0_OtherNode_click_time > 0: self.ttl_SPCM0_OtherNode.count(t_end_SPCM)
                    if SPCM1_OtherNode_click_time > 0: self.ttl_SPCM1_OtherNode.count(t_end_SPCM)

                if detector >= 0:
                    at_mu(mapping_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
                    ##todo: I am setting FORT drop in dma because FORT power barely drifts.

                    # delay(2 * us)
                    self.core_dma.playback_handle(parity_mapping_dma_handle)

                    # delay(5 * us)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    ##todo: I am setting FORT drop in dma because FORT power barely drifts.
                    # t_after_parity_mapping = now_mu()

                    self.core.break_realtime()

                    delay(5 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(5 * us)
                    atom_parity_shot(self)

                    ### actually this is not necessary;
                    ### Drain the herald channel's additional edges from the photon gate.
                    ### Exactly one detector clicked in this branch, so only that channel can hold leftovers.
                    if SPCM0_click_time > 0: self.ttl_SPCM0.count(t_end_SPCM)
                    if SPCM1_click_time > 0: self.ttl_SPCM1.count(t_end_SPCM)
                    if SPCM0_OtherNode_click_time > 0: self.ttl_SPCM0_OtherNode.count(t_end_SPCM)
                    if SPCM1_OtherNode_click_time > 0: self.ttl_SPCM1_OtherNode.count(t_end_SPCM)

                    self.core.break_realtime()
                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0 if detector == 0 else 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0 if detector == 1 else 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0 if detector == 2 else 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0 if detector == 3 else 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP

                    delay(1 * ms)

                    # t_delay_checker = (t_after_parity_mapping - t_after_OP_dma)
                    # print("t_delay_checker: ", t_delay_checker)

                    self.measurement += 1
                    break

            delay(20 * us)
            excitation_cycle += 1


        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        delay(1 * ms)
        second_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * ms)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        self.measurement -= 1
        end_measurement(self)

        delay(5 * ms)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    # delay(15 * ms)
    self.core.break_realtime()
    for i in range(self.n_measurements):
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        # self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        # self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton[i])
        self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(5 * ms)
