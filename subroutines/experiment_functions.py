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
def test_RF_pulse_experiment(self):
    """
    testing the effect of RF pulse with external RF antenna on scope
    """
    self.core.reset()

    self.measurement = 0
    # self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.set(frequency=2.0*MHz, amplitude=dB_to_V(-3.0))

    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    while self.measurement < self.n_measurements:
        delay(100*us)
        self.dds_MW_RF.sw.on()  ### turn on RF
        self.dds_FORT.sw.on()
        delay(1*us)
        self.dds_MW_RF.sw.off()  ### turn off RF
        self.dds_FORT.sw.off()

@kernel
def test_Repump_pulse_experiment(self):
    """
    testing the response of repump aom
    """
    self.core.reset()

    self.measurement = 0
    # self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))

    self.dds_AOM_A6.sw.on()
    delay(1*ms)
    while self.measurement < self.n_measurements:
        delay(1000*ms)
        self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        delay(1000*ms)
        self.ttl_repump_switch.on()  # turns the MOT RP AOM off

    self.dds_AOM_A6.sw.off()

@kernel
def test_BA_pulse_experiment(self):
    self.core.reset()

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(1 * s)
        self.dds_FORT.sw.on()
        delay(1*ms)

        self.core_dma.playback_handle(ba_dma_handle)

        delay(10*ms)

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

@kernel
def tune_coil_driver_experiment(self):
    self.core.reset()
    delay(1 * s)
    ### set coils for pumping
    for i in range(100000):
        ### Set the coils to MOT loading setting
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)

        delay(4*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(4 * ms)  # coil relaxation time

@kernel
def test_excitation_rise_time_experiment(self):

    """
    Simple code to test the rise/fall time of the excitation pulse.
    """
    self.core.reset()
    delay(1 * s)

    self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))

    # turning on excitation DDS and switch
    self.ttl_exc0_switch.off()  # EXC0 AOM ON
    self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet

    # delay(1* s)
    # for i in range(self.n_measurements):
    #
    #     t1 = now_mu()
    #     at_mu(t1 + 100)
    #     self.ttl_GRIN2_switch.off()        # turn on GRIN2 AOM
    #
    #     at_mu(t1 + 100+ 100)
    #     self.ttl_GRIN2_switch.on()          # turn off GRIN2 AOM
    #     # delay(4 * ms)  # coil relaxation time
    #     delay(1*us)

    delay(1* s)
    for i in range(10000000):
        self.ttl_GRIN2_switch.off()        # turn on GRIN2 AOM
        delay(30*ns)
        self.ttl_GRIN2_switch.on()          # turn off GRIN2 AOM
        # delay(4 * ms)  # coil relaxation time
        delay(100*us)

@kernel
def run_zotino_stabilization_Aqua(self):
    ####Start Kais's calibration code
    ref_voltage_V = 0.00000
    logger.debug('Hardware Setup')
    self.core.reset()
    self.zotino0.init()
    self.core.break_realtime()
    self.sampler1.init()
    self.core.break_realtime()
    delay(10 * ms)
    for q in [0,1,3,4,5,6,7]:
        delay(10 * ms)
        self.sampler1.set_gain_mu(q, 3)

    diff = [0.01] * 8
    MV = [0.0] * 8
    self.zotino0.init()
    delay(10 * ms)
    for m in range(8):
        delay(10 * ms)
        self.zotino0.write_dac(m, ref_voltage_V)
    self.zotino0.load()
    delay(10 * ms)
    for j in range(20):
        delay(10 * ms)
        n_channels = 8  # sets number of channels to read off
        smp = [0.0] * n_channels  # creates list of floating point variables
        self.sampler1.sample(smp)
        delay(10 * ms)
        # print (smp)
        for k in [0,1,3,4,5,6,7]:
            diff[k] = ref_voltage_V - smp[k]
            # print(diff[k])
        # print(k)

        delay(0.1)
        self.zotino0.init()
        self.core.break_realtime()
        delay(10 * ms)
        for n in [0,1,3,4,5,6,7]:
            MV[n] = MV[n] + diff[n]
            self.zotino0.write_offset(n, MV[n])
            # print(n)
        self.zotino0.load()
        delay(10 * ms)
    self.core.break_realtime()
    #####End Kais's calibration code

@kernel
def run_zotino_stabilization(self):
    """
    In the beginning of an experiment, measure the zotino offset and compensate it.

    :param self:
    :return:
    """

    ref_voltage_V = 5.00000
    self.core.reset()

    measurement_buf = np.array([0.0] * 8)

    measurement1 = 0.0  # 1
    measurement2 = 0.0  # 2
    measurement3 = 0.0  # 3
    measurement4 = 0.0  # 4
    measurement5 = 0.0
    measurement6 = 0.0
    measurement7 = 0.0
    measurement8 = 0.0

    self.zotino0.set_dac([ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V],
                         [8,9,10,11,12,13,14,15])

    delay(10*ms)

    avgs = 50
    for i in range(avgs):
        self.sampler2.sample(measurement_buf)

        delay(0.1 * ms)
        measurement1 += measurement_buf[0]
        measurement2 += measurement_buf[1]
        measurement3 += measurement_buf[2]
        measurement4 += measurement_buf[3]
        measurement5 += measurement_buf[4]
        measurement6 += measurement_buf[5]
        measurement7 += measurement_buf[6]
        measurement8 += measurement_buf[7]
        delay(0.1 * ms)

    measurement1 /= avgs
    measurement2 /= avgs
    measurement3 /= avgs
    measurement4 /= avgs
    measurement5 /= avgs
    measurement6 /= avgs
    measurement7 /= avgs
    measurement8 /= avgs

    self.set_dataset("zotino_test1_offset", measurement1 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test2_offset", measurement2 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test3_offset", measurement3 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test4_offset", measurement4 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test5_offset", measurement5 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test6_offset", measurement6 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test7_offset", measurement7 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test8_offset", measurement8 - ref_voltage_V, broadcast=True, persist=True)

    self.append_to_dataset("zotino_test1_offset_monitor", measurement1)
    self.append_to_dataset("zotino_test2_offset_monitor", measurement2)
    self.append_to_dataset("zotino_test3_offset_monitor", measurement3)
    self.append_to_dataset("zotino_test4_offset_monitor", measurement4)
    self.append_to_dataset("zotino_test5_offset_monitor", measurement5)
    self.append_to_dataset("zotino_test6_offset_monitor", measurement6)
    self.append_to_dataset("zotino_test7_offset_monitor", measurement7)
    self.append_to_dataset("zotino_test8_offset_monitor", measurement8)

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
        self.dds_FORT.sw.on()  ### turns FORT on
        delay(0.1*ms)

        record_FORT_MM_power(self)
        record_FORT_APD_power(self)

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

                ### restart the loading attempt with the (possibly) improved shims
                try_n = 0
                t_before_atom = now_mu()
                delay(0.1 * ms)
                continue

            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms) ### necessary to avoid underflow

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

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
    # delay(0.1 * ms)
    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

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
        self.stabilizer_FORT.run(setpoint_index=1)  # set to  science setpoint
        self.stabilizer_FORT.run(setpoint_index=2)  # the science holding setpoint

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
        if self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(self.t_PGC_after_loading)

    self.dds_cooling_DP.sw.off()

@kernel
def testing_record_chopped_FORT_and_TTL(self):
    ### Akbar 2025-11-25
    t_FORT_offset = self.core.seconds_to_mu(100 * ns)  ### delay wrt "t0" in each cycle
    t_FORT_OFF = self.core.seconds_to_mu(self.chop_RO_FORT_off) ### the length of the FORT OFF time

    t_cooling_offset = self.core.seconds_to_mu(self.chop_RO_pulse_offset) ### delay wrt "t0" in each cycle
    t_cooling_length = self.core.seconds_to_mu(self.chop_RO_pulse_length) ### duration of the cooling pulse

    t_gate_offset = self.core.seconds_to_mu(self.chop_RO_gate_offset) ### delay wrt t0 for gate to open
    t_gate_length = self.core.seconds_to_mu(self.chop_RO_gate_length) ### gate open time

    with self.core_dma.record("FORT_and_TTL_chopped"):
        for i in range(4000): ### number of cycles to repeat
            t0 = now_mu()

            at_mu(t0 + t_FORT_offset)
            self.dds_FORT.sw.off()
            at_mu(t0 + t_FORT_offset + t_FORT_OFF)
            self.dds_FORT.sw.on()

            at_mu(t0 + t_cooling_offset)
            self.dds_cooling_DP.sw.on()  ### Turn on cooling
            at_mu(t0 + t_cooling_offset + t_cooling_length)
            self.dds_cooling_DP.sw.off()  ### Turn off cooling

            at_mu(t0 + t_gate_offset)
            self.ttl0._set_sensitivity(1)
            self.ttl1._set_sensitivity(1)

            at_mu(t0 + t_gate_offset + t_gate_length)
            self.ttl0._set_sensitivity(0)
            self.ttl1._set_sensitivity(0)

            delay(0.5 * us)

@kernel
def testing_chopped_FORT_and_TTL_experiment(self):
    """
    xxx
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable
    delay(1 * ms)
    self.ttl_SPCM0._set_sensitivity(0)
    self.ttl_SPCM1._set_sensitivity(0)

    ### record and get handle
    testing_record_chopped_FORT_and_TTL(self)
    delay(10 * ms)
    test_chopped_handle = self.core_dma.get_handle("FORT_and_TTL_chopped")
    delay(20 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        self.dds_cooling_DP.sw.off()  ### Turn off cooling

        delay(1 * ms)
        self.dds_FORT.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        delay(1 * ms)
        self.ttl_SPCM0._set_sensitivity(0)
        self.ttl_SPCM1._set_sensitivity(0)

        delay(1 * ms)

        ### DMA playback
        self.core_dma.playback_handle(test_chopped_handle)

        delay(10 * ms)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(now_mu())
        self.SPCM1_RO1 = self.ttl_SPCM1.count(now_mu())
        self.AllSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)

        # delay(10*ms)
        self.core.break_realtime()

        end_measurement(self)

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
def test_shot(self):
    """
    Based on non-chopped atom readout. Can be used as a middle RO betweek RO1 and RO2 to opimize the readout settings.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP
    """
    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(5 * us)

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
    delay(5 * us)

    with parallel:
        self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_test_shot)

    self.SPCM0_test_RO = self.ttl_SPCM0_counter.fetch_count()
    self.SPCM1_test_RO = self.ttl_SPCM1_counter.fetch_count()
    self.SPCM0_OtherNode_test_RO = self.ttl_SPCM0_OtherNode_counter.fetch_count()
    self.SPCM1_OtherNode_test_RO = self.ttl_SPCM1_OtherNode_counter.fetch_count()
    self.AllSPCMs_test_RO = self.SPCM0_test_RO + self.SPCM1_test_RO + self.SPCM0_OtherNode_test_RO + self.SPCM1_OtherNode_test_RO

    delay(10 * us)
    self.dds_cooling_DP.sw.off()  ### turn off cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP
    delay(10 * us)

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
            self.ttl_pumping_repump_switch.on()
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
    if self.t_depumping + self.t_pumping > 1*ms:
        delay(1500 * us)  # we need extra slack

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light
    delay(10 * us)
    self.dds_AOM_A1.sw.off()
    self.dds_AOM_A2.sw.off()
    self.dds_AOM_A3.sw.off()
    self.dds_AOM_A4.sw.off()
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    delay(10 * us)

    if self.which_node == "alice":
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
        ### This is split into GRIN1 and GRIN2 in Alice
    else:
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(-5.0))
        ### This is just for GRIN1 in Bob
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))
    delay(1 * us)

    ### Tunring on pumping RP:
    self.ttl_pumping_repump_switch.off()
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
    delay(1 * ms)  # coil relaxation time

    self.core_dma.playback_handle(op_dma_handle)
    delay(self.t_depumping)
    self.dds_D1_pumping_DP.sw.off() ### turning off D1 DP
    self.ttl_pumping_repump_switch.on() ### turning off pumping RP

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
    delay(10*us)
    # self.dds_cooling_DP.sw.off()

@kernel
def CW_optical_pumping_node1(self):
    """
    - Includes turning on or off some AOMs, or setting their powers, setting the coils, etc.
    - D1 is on GRIN1
    - I am bypassing GRIN1 AOM to increase the D1 power at the atom to speed up O.P.
    """

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light
    delay(5 * us)

    ### Turning on fiber AOMs for delivery of the pumping repump
    # self.dds_AOM_A1.set(frequency=self.AOM_A1_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A3.set(frequency=self.AOM_A3_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))

    delay(5 * us)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    # ### so that D1 can pass
    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
    # delay(5 * us)
    # self.GRIN1and2_dds.sw.on()
    # self.ttl_GRIN1_switch.off()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.4 * ms)  # coil relaxation time

    ### Optical pumping phase
    self.ttl_pumping_repump_switch.off()
    self.dds_D1_pumping_DP.sw.on()

    delay(self.t_pumping)

    self.dds_D1_pumping_DP.sw.off()
    delay(self.t_depumping)
    self.ttl_pumping_repump_switch.on()

    delay(2*us)

    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

    delay(1 * us)

    # self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=self.stabilizer_AOM_A1.amplitude)
    # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=self.stabilizer_AOM_A2.amplitude)
    # self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)
    # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=self.stabilizer_AOM_A4.amplitude)
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    delay(5*us)
    # self.GRIN1and2_dds.sw.off()
    # self.ttl_GRIN1_switch.on()

@kernel
def CW_optical_pumping_node2(self):
    """
    *** previously named as "optical_pumping_GRIN1"
    *** OP with GRIN1
    optical pumping without chopping the FORT - swtching ON/OFF with TTL

    ** OP AOM is driven with external RF source.
    ** ONLY GRIN1 AOM is used to turn on/off the OP.

    Note: To avoid conflict with Node1 codes, I left the names of the dds channels that are now
    used for GRIN1 and GRIN2 dds.

    Name of the dds channels >>>>   what actually does here
    GRIN1and2_dds            >>>>   GRIN1 dds
    dds_D1_pumping_DP        >>>>   GRIN2 dds

    Eunji - for Node2

    """
    self.dds_cooling_DP.sw.off()  # no MOT cooling light
    self.ttl_repump_switch.on()   # no MOT RP AOM
    self.ttl_exc0_switch.on()     # no excitation
    self.ttl_D1_pumping.on()      # no D1
    delay(5 * us)

    ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))
    delay(5 * us)

    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.395 * ms)  # coil relaxation time

    ### so that D1 can pass
    self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
    self.GRIN1and2_dds.sw.on()  ## GRIN1 RF ON but not yet activated with external switch
    delay(5 * us)

    with parallel:
        self.ttl_pumping_repump_switch.off()  ## pumping repump ON
        self.ttl_D1_pumping.off()  ##turning D1 ON
        self.ttl_GRIN1_switch.off()  ## turning D1 ON

    delay(self.t_pumping)     ## pumping time

    with parallel:
        self.ttl_D1_pumping.on()  ## turning D1 OFF
        self.ttl_GRIN1_switch.on()  ## GRIN OFF

    delay(self.t_depumping)    ## depumping time

    self.ttl_pumping_repump_switch.on()  ## pumping repump OFF
    delay(2 * us)

    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1]) # back to science setpoint
    # delay(1*us)

    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    self.GRIN1and2_dds.sw.off()  ## GRIN1 RF OFF
    delay(1 * us)

    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
    delay(5 * us)

@kernel
def record_CW_optical_pumping_node1(self):
    with self.core_dma.record("CW_optical_pumping_node1"):
        ### Turn off unused light
        with parallel:
            self.ttl_repump_switch.on()
            self.ttl_exc0_switch.on()
            self.dds_cooling_DP.sw.off()

        ### Configure pumping-repump delivery AOMs
        ### These must remain sequential if they share the same Urukul SPI bus
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

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
        delay(0.4 * ms)

        ### Optical pumping
        with parallel:
            self.ttl_pumping_repump_switch.off()
            self.dds_D1_pumping_DP.sw.on()

        delay(self.t_pumping)

        self.dds_D1_pumping_DP.sw.off()

        delay(self.t_depumping)

        ### Turn off all pumping RF/light paths
        with parallel:
            self.ttl_pumping_repump_switch.on()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

@kernel
def record_CW_optical_pumping_node2(self):
    with self.core_dma.record("CW_optical_pumping_node2"):
        ### Turn off unused light
        with parallel:
            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            self.ttl_exc0_switch.on()
            self.ttl_D1_pumping.on()

        ### Configure pumping-repump delivery AOMs
        ### These must remain sequential if they share the same Urukul SPI bus
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

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
        delay(0.395 * ms)

        ### Configure GRIN1 DDS
        self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        self.GRIN1and2_dds.sw.on()
        delay(2*us)

        ### Optical pumping
        with parallel:
            self.ttl_pumping_repump_switch.off()
            self.ttl_D1_pumping.off()
            self.ttl_GRIN1_switch.off()

        delay(self.t_pumping)

        ### Turn off D1 pumping light
        with parallel:
            self.ttl_D1_pumping.on()
            self.ttl_GRIN1_switch.on()

        delay(self.t_depumping)

        ### Turn off all pumping RF/light paths
        with parallel:
            self.ttl_pumping_repump_switch.on()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            self.GRIN1and2_dds.sw.off()

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        delay(1.7 * us)  ### to match the timing with node1

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
    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))

    # UNBLOCKING GRIN1
    self.ttl_GRIN1_switch.off()
    # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(0.0))
    self.GRIN1and2_dds.sw.on()

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
    self.append_to_dataset("Magnetometer_Zero_X", MagnetometerZ * 350) ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_Zero_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_Zero_Z", MagnetometerY * 350)
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
    self.append_to_dataset("Magnetometer_OP_X", MagnetometerZ * 350)  ### 1V corresponds to 350 mG
    self.append_to_dataset("Magnetometer_OP_Y", MagnetometerX * 350)  ### sensor's X axis is coils' Y axis, and vice versa.
    self.append_to_dataset("Magnetometer_OP_Z", MagnetometerY * 350)

    #####################################  Measure in the MOT phase
    ### Set the coils to MOT loading setting
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
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
    #     measure_REPUMP(self)
    # elif self.which_node == "bob":
    #     measure_Magnetometer_Node2(self)

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

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

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

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100 * ms)

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

# @kernel
# def trap_frequency_experiment(self):
#     """
#     For spectroscopy of the trap vibrational frequencies with the Rigol D11022Z function generator.
#
#     One way to use this is with GeneralVariableScan and scan_variable1 = f_Rigol_modulation,
#     where f_Rigol_modulation must be within the boundaries Rigol_carrier_frequency +/- Rigol_FM_deviation
#     specified in ExpermientVariables (which in turn must be set on the Rigol D11022Z).
#
#     :param self: an experiment instance.
#     :return:
#     """
#
#     set_RigolDG1022Z(frequency=self.Rigol_carrier_frequency,
#                      vpp=self.Rigol_V_pp,
#                      vdc=self.Rigol_V_DC)
#
#     self.core.reset()
#
#     self.SPCM0_RO1 = 0
#     self.SPCM0_RO2 = 0
#
#     self.require_D1_lock_to_advance = False # override experiment variable
#
#     # self.set_dataset(self.SPCM0_rate_dataset,
#     #                  [0.0],
#     #                  broadcast=True)
#
#     self.measurement = 0
#     while self.measurement < self.n_measurements:
#
#         #TODO: just set the rigol frequency using pyvisa
#
#         if self.enable_laser_feedback:
#             run_feedback_and_record_FORT_MM_power(self)
#
#         load_MOT_and_FORT(self)
#
#         delay(0.1*ms)
#         self.zotino0.set_dac(
#             [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
#             channels=self.coil_channels)
#
#         # set the FORT AOM to the science setting. this is only valid if we have run
#         # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT
#         self.dds_FORT.set(frequency=self.f_FORT,
#                                 amplitude=self.stabilizer_FORT.amplitudes[1])
#
#         # set the cooling DP AOM to the readout settings
#         self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
#                                 amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)
#
#         if not self.no_first_shot:
#             first_shot(self)
#
#             if self.t_recooling_after_first_shot > 0:
#                 recooling_after_first_shot(self)
#
#         ##############################################################################
#         # modulate the FORT
#         ##############################################################################
#
#         self.FORT_mod_switch.on()  # toggle the modulation to the VCA
#         delay(10 * ms)
#         self.FORT_mod_switch.off()
#         delay(1 * ms)
#
#         second_shot(self)
#
#         end_measurement(self)
#
#     self.dds_FORT.sw.off()

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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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

        ### with cw pumping:
        if self.t_pumping > 0.0:
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)

        ##todo: this is hardcoded delay to acount for coil drift- set it smaller so that mapping is tuned correclty
        delay(10*us)

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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

    self.core.break_realtime()

    # delay(1 * ms)
    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

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

        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom_recycle(self)
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

        ### with cw pumping:
        if self.t_pumping > 0.0:
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
            delay(5.5 * us) ##todo: this is hardcoded delay to acount for coil drift- set it smaller so that mapping is tuned correctly

        ############################ microwave phase to transfer population from F=1,mF=0 to F=2,mF=1

        ### Changing the bias field for microwave
        # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
        #                       self.AX_volts_microwave, self.AY_volts_microwave],
        #                      channels=self.coil_channels)

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

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2*us)

        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(2 * us)

        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
        delay(2 * us)

        ##### Mapping from |2,1> to |1,1>
        if self.t_microwave_11_pulse > 0.0:
            self.ttl_microwave_switch.off()
            delay(self.t_microwave_11_pulse)
            self.ttl_microwave_switch.on()

        delay(5 * us)
        # FORT_ramp2_smoothstep(self, direction="up")
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
        delay(2 * us)

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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)
    # record_chopped_optical_pumping(self)
    # delay(100 * ms)
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
        #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
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
        #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

    self.core.break_realtime()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")

    self.core.break_realtime()

    delay(1 * ms)

    self.dds_microwaves.sw.on()  ### turns on the DDS not the switches.

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        self.ttl_exc0_switch.on()  # turns off the excitation

        # load_MOT_and_FORT_until_atom(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
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
            if self.which_node == "alice":
                # CW_optical_pumping_node1(self)
                self.core_dma.playback_handle(CW_OP_node1_handle)
            else:
                # CW_optical_pumping_node2(self)
                self.core_dma.playback_handle(CW_OP_node2_handle)
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
def single_photon_experiment_3_atom_loading_advance_node2_AllSPCMs(self):
    """

    This is based on load_MOT_and_FORT_until_atom_and_recycle.
    Does not check for atom after each excitation attempt:

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

    AllSPCMs_RO_atom_check_array = [0]

    delay(100 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles / self.atom_check_every_n)
        tStamps_t1 = [0.0] * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        delay(100 * ms)  ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)
        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### 1) use ttl_D1_pumping/ttl_excitation to send D1/Exc to GRIN1/GRIN2
        ### 2) use ttl_GRIN1_switch/ttl_GRIN2_switch to send D1/EXC to atoms

        self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5 * us)
        self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 1  ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            delay(10 * ms)

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                CW_optical_pumping_node2(self)
                delay(10 * us)

                # ############ microwave phase - ONLY USED FOR VERIFYING OP.
                # if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                #     self.ttl_microwave_switch.off()  # todo: switching on with external RF switch creates a lag.
                #     delay(1 * us)
                #     t1 = now_mu()
                #     self.dds_microwaves.sw.on()
                #
                #     at_mu(t1 + int(self.t_microwave_pulse / ns))
                #     self.dds_microwaves.sw.off()
                #     self.ttl_microwave_switch.on()
                #
                #     delay(0.1 * ms)
                #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                # ############ blow-away phase - push out atoms in F=2 only
                # if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                #     self.core_dma.playback_handle(ba_dma_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(5 * us)
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                # # ############## previous code
                # # FORT_offset = 145
                # # FORT_offset = 90   ## as of 2026/03/08 Excelitas SPCMs
                # FORT_offset = 280   ## as of 2026/05/27 same Excelitas SPCMs; drifted.. why!
                #
                # t1 = now_mu()
                # at_mu(t1 + int(self.t_excitation_offset_mu))
                # self.ttl_GRIN2_switch.off()  # turns on excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                # self.ttl_GRIN2_switch.on()  # turns off excitation
                #
                # #todo: Turn FORT OFF earlier than 100ns. fall time is ~100ns
                # # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 145)  # turn off the FORT 145ns earlier than EXC pulse
                #
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # self.dds_FORT.sw.off()  ### turns FORT off
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset + int(self.t_photon_collection_time / ns)) # turning FORT ON 50ns later than gate end
                # self.dds_FORT.sw.on()  ### turns FORT on
                #
                # ## open SPCM gate 50ns earlier
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.gate_start_offset_mu) -50)  # turn the gate on 50ns earlier than EXC pulse
                # with parallel:
                #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                #
                # # ############## newer code 2026-05-28

                ###todo: 2026/05/28 Note - Eunji
                ### To make FORT, EXC, gate to happen at the same time, I need the offset to be:
                ### FORT at t_FORT
                ### EXC  at t_FORT - 280ns; t_EXC_hardware_offset_mu = -280
                ### GATE at t_FORT + 465ns;  t_GATE_hardware_offset_mu  = +465
                ### These variables do not change so I have it hard coded here.
                ### ### CASE - tested with this single photon code
                ### **** To turn EXC pulse on 50ns after FORT is turned off, I need to set
                ###      t_excitation_offset_mu = -280 + 50 = -230
                ### **** To turn GATE pulse on 25ns after FORT is turned off, I need to set
                ###      gate_start_offset_mu = +465 + 25 = 490 which is 50ns before excitation pulse
                ### **** To turn FORT ON after gate pulse is done,
                ###      Turn off at t_FORT t_photon_collection_time + 50(additional delay to account for gate delay)
                ############

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on 50ns after photon collection ended

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation after t_excitation_offset_mu

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))  ## turns on gate after gate_start_offset_mu
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1


                # at_mu(t1 + 30000)
                tStamps_t1[
                    excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
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
                self.dds_AOM_A4.sw.on()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

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
                self.dds_AOM_A4.sw.on()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())

                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
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

            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()

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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle / self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1) * self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_Photon_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_Photon_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_Photon_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_Photon_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(15 * ms)

@kernel
def single_photon_experiment_3_atom_loading_advance_AllSPCM(self):
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

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)

        first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        # ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        # self.GRIN1and2_dds.sw.on()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            # ############################### chopped optical pumping phase - pumps atoms into F=1,m_F=0
            # if self.t_pumping > 0.0:
            #
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
            #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0)) ### set to 5V for optical pumping
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            #     delay(1 * us)
            #
            #     ### Tunring on pumping RP:
            #     self.ttl_pumping_repump_switch.off()
            #     self.dds_AOM_A5.sw.on()
            #     self.dds_AOM_A6.sw.on()
            #
            #     # delay(1 * ms)
            #
            #     self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1
            #     delay(10 * us)
            #
            #     self.core_dma.playback_handle(op_dma_handle)
            #     delay(self.t_depumping)
            #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            #     self.ttl_pumping_repump_switch.on()  ### turning off pumping RP
            #
            #     delay(2 * us)
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(100 * us)
            #
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            #     # delay(1 * ms)
            #
            #     self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1
            #     delay(10 * us)

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)

                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            # delay(10*us)

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)

            # self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            # delay(2 * ms)
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                delay(5 * us)

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0 events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1 events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                # at_mu(t1 + 30000)
                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            # self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
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
                # delay(0.4*ms)

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(1 * us)

                # with parallel:
                #     self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                #     self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                #
                # SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                # SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                # AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                # AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

                # #### these added on 2025-09-09 based on Eunji's comment. I (AkS) did not have these in my photon experiments:
                # self.dds_cooling_DP.sw.off()
                # self.ttl_repump_switch.on()
                # delay(1 * us)
                # self.dds_AOM_A1.sw.off()
                # self.dds_AOM_A2.sw.off()
                # self.dds_AOM_A3.sw.off()
                # self.dds_AOM_A4.sw.off()
                # self.dds_AOM_A5.sw.off()
                # self.dds_AOM_A6.sw.off()
                # delay(1 * us)
                # ##############################

            delay(10 * us)

        # delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

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

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_3_atom_heat_test(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    This is similar to single_photon_experiment_3 that I used for photon generation, but for the purpose of optimizing the cycle.
    Here I am going to see how many excitation + OP cycle I can run (max_excitation_cycles) on the atom before loosing the atom
    without recooling, or keep max_excitation_cycles fixed, and change t_recooling and see how that affects the retention from
    the 2nd readout (RO2).

    Akbar 2026-03-18
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    # record_chopped_optical_pumping(self)
    # delay(100 * ms)
    #
    # if self.verify_OP_in_photon_experiment:
    #     record_chopped_blow_away(self)
    #     ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    #
    #     self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    #     delay(10 * ms)
    #     self.dds_microwaves.sw.on()
    #     delay(100 * ms)
    #
    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

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
        # tStamps_t1 = [0.0] * (self.max_excitation_cycles * self.n_excitation_attempts)
        # SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        # SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
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
        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 1  ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):

            self.core.break_realtime()

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            # ############################### copped optical pumping phase - pumps atoms into F=1,m_F=0
            # if self.t_pumping > 0.0:
            #
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
            #     self.ttl_pumping_repump_switch.off()
            #     self.dds_AOM_A5.sw.on()
            #     self.dds_AOM_A6.sw.on()
            #
            #     delay(1 * us)
            #
            #     self.ttl_GRIN1_switch.off()  ### turn ON GRIN1 AOM for D1
            #     delay(10 * us)
            #
            #     self.core_dma.playback_handle(op_dma_handle)
            #     delay(self.t_depumping)
            #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            #     self.ttl_pumping_repump_switch.on()  ### turning off pumping RP
            #
            #     delay(2 * us)
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(100 * us)
            #
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            #     # delay(1 * us)
            #
            #     self.ttl_GRIN1_switch.on()  ### turn OFF GRIN1 AOM
            #     delay(10 * us)

            ### with cw pumping:
            if self.t_pumping > 0.0:
                delay(2 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(2 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(100 * us)

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

                # at_mu(t1 + int(self.gate_start_offset_mu))
                # with parallel:
                #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                #
                # ### timestamping SPCM0 events
                # while SPCM0_click_counter < max_clicks:
                #     SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                #     if SPCM0_click_time == -1.0:
                #         break
                #     # SPCM0_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                #     #     SPCM0_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)
                #     SPCM0_click_counter += 1
                #
                # ### timestamping SPCM1 events
                # while SPCM1_click_counter < max_clicks:
                #     SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                #     if SPCM1_click_time == -1.0:
                #         break
                #     # SPCM1_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                #     #     SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                #     SPCM1_click_counter += 1

                # at_mu(t1 + 30000)
                # tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization,
                     self.AX_volts_PGC_optimization, self.AY_volts_PGC_optimization],
                    channels=self.coil_channels)

                ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)
                delay(1 * ms)  ### coils relaxation time

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        # delay(0.1 * ms)
        self.core.break_realtime()

        end_measurement(self)

        delay(5 * ms)

        delay(1 * ms)
        # for i in range((excitation_cycle + 1) * self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle+1)

        delay(1 * ms)


    delay(15 * ms)

@kernel
def single_photon_experiment_4_atom_loading_advance(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    This is similar to single_photon_experiment_3_atom_loading_advance but with modified OP to speed up the rate.
    If max_excitation_cycles is large (like 3000), the timestamp arrays get too large and artiq may freeze
    randomly after sometime. So, keep max_excitation_cycles <= 2000 and n_measurements <=100.

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

    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

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

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
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

        self.ttl_GRIN2_switch.on()  # turns off excitation

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light

        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()

        self.ttl_exc0_switch.off()  # turns on the excitation0 AOM


        ################ Preparing for fast O.P.
        self.ttl_repump_switch.on()  # turns off the MOT RP AOM
        self.dds_cooling_DP.sw.off()  # no cooling light

        ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        delay(5 * us)
        # self.dds_AOM_A5.sw.on()
        # self.dds_AOM_A6.sw.on()

        ### set coils for pumping
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time
        ##########################################


        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                ### Optical pumping phase
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                self.ttl_pumping_repump_switch.off()
                self.dds_D1_pumping_DP.sw.on()

                delay(self.t_pumping)

                self.dds_D1_pumping_DP.sw.off()
                delay(self.t_depumping)
                self.ttl_pumping_repump_switch.on()

                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                delay(1 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)


            # delay(2 * ms)
            self.core.break_realtime()

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

            # delay(20 * us)
            # self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

                delay(0.4 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(5 * us)
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                delay(1 * us)

                ### to be prepared for OP
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                    delay(5*us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
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
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(5*us)
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break


                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

                delay(20*us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

            delay(10 * us)

        # delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

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

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_5_atom_loading_advance_AllSPCM(self):
    """
    ### this is specifically for node1

    This is similar to single_photon_experiment_4_atom_loading_advance but I have reduced the size of the timestamp array
    to avoid freezing artiq. For this:
        1- I am not looking for 2nd clicks in each SPM window. Effectively, this is like max_clicks = 1.
        2- I am not appending the timestamps when a photon is not registered.
        3- I initialize and use a shorter array than max_excitation_cycles because we detect 5% of the times on average.


    IMPORTANT: for this function to work, you need to modify BaseExperiment to be like the following:
        # self.experiment.set_dataset("SPCM0_SinglePhoton_tStamps", [[0.0,0.0]], broadcast=True)
        # self.experiment.set_dataset("SPCM1_SinglePhoton_tStamps", [[0.0,0.0]], broadcast=True)
        self.experiment.set_dataset("SPCM0_SinglePhoton_tStamps", [0.0], broadcast=True)
        self.experiment.set_dataset("SPCM1_SinglePhoton_tStamps", [0.0], broadcast=True)

    """

    self.core.reset()
    delay(1 * ms)

    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

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

        self.ttl_GRIN2_switch.on()  # turns off excitation

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light

        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()

        # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

        ################ Preparing for fast O.P.
        self.ttl_repump_switch.on()  # turns off the MOT RP AOM
        self.dds_cooling_DP.sw.off()  # no cooling light

        ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
        delay(5 * us)
        # self.dds_AOM_A5.sw.on()
        # self.dds_AOM_A6.sw.on()

        ### set coils for pumping
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time
        ##########################################


        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                ### Optical pumping phase
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                self.ttl_pumping_repump_switch.off()
                self.dds_D1_pumping_DP.sw.on()

                delay(self.t_pumping)

                self.dds_D1_pumping_DP.sw.off()
                delay(self.t_depumping)
                self.ttl_pumping_repump_switch.on()

                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                delay(1 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)


            # delay(2 * ms)
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


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

                delay(30 * us)  ### 20us is not enough

            # delay(20 * us)
            # self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

                delay(0.4 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(5 * us)
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                delay(1 * us)

                ### to be prepared for OP
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                    delay(5*us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
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
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(5*us)
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

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
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break


                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

                delay(20*us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

            delay(10 * us)

        # delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

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


        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_6_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_5_atom_loading_advance. But works for any nodes.
    """

    self.core.reset()
    delay(1 * ms)

    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

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

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)


            # delay(2 * ms)
            # self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                delay(5 * us)

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


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

                delay(30 * us)  ### 20us is not enough

            # delay(20 * us)
            # self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

                delay(0.4 * ms)  ### coils relaxation time

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
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
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

            delay(10 * us)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
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

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_7_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_6_atom_loading_advance.
    But uses DMA for
    * OP
    * recooling
    """

    self.core.reset()
    delay(1 * ms)

    #### recording DMA
    record_CW_optical_pumping_node2(self)
    record_recooling(self)

    #### getting DMA handles
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
    recooling_dma_handle = self.core_dma.get_handle("recooling")

    delay(1 * ms)
    self.core.break_realtime()
    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

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

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    delay(10 * us)
                    CW_optical_pumping_node1(self)
                    delay(10 * us)
                else:
                    self.core_dma.playback_handle(CW_OP_node2_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                delay(5 * us)

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


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

                delay(10 * us)

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.core_dma.playback_handle(recooling_dma_handle)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(0.4*ms)

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

            delay(10 * us)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
        delay(0.1 * ms)

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        delay(5 * us)

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
def single_photon_experiment_8_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_7_atom_loading_advance. with DMA in excitation phase
    But uses DMA for
    * OP
    * recooling
    * EXC
    """

    self.core.reset()
    delay(1 * ms)

    #### recording DMA
    record_CW_optical_pumping_node2(self)
    record_recooling(self)
    record_excitation_and_photon_collection(self)

    #### getting DMA handles
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
    recooling_dma_handle = self.core_dma.get_handle("recooling")
    excitation_dma_handle = self.core_dma.get_handle("excitation_and_photon_collection")

    delay(1 * ms)
    self.core.break_realtime()
    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

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

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    delay(10 * us)
                    CW_optical_pumping_node1(self)
                    delay(10 * us)
                else:
                    self.core_dma.playback_handle(CW_OP_node2_handle)

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

                delay(5 * us)

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.core_dma.playback_handle(recooling_dma_handle)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                # self.zotino0.set_dac(
                #     [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                #     channels=self.coil_channels)
                #
                # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                # delay(0.4*ms)
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

            delay(5 * us)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
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
def atom_photon_parity_2_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Used to get nice parity oscillation before state rotation.

    A simple parity oscillation experiment. To speed up the experiment, I reuse the atom and
    do OP after say 5 excitation attempts and try excitation again. I repeat the loop as long as
    there is an atom measured every atom_check_every_n.
    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection of a photon (either from SPCM0 or 1) proceeds with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- Blow away F=2 and measure retention
    7- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    # AllSPCMs_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        # AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        # SPCM0_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        # SPCM1_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)

        self.core.break_realtime()
        # delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(1000 * us)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

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
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off() ### alowing D1 to pass
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on() ### blocking D1
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
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_11_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(200 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    delay(20 * us)

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
                delay(1 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)


            ############################# readout every self.atom_check_every_n to see if the atom survived
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(1 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
                delay(10 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

                delay(10*us)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
        # for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
        #     self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        # delay(1 * ms)
        # for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(15 * ms)

@kernel
def atom_photon_parity_3_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Similar to parity_2_experiment, but with MW+RF pi/2 pulse to rotate to x-basis.

    To speed up the experiment, I reuse the atom and
    do OP after say 5 excitation attempts and try excitation again. I repeat the loop as long as
    there is an atom measured every atom_check_every_n.
    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    # self.dds_MW_RF.set_phase_mode(PHASE_MODE_TRACKING)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    delay(1 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        self.core.break_realtime()
        # delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(1000 * us)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

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
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

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
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
            delay(2*us)

            self.n_excitation_attempts = 1

            for excitation_attempt in range(self.n_excitation_attempts):
                delay(20 * us)
                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(20 * us)
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

                delay(2 * us) ### 1us is not enough.
                self.ttl_microwave_switch.off()
                delay(self.t_microwave_11_pulse)
                self.ttl_microwave_switch.on()

                delay(2*us)

                if self.t_MW_RF_pulse > 0:
                    self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    delay(5 * us)

                    with parallel:
                        self.ttl_microwave_switch.off()  ### turn on MW
                        self.dds_MW_RF.sw.on()  ### turn on RF

                    delay(self.t_MW_RF_pulse)

                    with parallel:
                        self.ttl_microwave_switch.on()  ### turn off MW
                        self.dds_MW_RF.sw.off()  ### turn off RF

                delay(2 * us)

                if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(20 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    delay(20 * us)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(100*us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(1 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)


            ############################# readout every self.atom_check_every_n to see if the atom survived
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(1 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
                delay(0.1 * ms)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

                delay(1*ms)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
        # for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
        #     self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        # delay(1 * ms)
        # for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def atom_photon_parity_4_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Atom-photon parity oscillation in X-basis. Compared to #3, here I am timetagging the photons and starting the MW pulses
    at a fixed time relative to photon detection.

    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    self.set_dataset("Slack_tracking", [0.0], broadcast=True)
    delay(1*ms)

    record_chopped_optical_pumping(self)
    delay(200*ms)

    record_chopped_blow_away(self)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switch
    self.ttl_microwave_switch.on() ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1] * self.n_measurements
    angle_780_QWP = [-1] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
        self.ttl_microwave_switch.on()  ### close the switch

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(10*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

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
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(10*us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)

                if SPCM0_click_time>0 and SPCM1_click_time<0:
                #     at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                #
                #     self.ttl_microwave_switch.off()
                #     delay(self.t_microwave_11_pulse)
                #     self.ttl_microwave_switch.on()
                #
                #     delay(2 * us)
                #
                #     if self.t_MW_RF_pulse > 0:
                #         self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                #         delay(5 * us)
                #
                #         with parallel:
                #             self.ttl_microwave_switch.off()  ### turn on MW
                #             self.dds_MW_RF.sw.on()  ### turn on RF
                #
                #         delay(self.t_MW_RF_pulse)
                #
                #         with parallel:
                #             self.ttl_microwave_switch.on()  ### turn off MW
                #             self.dds_MW_RF.sw.off()  ### turn off RF

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase = 0.0)

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase = 0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)



                    # self.core.break_realtime()
                    # self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    # self.append_to_dataset('SPCM0_SinglePhoton', 1)
                    # self.append_to_dataset('SPCM1_SinglePhoton', 0)
                    # self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    # self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    # delay(20 * us)
                    # self.core.break_realtime()


                    delay(1*ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time<0 and SPCM1_click_time>0:
                    # at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    #
                    # self.ttl_microwave_switch.off()
                    # delay(self.t_microwave_11_pulse)
                    # self.ttl_microwave_switch.on()
                    #
                    # delay(2 * us)
                    #
                    # if self.t_MW_RF_pulse > 0:
                    #     self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    #     delay(5 * us)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.off()  ### turn on MW
                    #         self.dds_MW_RF.sw.on()  ### turn on RF
                    #
                    #     delay(self.t_MW_RF_pulse)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.on()  ### turn off MW
                    #         self.dds_MW_RF.sw.off()  ### turn off RF

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)



                    # self.core.break_realtime()
                    # self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    # self.append_to_dataset('SPCM0_SinglePhoton', 0)
                    # self.append_to_dataset('SPCM1_SinglePhoton', 1)
                    # self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    # self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    # delay(20 * us)
                    # self.core.break_realtime()

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(20*us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            ####################################### atom check
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(1 * ms)

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.off()
            delay(1 * us)
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
            delay(0.1 * ms)

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

            SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

            delay(1 * ms)

            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            delay(1 * us)
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(1 * us)

            ### stopping the excitation cycle after the atom is lost
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            else:
                atom_loaded = False

            ################################### atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(1 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

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

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    # delay(15 * ms)
    self.core.break_realtime()
    for i in range(self.n_measurements):
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_5_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Basically the same as parity_4_experiment with minor differences/improvements

    Atom-photon parity oscillation in X-basis. Compared to #3, here I am timetagging the photons and starting the MW pulses
    at a fixed time relative to photon detection.

    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment


    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_optical_pumping(self)
    delay(100 * ms)

    record_chopped_blow_away(self)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1] * self.n_measurements
    angle_780_QWP = [-1] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
        self.ttl_microwave_switch.on()  ### close the switch

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment freezes and does not advance after a while!!
            ### The problem is in get_handle in the function that messes up with RTIO timeline in the loop.
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)
            #     self.GRIN1and2_dds.sw.on() ### turning back on for excitation after chopped_optical_pumping

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(10 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(10 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(10 * us)

                ### Tunring on pumping RP:
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(10 * us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)

                if SPCM0_click_time > 0 and SPCM1_click_time < 0:
                    #     at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    #
                    #     self.ttl_microwave_switch.off()
                    #     delay(self.t_microwave_11_pulse)
                    #     self.ttl_microwave_switch.on()
                    #
                    #     delay(2 * us)
                    #
                    #     if self.t_MW_RF_pulse > 0:
                    #         self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    #         delay(5 * us)
                    #
                    #         with parallel:
                    #             self.ttl_microwave_switch.off()  ### turn on MW
                    #             self.dds_MW_RF.sw.on()  ### turn on RF
                    #
                    #         delay(self.t_MW_RF_pulse)
                    #
                    #         with parallel:
                    #             self.ttl_microwave_switch.on()  ### turn off MW
                    #             self.dds_MW_RF.sw.off()  ### turn off RF

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0:
                    # at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    #
                    # self.ttl_microwave_switch.off()
                    # delay(self.t_microwave_11_pulse)
                    # self.ttl_microwave_switch.on()
                    #
                    # delay(2 * us)
                    #
                    # if self.t_MW_RF_pulse > 0:
                    #     self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    #     delay(5 * us)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.off()  ### turn on MW
                    #         self.dds_MW_RF.sw.on()  ### turn on RF
                    #
                    #     delay(self.t_MW_RF_pulse)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.on()  ### turn off MW
                    #         self.dds_MW_RF.sw.off()  ### turn off RF

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            ####################################### atom check
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(1 * ms)

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.off()
            delay(10 * us)
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            if not self.PGC_and_RO_with_on_chip_beams:
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
            delay(10 * us)

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

            SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

            delay(1 * ms)

            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            delay(10 * us)
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(10 * us)

            ### stopping the excitation cycle after the atom is lost
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            else:
                atom_loaded = False

            ################################### atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(1 * ms)  ### coils relaxation time

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

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

            excitation_cycle += 1

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1 * ms)
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

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    # delay(15 * ms)
    self.core.break_realtime()
    for i in range(self.n_measurements):
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_6_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Trying to reduce the delay after photon detection. Do we need to set all the phases to 0.0?

    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    # record_chopped_optical_pumping(self)
    # delay(100 * ms)

    record_chopped_blow_away(self)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    # delay(10 * ms)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1] * self.n_measurements
    angle_780_QWP = [-1] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
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
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            # ############################### chopped optical pumping phase - pumps atoms into F=1,m_F=0
            # ### strange that with chopped_optical_pumping function, the experiment freezes and does not advance after a while!!
            # ### The problem is in get_handle in the function that messes up with RTIO timeline in the loop.
            # # if self.t_pumping > 0.0:
            # #     chopped_optical_pumping(self)
            # #     delay(1 * ms)
            # #     self.GRIN1and2_dds.sw.on() ### turning back on for excitation after chopped_optical_pumping
            #
            # if self.t_pumping > 0.0:
            #     self.ttl_repump_switch.on()  # turns off the MOT RP AOM
            #     self.ttl_exc0_switch.on()  # turns off the excitation
            #     self.dds_cooling_DP.sw.off()  # no cooling light
            #     delay(10 * us)
            #     self.dds_AOM_A1.sw.off()
            #     self.dds_AOM_A2.sw.off()
            #     self.dds_AOM_A3.sw.off()
            #     self.dds_AOM_A4.sw.off()
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(10 * us)
            #
            #     ### set coils for pumping
            #     self.zotino0.set_dac(
            #         [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            #         channels=self.coil_channels)
            #     delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.
            #
            #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            #     delay(10 * us)
            #
            #     ### Tunring on pumping RP:
            #     self.ttl_pumping_repump_switch.off()
            #     self.dds_AOM_A5.sw.on()
            #     self.dds_AOM_A6.sw.on()
            #
            #     delay(1 * ms)
            #
            #     self.ttl_GRIN1_switch.off()
            #     delay(10 * us)
            #
            #     self.core_dma.playback_handle(op_dma_handle)
            #     delay(self.t_depumping)
            #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
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
            #     self.ttl_GRIN1_switch.on()
            #     delay(10 * us)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            # FORT_ramp2_smoothstep(self, direction="down")
            # delay(10 * us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(20 * us)
                # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation
                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                at_mu(t1 + int(self.gate_start_offset_mu))

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)

                if SPCM0_click_time > 0 and SPCM1_click_time < 0:
                    #     at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    #
                    #     self.ttl_microwave_switch.off()
                    #     delay(self.t_microwave_11_pulse)
                    #     self.ttl_microwave_switch.on()
                    #
                    #     delay(2 * us)
                    #
                    #     if self.t_MW_RF_pulse > 0:
                    #         self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    #         delay(5 * us)
                    #
                    #         with parallel:
                    #             self.ttl_microwave_switch.off()  ### turn on MW
                    #             self.dds_MW_RF.sw.on()  ### turn on RF
                    #
                    #         delay(self.t_MW_RF_pulse)
                    #
                    #         with parallel:
                    #             self.ttl_microwave_switch.on()  ### turn off MW
                    #             self.dds_MW_RF.sw.off()  ### turn off RF

                    # at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(20 * us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    ################################### atom cooling phase with PGC settings
                    if self.t_recooling > 0:
                        self.zotino0.set_dac(
                            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                            channels=self.coil_channels)

                        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                        delay(1 * ms)  ### coils relaxation time

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
                        else:
                            self.dds_AOM_A5.sw.off()
                            self.dds_AOM_A6.sw.off()

                        delay(self.t_recooling)

                        self.dds_cooling_DP.sw.off()
                        self.ttl_repump_switch.on()
                        delay(1 * us)
                        self.dds_AOM_A1.sw.off()
                        self.dds_AOM_A2.sw.off()
                        self.dds_AOM_A3.sw.off()
                        self.dds_AOM_A4.sw.off()
                        self.dds_AOM_A5.sw.off()
                        self.dds_AOM_A6.sw.off()

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0:
                    # at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    #
                    # self.ttl_microwave_switch.off()
                    # delay(self.t_microwave_11_pulse)
                    # self.ttl_microwave_switch.on()
                    #
                    # delay(2 * us)
                    #
                    # if self.t_MW_RF_pulse > 0:
                    #     self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    #     delay(5 * us)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.off()  ### turn on MW
                    #         self.dds_MW_RF.sw.on()  ### turn on RF
                    #
                    #     delay(self.t_MW_RF_pulse)
                    #
                    #     with parallel:
                    #         self.ttl_microwave_switch.on()  ### turn off MW
                    #         self.dds_MW_RF.sw.off()  ### turn off RF

                    # at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(20*us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    ################################### atom cooling phase with PGC settings
                    if self.t_recooling > 0:
                        self.zotino0.set_dac(
                            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                            channels=self.coil_channels)

                        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                        delay(1 * ms)  ### coils relaxation time

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
                        else:
                            self.dds_AOM_A5.sw.off()
                            self.dds_AOM_A6.sw.off()

                        delay(self.t_recooling)

                        self.dds_cooling_DP.sw.off()
                        self.ttl_repump_switch.on()
                        delay(1 * us)
                        self.dds_AOM_A1.sw.off()
                        self.dds_AOM_A2.sw.off()
                        self.dds_AOM_A3.sw.off()
                        self.dds_AOM_A4.sw.off()
                        self.dds_AOM_A5.sw.off()
                        self.dds_AOM_A6.sw.off()

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            # FORT_ramp2_smoothstep(self, direction="up")

            ####################################### atom check
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(1 * ms)

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.off()
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

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

            SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

            delay(10 * us)

            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            delay(10 * us)
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(10 * us)

            ### stopping the excitation cycle after the atom is lost
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            else:
                atom_loaded = False

            ################################### atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(1 * ms)  ### coils relaxation time

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

            excitation_cycle += 1

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1 * ms)
        # self.GRIN1and2_dds.sw.off()

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
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_6_node2_AllSPCM_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Trying to reduce the delay after photon detection. Do we need to set all the phases to 0.0?

    """

    self.core.reset()
    delay(100 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_blow_away(self)
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    # ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    # delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.off()  ### turns off the DDS
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.off()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    SPCM0_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1.0] * self.n_measurements
    angle_780_QWP = [-1.0] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()
        delay(100 * ms)

        self.ttl_exc0_switch.on()  # turns off the excitation
        self.dds_microwaves.sw.off()  ### turns off the DDS not the switch
        self.ttl_microwave_switch.on()  ### close the switch
        delay(10 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        ### 1) use ttl_D1_pumping/ttl_excitation to send D1/Exc to GRIN1/GRIN2
        ### 2) use ttl_GRIN1_switch/ttl_GRIN2_switch to send D1/EXC to atoms

        # self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        # self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        # self.GRIN1and2_dds.sw.on()      # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        # self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        # delay(10*ms)
        self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5*us)
        self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()      # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)
        self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:
            delay(1 *ms)
            # delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment freezes and does not advance after a while!!
            ### The problem is in get_handle in the function that messes up with RTIO timeline in the loop.
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)
            #     self.GRIN1and2_dds.sw.on() ### turning back on for excitation after chopped_optical_pumping

            if self.t_pumping > 0.0:
                CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(5*us)
            self.core.break_realtime()

            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(100 * us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()
                # todo: node2- change this to turning on the ttl;
                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                # #### node2
                # FORT_offset = 90   ## as of 2026/03/08 Excelitas SPCMs
                #
                # t1 = now_mu()
                # at_mu(t1 + int(self.t_excitation_offset_mu))
                # self.ttl_GRIN2_switch.off()  # turns on excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                # self.ttl_GRIN2_switch.on()  # turns off excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # self.dds_FORT.sw.off()  ### turns FORT off
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset + int(self.t_photon_collection_time / ns))
                # self.dds_FORT.sw.on()  ### turns FORT on
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.gate_start_offset_mu))  # turn the gate on simultaneously with the EXC pulse
                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on 50ns after photon collection ended

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation after t_excitation_offset_mu

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))  ## turns on gate after gate_start_offset_mu

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)

                # self.ttl_microwave_switch.off()
                # delay(5*us)
                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()   #### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()  #### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF


                    delay(1*ms)
                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(100 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()   #### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()  #### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF


                    delay(1*ms)
                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(100 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:

                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds,
                                                amplitude=dB_to_V(self.p_microwaves))

                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(
                            self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                # print("what happens here")
            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            ####################################### atom check
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(1 * ms)

            self.dds_cooling_DP.sw.on()
            self.ttl_repump_switch.off()
            delay(10 * us)
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            if not self.PGC_and_RO_with_on_chip_beams:
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
            delay(10 * us)

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)

            # SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            # BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
            AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                      self.ttl_SPCM1_counter.fetch_count() + \
                                      self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                      self.ttl_SPCM1_OtherNode_counter.fetch_count())


            delay(1 * ms)

            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            delay(10 * us)
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            delay(10 * us)

            ### stopping the excitation cycle after the atom is lost
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True

            else:
                atom_loaded = False

            ################################### atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(1 * ms)  ### coils relaxation time

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

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

            excitation_cycle += 1

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(10 * ms)
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

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(50 * ms)
    self.core.break_realtime()
    for i in range(self.n_measurements):
        delay(10*ms)
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton[i])
        self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])

    delay(50 * ms)

@kernel
def atom_photon_parity_7_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    In parity_6_experiment, we lose most of the atoms in blowaway phase. In this experiment
    I am trying to debug this. This experiment does not care about the photon detection
    to run faster.

    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1] * self.n_measurements
    angle_780_QWP = [-1] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        delay(1 * ms)
        self.ttl_microwave_switch.on()  ### close the switch
        delay(20 * us)
        self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

        first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation
        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # for excitation_attempt in range(self.n_excitation_attempts):
            slack = now_mu() - self.core.get_rtio_counter_mu()
            if slack < 1e5:
                # self.print_async("slack added in measurement:", self.measurement)
                self.core.break_realtime()

            self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
            delay(5 * us)
            # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            t1 = now_mu()

            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation
            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
            self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

            at_mu(t1 + int(self.gate_start_offset_mu))

            with parallel:
                t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

            SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
            SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)


            delay(self.t_start_MW_mapping_mu*ns)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

            # delay(self.t_start_MW_mapping_mu*ns + 1000*ns)
            # self.ttl_microwave_switch.off()
            # delay(self.t_start_MW_mapping_mu*ns + 1000*ns + self.t_microwave_11_pulse)
            # self.ttl_microwave_switch.on()

            ############################ blow-away phase - push out atoms in F=2 only
            delay(20 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            delay(10 * us)
            self.core_dma.playback_handle(ba_dma_handle)

            delay(10 * us)
            atom_parity_shot(self)

            delay(1 * ms)
            AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
            SPCM0_SinglePhoton[self.measurement] = 1.0
            SPCM1_SinglePhoton[self.measurement] = 1.0
            angle_780_HWP[self.measurement] = self.target_780_HWP
            angle_780_QWP[self.measurement] = self.target_780_QWP
            delay(1 * ms)

            self.measurement += 1

            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            FORT_ramp2_smoothstep(self, direction="up")

            ### stopping the excitation cycle after the atom is lost
            if self.AllSPCMs_parity_RO / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                # atom_loaded = True
                atom_loaded = False

            else:
                atom_loaded = False

            # ################################### atom cooling phase with PGC settings
            # if self.t_recooling > 0:
            #     self.zotino0.set_dac(
            #         [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            #         channels=self.coil_channels)
            #
            #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
            #     delay(1 * ms)  ### coils relaxation time
            #
            #     self.dds_cooling_DP.sw.on()
            #     self.ttl_repump_switch.off()
            #     delay(1 * us)
            #     self.dds_AOM_A1.sw.on()
            #     self.dds_AOM_A2.sw.on()
            #     self.dds_AOM_A3.sw.on()
            #     self.dds_AOM_A4.sw.on()
            #     if not self.PGC_and_RO_with_on_chip_beams:
            #         self.dds_AOM_A5.sw.on()
            #         self.dds_AOM_A6.sw.on()
            #     else:
            #         self.dds_AOM_A5.sw.off()
            #         self.dds_AOM_A6.sw.off()
            #
            #     delay(self.t_recooling)
            #
            #     self.dds_cooling_DP.sw.off()
            #     self.ttl_repump_switch.on()
            #     delay(1 * us)
            #     self.dds_AOM_A1.sw.off()
            #     self.dds_AOM_A2.sw.off()
            #     self.dds_AOM_A3.sw.off()
            #     self.dds_AOM_A4.sw.off()
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(1 * us)

            excitation_cycle += 1

        # delay(20 * us)
        # self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1 * ms)
        # self.GRIN1and2_dds.sw.off()

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
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_8_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    this is yet another debugging experiment similar to atom_photon_parity_7_experiment, but based on Rabi_2
    with added excitation. The goal is to see why blowaway removes the atoms from the trap after excitation.
    All atoms should be in f=1 manifold and should survive the blowaway.
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

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    self.core.break_realtime()

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    # delay(1 * ms)
    # move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    # move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    # delay(10 * ms)
    # self.core.reset()

    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.ttl_GRIN2_switch.on()  # turns off excitation switch
    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
    delay(5 * us)
    self.GRIN1and2_dds.sw.on()

    self.measurement = 0
    while self.measurement < self.n_measurements:
        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

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

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay (10 * us)
            if self.which_node == "alice":
                CW_optical_pumping_node1(self)
            else:
                CW_optical_pumping_node2(self)
            delay(10*us)

        ############################
        # excitation phase
        ############################

        self.GRIN1and2_dds.sw.on() ### require this here to turn on the channel. It will be off without this.

        self.ttl_GRIN2_switch.on()
        delay(5*us)

        # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
        # delay(1*us)
        # delay(self.t_delay_between_shots)
        # self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

        self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

        self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
        delay(5 * us)
        t1 = now_mu()

        self.dds_FORT.sw.off()  ### turns FORT off

        at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
        self.dds_FORT.sw.on()  ### turns FORT on

        at_mu(t1 + int(self.t_excitation_offset_mu))
        self.ttl_GRIN2_switch.off()  # turns on excitation

        at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
        self.ttl_GRIN2_switch.on()  # turns off excitation
        at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
        self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

        ############################
        # microwave phase
        ############################
        delay(10*us)
        if self.t_microwave_pulse > 0.0:
            # ### Changing the bias field
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            # FORT_ramp2_smoothstep(self, direction="down")
            delay(1 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(20 * us)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            # FORT_ramp2_smoothstep(self, direction="up")
            delay(10 * us)

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

    delay(1 * ms)
    self.dds_microwaves.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def atom_photon_parity_9_AllSPCM_experiment(self):
    """
    *** compatible with either node1 or node2

    After optimizing the sequence using parity_7 and parity_8, this experiment is an improved version of
    parity_6 experiment. I have tried to minimize atom loss during the process.

    In this experiment, I don't use while atom_loaded. This is already being checked in atom loading function.
    For some reason, this reduces atom loss.

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

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

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

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

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

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

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

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                # delay(5 * us)
                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                delay(5 * us)

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns) + 1000)
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                at_mu(t1 + int(self.gate_start_offset_mu))

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)

                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10 * us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10*us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10 * us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10*us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break


            delay(20 * us)
            excitation_cycle += 1

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
    delay(50 * ms)


@kernel
def atom_photon_parity_10_AllSPCM_experiment(self):
    """
    *** compatible with either node1 or node2
    *** DMA

    After optimizing the sequence using parity_7 and parity_8, this experiment is an improved version of
    parity_6 experiment. I have tried to minimize atom loss during the process.

    In this experiment, I don't use while atom_loaded. This is already being checked in atom loading function.
    For some reason, this reduces atom loss.

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
    record_CW_optical_pumping_node2(self)
    # record_recooling(self)
    record_excitation_and_photon_collection(self)
    record_parity_mapping(self)

    #### getting DMA handles
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
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

        load_until_atom_smooth_FORT_recycle(self)
        # load_MOT_and_FORT_until_atom_recycle(self)

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

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                    delay(10 * us)
                else:
                    # CW_optical_pumping_node2(self)
                    self.core_dma.playback_handle(CW_OP_node2_handle)
            t_after_OP_dma = now_mu()

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

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

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
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

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

#
# ###############################################################################
# # 4. TWO NODES EXPERIMENT FUNCTIONS WITH TTL SYNCHRONIZATION
# # These are functions that are used for two nodes experiment with TTL
# # symnchronization. Working as of 20260810. Soon will be deprecated after
# # upgrading to master-satellite configuration. - Eunji
# ###############################################################################

@kernel
def two_nodes_synchronization(self):
    """
    1. Initial coarse synchronization(Two-way communication)
    - ttl_node1_output1 and ttl_node2_output1 turned ON
    - IF both nodes receives this signal, move on to finer synchronization
    2. Fine synchronization(One-way communication)
    - ttl_node1_output2 is pulsed @ t_node1_ref_0
    - node2 receives this siganl and time tags @ t_node2_ref_0
    3. Sets each node's time cursor to 150us after this.
    """
    # ############################################################
    # # Starting synchronization
    # # Signal to the other node.
    # ############################################################
    other_node_ready = False
    readout = 0
    # sync_starting_at = now_mu()

    # self.core.break_realtime()

    delay(10 * us)

    if self.which_node == "alice":
        self.ttl_node1_output1.on()
    else:
        self.ttl_node2_output1.on()

    delay(1 * us)

    ############################################################
    # 1. Initial coarse synchronization
    # Wait until the other node is also ready.
    # Check every wait_time
    ############################################################
    while not other_node_ready:
        wait_time = 4*us
        delay(wait_time)
        if self.which_node == "alice":
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node2_input1.sample_get())
            ## this can be replaced with sample_get_nonrt; but it has core.break_realtime() so needs to be careful

            if readout == 1:
                other_node_ready = True
            else:
                other_node_ready = False

        else:
            self.ttl_node1_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node1_input1.sample_get())

            if readout == 1:
                other_node_ready = True
            else:
                other_node_ready = False

    delay(10 * us)

    ############################################################
    # 2. Fine synchronization
    # Synchronizing with time tag
    ############################################################
    #todo: check if these timings can be reduced or needs to be increased
    #      now it takes 221us for the entire syncrhonization to be done;
    ##### Syncing by sending a TTL pulse from node1 to node2
    t_node1_ref_0 = -1
    t_node2_ref_0 = -1

    if self.which_node == "alice":
        delay(50*us)
        # delay(15 * us)
        t_node1_ref_0 = now_mu()
        self.ttl_node1_output2.pulse(1*us)
    else:
        t_gate_end = self.ttl_node1_input2.gate_rising(150*us)
        # t_gate_end = self.ttl_node1_input2.gate_rising(30 * us)
        t_node2_ref_0 = self.ttl_node1_input2.timestamp_mu(t_gate_end)

        if t_node2_ref_0 == -1:
            self.print_async("Node2 did not detect pulse from Node1")

    delay(10 * us)
    t_node1_ref_1 = t_node1_ref_0 + self.core.seconds_to_mu(200 * us) + 200
    t_node2_ref_1 = t_node2_ref_0 + self.core.seconds_to_mu(200 * us)

    if self.which_node == "alice":
        at_mu(t_node1_ref_1)
        self.ttl_node1_output1.off()
        # sync_time_took = int((t_node1_ref_1 - sync_starting_at) / 1e3)
    else:
        at_mu(t_node2_ref_1)
        self.ttl_node2_output1.off()
        # sync_time_took = int((t_node2_ref_1 - sync_starting_at) / 1e3)

    delay(10 * us)

    # sync_time_took = int((t_node1_ref_1 - sync_starting_at) / 1e3)
    # self.print_async("sync_time_took :", sync_time_took, " us")
    # self.append_to_dataset("sync_time_took", sync_time_took)

@kernel
def two_nodes_synchronization2(self) -> TBool:
    t_node1_ref_0 = -1
    t_node2_ref_0 = -1

    poll_period = 1 * ms
    timeout_cycles = 10000     # 10 s if poll_period = 1 ms and timeout_cycles = 10000
    ack_delay_mu = self.core.seconds_to_mu(50 * us)
    common_delay_mu = self.core.seconds_to_mu(2 * ms)

    alice_extra_offset_mu = 200  ### Calibrated this on the scope.

    self.core.break_realtime()
    delay(100 * us)

    self.print_async("******* Syncing the nodes **********")

    if self.which_node == "bob":

        #### Open Bob's input gate first.
        t_gate_start = now_mu() + self.core.seconds_to_mu(100 * us)
        at_mu(t_gate_start)
        t_gate_end = self.ttl_node1_input1.gate_rising(10.01 * s)

        ### Now assert READY while the gate is open.
        ### This is scheduled after gate opening but before calling timestamp_mu().
        at_mu(t_gate_start + self.core.seconds_to_mu(10 * us))
        self.ttl_node2_output1.on()

        ### Now wait for Alice's reference pulse.
        t_node2_ref_0 = self.ttl_node1_input1.timestamp_mu(t_gate_end)

        if t_node2_ref_0 < 0:
            self.core.break_realtime()
            self.ttl_node2_output1.off()
            return False

        ### Drop READY after successful sync.
        self.core.break_realtime()
        self.ttl_node2_output1.off()

        ### Put Bob at the common future reference.
        at_mu(t_node2_ref_0 + common_delay_mu)
        return True

    else:
        ### Alice polls Bob's READY level.
        ready = False
        for i in range(timeout_cycles):
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            if int(self.ttl_node2_input1.sample_get()):
                ready = True
                break
            delay(poll_period)

        if not ready:
            self.core.break_realtime()
            return False

        ### Alice sends the reference pulse at a deterministic future time.
        t_node1_ref_0 = now_mu() + ack_delay_mu
        at_mu(t_node1_ref_0)
        self.ttl_node1_output1.pulse(500 * ns)

        ### Put Alice at the corresponding future reference.
        at_mu(t_node1_ref_0 + common_delay_mu + alice_extra_offset_mu)
        return True

@kernel
def two_nodes_synchronization_trial(self) -> TBool:
    """
    Faisal's code based on two_nodes_synchronization2
    """
    delay(1*ms)

    t_node1_ref_0 = np.int64(-1)
    t_node2_ref_0 = np.int64(-1)

    poll_period = 1 * ms
    timeout_cycles = 10000     # 10 s if poll_period = 1 ms and timeout_cycles = 10000
    ack_delay_mu = self.core.seconds_to_mu(50 * us)
    common_delay_mu = self.core.seconds_to_mu(2 * ms)

    alice_extra_offset_mu = np.int64(200)  ### Calibrated this on the scope.

    self.core.break_realtime()

    delay(100 * us)
    self.print_async("******* Syncing the nodes **********")

    if self.which_node == "bob":

        #### Open Bob's input gate first.
        t_gate_start = now_mu() + self.core.seconds_to_mu(100 * us)
        gate_duration_mu = self.core.seconds_to_mu(10*s)
        t_gate_end = t_gate_start + gate_duration_mu


        at_mu(t_gate_start)
        with parallel:
            self.ttl_node1_input1.gate_rising_mu(gate_duration_mu)
            with sequential:
                delay(10*us)
                self.ttl_node2_output1.on()


        ### Now wait for Alice's reference pulse.
        t_node2_ref_0 = self.ttl_node1_input1.timestamp_mu(t_gate_end)

        if t_node2_ref_0 < 0:
            self.core.break_realtime()
            self.ttl_node2_output1.off()
            return False

        ### Drop READY after successful sync.
        self.core.break_realtime()
        self.ttl_node2_output1.off()

        ### Put Bob at the common future reference.
        at_mu(t_node2_ref_0 + common_delay_mu)
        return True

    else:

        self.ttl_node2_input1.count(now_mu())

        ### Alice polls Bob's READY level.
        ready = False
        for i in range(timeout_cycles):
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            if int(self.ttl_node2_input1.sample_get()):
                ready = True
                break
            delay(poll_period)

        if not ready:
            self.core.break_realtime()
            return False

        ### Alice sends the reference pulse at a deterministic future time.
        t_node1_ref_0 = now_mu() + ack_delay_mu
        at_mu(t_node1_ref_0)
        self.ttl_node1_output1.pulse(500 * ns)

        ### Put Alice at the corresponding future reference.
        at_mu(t_node1_ref_0 + common_delay_mu + alice_extra_offset_mu)
        return True

@kernel
def load_atom_in_both_nodes_individually(self):
    """
    Eunji - this is assuming that both nodes can try atom loading/readout simultaneously without distruption
    - this maybe useful when we integrate a fast shutter most likeley a MEMS switch.
    - but we need to think more carefully how much loss it has with coupling/mode matching/pol dependence
    - if we need to use MEMS switch, it might need to be spliced.

    Load/recycle atoms until both nodes have atoms simultaneously.

    Assumption:
        This node's TTL output is a live atom-status indicator.

            HIGH = this node currently has an atom
            LOW  = this node currently does not have an atom

        The other node's TTL input is used to determine whether the other node
        currently has an atom.

    This function exits only when:
        this node has atom AND the other node has atom.

    Datasets used:
        time_without_atom_other_node:
            Live-updated time since the other node was last seen without atom.

        Atom_loading_time_other_node:
            Appended once when the other node becomes ready.
    """

    atom_loaded = False
    other_node_has_atom = False
    both_nodes_have_atom = False
    readout = 0

    ############################################################
    # Timer for how long the OTHER node has been without atom.
    ############################################################
    t_before_other_atom = now_mu()
    t_after_other_atom = now_mu()
    t_no_other_atom = now_mu()

    time_without_other_atom = 0.0
    self.atom_loading_time_other_node = 0.0
    other_node_loading_time_recorded = False
    other_node_was_without_atom = False

    delay(100 * us)

    ############################################################
    # Initial check for THIS node using latest RO2 information.
    ############################################################
    if self.measurement > 0:
        if self.AllSPCMs_RO2 / self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True
            self.print_async("This node starts with atom.")
        else:
            atom_loaded = False
            self.print_async("This node starts without atom.")
    else:
        atom_loaded = False
        self.print_async("No previous measurement. This node needs loading.")

    ############################################################
    # Update this node TTL according to the actual initial status.
    ############################################################
    if atom_loaded:
        if self.which_node == "alice":
            self.ttl_node1_output1.on()
            self.print_async("Node1 atom TTL is ON.")
        else:
            self.ttl_node2_output1.on()
            self.print_async("Node2 atom TTL is ON.")
    else:
        if self.which_node == "alice":
            self.ttl_node1_output1.off()
            self.print_async("Node1 atom TTL is OFF.")
        else:
            self.ttl_node2_output1.off()
            self.print_async("Node2 atom TTL is OFF.")

    delay(100 * us)

    ############################################################
    # Main loop.
    # Only exits when this node has atom AND other node has atom.
    ############################################################
    while not (atom_loaded and both_nodes_have_atom):

        ########################################################
        # If this node does not have atom, load/reload it.
        ########################################################
        if not atom_loaded:
            self.print_async("Loading atom in this node...")

            load_MOT_and_FORT_until_atom_recycle(self)

            atom_loaded = True
            self.print_async("Atom loaded in this node.")

            if self.which_node == "alice":
                self.ttl_node1_output1.on()
                self.print_async("Node1 atom TTL is ON after loading.")
            else:
                self.ttl_node2_output1.on()
                self.print_async("Node2 atom TTL is ON after loading.")

            delay(100 * us)

        ########################################################
        # Check whether the OTHER node has atom.
        ########################################################
        if self.which_node == "alice":
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node2_input1.sample_get())

            if readout == 1:
                other_node_has_atom = True
                self.print_async("Node2 has atom.")
            else:
                other_node_has_atom = False
                self.print_async("Node2 does NOT have atom.")

        else:
            self.ttl_node1_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node1_input1.sample_get())

            if readout == 1:
                other_node_has_atom = True
                self.print_async("Node1 has atom.")
            else:
                other_node_has_atom = False
                self.print_async("Node1 does NOT have atom.")

        ########################################################
        # Track how long the OTHER node has been without atom.
        ########################################################
        if other_node_has_atom:
            self.set_dataset("time_without_atom_other_node", 0.0, broadcast=True)

            if other_node_was_without_atom and not other_node_loading_time_recorded:
                t_after_other_atom = now_mu()

                self.atom_loading_time_other_node = self.core.mu_to_seconds(
                    t_after_other_atom - t_before_other_atom
                )

                self.append_to_dataset(
                    "Atom_loading_time_other_node",
                    self.atom_loading_time_other_node
                )

                other_node_loading_time_recorded = True
                other_node_was_without_atom = False

                self.print_async("Other node loading time recorded.")

        else:
            ####################################################
            # If this is the first time seeing the other node
            # without atom, start the timer.
            ####################################################
            if not other_node_was_without_atom:
                t_before_other_atom = now_mu()
                other_node_was_without_atom = True
                other_node_loading_time_recorded = False

            delay(0.1 * ms)

            t_no_other_atom = now_mu()

            time_without_other_atom = self.core.mu_to_seconds(
                t_no_other_atom - t_before_other_atom
            )

            self.set_dataset(
                "time_without_atom_other_node",
                time_without_other_atom,
                broadcast=True
            )

        ########################################################
        # If both nodes have atoms, exit the loop.
        ########################################################
        if atom_loaded and other_node_has_atom:
            both_nodes_have_atom = True
            self.print_async("Both nodes have atoms. Exiting loading loop.")
            break
        else:
            both_nodes_have_atom = False

        ########################################################
        # Other node does not have atom yet.
        # Wait 200 ms, then check whether THIS node still has atom.
        ########################################################
        delay(200 * ms)

        ########################################################
        # Read out this node's current atom using second_shot.
        ########################################################
        second_shot(self)

        ########################################################
        # Check if this node still has atom.
        ########################################################
        if self.AllSPCMs_RO2 / self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True
            self.print_async("This node still has atom while waiting.")

            ####################################################
            # Keep this node TTL HIGH.
            ####################################################
            if self.which_node == "alice":
                self.ttl_node1_output1.on()
            else:
                self.ttl_node2_output1.on()

            ####################################################
            # Recool after the readout.
            ####################################################
            recooling_after_first_shot(self)

        else:
            atom_loaded = False
            both_nodes_have_atom = False
            self.print_async("This node lost atom while waiting.")

            ####################################################
            # This node no longer has atom, so set TTL LOW.
            ####################################################
            if self.which_node == "alice":
                self.ttl_node1_output1.off()
                self.print_async("Node1 atom TTL is OFF.")
            else:
                self.ttl_node2_output1.off()
                self.print_async("Node2 atom TTL is OFF.")

            ####################################################
            # Important:
            # Do NOT reset time_without_atom_other_node here.
            # That timer belongs to the other node's atom status,
            # not this node's atom status.
            ####################################################

            delay(100 * us)

    ############################################################
    # Function exits only here, when both nodes have atoms.
    ############################################################
    self.print_async("Confirmed: both nodes have atoms.")
    delay(100 * us)

@kernel
def load_until_atom_in_both_nodes_recycle(self):
    """
    Two-node combined-loading version.

    * based on load_MOT_and_FORT_until_atom_recycle
    * with synchronization added in the while loop
    * Checks for atoms in both nodes using two_atom_threshold_in_loading

    """
    ### First check if there is already an atom in the FORT based on RO2
    ### This will miss the first measurement of each iteration (except the first iteration)
    ### and results in the first RO1 to be less than threshold. But it is OK.
    delay(100 * us)
    if self.AllSPCMs_RO2/self.t_SPCM_second_shot > self.two_atom_threshold:
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
            delay(1 * ms)
            Sampler_test(self)
            delay(1 * ms)
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
        if self.which_node == "alice":
            self.zotino0.set_dac([3.5], self.UV_trig_channel)
        delay(1 * ms)

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        AllSPCMs_atom_check_loaded = 0  ### for initilization
        AllSPCMs_atom_check_not_loaded = 0

        shim_tune_runs = 0

        n_synchronization_in_loading = 0
        while True:
            two_nodes_synchronization(self)  ### Syncing the two nodes.

            n_synchronization_in_loading += 1

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

                if AllSPCMs_atom_check / atom_check_time > self.two_atom_threshold_for_loading:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True
                    AllSPCMs_atom_check_loaded = AllSPCMs_atom_check

            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True)  ### resetting time_without_atom when we load an atom
                t_after_atom = now_mu()

                ### just to check the histogram during atom loading to find a good two_atom_threshold_for_loading
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_loaded)
                self.append_to_dataset("AllSPCMs_atom_check_in_loading", AllSPCMs_atom_check_not_loaded)
                delay(1 * ms)

                # self.print_async("n_synchronization_in_loading: ", n_synchronization_in_loading)
                break  ### Exit the outer loop if an atom is loaded

            #### time_without_atom shows how long is passed from the previous atom loading. Calculated only when try_n > max_tries
            delay(0.1 * ms)
            t_no_atom = now_mu()
            time_without_atom = self.core.mu_to_seconds(t_no_atom - t_before_atom)
            self.set_dataset("time_without_atom", time_without_atom, broadcast=True)

            # ### If max_tries reached and atom loading is too bad, run shim tuning.
            # t_shim_tune_trigger = 10.0  # seconds
            # max_shim_tune_runs = 3  # prevents infinite tuning loop
            #
            # if self.tune_shims_when_loading_is_bad and (time_without_atom > t_shim_tune_trigger) and (shim_tune_runs < max_shim_tune_runs):
            #     self.print_async("Atom loading is bad. Tuning X and Y shims.")
            #     tune_shims_for_atom_loading(self)
            #     shim_tune_runs += 1
            #
            #     ### restart the loading attempt with the (possibly) improved shims
            #     try_n = 0
            #     t_before_atom = now_mu()
            #     delay(0.1 * ms)
            #     continue

            ### If max_tries reached and still no atom, run feedback
            if self.enable_laser_feedback:
                delay(0.1 * ms)  ### necessary to avoid underflow
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                run_feedback_and_record_FORT_MM_power(self, record_power=False)
                self.n_feedback_per_iteration += 1
                # bug -- microwave dds and FORT are off after AOM feedback; not clear why yet. for now, just turn them back on
                self.dds_microwaves.sw.on()
                self.dds_FORT.sw.on()
                delay(0.1 * ms)

                try_n = 0

        delay(1 * ms)
        if self.which_node == "alice":
            self.zotino0.set_dac([0.0], self.UV_trig_channel)
        delay(100 * us)

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

        ### saving the atom loading time for each loaded atom.
        self.atom_loading_time = self.core.mu_to_seconds(t_after_atom - t_before_atom)
        self.append_to_dataset("Atom_loading_time", self.atom_loading_time)
        delay(1 * ms)
        self.append_to_dataset("atom_loading_wall_clock", now_mu())  ### just to plot Atom_loading_time vs actual time in analysis
        self.n_atom_loaded_per_iteration += 1
        delay(10 * ms)

@kernel
def two_node_alternating_shot(self):
    """
    Alternating two-node readout. -Eunji as of 2026.07.14

    Alice windows: even windows
    Bob windows:   odd windows

    Each node gets:
        n_readout_windows_per_node windows
        1 ms counter gate per window

    Timing:
        Alice window starts at 0 ms
        Bob window starts 2.3 ms later
        Alice next window starts 4.6 ms later...

    ** Note that unlike other readout functions, I only save AllSPCM data rather than
        saving each dataset separately for each SPCM. Note that at some point if you
        want to implement this, alice and bob should have separate dataset.
        for example:
            SPCM0_alternating_RO_alice, SPCM0_alternating_RO_bob,
            SPCM1_alternating_RO_alice, SPCM1_alternating_RO_bob,
            SPCM0_OtherNode_alternating_RO_alice, SPCM0_OtherNode_alternating_RO_bob,
            SPCM1_OtherNode_alternating_RO_alice, SPCM1_OtherNode_alternating_RO_bob,
    """

    self.AllSPCMs_alternating_RO_alice = 0
    self.AllSPCMs_alternating_RO_bob = 0

    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(1 * ms)  ## coils relaxation time

    # Set cooling DP AOM to readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

    # Set FORT AOM to science settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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

    n_readout_windows_per_node = 10

    t_readout_window = 1.0 * ms

    self.ttl_repump_switch.on()  # turn off MOT RP
    self.dds_cooling_DP.sw.off()

    # 2*n_readout_windows_per_node because even windows are Alice, odd windows are Bob
    for i in range(2 * n_readout_windows_per_node):

        alice_window = (i % 2 == 0)

        # Only the active node turns on readout light.
        if alice_window:
            if self.which_node == "alice":
                self.ttl_repump_switch.off()   # turn on MOT RP
                self.dds_cooling_DP.sw.on()    # turn on cooling/readout=
        else:
            if self.which_node == "bob":
                self.ttl_repump_switch.off()   # turn on MOT RP
                self.dds_cooling_DP.sw.on()    # turn on cooling/readout

        # waiting for RO beams to turn on
        delay(0.1 * ms)

        # Open all counters during this readout window
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(t_readout_window)
            self.ttl_SPCM1_counter.gate_rising(t_readout_window)
            self.ttl_SPCM0_OtherNode_counter.gate_rising(t_readout_window)
            self.ttl_SPCM1_OtherNode_counter.gate_rising(t_readout_window)

        delay(0.1 * ms)
        self.ttl_repump_switch.on()  # turn off MOT RP
        self.dds_cooling_DP.sw.off()

        SPCM0_alternating_RO = self.ttl_SPCM0_counter.fetch_count()
        SPCM1_alternating_RO = self.ttl_SPCM1_counter.fetch_count()
        SPCM0_OtherNode_alternating_RO = self.ttl_SPCM0_OtherNode_counter.fetch_count()
        SPCM1_OtherNode_alternating_RO = self.ttl_SPCM1_OtherNode_counter.fetch_count()

        window_count = (SPCM0_alternating_RO + SPCM1_alternating_RO
                        + SPCM0_OtherNode_alternating_RO + SPCM1_OtherNode_alternating_RO)

        if alice_window:
            self.AllSPCMs_alternating_RO_alice += window_count
        else:
            self.AllSPCMs_alternating_RO_bob += window_count

        delay(0.1*ms)

@kernel
def Two_nodes_atom_loading_experiment(self):
    """
    Simple atom loading experiment in both nodes
    - checking for atoms in both nodes at the same time
    - based on atom_loading_2_experiment

    Sequence as of 2026.07.13
    1. Load atoms in both nodes simultaneously
    2. First shot
    3. FORT drop if > 0
    4. testing individual node shot - Alternating RO
    5. Second shot
    6. end_measurement
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        run_feedback_and_record_FORT_MM_power(self)

    ############################################################
    # Turn OFF this node's TTL to be ready for synchronization.
    ############################################################
    if self.which_node == "alice":
        self.ttl_node1_output1.off()
    else:
        self.ttl_node2_output1.off()
    delay(50 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        two_nodes_synchronization(self)

        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

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

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)

        two_nodes_synchronization(self)

        delay(1*ms)
        two_node_alternating_shot(self)
        delay(1 * ms)

        two_nodes_synchronization(self)

        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def Two_node_single_photon_experiment(self):
    """
    **** This is based on single_photon_experiment_3_atom_loading_advance
    **** Note that 5 is the latest version as of 2026-06-19
    **** But version 5 has not yet been tested in node2 so I will be using version 3

    **** This experiment is used to match the timing between the single photon arrival
         time from each node. Here, to do so inidividually, I physically block the excitation
         light of the other node while collecting the photon arrival histogram from one.

    **** It starts by loading atoms in both nodes, and then does all the excitation cycles
         in both nodes. This is to make sure that two nodes are running the same sequence
         and being synchronized to the next step correctly.

    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

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

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            # self.core.break_realtime()  ## Commenting this out to avoid timing error

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            # ############################### chopped optical pumping phase - pumps atoms into F=1,m_F=0
            # if self.t_pumping > 0.0:
            #
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
            #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0)) ### set to 5V for optical pumping
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            #     delay(1 * us)
            #
            #     ### Tunring on pumping RP:
            #     self.ttl_pumping_repump_switch.off()
            #     self.dds_AOM_A5.sw.on()
            #     self.dds_AOM_A6.sw.on()
            #
            #     # delay(1 * ms)
            #
            #     self.ttl_GRIN1_switch.off() ### was used when D1 was on GRIN1
            #     delay(10 * us)
            #
            #     self.core_dma.playback_handle(op_dma_handle)
            #     delay(self.t_depumping)
            #     self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
            #     self.ttl_pumping_repump_switch.on()  ### turning off pumping RP
            #
            #     delay(2 * us)
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(100 * us)
            #
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
            #     # delay(1 * ms)
            #
            #     self.ttl_GRIN1_switch.on() ### was used when D1 was on GRIN1
            #     delay(10 * us)

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            delay(1*ms)
            if self.t_pumping > 0.0:
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10*us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
            # delay(10*us)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(5 * us)

            self.core.break_realtime()

            two_nodes_synchronization(self)

            if self.which_node == "bob":
                delay_mu(int(self.t_delay_in_bob_mu))

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
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                # at_mu(t1 + 30000)
                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            #todo: check Node2 + timing
            self.core.break_realtime()
            two_nodes_synchronization(self)

            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

            self.core.break_realtime()
            two_nodes_synchronization(self)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                # delay(0.4*ms)

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())


                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.two_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

                # #### these added on 2025-09-09 based on Eunji's comment. I (AkS) did not have these in my photon experiments:
                # self.dds_cooling_DP.sw.off()
                # self.ttl_repump_switch.on()
                # delay(1 * us)
                # self.dds_AOM_A1.sw.off()
                # self.dds_AOM_A2.sw.off()
                # self.dds_AOM_A3.sw.off()
                # self.dds_AOM_A4.sw.off()
                # self.dds_AOM_A5.sw.off()
                # self.dds_AOM_A6.sw.off()
                # delay(1 * us)
                # ##############################

            self.core.break_realtime()
            two_nodes_synchronization(self)
            delay(10 * us)

        delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.1 * ms)

        two_nodes_synchronization(self)

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

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def Two_node_single_photon_2_experiment(self):
    """
    **** This is for Hong-Ou_Mandel Experiment
    **** This is based on Two_node_single_photon_experiment

    **** In HOM Measurement, we want to vary the excitation timing of one node to see
         see the HOM dip. However while we do that, we might accidentally keep the
         FORT from the other node on while one node is collecting photons, making
         false count measurements. Thus, here I added some timing changes in
         excitation cycle loop.

    ****
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    self.core.break_realtime()

    two_nodes_synchronization(self)
    delay(10*us)

    #### So that the other node is not turning the FORT on while photon collection in one node;
    t_diff_mu = int(self.t_delay_in_bob_mu - 189)

    if self.which_node == "alice":
        # t_diff = int(self.t_delay_in_bob_mu - 189)
        # t_FORT_OFF_mu = 0

        if t_diff_mu < 0:
            t_gate_mu = t_diff_mu
            t_photon_add_mu = -1 * t_diff_mu
            t_FORT_OFF_mu = t_diff_mu
        else:
            t_gate_mu = 0
            t_photon_add_mu = t_diff_mu
            t_FORT_OFF_mu = 0

        # t_photon = abs(t_diff)
        # self.t_photon_collection_time = self.t_photon_collection_time + int(t_photon_add) * ns

    else:
        # t_gate_mu = 0
        # t_photon = 0
        if t_diff_mu < 0:
            t_gate_mu = 0
            t_photon_add_mu = -1 * t_diff_mu
            t_FORT_OFF_mu = 0
        else:
            t_gate_mu = -1 * t_diff_mu
            t_photon_add_mu = t_diff_mu
            t_FORT_OFF_mu = -1 * t_diff_mu

    # t_diff_mu = 0
    # t_gate_mu = 0
    # t_photon_add_mu = 0
    # t_FORT_OFF_mu = 0

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        delay(10*us)
        two_nodes_synchronization(self)

        ##todo: add synchronization here also??
        ##      also, add how much time it took for synchronization.
        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

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

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            delay(10*us)
            if self.t_pumping > 0.0:
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10*us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
            # delay(10*us)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(5 * us)

            self.core.break_realtime()

            two_nodes_synchronization(self)

            for excitation_attempt in range(self.n_excitation_attempts):

                if self.which_node == "bob":
                    delay_mu(int(self.t_delay_in_bob_mu))

                t0 = now_mu()
                t1 = int(t0 + 150)

                at_mu(t1 - int(t_FORT_OFF_mu))
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns) + int(t_photon_add_mu))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu) + int(t_gate_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(10*us)

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings

            if self.t_recooling > 0 and (excitation_cycle + 1) % self.recool_every_n_OP == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                delay(1 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(1 * us)

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(1*us)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                # delay(0.4*ms)

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
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())


                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.two_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(1*us)
            delay(10 * us)

        # delay(1 * ms)
        # self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.01 * ms)

        # two_nodes_synchronization(self)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        self.core.break_realtime()

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def Two_node_single_photon_2_optimization_experiment(self):
    """
    **** This is for Hong-Ou_Mandel Experiment
    **** This is based on Two_node_single_photon_experiment

    **** In HOM Measurement, we want to vary the excitation timing of one node to see
         see the HOM dip. However while we do that, we might accidentally keep the
         FORT from the other node on while one node is collecting photons, making
         false count measurements. Thus, here I added some timing changes in
         excitation cycle loop.

    ****
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    self.core.break_realtime()

    t_diff_mu = 0
    t_gate_mu = 0
    t_photon_add_mu = 0
    t_FORT_OFF_mu = 0

    while self.measurement < self.n_measurements:
        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

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

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
        delay(10*us)
        if self.t_pumping > 0.0:
            if self.which_node == "alice":
                CW_optical_pumping_node1(self)
            else:
                CW_optical_pumping_node2(self)
            delay(10*us)

        ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
        # delay(10*us)

        self.ttl_exc0_switch.off() # turns on the excitation0 AOM
        delay(5 * us)

        self.core.break_realtime()

        two_nodes_synchronization(self)

        for excitation_attempt in range(self.n_excitation_attempts):

            if self.which_node == "bob":
                delay_mu(int(self.t_delay_in_bob_mu))

            t1 = now_mu()

            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation


            self.core.break_realtime()
            two_nodes_synchronization(self)

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.1 * ms)

        two_nodes_synchronization(self)

        if self.t_microwave_pulse > 0.0:
            # ### Changing the bias field
            # if self.which_node == "alice":
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            # delay(1 * ms)

            # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
            # delay(0.1 * ms)
            # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(2 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            # FORT_ramp2_smoothstep(self, direction="up")
            delay(2 * us)


        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        two_nodes_synchronization(self)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()