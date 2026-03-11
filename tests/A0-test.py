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
    delay(100 * us)
    if self.measurement > 0:
        if self.BothSPCMs_RO2/self.t_SPCM_second_shot > self.single_atom_threshold:
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

        BothSPCMs_atom_check_loaded = 0  ### for initilization
        BothSPCMs_atom_check_not_loaded = 0

        while True:
            while not atom_loaded and try_n < max_tries:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(atom_check_time)
                    self.ttl_SPCM1_counter.gate_rising(atom_check_time)

                BothSPCMs_atom_check = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)

                try_n += 1

                ### To save only one photon counts of unloaded case for each loaded atom. Otherwise, the unloaded counts
                ### would overwhelm the dataset.
                if try_n == 1:
                    BothSPCMs_atom_check_not_loaded = BothSPCMs_atom_check

                if BothSPCMs_atom_check / atom_check_time > self.single_atom_threshold_for_loading:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True
                    BothSPCMs_atom_check_loaded = BothSPCMs_atom_check


            if atom_loaded:
                self.set_dataset("time_without_atom", 0.0, broadcast=True) ### resetting time_without_atom when we load an atom
                t_after_atom = now_mu()

                ### just to check the histogram during atom loading to find a good single_atom_threshold_for_loading
                self.append_to_dataset("BothSPCMs_atom_check_in_loading", BothSPCMs_atom_check_loaded)
                self.append_to_dataset("BothSPCMs_atom_check_in_loading", BothSPCMs_atom_check_not_loaded)
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

                # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
                # delay(0.1 * ms)
                # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

                ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
                ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(0.1 * ms)
                self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
                self.laser_stabilizer.run()
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
        if self.which_node == 'alice':
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            FORT_ramp1_smoothstep(self, direction="down")
        elif self.which_node == 'bob':
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

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
    ### set the cooling DP AOM to the readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                            amplitude=self.ampl_cooling_DP_RO)

    ### set the FORT AOM to the science settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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
        self.BothSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)
        delay(0.1 * ms)
        self.dds_cooling_DP.sw.off()  ### turn off cooling
        self.ttl_repump_switch.on()  ### turn off MOT RP
        delay(10 * us)

    else:
        self.ttl_repump_switch.off()  ### turn on MOT RP
        self.dds_cooling_DP.sw.on()  ### Turn on cooling
        delay(0.1 * ms)
        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
        # delay(0.1 * ms)
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)
        # delay(1 * ms)
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
        self.BothSPCMs_RO2 = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)
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

        self.SPCM0_RO2 = self.ttl_SPCM0_counter.fetch_count()
        self.SPCM1_RO2 = self.ttl_SPCM1_counter.fetch_count()
        self.BothSPCMs_RO2 = int((self.SPCM0_RO2 + self.SPCM1_RO2) / 2)

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

    if self.which_node == "alice":
        # delay(1 * ms)
        # measure_GRIN1(self)
        # delay(1 * ms)
        # measure_PUMPING_REPUMP(self)
        # delay(1 * ms)
        # measure_Magnetometer(self)
        # delay(1*ms)
        # Sampler0_test(self)
        # delay(1*ms)
        # measure_coil_driver(self)
        # delay(1*ms)
        # measure_MOT_end(self)
        # delay(1*ms)
        measure_REPUMP(self)
    else:
        """
        test used for monitring MOT power
        MOT_end_monitor1 defined

        This is in end_measurement

        AOM1: Sampler0, 7
        
        AOM1 test: Sampler2, 1
        
        """
        ao_s1 = 7
        ao_s1_test = 1

        # avgs = 50
        #
        # self.dds_FORT.sw.off()
        # self.ttl_repump_switch.on()  # turns the RP AOM off
        # self.dds_cooling_DP.sw.on()
        #
        # delay(0.1 * ms)
        #
        # ### MOT1 & MOT2 & MOT5
        # measurement_buf = np.array([0.0] * 8)
        # measurement1 = 0.0  # MOT1
        # measurement2 = 0.0  # MOT1
        #
        # self.dds_AOM_A1.sw.on()
        #
        # delay(0.1 * ms)
        #
        # for i in range(avgs):
        #     self.sampler0.sample(measurement_buf)
        #     delay(0.1 * ms)
        #     measurement1 += measurement_buf[ao_s1]  # MOT1
        #
        #     self.sampler2.sample(measurement_buf)
        #     delay(0.1 * ms)
        #     measurement1 += measurement_buf[ao_s1_test]  # MOT1 test
        #
        #
        # measurement1 /= avgs
        # measurement2 /= avgs
        #
        # self.append_to_dataset("MOT1_end_monitor", measurement1)
        # self.append_to_dataset("MOT2_end_monitor", measurement2)

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
        if not in_health_check:  ## advance and in_health_check are different type so can't be mixed.
            self.append_to_dataset('SPCM0_RO1', self.SPCM0_RO1)
            self.append_to_dataset('SPCM1_RO1', self.SPCM1_RO1)
            self.append_to_dataset('SPCM0_RO2', self.SPCM0_RO2)
            self.append_to_dataset('SPCM1_RO2', self.SPCM1_RO2)
            self.append_to_dataset('BothSPCMs_RO1', self.BothSPCMs_RO1)
            self.append_to_dataset('BothSPCMs_RO2', self.BothSPCMs_RO2)
            self.append_to_dataset('SPCM0_test_RO', self.SPCM0_test_RO)
            delay(1 * ms)
        else:
            self.append_to_dataset('SPCM0_RO1_in_health_check', self.SPCM0_RO1)
            self.append_to_dataset('SPCM1_RO1_in_health_check', self.SPCM1_RO1)
            self.append_to_dataset('SPCM0_RO2_in_health_check', self.SPCM0_RO2)
            self.append_to_dataset('SPCM1_RO2_in_health_check', self.SPCM1_RO2)
            self.append_to_dataset('BothSPCMs_RO1_in_health_check', self.BothSPCMs_RO1)
            self.append_to_dataset('BothSPCMs_RO2_in_health_check', self.BothSPCMs_RO2)

@kernel
def FORT_ramp1_smoothstep(self, direction="down"):
    """
    For ramping FORT from loading setpoint to science and vice versa.
    Smoothly ramp FORT power using a quintic smoothstep profile. If t_FORT_ramp is too short (<1ms), it uses less
    number of steps to avoid Underflow errors. This can handle any t_FORT_ramp, from 1us to 10ms, for example.

    direction: "down" or "up"
    """
    # assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"
    #
    # p_high = self.stabilizer_FORT.amplitudes[0]
    # p_low = self.stabilizer_FORT.amplitudes[1]
    # n_steps = 100
    # step_delay = self.t_FORT_ramp / n_steps
    #
    # for i in range(n_steps):
    #     x = i / (n_steps - 1)  # normalized ramp position in [0,1]
    #     smoothstep = 6 * x ** 5 - 15 * x ** 4 + 10 * x ** 3
    #
    #     if direction == "down":
    #         p_FORT = p_high - smoothstep * (p_high - p_low)
    #     else:  # direction == "up"
    #         p_FORT = p_low + smoothstep * (p_high - p_low)
    #
    #     delay(step_delay / 2)
    #     self.dds_FORT.set(frequency=self.f_FORT, amplitude=p_FORT)
    #     delay(step_delay / 2)

    assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"

    p_high = self.stabilizer_FORT.amplitudes[0]
    p_low = self.stabilizer_FORT.amplitudes[1]
    n_steps_max = 100
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

@kernel
def FORT_ramp2_smoothstep(self, direction="down"):
    """
    For ramping FORT from science setpoint to holding (microwave) and vice versa.
    Smoothly ramp FORT power using a quintic smoothstep profile. If t_FORT_ramp is too short (<1ms), it uses less
    number of steps to avoid Underflow errors. This can handle any t_FORT_ramp, from 1us to 10ms, for example.

    direction: "down" or "up"
    """

    assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"

    p_high = self.stabilizer_FORT.amplitudes[1]
    p_low = self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1]
    n_steps_max = 100
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

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        load_until_atom_smooth_FORT_recycle(self)

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

