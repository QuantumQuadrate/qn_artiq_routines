from artiq.experiment import *
import logging
import numpy as np
from abc import ABC

import os, sys
cwd = os.getcwd() + "\\"

sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.conversions import dB_to_V_kernel

# todo: think about how to implement the experiment functions in a wrapper.
#  need to pass the experiments by reference.
"""
Table of contents:
1. Subroutine functions
2. Experiment functions
"""

###############################################################################
# 1. SUBROUTINE FUNCTIONS
# These are functions that get used within various experiments. they are not
# intended to be run in a standalone fashion, e.g. from GeneralVariableScan.
# Consequently, note that the name should not end in "experiment"
###############################################################################

@kernel
def load_MOT_and_FORT(self):
    """
    The FORT loading sequence, from loading a MOT to hopefully trapping one atom

    Loads the MOT then the turns on the FORT and waits a fixed amount of time
    to load an atom. Optionally, the FORT can be turned on right away if
    FORT_on_at_start is True.

    :param self: the experiment instance
    :return:
    """

    self.dds_FORT.sw.on()

    if not self.FORT_on_at_MOT_start:
        self.dds_FORT.sw.off()
    else:
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    # self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

    self.dds_cooling_DP.sw.on()

    delay(1*ms)

    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(1*ms) # if this delay is not here, the following line setting the dac doesn't execute

    # Turn on the MOT coils and cooling light - this is being ignored
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    self.ttl_UV.pulse(self.t_UV_pulse)

    # wait for the MOT to load
    delay(self.t_MOT_loading - self.t_MOT_phase2)

    if self.t_MOT_phase2 > 0:

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT_phase2, self.AZ_top_volts_MOT_phase2, self.AX_volts_MOT_phase2,
             self.AY_volts_MOT_phase2],
            channels=self.coil_channels)

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT_phase2,
                                amplitude=self.ampl_cooling_DP_MOT) # todo: make a variable for phase 2
        delay(self.t_MOT_phase2)

    # turn on the dipole trap and wait to load atoms
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    if not self.FORT_on_at_MOT_start:
        delay_mu(self.t_FORT_loading_mu)

    self.stabilizer_FORT.run(setpoint_index=1) # the science setpoint

    self.dds_cooling_DP.sw.off()
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate
    self.dds_cooling_DP.sw.on()

    if self.do_PGC_in_MOT and self.t_PGC_in_MOT > 0:

        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        self.zotino0.set_dac([self.AZ_bottom_volts_PGC, self.AZ_top_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                             channels=self.coil_channels)

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_PGC)
        delay(self.t_PGC_in_MOT)

    self.dds_cooling_DP.sw.off()

@kernel
def record_chopped_blow_away(self):
    """

    :param self:
    :return:
    """

    # todo: change OP -> BA

    n_chop_cycles = int(self.t_blowaway/self.t_BA_chop_period + 0.5)
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

    # set coils for blowaway
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_blowaway, self.AZ_top_volts_blowaway,
         self.AX_volts_blowaway, self.AY_volts_blowaway],
        channels=self.coil_channels)
    delay(0.1 * ms)

    with sequential:

        self.dds_cooling_DP.set(
            frequency=self.f_cooling_DP_blowaway,
            amplitude=self.ampl_cooling_DP_MOT)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()

        if self.blowaway_light_off:
            self.dds_cooling_DP.sw.off()
            self.dds_AOM_A6.sw.off()
        else:
            # just turn the AOM up all the way. as long as we're 'saturating' the blowaway, it's okay if this doesn't
            # always give the same optical power
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,
                                amplitude=dB_to_V_kernel(-7.0))
            self.dds_AOM_A6.sw.on()
            self.dds_cooling_DP.sw.on()

    self.core_dma.playback_handle(ba_dma_handle)

    # reset AOM RF powers
    self.dds_cooling_DP.sw.off()
    self.dds_cooling_DP.set(
        frequency=self.f_cooling_DP_RO,
        amplitude=self.ampl_cooling_DP_MOT)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,
                        amplitude=self.stabilizer_AOM_A6.amplitude)
    delay(0.1*ms)
    self.dds_AOM_A1.sw.on()
    self.dds_AOM_A2.sw.on()
    self.dds_AOM_A3.sw.on()
    self.dds_AOM_A4.sw.on()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()
    self.ttl_repump_switch.off()  # turns on the RP AOM

@kernel
def record_chopped_optical_pumping(self):
    """

    :param self:
    :return:
    """

    n_chop_cycles = int(self.t_pumping/self.t_OP_chop_period + 0.5)
    assert n_chop_cycles >= 1, "t_pumping should be > t_OP_chop_period"

    # todo: use duty cycle ExperimentVariables
    OP_pulse = self.t_OP_chop_period * 0.3
    FORT_pulse = self.t_OP_chop_period - OP_pulse

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
        OP_on_mu = self.core.seconds_to_mu(0.5 * us)

        if not (self.pumping_light_off or self.D1_off_in_OP_phase):
            for i in range(n_chop_cycles):
                at_mu(start+i*period_mu+FORT_on_mu)
                self.dds_FORT.sw.off()
                delay_mu(OP_pulse_length_mu)
                self.dds_FORT.sw.on()
                at_mu(start+i*period_mu+OP_on_mu)
                self.dds_D1_pumping_SP.sw.on()
                delay_mu(OP_pulse_length_mu)
                self.dds_D1_pumping_SP.sw.off()
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
                    self.dds_D1_pumping_SP.sw.on()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_D1_pumping_SP.sw.off()
            else:
                for i in range(n_depump_chop_cycles):
                    at_mu(start + i * period_mu + FORT_on_mu)
                    self.dds_FORT.sw.off()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_FORT.sw.on()

@kernel
def chopped_optical_pumping(self):
    """

    :param self:
    :return:
    """

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    if self.t_depumping + self.t_pumping > 3*ms:
        delay(2000 * us)  # we need extra slack
    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.dds_cooling_DP.sw.off()  # no cooling light

    # ramp up the fiber AOMs to maximize the amount of pumping repump we get
    if not self.pumping_light_off:
        self.dds_pumping_repump.sw.on()
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V_kernel(-7.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V_kernel(-7.0))
    delay(1*us)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    # set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.4 * ms)  # coil relaxation time

    with sequential:

        delay(1*us)

        self.core_dma.playback_handle(op_dma_handle)
        delay(self.t_depumping)

        self.dds_D1_pumping_SP.sw.off()
        self.dds_pumping_repump.sw.off()

        # reset MOT power
        self.dds_cooling_DP.sw.off()
        self.dds_cooling_DP.set(
            frequency=self.f_cooling_DP_RO,
            amplitude=self.ampl_cooling_DP_MOT)

        # reset the fiber AOM amplitudes
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

@kernel
def end_measurement(self):
    """
    End the measurement by setting datasets and deciding whether to increment the measuement index
    :param self:
    :return measurement: TInt32, the measurement index
    """

    # update the datasets
    if not self.no_first_shot:
        self.append_to_dataset('photocounts_current_iteration', self.counts)
        self.counts_list[self.measurement] = self.counts

    # update the datasets
    self.set_dataset(self.measurements_progress, 100*self.measurement/self.n_measurements, broadcast=True)
    self.append_to_dataset('photocounts2_current_iteration', self.counts2)
    self.counts2_list[self.measurement] = self.counts2

    advance = 1
    if self.__class__.__name__ != 'ExperimentCycler':
        if self.require_atom_loading_to_advance:
            if not self.counts/self.t_SPCM_first_shot > self.single_atom_counts_per_s:
                advance *= 0
        if self.require_D1_lock_to_advance:
            self.ttl_D1_lock_monitor.sample_input()
            delay(0.1 * ms)
            advance *= int(1 - self.ttl_D1_lock_monitor.sample_get())
            if not advance:
                logging.warning("D1 laser not locked")

    if advance:
        self.measurement += 1
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', self.counts)
        self.append_to_dataset('photocounts2', self.counts2)

###############################################################################
# 2. EXPERIMENT FUNCTIONS
# These are the experiments we run, and the name of each should end with
# experiment in order to have it show up in GeneralVariableScans
###############################################################################

@kernel
def test_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """
    x = self.t_blowaway
    self.append_to_dataset('test_dataset', x)
    self.print_async(x)

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

    self.counts = 0
    self.counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the science setting. this is only valid if we have run
        # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT

        self.dds_FORT.set(frequency=self.f_FORT,
                                amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            self.counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)

        # take the second shot
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        self.counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def optical_pumping_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.counts = 0
    self.counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)
    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            self.counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()
        delay(1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM off

        # # set the FORT to the holding setting, i.e. for doing nothing
        # self.dds_FORT.set(frequency=self.f_FORT,
        #                   amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_holding)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            self.ttl7.pulse(10 * us)  # in case we want to look at signals on an oscilloscope
            chopped_optical_pumping(self)

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # take the second shot
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM on
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        self.counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def microwave_Rabi_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    self.counts = 0
    self.counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)
    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            self.counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()
        delay(1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM off

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)

        ############################
        # microwave phase
        ############################

        if self.t_microwave_pulse > 0.0:
            self.ttl_repump_switch.on()  # turns off the RP AOM

            # todo: set coils for microwaves. good for diagnostics-- we can use this phase to zero the B-field
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP,
                 self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(0.1*ms)
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, self.AZ_top_volts_microwave,
            #      self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)
            # delay(0.1*ms)

            self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

            self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=self.ampl_microwaves)
            self.dds_microwaves.sw.on()
            self.ttl_microwave_switch.off()

            delay(self.t_microwave_pulse)

            self.dds_microwaves.sw.off()
            self.ttl_microwave_switch.on()

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # take the second shot
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM on
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        self.counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)

    self.dds_FORT.sw.off()

@kernel
def single_photon_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    self.counts = 0
    self.counts2 = 0
    excitation_counts = 0
    excitation_counts_array = [0] # overwritten below but initialized here so it is always initialized

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    record_chopped_optical_pumping(self)
    delay(100*ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        excitation_counts_array = [0] * self.n_excitation_cycles

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs
        delay(10*ms)
        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        if not self.no_first_shot:
            self.ttl_SPCM_gate.off()
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            self.counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()
        delay(1 * ms)
        self.ttl_repump_switch.on()  # turns the RP AOM off

        # # set the FORT to the holding setting, i.e. for doing nothing
        # self.dds_FORT.set(frequency=self.f_FORT,
        #                   amplitude=self.stabilizer_FORT.amplitude)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        self.ttl7.pulse(10 * us)  # in case we want to look at signals on an oscilloscope

        excitation_counts = 0
        for excitaton_cycle in range(self.n_excitation_cycles):

            delay(0.5*ms)
            self.ttl_SPCM_gate.on() # blocks the SPCM output

            ############################
            # optical pumping phase - pumps atoms into F=1,m_F=0
            ############################

            chopped_optical_pumping(self)

            ############################
            # excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            ############################

            now = now_mu()

            at_mu(now+1)
            with parallel:
                self.dds_FORT.sw.off()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
            at_mu(now+10)
            self.ttl_repump_switch.off()  # repump AOM is on for excitation
            mu_offset = 800 # accounts for various latencies
            at_mu(now + mu_offset - 200) # make sure stuff is off, no more Raman photons from FORT
            # self.ttl_repump_switch.off()  # repump AOM is on for excitation
            at_mu(now + mu_offset+100+self.gate_start_offset_mu) # allow for repump rise time and FORT after-pulsing
            t_collect = now_mu()
            t_gate_end = self.ttl0.gate_rising(self.n_excitation_attempts * (self.t_excitation_pulse + 100 * ns))
            t_excite = now_mu()
            pulses_over_mu = 0
            for attempt in range(self.n_excitation_attempts):
                at_mu(now + mu_offset + 201 + int(attempt * (self.t_excitation_pulse / ns + 100))
                      + self.gate_start_offset_mu)
                self.dds_excitation.sw.pulse(self.t_excitation_pulse)
                at_mu(now + mu_offset + 741 + int(attempt * (self.t_excitation_pulse / ns + 100) -
                                                  0.1*self.t_excitation_pulse / ns) +self.gate_start_offset_mu)
                # fast switch to gate SPCM output
                self.ttl_SPCM_gate.off()
                delay(0.2*self.t_excitation_pulse+100*ns + self.gate_switch_offset)
                self.ttl_SPCM_gate.on()

                pulses_over_mu = now_mu()
            self.ttl_repump_switch.on()
            at_mu(pulses_over_mu - 200) # fudge factor
            self.ttl_SPCM_gate.on() # TTL high turns switch off, i.e. signal blocked
            self.dds_FORT.sw.on()
            excitation_counts = self.ttl0.count(
                t_gate_end)  # this is the number of clicks we got over n_excitation attempts
            excitation_counts_array[excitaton_cycle] = excitation_counts
            delay(0.1*ms) # ttl count consumes all the RTIO slack.

        delay(1 * ms)
        # self.ttl_SPCM_gate.on()  # enables the SPCM

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitudes[1])

        # take the second shot
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)
        self.ttl_SPCM_gate.off()
        self.ttl_repump_switch.off()  # turns the RP AOM on
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        self.counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)
        for val in excitation_counts_array:
            self.append_to_dataset('excitation_counts', val)

        delay(10*ms)

    self.dds_FORT.sw.off()