from artiq.experiment import *
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
        # FORT is on to thermalize but frequency shifted to not couple to the fiber
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
    else:
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

    self.dds_cooling_DP.sw.on()

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

    self.dds_cooling_DP.sw.off()
    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate
    self.dds_cooling_DP.sw.on()

    if self.do_PGC_in_MOT and self.t_PGC_in_MOT > 0:

        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_PGC)

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

    n_chop_cycles = int(self.t_blowaway /self.t_BA_chop_period)
    assert n_chop_cycles > 1, "t_blowaway should be > t_BA_chop_period"

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
    self.dds_FORT.set(
        frequency=self.f_FORT,
        amplitude=self.stabilizer_FORT.amplitude*self.p_FORT_RO)
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
def blow_away(self):
    self.ttl_repump_switch.on()  # turns off the RP AOM

    # set coils for blowaway
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_blowaway, self.AZ_top_volts_blowaway,
         self.AX_volts_blowaway, self.AY_volts_blowaway],
        channels=self.coil_channels)
    delay(0.1 * ms)

    with sequential:

        # lower FORT power
        self.dds_FORT.set(
            frequency=self.f_FORT,
            amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_blowaway)

        self.dds_cooling_DP.set(
            frequency=self.f_cooling_DP_blowaway,  # resonant_2_to_3,
            amplitude=self.ampl_cooling_DP_MOT)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()

        if self.blowaway_light_off:
            self.dds_AOM_A6.sw.off()
        else:
            # just turn the AOM up all the way. as long as we're 'saturating' the blowaway, it's okay if this doesn't
            # always give the same optical power
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,
                                amplitude=dB_to_V_kernel(-7.0))
            self.dds_cooling_DP.sw.on()

    delay(self.t_blowaway)

    # reset AOM RF powers
    self.dds_cooling_DP.sw.off()
    self.dds_cooling_DP.set(
        frequency=self.f_cooling_DP_RO,
        amplitude=self.ampl_cooling_DP_MOT)
    self.dds_FORT.set(
        frequency=self.f_FORT,
        amplitude=self.stabilizer_FORT.amplitude*self.p_FORT_RO)
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

    n_chop_cycles = int(self.t_pumping /self.t_OP_chop_period)
    assert n_chop_cycles > 1, "t_pumping should be > t_OP_chop_period"

    OP_pulse = self.t_OP_chop_period * 0.35
    FORT_pulse = self.t_OP_chop_period - OP_pulse

    self.core.reset()

    with self.core_dma.record("chopped_optical_pumping"):

        start = now_mu()
        period_mu = self.core.seconds_to_mu(self.t_OP_chop_period)
        OP_pulse_length_mu = self.core.seconds_to_mu(OP_pulse)
        OP_on_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_pulse_length_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_on_mu = self.core.seconds_to_mu(OP_pulse)

        self.dds_FORT.sw.off()
        delay_mu(OP_pulse_length_mu)

        if not self.pumping_light_off:
            for i in range(n_chop_cycles):
                at_mu(start+i*period_mu+FORT_on_mu)
                self.dds_FORT.sw.on()
                delay_mu(FORT_pulse_length_mu)
                self.dds_FORT.sw.off()
                at_mu(start+i*period_mu+OP_on_mu)
                self.dds_D1_pumping_SP.sw.on()
                delay_mu(OP_pulse_length_mu)
                self.dds_D1_pumping_SP.sw.off()
            self.dds_FORT.sw.on()
        else:
            for i in range(n_chop_cycles):
                at_mu(start + i * period_mu + FORT_on_mu)
                self.dds_FORT.sw.on()
                delay_mu(FORT_pulse_length_mu)
                self.dds_FORT.sw.off()
                at_mu(start + i * period_mu + OP_on_mu)
                # self.dds_D1_pumping_SP.sw.on()
                delay_mu(OP_pulse_length_mu)
                # self.dds_D1_pumping_SP.sw.off()
            self.dds_FORT.sw.on()


@kernel
def chopped_optical_pumping(self):
    """

    :param self:
    :return:
    """

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.dds_cooling_DP.sw.off()  # no cooling light

    # ramp up the fiber AOMs to maximize the amount of pumping repump we get
    if not self.pumping_light_off:
        self.dds_pumping_repump.sw.on()
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V_kernel(-7.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V_kernel(-7.0))

    # set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.1 * ms)  # coil relaxation time

    with sequential:

        # adjust FORT power
        self.dds_FORT.set(
            frequency=self.f_FORT,
            amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_OP)

        delay(1*us)

        self.core_dma.playback_handle(op_dma_handle)

        self.dds_D1_pumping_SP.sw.off()
        self.dds_pumping_repump.sw.off()

        self.dds_FORT.set(
            frequency=self.f_FORT,
            amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # reset MOT power
        self.dds_cooling_DP.sw.off()
        self.dds_cooling_DP.set(
            frequency=self.f_cooling_DP_RO,
            amplitude=self.ampl_cooling_DP_MOT)

        # reset the fiber AOM amplitudes
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

def optical_pumping(self):
    """

    :param self:
    :return:
    """

    self.ttl_repump_switch.on()  # turns off the RP AOM
    self.dds_cooling_DP.sw.off()  # no cooling light

    # set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(0.1 * ms)  # coil relaxation time

    with sequential:

        # lower FORT power
        self.dds_FORT.set(
            frequency=self.f_FORT,
            amplitude=self.stabilizer_FORT.amplitude * p_FORT_OP)

        if self.pumping_light_off:
            self.dds_D1_pumping_SP.sw.off()
            self.dds_pumping_repump.sw.off()
        else:
            self.dds_D1_pumping_SP.sw.on()
            self.dds_pumping_repump.sw.on()

        delay(self.t_pumping)

        self.dds_D1_pumping_SP.sw.off()
        self.dds_pumping_repump.sw.off()

        # reset MOT power
        self.dds_cooling_DP.sw.off()
        self.dds_cooling_DP.set(
            frequency=self.f_cooling_DP_RO,
            amplitude=self.ampl_cooling_DP_MOT)


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

    counts = 0
    counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    for measurement in range(self.n_measurements):

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                                amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

        # set the FORT to the holding setting, i.e. for doing nothing
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_holding)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # take the second shot
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)
            self.counts_list[measurement] = counts

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)
        self.counts2_list[measurement] = counts2

    # effectively turn the FORT AOM off
    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

@kernel
def optical_pumping_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    counts = 0
    counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    if self.t_pumping > 0.0:
        record_chopped_optical_pumping(self)
        delay(100*ms)
    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100*ms)

    for measurement in range(self.n_measurements):

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()
        delay(1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM off

        # set the FORT to the holding setting, i.e. for doing nothing
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_holding)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)

        # # raise FORT power back to holding level
        # self.dds_FORT.set(
        #     frequency=self.f_FORT,
        #     amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_holding)

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:
            # blow_away(self)
            chopped_blow_away(self)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # take the second shot
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)
        self.ttl_repump_switch.off()  # turns the RP AOM on
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)
            self.counts_list[measurement] = counts

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)
        self.counts2_list[measurement] = counts2

    # effectively turn the FORT AOM off
    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

@kernel
def microwave_Rabi_experiment(self):
    """
    atom loading experiment with a microwave pulse of duration t_microwave_pulse
    between the readouts

    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.reset()

    counts = 0
    counts2 = 0

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

    self.dds_D1_pumping_SP.sw.off()

    for measurement in range(self.n_measurements):

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT(self)

        delay(0.1 * ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT * self.p_cooling_DP_RO)

        if not self.no_first_shot:
            self.ttl_repump_switch.off()  # turns the RP AOM on
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()
        delay(1*ms)
        self.ttl_repump_switch.off() # turns the RP AOM off

        # set the FORT to the holding setting, i.e. for doing nothing
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_holding)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################
        if self.t_pumping > 0.0:
            chopped_optical_pumping(self)

        # raise FORT power back to holding level
        self.dds_FORT.set(
            frequency=self.f_FORT,
            amplitude=self.stabilizer_FORT.amplitude*self.p_FORT_holding)

        ############################
        # microwave phase
        ############################

        if self.t_microwave_pulse > 0.0:
            self.ttl_repump_switch.on()  # turns off the RP AOM

            # set coils for microwaves. good for diagnostics-- we can use this phase to zero the B-field
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_microwave, self.AZ_top_volts_microwave,
                 self.AX_volts_microwave, self.AY_volts_microwave],
                channels=self.coil_channels)
            delay(0.1*ms)

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
            blow_away(self)
        delay(self.t_delay_between_shots - self.t_blowaway)

        # set the FORT AOM to the readout settings
        self.dds_FORT.set(frequency=self.f_FORT,
                          amplitude=self.stabilizer_FORT.amplitude * self.p_FORT_RO)

        # take the second shot
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)
        self.ttl_repump_switch.off() # turns the RP AOM on
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)
            self.counts_list[measurement] = counts

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)
        self.counts2_list[measurement] = counts2

    # effectively turn the FORT AOM off
    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

# @kernel
# def microwave_experiment(self):
#     """
#     atom loading experiment with a microwave pulse of duration t_microwave_pulse
#     between the readouts
#
#     self is the experiment instance to which ExperimentVariables are bound
#     """
#
#     self.core.break_realtime()
#
#     # initialize these or ARTIQ complains
#     counts = 0
#     counts2 = 0
#
#     for measurement in range(self.n_measurements):
#
#         if self.enable_laser_feedback:
#             if measurement % 10 == 0:
#                 self.laser_stabilizer.run() # this tunes the MOT and FORT AOMs
#
#         ############################
#         # load the MOT and FORT
#         ############################
#
#         load_MOT_and_FORT(self)
#
#         ############################
#         # take the first shot
#         ############################
#         if not self.no_first_shot:
#             self.dds_cooling_DP.sw.on()
#             t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
#             counts = self.ttl0.count(t_gate_end)
#             delay(1 * ms)
#             self.dds_cooling_DP.sw.off()
#
#         delay(3 * ms)
#
#         ############################
#         # optical pumping phase - pumps atoms into F=1,m_F=0
#         ############################
#
#         if self.t_pumping > 0.0:
#             self.ttl_repump_switch.on()  # turns off the RP AOM
#             self.dds_cooling_DP.sw.off()  # no MOT light
#
#             # set coils for pumping
#             self.zotino0.set_dac(
#                 [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
#                 channels=self.coil_channels)
#             delay(0.1 * ms) # coil relaxation time
#
#             with sequential:
#
#                 # lower FORT power
#                 self.dds_FORT.set(
#                     frequency=self.f_FORT,
#                     amplitude=self.ampl_FORT_OP)
#
#                 # this could be condensed but is left as is for clarity
#                 if self.pumping_light_off:
#                     self.dds_D1_pumping_SP.sw.off()
#                     self.dds_pumping_repump.sw.off()
#                 elif self.control_experiment and measurement % 2 == 0:
#                     self.dds_D1_pumping_SP.sw.off()
#                     self.dds_pumping_repump.sw.off()
#                 else:
#                     self.dds_D1_pumping_SP.sw.on()
#                     self.dds_pumping_repump.sw.on()
#
#                 delay(self.t_pumping)
#
#                 self.dds_D1_pumping_SP.sw.off()
#                 self.dds_pumping_repump.sw.off()
#
#                 # reset MOT power
#                 self.dds_cooling_DP.sw.off()
#                 self.dds_cooling_DP.set(
#                     frequency=self.f_cooling_DP_RO,
#                     amplitude=self.ampl_cooling_DP_MOT)
#
#         # raise FORT power back to default
#         self.dds_FORT.set(
#             frequency=self.f_FORT,
#             amplitude=self.stabilizer_FORT.amplitude)
#
#         ############################
#         # microwave phase
#         ############################
#
#         self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=self.ampl_microwaves)
#         self.dds_microwaves.sw.on()
#         self.ttl_microwave_switch.off()
#
#         ############################
#         # blow-away phase - push out atoms in F=2 only
#         ############################
#
#         if self.t_blowaway > 0.0:
#
#             self.ttl_repump_switch.on()  # turns off the RP AOM
#
#             # set coils for readout # todo: update with blowaway specific values later
#             self.zotino0.set_dac(
#                 [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
#                 channels=self.coil_channels)
#             delay(0.1 * ms)
#
#             with sequential:
#
#                 # lower FORT power
#                 self.dds_FORT.set(
#                     frequency=self.f_FORT,
#                     amplitude=self.ampl_FORT_blowaway)
#
#                 # set cooling light to be resonant with a free-space atom
#                 self.dds_cooling_DP.set(
#                     frequency=self.f_cooling_DP_resonant_2_to_3,
#                     amplitude=self.ampl_cooling_DP_MOT)
#
#                 self.dds_AOM_A1.sw.off()
#                 self.dds_AOM_A2.sw.off()
#                 self.dds_AOM_A3.sw.off()
#                 self.dds_AOM_A4.sw.off()
#                 self.dds_AOM_A5.sw.off()
#
#                 if self.blowaway_light_off:
#                     self.dds_AOM_A6.sw.off()
#                 else:
#                     self.dds_cooling_DP.sw.on()
#
#             delay(self.t_blowaway)
#
#             # reset MOT power
#             self.dds_cooling_DP.sw.off()
#             self.dds_cooling_DP.set(
#                 frequency=self.f_cooling_DP_RO,
#                 amplitude=self.ampl_cooling_DP_MOT)
#
#             # reset FORT power
#             self.dds_FORT.set(
#                 frequency=self.f_FORT,
#                 amplitude=self.stabilizer_FORT.amplitude)
#
#             self.dds_AOM_A1.sw.on()
#             self.dds_AOM_A2.sw.on()
#             self.dds_AOM_A3.sw.on()
#             self.dds_AOM_A4.sw.on()
#             self.dds_AOM_A5.sw.on()
#             self.dds_AOM_A6.sw.on()
#             self.ttl_repump_switch.off()  # turns on the RP AOM
#
#         delay(0.1 * ms)
#
#         ############################
#         # take the second shot
#         ############################
#
#         # set coils for readout
#         self.zotino0.set_dac(
#             [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
#             channels=self.coil_channels)
#         delay(0.1 * ms)
#
#         self.dds_cooling_DP.sw.on()
#         t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
#         counts2 = self.ttl0.count(t_gate_end)
#         delay(1 * ms)
#         self.dds_cooling_DP.sw.off()
#
#         delay(2 * ms)
#
#         self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope
#
#         # update the datasets
#         if not self.no_first_shot:
#             self.append_to_dataset('photocounts', counts)
#             self.append_to_dataset('photocounts_current_iteration', counts)
#
#         # update the datasets
#         self.append_to_dataset('photocounts2', counts2)
#         self.append_to_dataset('photocounts2_current_iteration', counts2)
#
#         self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope
#
#     # effectively turn the FORT AOM off
#     self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
#
# @kernel
# def single_photons_experiment(self): # for debugging
#     """
#     self is the experiment instance to which ExperimentVariables are bound
#     """
#
#     self.core.break_realtime()
#
#     # initialize these or ARTIQ complains
#     counts = 0
#     counts2 = 0
#
#     for measurement in range(self.n_measurements):
#
#         # make sure any beams that can mess up the fW measurement (e.g. OP, excitation) are turned off
#         if self.enable_laser_feedback:
#             if measurement % 10 == 0:
#                 self.laser_stabilizer.run()
#                 delay(1 * ms)
#             self.dds_FORT.sw.on()
#             self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
#
#         self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope
#
#         ############################
#         # load the MOT
#         ############################
#         self.zotino0.set_dac(
#             [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
#             channels=self.coil_channels)
#         self.ttl_repump_switch.off()  # turns on the RP AOM
#         self.dds_cooling_DP.sw.on()
#         self.dds_AOM_A1.sw.on()
#         self.dds_AOM_A2.sw.on()
#         self.dds_AOM_A3.sw.on()
#         self.dds_AOM_A4.sw.on()
#         self.dds_AOM_A5.sw.on()
#         self.dds_AOM_A6.sw.on()
#
#         # wait for the MOT to load
#         delay_mu(self.t_MOT_loading_mu)
#
#         # load atom from a PGC phase
#         if self.do_PGC_in_MOT:
#             self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)
#             self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
#             delay(self.t_PGC_in_MOT)
#
#         # turn on the dipole trap and wait to load atoms
#         self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
#         delay_mu(self.t_FORT_loading_mu)
#
#         # turn off the coils
#         if not self.do_PGC_in_MOT:
#             self.zotino0.set_dac(
#                 [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
#                 channels=self.coil_channels)
#
#         delay(3 * ms)  # should wait for the MOT to dissipate
#
#         # set the cooling DP AOM to the readout settings
#         self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)
#
#         ############################
#         # take the first shot
#         ############################
#         if not self.no_first_shot:
#             self.dds_cooling_DP.sw.on()
#             t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
#             counts = self.ttl0.count(t_gate_end)
#             delay(1 * ms)
#             self.dds_cooling_DP.sw.off()
#
#         delay(3 * ms)
#
#         ############################
#         # optical pumping phase - pumps atoms into F=1,m_F=0
#         ############################
#
#         if self.t_pumping > 0.0:
#             self.ttl_repump_switch.on()  # turns off the RP AOM
#             self.dds_cooling_DP.sw.off()  # no MOT light
#
#             # set coils for pumping
#             self.zotino0.set_dac(
#                 [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
#                 channels=self.coil_channels)
#             delay(0.1 * ms) # coil relaxation time
#
#             with sequential:
#
#                 # lower FORT power
#                 self.dds_FORT.set(
#                     frequency=self.f_FORT,
#                     amplitude=self.ampl_FORT_OP)
#
#                 # this could be condensed but is left as is for clarity
#                 if self.pumping_light_off:
#                     self.dds_D1_pumping_SP.sw.off()
#                     self.dds_pumping_repump.sw.off()
#                 elif self.control_experiment and measurement % 2 == 0:
#                     self.dds_D1_pumping_SP.sw.off()
#                     self.dds_pumping_repump.sw.off()
#                 else:
#                     self.dds_D1_pumping_SP.sw.on()
#                     self.dds_pumping_repump.sw.on()
#
#                 delay(self.t_pumping)
#
#                 self.dds_D1_pumping_SP.sw.off()
#                 self.dds_pumping_repump.sw.off()
#
#                 # reset MOT power
#                 self.dds_cooling_DP.sw.off()
#                 self.dds_cooling_DP.set(
#                     frequency=self.f_cooling_DP_RO,
#                     amplitude=self.ampl_cooling_DP_MOT)
#
#         delay(0.1 * ms)
#
#         ###########################
#         # excitation and photon collection
#         ###########################
#
#         # try to mitigate any leakage
#         self.dds_AOM_A1.sw.off()
#         self.dds_AOM_A2.sw.off()
#         self.dds_AOM_A3.sw.off()
#         self.dds_AOM_A4.sw.off()
#         self.dds_AOM_A5.sw.off()
#         self.dds_AOM_A6.sw.off()
#
#         # turn off the FORT
#         self.dds_FORT.sw.off() # fully off, rather than frequency shift just in case
#
#         # excite the atom
#         self.dds_excitation.sw.pulse(30*ns)
#
#         delay(self.t_photon_collection_delay)
#         t_gate_end = self.ttl0.gate_rising(self.t_photon_collection_time)
#         counts2 = self.ttl0.count(t_gate_end)
#
#         delay(1 * ms) # ttl count consumes all the RTIO slack
#
#         self.dds_FORT.sw.on()
#
#         # turn AOMs back on
#         self.dds_AOM_A1.sw.on()
#         self.dds_AOM_A2.sw.on()
#         self.dds_AOM_A3.sw.on()
#         self.dds_AOM_A4.sw.on()
#         self.dds_AOM_A5.sw.on()
#         self.dds_AOM_A6.sw.on()
#
#         # effectively turn the FORT AOM off
#         self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
#         # set the cooling DP AOM to the MOT settings
#         self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
#
#         delay(2 * ms)
#
#         # update the datasets
#         if not self.no_first_shot:
#             self.append_to_dataset('photocounts', counts)
#             self.append_to_dataset('photocounts_current_iteration', counts)
#
#         # update the datasets
#         self.append_to_dataset('photocounts2', counts2)
#         self.append_to_dataset('photocounts2_current_iteration', counts2)