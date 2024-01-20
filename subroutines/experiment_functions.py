from artiq.experiment import *
import numpy as np
from abc import ABC

# todo: think about how to implement the experiment functions in a wrapper.
#  need to pass the experiments by reference.
"""
arguments for having a wrapper class for each experiment
 - dataset initialization will inevitably vary by experiment so we don't want to 
 hardcode that in the variable scan code
 - experiments often are set to end a certain way (e.g., leave the MOT on,
 but not the FORT), and the function that does this should be experiment specific
 - all of the usual OOP arguments
"""

###############################################################################
# SUBROUTINE FUNCTIONS
# These are functions
###############################################################################


@kernel
def load_MOT_and_FORT(self):
    """
    The FORT loading sequence, from loading a MOT to hopefully trapping one atom

    :param self: the experiment instance
    :return:
    """

    self.dds_FORT.sw.on()
    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)

    self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)
    # delay(2 * ms)
    self.dds_cooling_DP.sw.on()

    # wait for the MOT to load
    delay(self.t_MOT_loading - self.t_MOT_phase2)

    if self.t_MOT_phase2 > 0:
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT_phase2,
                                amplitude=self.ampl_cooling_DP_MOT) # todo: make a variable for phase 2
        delay(self.t_MOT_phase2)

    # todo: try loading from a PGC phase

    # turn on the dipole trap and wait to load atoms
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_loading)
    delay_mu(self.t_FORT_loading_mu)

    # turn off the coils
    self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                         channels=self.coil_channels)

    if self.do_PGC_in_MOT and self.t_PGC_in_MOT > 0:
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
        delay(self.t_PGC_in_MOT)

    delay(self.t_MOT_dissipation)  # should wait several ms for the MOT to dissipate

@kernel
def test_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """
    x = self.t_blowaway
    self.append_to_dataset('test_dataset', x)
    self.print_async(x)

@kernel
def atom_loading_experiment(self):
    """
    :param self: an experiment instance.
    :return:
    """

    self.core.reset()

    counts = 0
    counts2 = 0

    for measurement in range(self.n_measurements):

        if self.enable_laser_feedback:
            if measurement % 10 == 0:
                self.laser_stabilizer.run()
            self.fast_laser_stabilizer.run() # this controls the FORT AOM
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)


        load_MOT_and_FORT(self)

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

        if not self.no_first_shot:
            # take the first shot
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

        delay(self.t_delay_between_shots)

        # take the second shot
        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        counts2 = self.ttl0.count(t_gate_end)
        delay(1 * ms)

        # todo: check the FORT extinction ratio here
        # effectively turn the FORT AOM off
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)
        # set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(2 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)

@kernel
def optical_pumping_experiment(self):
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.break_realtime()

    # initialize these or ARTIQ complains
    counts = 0
    counts2 = 0

    for measurement in range(self.n_measurements):

        self.core.reset()
        delay(1 * ms)

        # make sure any beams that can mess up the fW measurement (e.g. OP, excitation) are turned off
        if self.enable_laser_feedback:
            if measurement % 10 == 0:
                self.laser_stabilizer.run()
            self.fast_laser_stabilizer.run()  # this controls the FORT AOM
            self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)

        ############################
        # load the MOT and FORT
        ############################

        load_MOT_and_FORT(self)

        ############################
        # take the first shot
        ############################
        if not self.no_first_shot:
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

        delay(3 * ms)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        if self.t_pumping > 0.0:
            self.ttl_repump_switch.on()  # turns off the RP AOM
            self.dds_cooling_DP.sw.off()  # no MOT light

            # set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(0.1 * ms) # coil relaxation time

            with sequential:

                # lower FORT power
                self.dds_FORT.set(
                    frequency=self.f_FORT,
                    amplitude=self.ampl_FORT_OP)

                # this could be condensed but is left as is for clarity
                if self.pumping_light_off:
                    self.dds_D1_pumping_SP.sw.off()
                    self.dds_pumping_repump.sw.off()
                elif self.control_experiment and measurement % 2 == 0:
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

        ############################
        # blow-away phase - push out atoms in F=2 only
        ############################

        if self.t_blowaway > 0.0:

            self.ttl_repump_switch.on()  # turns off the RP AOM

            # set coils for readout # todo: update with blowaway specific values later
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                channels=self.coil_channels)
            delay(0.1 * ms)

            with sequential:

                # lower FORT power
                self.dds_FORT.set(
                    frequency=self.f_FORT,
                    amplitude=self.ampl_FORT_blowaway)

                # set cooling light to be resonant with a free-space atom
                self.dds_cooling_DP.set(
                    frequency=self.f_cooling_DP_resonant_2_to_3,
                    amplitude=self.ampl_cooling_DP_MOT)

                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()

                if self.blowaway_light_off:
                    self.dds_AOM_A6.sw.off()
                else:
                    self.dds_cooling_DP.sw.on()

            delay(self.t_blowaway)

            # reset MOT power
            self.dds_cooling_DP.sw.off()
            self.dds_cooling_DP.set(
                frequency=self.f_cooling_DP_RO,
                amplitude=self.ampl_cooling_DP_MOT)

            # reset FORT power
            self.dds_FORT.set(
                frequency=self.f_FORT,
                amplitude=self.ampl_FORT_loading)

            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()
            self.ttl_repump_switch.off()  # turns on the RP AOM

        delay(0.1 * ms)

        ############################
        # take the second shot
        ############################

        # set coils for readout
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
            channels=self.coil_channels)
        delay(0.1 * ms)

        self.dds_cooling_DP.sw.on()
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
        counts2 = self.ttl0.count(t_gate_end)
        delay(1 * ms)
        self.dds_cooling_DP.sw.off()

        # effectively turn the FORT AOM off
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)
        # set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(2 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)

@kernel
def single_photons_experiment(self): # for debugging
    """
    self is the experiment instance to which ExperimentVariables are bound
    """

    self.core.break_realtime()

    # initialize these or ARTIQ complains
    counts = 0
    counts2 = 0

    for measurement in range(self.n_measurements):

        # make sure any beams that can mess up the fW measurement (e.g. OP, excitation) are turned off
        if self.enable_laser_feedback:
            if measurement % 10 == 0:
                self.laser_stabilizer.run()
                delay(1 * ms)
            self.dds_FORT.sw.on()
            self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)

        self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

        ############################
        # load the MOT
        ############################
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
            channels=self.coil_channels)
        self.ttl_repump_switch.off()  # turns on the RP AOM
        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        # wait for the MOT to load
        delay_mu(self.t_MOT_loading_mu)

        # load atom from a PGC phase
        if self.do_PGC_in_MOT:
            self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
            delay(self.t_PGC_in_MOT)

        # turn on the dipole trap and wait to load atoms
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_loading)
        delay_mu(self.t_FORT_loading_mu)

        # turn off the coils
        if not self.do_PGC_in_MOT:
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                channels=self.coil_channels)

        delay(3 * ms)  # should wait for the MOT to dissipate

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

        ############################
        # take the first shot
        ############################
        if not self.no_first_shot:
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

        delay(3 * ms)

        ############################
        # optical pumping phase - pumps atoms into F=1,m_F=0
        ############################

        if self.t_pumping > 0.0:
            self.ttl_repump_switch.on()  # turns off the RP AOM
            self.dds_cooling_DP.sw.off()  # no MOT light

            # set coils for pumping
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(0.1 * ms) # coil relaxation time

            with sequential:

                # lower FORT power
                self.dds_FORT.set(
                    frequency=self.f_FORT,
                    amplitude=self.ampl_FORT_OP)

                # this could be condensed but is left as is for clarity
                if self.pumping_light_off:
                    self.dds_D1_pumping_SP.sw.off()
                    self.dds_pumping_repump.sw.off()
                elif self.control_experiment and measurement % 2 == 0:
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

        delay(0.1 * ms)

        ###########################
        # excitation and photon collection
        ###########################

        # try to mitigate any leakage
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        # turn off the FORT
        self.dds_FORT.sw.off() # fully off, rather than frequency shift just in case

        # excite the atom
        self.dds_excitation.sw.pulse(30*ns)

        delay(self.t_photon_collection_delay)
        t_gate_end = self.ttl0.gate_rising(self.t_photon_collection_time)
        counts2 = self.ttl0.count(t_gate_end)

        delay(1 * ms) # ttl count consumes all the RTIO slack

        self.dds_FORT.sw.on()

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        # effectively turn the FORT AOM off
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)
        # set the cooling DP AOM to the MOT settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(2 * ms)

        # update the datasets
        if not self.no_first_shot:
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)

        # update the datasets
        self.append_to_dataset('photocounts2', counts2)
        self.append_to_dataset('photocounts2_current_iteration', counts2)