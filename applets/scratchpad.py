@kernel
def single_photon_experiment_modified_for_readout_every_cyle(self):
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
    self.counts = 0
    self.counts2 = 0
    excitation_counts = 0
    excitation_counts1 = 0
    excitation_counts_array = [0]

    if self.record_every_shot:
        self.counts_by_cycle = [0] * (self.n_excitation_cycles + 2)

        #rtio_log("2nd_shot_block", 0) # todo: delete. for debugging.

    self.set_dataset(self.count_rate_dataset,
                     [0.0],
                     broadcast=True)

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

    self.measurement = 0
    while self.measurement < self.n_measurements:

        excitation_counts_array = [0] * self.n_excitation_cycles
        excitation_counts_array1 = [0] * self.n_excitation_cycles

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs

            # bug -- microwave dds is off after AOM feedback; not clear why yet. for now, just turn it back on
            if self.verify_OP_in_photon_experiment:
                self.dds_microwaves.sw.on()

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
            first_shot(self)
            if self.record_every_shot:
                self.counts_by_cycle[0] = self.counts

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

        delay(1*us)

        self.ttl_SPCM_gate.on()  # blocks the SPCM output - this is related to the atom readouts undercounting
        loop_start_mu = now_mu()

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)
        delay(0.4 * ms)  # coil relaxation time

        # this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # use ttl_excitation to swith on/off D1 or Exc light
        self.dds_excitation.sw.on()

        # don't use gate_rising. set the sensitivity, do the pumping and excitation sequence, and count the photons
        # after excitation attempts. we'll block the counter channel with an external switch so we don't get any clicks
        # during the OP phases.
        self.ttl_SPCM0._set_sensitivity(1)
        self.ttl_SPCM1._set_sensitivity(1)
        for excitaton_cycle in range(self.n_excitation_cycles): # todo: revert later

            delay(0.5*ms)

            # low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            # todo: make sure this is consistent with any updates in chopped_optical_pumping function
            ############################
            # optical pumping phase - pumps atoms into F=1,m_F=0
            ############################

            # todo: D1 feedback
            if self.t_pumping > 0.0:
                self.dds_excitation.set(frequency=self.f_excitation, amplitude=dB_to_V(0.0))
                # self.dds_excitation.set(frequency=self.f_excitation, amplitude=dB_to_V(1.0))
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM

                if not self.pumping_light_off:
                    self.dds_pumping_repump.sw.on()

                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(.1 * ms) # maybe this can be even shorter

                self.ttl_excitation_switch.off()

                with sequential:

                    delay(1 * us)

                    self.core_dma.playback_handle(op_dma_handle)
                    delay(self.t_depumping)

                    self.dds_D1_pumping_DP.sw.off()
                    self.dds_pumping_repump.sw.off() # turn the repump back on

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                self.ttl_excitation_switch.on()

                ############################
                # microwave phase - ONLY USED FOR VERIFYING OP.
                ############################

                if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_pulse)
                    self.ttl_microwave_switch.on()
                    delay(0.1*ms)

                ############################
                # blow-away phase - push out atoms in F=2 only
                ############################

                if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                    chopped_blow_away(self)

            ############################
            # excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            ############################
            self.dds_excitation.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.ttl_repump_switch.off()  # repump AOM is on for excitation

            now = now_mu()

            # messy and confusing. todo: try to improve this
            # at_mu(now+10)
            # self.ttl_repump_switch.off()  # repump AOM is on for excitation
            # mu_offset = 800 # accounts for various latencies
            # at_mu(now + mu_offset - 200) # make sure stuff is off, no more Raman photons from FORT
            # at_mu(now + mu_offset+100+self.gate_start_offset_mu) # allow for repump rise time and FORT after-pulsing
            # t_collect = now_mu()
            # t_excite = now_mu()
            # pulses_over_mu = 0
            # for attempt in range(self.n_excitation_attempts):
            #     at_mu(now + mu_offset + 201 + int(attempt * (self.t_excitation_pulse / ns + 100))
            #           + self.gate_start_offset_mu)
            #     self.dds_excitation.sw.pulse(self.t_excitation_pulse)
            #     at_mu(now + mu_offset + 741 + int(attempt * (self.t_excitation_pulse / ns + 100) -
            #                                       0.1*self.t_excitation_pulse / ns) +self.gate_start_offset_mu)
            #
            #     # fast switch to gate SPCM output - why am I using an external switch
            #     # instead of the gate input on the SPCM?
            #     self.ttl_SPCM_gate.off()
            #     delay(0.2*self.t_excitation_pulse+100*ns + self.gate_switch_offset)
            #     self.ttl_SPCM_gate.on()
            #
            #     pulses_over_mu = now_mu() # overwrite each loop iteration

            self.dds_FORT.sw.off()
            at_mu(now+150)
            self.ttl_excitation_switch.off()
            at_mu(now + 150 + self.gate_start_offset_mu)
            self.ttl_SPCM_gate.off()
            at_mu(now + 150 + int(self.t_excitation_pulse/ns))
            self.ttl_excitation_switch.on()

            at_mu(now + 150 + int(self.t_photon_collection_time / ns + self.gate_start_offset_mu))

            self.ttl_SPCM_gate.on()

            at_mu(now + 150 + int(self.t_photon_collection_time / ns))

            self.dds_FORT.sw.on()
            pulses_over_mu = now_mu()

            delay(1*us)
            self.ttl_repump_switch.on() # block MOT repump (and therefore also the excitation light)



            # todo: terrible way of doing this. set the sensitivity of the gate at the beginning of the loop.
            #  it will only register events when the SPCM switch lets events through.
            excitation_counts = self.ttl_SPCM0.count(
                pulses_over_mu)  # this is the number of clicks we got over n_excitation attempts
            excitation_counts1 = self.ttl_SPCM1.count(
                pulses_over_mu)  # this is the number of clicks we got over n_excitation attempts
            excitation_counts_array[excitaton_cycle] = excitation_counts
            excitation_counts_array1[excitaton_cycle] = excitation_counts1
            delay(0.1*ms) # ttl count consumes all the RTIO slack.
            # loop_over_mu = now_mu()

            # todo: delete
            # self.print_async("bottom of excitation loop",now_mu() - loop_start_mu)

            ############################
            # recooling phase
            ############################

            # # todo: use a specific detuning for this stage?
            delay(1*ms)
            if self.t_recooling > 0:

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                    channels=self.coil_channels)

                delay(0.4*ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                if not self.record_every_shot:
                    delay(self.t_recooling)
                else:
                    # recool_and_shot(self)
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.counts_by_cycle[excitaton_cycle + 1] = self.ttl_SPCM0_counter.fetch_count()
                    delay(0.1 * ms)

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
                    [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(0.4 * ms)  # coil relaxation time


        self.ttl_SPCM0._set_sensitivity(0) # close the gating window on the TTL channel
        self.dds_excitation.sw.off()

        delay(1*ms)

        # turn AOMs back on
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        # todo: delete
        # self.core.wait_until_mu(loop_over_mu+1000)
        # at_mu(loop_over_mu+1000)
        # self.print_async("out of the loop",now_mu() - loop_start_mu)

        delay(1*ms)

        #rtio_log("2nd_shot_block",1)
        # self.print_async("second readout",now_mu() - loop_start_mu) # todo: delete
        with sequential:

            self.ttl_SPCM_gate.off() # enables the SPCM
            self.ttl_repump_switch.off()  # turns the RP AOM on

            # take the second shot
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                channels=self.coil_channels)
            delay(0.1 * ms)


            second_shot(self)
            if self.record_every_shot:
                self.counts_by_cycle[self.n_excitation_cycles + 1] = self.counts2
        #rtio_log("2nd_shot_block",0)

        end_measurement(self)

        for val in excitation_counts_array:
            self.append_to_dataset('excitation_counts', val)
        for val in excitation_counts_array1:
            self.append_to_dataset('excitation_counts1', val)

        if self.record_every_shot:
            self.counts_by_measurement[self.measurement] = self.counts_by_cycle
            self.append_to_dataset('counts_by_iteration', self.counts_by_measurement)
        delay(10*ms)

    self.dds_FORT.sw.off()