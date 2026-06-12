@kernel
def atom_photon_parity_6_experiment(self):
    """
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
    delay(100 * ms)

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
                CW_optical_pumping_node1(self)
                delay(10 * us)

            ### Changing the bias field.
            self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
                                  self.AX_volts_microwave, self.AY_volts_microwave],
                                 channels=self.coil_channels)
            delay(1 * ms)

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
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])

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
                    chopped_blow_away(self)

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
                    AllSPCMs_parity_RO[self.measurement] = self.BothSPCMs_parity_RO
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
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])

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
                    chopped_blow_away(self)

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
                    BothSPCMs_parity_RO[self.measurement] = self.BothSPCMs_parity_RO
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
            BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

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
            if BothSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
        self.append_to_dataset('BothSPCMs_parity_RO', BothSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_9_experiment(self):
    """
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
    delay(100 * ms)

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
        self.stabilizer_FORT.run(setpoint_index=1)  # the science setpoint
        self.laser_stabilizer.run()
        self.dds_microwaves.sw.on()

    self.measurement = 0  # advances in end_measurement

    BothSPCMs_parity_RO = [-1] * self.n_measurements
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

        if self.BothSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                CW_optical_pumping_node1(self)
                delay(10 * us)

            ### Changing the bias field.
            self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
                                  self.AX_volts_microwave, self.AY_volts_microwave],
                                 channels=self.coil_channels)
            delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])

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
                    chopped_blow_away(self)

                    self.ttl_SPCM0_logic.pulse(1 * us)
                    delay(1 * us)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    BothSPCMs_parity_RO[self.measurement] = self.BothSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0:
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])

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
                    chopped_blow_away(self)

                    self.ttl_SPCM0_logic.pulse(1 * us)
                    delay(1 * us)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    BothSPCMs_parity_RO[self.measurement] = self.BothSPCMs_parity_RO
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
            excitation_cycle += 1

        delay(1 * ms)
        second_shot(self)

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
        self.append_to_dataset('BothSPCMs_parity_RO', BothSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)
