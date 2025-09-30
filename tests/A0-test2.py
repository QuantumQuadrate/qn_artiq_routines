atom_photon_parity_3_experiment:
self.measurement = 0  # advances in end_measurement
    while self.measurement < self.n_measurements:
        self.core.break_realtime()
        load_until_atom_smooth_FORT_recycle(self)
        first_shot(self)

        excitation_cycle = 0  ### just for initialization.
        if self.BothSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:
            delay(1000 * us)
            ############ optical pumping phase
            ############################## excitation phase
            for excitation_attempt in range(self.n_excitation_attempts):

                self.ttl_GRIN2_switch.off()  # turns on excitation

                ######### Using the edge_counter (works well):
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
                SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

                delay(2 * us)

                microwave_11_pulse
                MW_RF_pulse

                if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                    chopped_blow_away(self)
                    atom_parity_shot(self)
                    self.append_to_dataset('BothSPCMs_parity_RO', self.BothSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)

                    self.measurement += 1
                    break

            if self.measurement == self.n_measurements:
                break

            ############################# readout every self.atom_check_every_n to see if the atom survived
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(0.4 * ms)

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
                BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

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
                if BothSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

                else:
                    atom_loaded = False
            excitation_cycle +=1

        second_shot(self)
        self.measurement -= 1
        end_measurement(self)