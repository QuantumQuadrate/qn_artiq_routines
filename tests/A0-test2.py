@kernel
def record_chopped_RO(self):
    ### Akbar 2025-11-25
    t_FORT_ON = self.core.seconds_to_mu(1*us) ### the length of the FORT pulse
    t_RO = self.core.seconds_to_mu(0.8*us) ### duration of RO
    t_FORT_delay = self.core.seconds_to_mu(100*ns) ### delay wrt "t0" in each cycle
    t_RO_delay = self.core.seconds_to_mu(0*ns)

    with self.core_dma.record("chopped_RO"):
        for i in range(4000): ### number of cycles to repeat
            t0 = now_mu()

            at_mu(t0+t_FORT_delay)
            self.dds_FORT.sw.on()
            at_mu(t0 + t_FORT_delay + t_FORT_ON)
            self.dds_FORT.sw.off()

            at_mu(t0 + t_RO_delay)
            # self.ttl_GRIN1_switch.on() ### GRIN1_switch is simply used with set_sensitivity to represent that on scope.
            self.ttl0._set_sensitivity(1)
            self.ttl1._set_sensitivity(1)

            at_mu(t0 + t_RO_delay + t_RO)
            # self.ttl_GRIN1_switch.off()
            self.ttl0._set_sensitivity(0)
            self.ttl1._set_sensitivity(0)

            delay(1 * us)


@kernel
def first_shot(self):

    if self.use_chopped_readout:
        chopped_RO_handle = self.core_dma.get_handle("chopped_RO")
        delay(10 * ms)
        ### DMA playback
        self.core_dma.playback_handle(chopped_RO_handle)
        delay(10 * ms)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(now_mu())
        self.SPCM1_RO1 = self.ttl_SPCM1.count(now_mu())
        self.BothSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)
        delay(0.1 * ms)
        # self.dds_cooling_DP.sw.off()  ### turn off cooling
        # self.ttl_repump_switch.on()  ### turn off MOT RP
        # delay(10 * us)

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
def atom_loading_experiment(self):

    if self.use_chopped_readout:
        ### Akbar 2025-11-25
        ### record and get handle
        record_chopped_RO(self)
        delay(10 * ms)

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

        end_measurement(self)



