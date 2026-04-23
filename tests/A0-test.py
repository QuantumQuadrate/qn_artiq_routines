
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

        BothSPCMs_atom_check_loaded = 0  ### for initilization
        BothSPCMs_atom_check_not_loaded = 0

        shim_tune_runs = 0

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
def record_chopped_blow_away(self):
    """

    :param self:
    :return:
    """

    # todo: change OP -> BA

    n_chop_cycles = int(self.t_blowaway/self.t_BA_chop_period + 0.5)
    # self.print_async("blowaway cycles:", n_chop_cycles)

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

    delay(0.1 * ms)

    # set coils for blowaway
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_blowaway, -self.AZ_bottom_volts_blowaway,
         self.AX_volts_blowaway, self.AY_volts_blowaway],
        channels=self.coil_channels)
    delay(0.3 * ms)

    with sequential:

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_blowaway, amplitude=self.ampl_cooling_DP_MOT)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        delay(0.1 * ms)

        if self.blowaway_light_off:
            self.dds_cooling_DP.sw.off()
            self.dds_AOM_A6.sw.off()
        else:
            # just turn the AOM up all the way. as long as we're 'saturating' the blowaway, it's okay if this doesn't
            # always give the same optical power
            self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-7.0))
            self.dds_AOM_A6.sw.on()
            self.dds_cooling_DP.sw.on()

    self.core_dma.playback_handle(ba_dma_handle)


    self.dds_cooling_DP.sw.off() ### Turns off cooling DP
    self.ttl_repump_switch.on() ### turns off MOT RP
    self.dds_AOM_A6.sw.off() ### turns off fiber AOM6

    delay(0.1 * ms)

    ### reset AOM RF powers
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
    delay(0.1*ms)

@kernel
def CW_optical_pumping_node1(self):
    """
    - Includes turning on or off some AOMs, or setting their powers, setting the coils, etc.
    - D1 is on GRIN1
    - I am bypassing GRIN1 AOM to increase the D1 power at the atom to speed up O.P.

    Akbar 2026-03-27

    """

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light

    ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    delay(5 * us)
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

    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    delay(5*us)
    # self.GRIN1and2_dds.sw.off()
    # self.ttl_GRIN1_switch.on()

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
        # Sampler_test(self)
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

    delay(1 * ms)

    advance = 1
    if self.__class__.__name__ != 'ExperimentCycler':
        if self.require_atom_loading_to_advance:
            if not self.SPCM0_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
                advance = 0

        if self.require_D1_lock_to_advance:
            t1_D1_checked = now_mu()
            while True:
                self.ttl_D1_lock_monitor.sample_input()
                delay(0.1 * ms)
                laser_locked = int(1 - self.ttl_D1_lock_monitor.sample_get())  ## this is 1 when the laser is locked, it is 0 otherwise.

                if laser_locked:
                    advance = 1
                    self.set_dataset("time_without_D1", 0.0, broadcast=True)  ### resetting time_without_D1 when the laser is locked
                    delay(10 * us)
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
    For ramping FORT from science setpoint to holding (microwave) and vice versa.
    Smoothly ramp FORT power using a quintic smoothstep profile. If t_FORT_ramp is too short (<1ms), it uses less
    number of steps to avoid Underflow errors. This can handle any t_FORT_ramp, from 1us to 10ms, for example.

    direction: "down" or "up"
    """

    assert (direction == "down" or direction == "up"), "Direction must be 'down' or 'up'"

    p_high = self.stabilizer_FORT.amplitudes[1]
    p_low = self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1]
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

@kernel
def single_photon_experiment_3_atom_loading_advance(self):
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

    BothSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

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

        BothSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
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

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

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
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light

                ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                delay(5 * us)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(5 * us)
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

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

                delay(2 * us)

                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                delay(1 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

                delay(5 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(10*us)

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
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

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                BothSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = BothSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if BothSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
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
        for val in BothSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('BothSPCMs_RO_atom_check', val)

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
def single_photon_experiment_4_atom_loading_advance(self):
    """
    This is similar to single_photon_experiment_3_atom_loading_advance but with modified OP to speed up the rate.

    Has bugs. Runs but with strange timing. Look at the atom loading time or BothSPCMs atom check in loading applet.





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

    BothSPCMs_RO_atom_check_array = [0]

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

        BothSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
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
            # self.core.break_realtime()

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
                BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                BothSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = BothSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if BothSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
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
        for val in BothSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('BothSPCMs_RO_atom_check', val)

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