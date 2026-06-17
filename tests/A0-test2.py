@kernel
def two_nodes_synchronization2(self):

    ##### Syncing by sending a TTL pulse from node1 to node2
    t_node1_ref_0 = -1
    t_node2_ref_0 = -1

    if self.which_node == "bob":
        t0 = now_mu()
        at_mu(t0 + 100)
        t_gate_end = self.ttl_node1_input1.gate_rising(10 * s)
        t_node2_ref_0 = self.ttl_node1_input1.timestamp_mu(t_gate_end)

        at_mu(t0 + 10000)
        self.ttl_node2_output1.on()

    else:
        self.ttl_node2_input1.sample_input()
        delay(1 * us)
        ttl_node2_input1_signal = int(self.ttl_node2_input1.sample_get())
        if ttl_node2_input1_signal:
            t_node1_ref_0 = now_mu()
            self.ttl_node1_output1.pulse(500 * ns)

    delay(1 * ms)
    t_node1_ref_1 = t_node1_ref_0 + self.core.seconds_to_mu(2 * ms) + 200
    t_node2_ref_1 = t_node2_ref_0 + self.core.seconds_to_mu(2 * ms)

    if self.which_node == "alice":
        at_mu(t_node1_ref_1)
    else:
        at_mu(t_node2_ref_1)

@kernel
def load_until_atom_in_both_nodes_recycle2(self):
    """
    Two-node combined-loading version. in progress ...

    todo: This technically does not do recycle - add this later

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

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        AllSPCMs_atom_check_loaded = 0  ### for initilization
        AllSPCMs_atom_check_not_loaded = 0

        shim_tune_runs = 0

        while True:
            two_nodes_synchronization2(self) ### Syncing the two nodes. Todo: to be written.
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
def Two_nodes_atom_loading_3_experiment(self):
    """
    Simple atom loading experiment in both nodes

    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        run_feedback_and_record_FORT_MM_power(self)

    ### For testing only. Todo: remove later
    # delay(10*ms)
    # self.ttl_pumping_repump_switch.off()
    # load_until_atom_in_both_nodes_together_recycle(self)

    ############################################################
    # Turn OFF this node's TTL to be ready for synchronization.
    ############################################################
    if self.which_node == "alice":
        self.ttl_node1_output1.off()
    else:
        self.ttl_node2_output1.off()
    delay(1 * ms)


    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        load_until_atom_in_both_nodes_recycle2(self)

        # Turn OFF this node's TTL to be ready for synchronization.
        if self.which_node == "alice":
            self.ttl_node1_output1.off()
        else:
            self.ttl_node2_output1.off()
        delay(1 * ms)



        ################### Syncing by sending a TTL pulse from node1 to node2
        two_nodes_synchronization2(self) ### Syncing the two nodes. Todo: to be written.
        ##########################################################




        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)



        ################### Syncing by sending a TTL pulse from node1 to node2
        two_nodes_synchronization2(self)  ### Syncing the two nodes. Todo: to be written.
        ##########################################################



        second_shot(self)

        end_measurement(self)