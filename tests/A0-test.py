@kernel
def atom_photon_parity_7_experiment(self):
    """
    In parity_6_experiment, we lose most of the atoms in blowaway phase, even if the atoms are in f=1 states. In this experiment
    I am trying to debug this.

    """

    self.measurement = 0  # advances in end_measurement

    BothSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1] * self.n_measurements
    angle_780_QWP = [-1] * self.n_measurements

    while self.measurement < self.n_measurements:

        load_until_atom_smooth_FORT_recycle(self)

        first_shot(self)

        excitation_cycle = 0  ### just for initialization.

        if self.BothSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            CW_optical_pumping_node1(self)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            self.dds_FORT.sw.off()  ### turns FORT off
            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            chopped_blow_away(self)

            atom_parity_shot(self)

            BothSPCMs_parity_RO[self.measurement] = self.BothSPCMs_parity_RO

            self.measurement += 1

            if self.measurement == self.n_measurements:
                break

            ####################################### atom check
            if BothSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                atom_loaded = True
            else:
                atom_loaded = False

            ################################### atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.dds_cooling_DP.sw.on()

            excitation_cycle += 1

        second_shot(self)

        self.measurement -= 1
        end_measurement(self)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

    for i in range(self.n_measurements):
        self.append_to_dataset('BothSPCMs_parity_RO', BothSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
