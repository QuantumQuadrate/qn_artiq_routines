"""
Dummy experiment to be deleted
"""
@kernel
def record_chopped_blow_away(self):
    n_chop_cycles = int(self.t_blowaway/self.t_BA_chop_period + 0.5)
    self.core.reset()

    with self.core_dma.record("chopped_blow_away"):
        start = now_mu()
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

    self.core_dma.playback_handle(ba_dma_handle)

@kernel
def microwave_Rabi_2_experiment(self):

    if self.t_blowaway > 0.0:
        record_chopped_blow_away(self)
        delay(100 * ms)


    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.t_blowaway > 0.0:
            chopped_blow_away(self)

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        second_shot(self)




