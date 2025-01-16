"""
    This is a scratchpad just for keeping temporary stuff

    single photon experiment ex

"""


# todo: Before first_shot:
#       ttl_SPCM_gate.off() at initialize_hardware() in BaseExperiment

# first_shot starts with SPCM gate open

def first_shot(self):
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.counts = self.ttl_SPCM0.count(t_gate_end)


def first_shot_edge_counter(self):
    #
    t_gate_end = self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_first_shot)
    self.counts = self.ttl_SPCM0_counter.fetch_count()


# todo:   Before Excitation loop:
#         ttl_SPCM_gate.on()
#         self.ttl_SPCM0._set_sensitivity(1)
#         self.ttl_SPCM1._set_sensitivity(1)


def excitation_phase_in_loop(self):
    now = now_mu()

    at_mu(now + 150 + self.gate_start_offset_mu)
    self.ttl_SPCM_gate.off()

    at_mu(now + 150 + int(self.t_photon_collection_time / ns + self.gate_start_offset_mu))
    self.ttl_SPCM_gate.on()

    at_mu(now + 150 + int(self.t_photon_collection_time / ns))
    self.dds_FORT.sw.on()
    pulses_over_mu = now_mu()

    # this is the number of clicks we got over n_excitation attempts
    SPCM0_SinglePhoton = self.ttl_SPCM0.count(pulses_over_mu)
    SPCM1_SinglePhoton = self.ttl_SPCM1.count(pulses_over_mu)

    SPCM0_SinglePhoton_array[excitaton_cycle] = SPCM0_SinglePhoton
    SPCM1_SinglePhoton_array[excitaton_cycle] = SPCM1_SinglePhoton

# todo:     After Excitation phase_in loop, Recooling & Readout starts


def readout_every_cycle(self):
    t_SPCM_recool_and_shot_mu = 20000000

    now = now_mu()
    self.ttl_SPCM1._set_sensitivity(0)  # closing the gating window for SPCM1
    self.ttl_SPCM_gate.off()  # gate turned on

    at_mu(now + t_SPCM_recool_and_shot_mu)
    self.ttl_SPCM_gate.on()  # gate turned off

    after_shot = now_mu()

    every_shot_count = self.ttl_SPCM0.count(after_shot)
    self.ttl_SPCM1._set_sensitivity(1)  # opening the gating window for SPCM1

    readout_counts_array[excitaton_cycle] = every_shot_count



# todo: After Exciation Loop:
#         self.ttl_SPCM0._set_sensitivity(0) # close the gating window on the TTL channel
#         self.ttl_SPCM1._set_sensitivity(0)  # close the gating window on the TTL channel