"""
Legacy functions we used in the past.

They were moved out of experiment_functions.py on 2026-10-20 (the six
optical-pumping helpers at the start of section 2 on 2026-10-23) and are kept
here only in case we need them again or want to look them up. Nothing imports
this file and it has no imports of its own, so none of these functions shows
up in a GeneralVariableScan menu or can run as it is; the master's repository
scan reports an import error for this file, which is expected. To use one
again, move it back into experiment_functions.py.
"""

###############################################################################
# 1. Test FUNCTIONS
###############################################################################


@kernel
def wait_for_ac_trigger(self, phase_delay_mu: TInt64) -> TInt64:
    """
    Wait for a 60 Hz zero-crossing rising edge and return a future
    RTIO timestamp having the requested phase delay.

    * ttl_node1_input2 is temporarily being used for AC trigger input

    ### add this line where you want to trigger to AC line
        ### 60Hz triggering
        phase_delay_mu = self.core.seconds_to_mu(self.dummy_variable)
        t_sequence_mu = wait_for_ac_trigger(self, phase_delay_mu)
        at_mu(t_sequence_mu)
    """

    T60_mu = self.core.seconds_to_mu(1.0 / 60.0) ##16.6ms
    margin_mu = self.core.seconds_to_mu(500 * us)

    self.core.break_realtime()

    t0 = np.int64(-1)

    # Use a gate shorter than one 60 Hz period, so that at most
    # one rising edge can be recorded in each gate.
    while t0 < 0:
        gate_end_mu = self.ttl_node1_input2.gate_rising(10 * ms)
        t0 = self.ttl_node1_input2.timestamp_mu(gate_end_mu)

        if t0 < 0:
            # The edge occurred outside this 10 ms window.
            # Restore RTIO slack and open another gate.
            self.core.break_realtime()

    # Use the next equivalent AC cycle, not the detected edge itself.
    t_start_mu = t0 + T60_mu + phase_delay_mu

    # Safety check: guarantee that the event remains sufficiently
    # far ahead of the hardware RTIO counter.
    while t_start_mu <= self.core.get_rtio_counter_mu() + margin_mu:
        t_start_mu += T60_mu

    return t_start_mu

@kernel
def test_RF_pulse_experiment(self):
    """
    testing the effect of RF pulse with external RF antenna on scope
    """
    self.core.reset()

    self.measurement = 0
    # self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.set(frequency=2.0*MHz, amplitude=dB_to_V(-3.0))

    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    while self.measurement < self.n_measurements:
        delay(100*us)
        self.dds_MW_RF.sw.on()  ### turn on RF
        self.dds_FORT.sw.on()
        delay(1*us)
        self.dds_MW_RF.sw.off()  ### turn off RF
        self.dds_FORT.sw.off()

@kernel
def test_Repump_pulse_experiment(self):
    """
    testing the response of repump aom
    """
    self.core.reset()

    self.measurement = 0
    # self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))

    self.dds_AOM_A6.sw.on()
    delay(1*ms)
    while self.measurement < self.n_measurements:
        delay(1000*ms)
        self.ttl_repump_switch.off()  # turns the MOT RP AOM on
        delay(1000*ms)
        self.ttl_repump_switch.on()  # turns the MOT RP AOM off

    self.dds_AOM_A6.sw.off()

@kernel
def test_BA_pulse_experiment(self):
    self.core.reset()

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(1 * s)
        self.dds_FORT.sw.on()
        delay(1*ms)

        self.core_dma.playback_handle(ba_dma_handle)

        delay(10*ms)

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

@kernel
def tune_coil_driver_experiment(self):
    self.core.reset()
    delay(1 * s)
    ### set coils for pumping
    for i in range(100000):
        # ### Set the coils to MOT loading setting
        # self.zotino0.set_dac(
        #     [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        #     channels=self.coil_channels)

        # self.zotino0.set_dac([0.0, 0.0, 0.0, 3.0], channels=self.coil_channels)

        self.zotino0.set_dac(
            [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
            channels=self.coil_channels)

        delay(10 * ms)

        # self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        delay(10 * ms)  # coil relaxation time

@kernel
def test_excitation_rise_time_experiment(self):

    """
    Simple code to test the rise/fall time of the excitation pulse.
    """
    self.core.reset()
    delay(1 * s)

    self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))

    # turning on excitation DDS and switch
    self.ttl_exc0_switch.off()  # EXC0 AOM ON
    self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet

    # delay(1* s)
    # for i in range(self.n_measurements):
    #
    #     t1 = now_mu()
    #     at_mu(t1 + 100)
    #     self.ttl_GRIN2_switch.off()        # turn on GRIN2 AOM
    #
    #     at_mu(t1 + 100+ 100)
    #     self.ttl_GRIN2_switch.on()          # turn off GRIN2 AOM
    #     # delay(4 * ms)  # coil relaxation time
    #     delay(1*us)

    delay(1* s)
    for i in range(10000000):
        self.ttl_GRIN2_switch.off()        # turn on GRIN2 AOM
        delay(30*ns)
        self.ttl_GRIN2_switch.on()          # turn off GRIN2 AOM
        # delay(4 * ms)  # coil relaxation time
        delay(100*us)

@kernel
def run_zotino_stabilization_Aqua(self):
    ####Start Kais's calibration code
    ref_voltage_V = 0.00000
    logger.debug('Hardware Setup')
    self.core.reset()
    self.zotino0.init()
    self.core.break_realtime()
    self.sampler1.init()
    self.core.break_realtime()
    delay(10 * ms)
    for q in [0,1,3,4,5,6,7]:
        delay(10 * ms)
        self.sampler1.set_gain_mu(q, 3)

    diff = [0.01] * 8
    MV = [0.0] * 8
    self.zotino0.init()
    delay(10 * ms)
    for m in range(8):
        delay(10 * ms)
        self.zotino0.write_dac(m, ref_voltage_V)
    self.zotino0.load()
    delay(10 * ms)
    for j in range(20):
        delay(10 * ms)
        n_channels = 8  # sets number of channels to read off
        smp = [0.0] * n_channels  # creates list of floating point variables
        self.sampler1.sample(smp)
        delay(10 * ms)
        # print (smp)
        for k in [0,1,3,4,5,6,7]:
            diff[k] = ref_voltage_V - smp[k]
            # print(diff[k])
        # print(k)

        delay(0.1)
        self.zotino0.init()
        self.core.break_realtime()
        delay(10 * ms)
        for n in [0,1,3,4,5,6,7]:
            MV[n] = MV[n] + diff[n]
            self.zotino0.write_offset(n, MV[n])
            # print(n)
        self.zotino0.load()
        delay(10 * ms)
    self.core.break_realtime()
    #####End Kais's calibration code

@kernel
def run_zotino_stabilization(self):
    """
    In the beginning of an experiment, measure the zotino offset and compensate it.

    :param self:
    :return:
    """

    ref_voltage_V = 5.00000
    self.core.reset()

    measurement_buf = np.array([0.0] * 8)

    measurement1 = 0.0  # 1
    measurement2 = 0.0  # 2
    measurement3 = 0.0  # 3
    measurement4 = 0.0  # 4
    measurement5 = 0.0
    measurement6 = 0.0
    measurement7 = 0.0
    measurement8 = 0.0

    self.zotino0.set_dac([ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V,
                          ref_voltage_V],
                         [8,9,10,11,12,13,14,15])

    delay(10*ms)

    avgs = 50
    for i in range(avgs):
        self.sampler2.sample(measurement_buf)

        delay(0.1 * ms)
        measurement1 += measurement_buf[0]
        measurement2 += measurement_buf[1]
        measurement3 += measurement_buf[2]
        measurement4 += measurement_buf[3]
        measurement5 += measurement_buf[4]
        measurement6 += measurement_buf[5]
        measurement7 += measurement_buf[6]
        measurement8 += measurement_buf[7]
        delay(0.1 * ms)

    measurement1 /= avgs
    measurement2 /= avgs
    measurement3 /= avgs
    measurement4 /= avgs
    measurement5 /= avgs
    measurement6 /= avgs
    measurement7 /= avgs
    measurement8 /= avgs

    self.set_dataset("zotino_test1_offset", measurement1 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test2_offset", measurement2 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test3_offset", measurement3 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test4_offset", measurement4 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test5_offset", measurement5 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test6_offset", measurement6 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test7_offset", measurement7 - ref_voltage_V, broadcast=True, persist=True)
    self.set_dataset("zotino_test8_offset", measurement8 - ref_voltage_V, broadcast=True, persist=True)

    self.append_to_dataset("zotino_test1_offset_monitor", measurement1)
    self.append_to_dataset("zotino_test2_offset_monitor", measurement2)
    self.append_to_dataset("zotino_test3_offset_monitor", measurement3)
    self.append_to_dataset("zotino_test4_offset_monitor", measurement4)
    self.append_to_dataset("zotino_test5_offset_monitor", measurement5)
    self.append_to_dataset("zotino_test6_offset_monitor", measurement6)
    self.append_to_dataset("zotino_test7_offset_monitor", measurement7)
    self.append_to_dataset("zotino_test8_offset_monitor", measurement8)

@kernel
def rotator_test_experiment(self):
    """
    ratator_test function to record the Sampler value while rotating the waveplate

    :param self:
    :return:
    """

    self.core.reset()
    delay(2 * s)

    # # GVS - set target_hwp_deg as the scan_variable.
    # move_to_target_deg(self, name="780_HWP", target_deg=self.target_hwp_deg)
    # wait_move(self, "780_HWP")
    # hwp780_pos = get_rotator_position(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_pos / self.deg_to_pos, ' deg')
    # delay(1*s)
    # record_PDA_power(self)

    move_to_target_deg(self, name="780_HWP", target_deg=10)
    hwp780_deg = get_rotator_deg(self, '780_HWP')
    delay(1 * s)
    self.print_async('hwp780 at ', hwp780_deg, ' deg')

    move_to_target_deg(self, name="780_HWP", target_deg=10)
    hwp780_deg = get_rotator_deg(self, '780_HWP')
    delay(1 * s)
    self.print_async('hwp780 at ', hwp780_deg, ' deg')

    # # delay(2*s)
    # # move_to_target_deg(self, name="780_HWP", target_deg=350)
    #
    # go_to_home(self, '780_HWP')
    #
    # hwp780_deg = get_rotator_deg(self, '780_HWP')
    # delay(1*s)
    # self.print_async('hwp780 at ', hwp780_deg, ' deg')

@kernel
def testing_record_chopped_FORT_and_TTL(self):
    ### Akbar 2025-11-25
    t_FORT_offset = self.core.seconds_to_mu(100 * ns)  ### delay wrt "t0" in each cycle
    t_FORT_OFF = self.core.seconds_to_mu(self.chop_RO_FORT_off) ### the length of the FORT OFF time

    t_cooling_offset = self.core.seconds_to_mu(self.chop_RO_pulse_offset) ### delay wrt "t0" in each cycle
    t_cooling_length = self.core.seconds_to_mu(self.chop_RO_pulse_length) ### duration of the cooling pulse

    t_gate_offset = self.core.seconds_to_mu(self.chop_RO_gate_offset) ### delay wrt t0 for gate to open
    t_gate_length = self.core.seconds_to_mu(self.chop_RO_gate_length) ### gate open time

    with self.core_dma.record("FORT_and_TTL_chopped"):
        for i in range(4000): ### number of cycles to repeat
            t0 = now_mu()

            at_mu(t0 + t_FORT_offset)
            self.dds_FORT.sw.off()
            at_mu(t0 + t_FORT_offset + t_FORT_OFF)
            self.dds_FORT.sw.on()

            at_mu(t0 + t_cooling_offset)
            self.dds_cooling_DP.sw.on()  ### Turn on cooling
            at_mu(t0 + t_cooling_offset + t_cooling_length)
            self.dds_cooling_DP.sw.off()  ### Turn off cooling

            at_mu(t0 + t_gate_offset)
            self.ttl0._set_sensitivity(1)
            self.ttl1._set_sensitivity(1)

            at_mu(t0 + t_gate_offset + t_gate_length)
            self.ttl0._set_sensitivity(0)
            self.ttl1._set_sensitivity(0)

            delay(0.5 * us)

@kernel
def testing_chopped_FORT_and_TTL_experiment(self):
    """
    xxx
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False # override experiment variable
    delay(1 * ms)
    self.ttl_SPCM0._set_sensitivity(0)
    self.ttl_SPCM1._set_sensitivity(0)

    ### record and get handle
    testing_record_chopped_FORT_and_TTL(self)
    delay(10 * ms)
    test_chopped_handle = self.core_dma.get_handle("FORT_and_TTL_chopped")
    delay(20 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        self.dds_cooling_DP.sw.off()  ### Turn off cooling

        delay(1 * ms)
        self.dds_FORT.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        delay(1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        delay(1 * ms)
        self.ttl_SPCM0._set_sensitivity(0)
        self.ttl_SPCM1._set_sensitivity(0)

        delay(1 * ms)

        ### DMA playback
        self.core_dma.playback_handle(test_chopped_handle)

        delay(10 * ms)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(now_mu())
        self.SPCM1_RO1 = self.ttl_SPCM1.count(now_mu())
        self.AllSPCMs_RO1 = int((self.SPCM0_RO1 + self.SPCM1_RO1) / 2)

        # delay(10*ms)
        self.core.break_realtime()

        end_measurement(self)


###############################################################################
# 1. ATOM LOADING FUNCTIONS
###############################################################################

@kernel
def load_MOT_and_FORT_for_Luca_scattering_measurement(self):
    """
    Modified version of load_MOT_and_FORT for imaging scattering in the chamber with the Luca
    and with an APD. Assumes the APD is connected to the Sampler card and channel specified in
    the code below.

    I typically run this with the MonitorFORTWithLuca experiment.

    The following shots are taken with the Luca:
    1. FORT scattering after feedback to loading setpoint
    2. FORT + MOT beam scattering

    :param self: the experiment instance
    :return:
    """

    self.dds_FORT.sw.on()

    if self.MOT_repump_off or self.MOT_light_off: # this is useful so that the SPCM0_RO1 dataset is all background
        self.ttl_repump_switch.on()  # turns the RP AOM off
    if self.MOT_light_off:
        self.dds_cooling_DP.sw.off()

    if not self.FORT_on_at_MOT_start:
        self.dds_FORT.sw.off()
    else:
        self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    # record FORT scattering with Luca and record Raman scattering from SM fiber
    self.ttl_Luca_trigger.pulse(5 * ms) # FORT loading scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_loading = self.ttl_SPCM0.count(t_gate_end)
    delay(10*us)

    self.APD_FORT_volts_loading = 0.0
    for i in range(self.APD_averages):
        self.sampler1.sample(self.APD_buffer)
        self.APD_FORT_volts_loading += self.APD_buffer[7]
        delay(1*ms)
    self.APD_FORT_volts_loading /= self.APD_averages

    delay(1 * ms)

    # set the cooling DP AOM to the MOT settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

    # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
    # delay(0.1 * ms)
    # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

    if not self.MOT_light_off:
        self.dds_cooling_DP.sw.on()

        delay(1*ms)

        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

    delay(1*ms) # if this delay is not here, the following line setting the dac doesn't execute

    # Turn on the MOT coils and cooling light
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
        channels=self.coil_channels)

    # wait for the MOT to load
    delay(self.t_MOT_loading/2)
    self.ttl_Luca_trigger.pulse(5 * ms) # total scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_and_MOT = self.ttl_SPCM0.count(t_gate_end)
    delay(1 * ms)
    delay(self.t_MOT_loading/2)

    # turn on the dipole trap and wait to load atoms if it is not already on
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

    if not self.FORT_on_at_MOT_start:
        delay_mu(self.t_FORT_loading_mu)
    delay(1*ms)
    if not self.no_feedback:
        self.stabilizer_FORT.run(setpoint_index=1)  # set to  science setpoint
        self.stabilizer_FORT.run(setpoint_index=2)  # the science holding setpoint

    self.dds_cooling_DP.sw.off()
    delay(1*ms)
    # record FORT scattering with Luca and record Raman scattering from SM fiber
    self.ttl_Luca_trigger.pulse(5 * ms)  # FORT loading scattering shot
    t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
    self.SPCM0_FORT_science = self.ttl_SPCM0.count(t_gate_end)
    delay(10*us)

    self.APD_FORT_volts_science = 0.0
    for i in range(self.APD_averages):
        self.sampler1.sample(self.APD_buffer)
        self.APD_FORT_volts_science += self.APD_buffer[7]
        delay(1*ms)

    self.APD_FORT_volts_science /= self.APD_averages
    delay(100*ms)  # this is 100 ms to ensure the Luca takes the next shot

    if not self.MOT_light_off:
        self.dds_cooling_DP.sw.on()

    if self.do_PGC_after_loading and self.t_PGC_after_loading > 0:

        if not self.no_feedback:
            self.dds_FORT.set(frequency=self.f_FORT,
                              amplitude=self.stabilizer_FORT.amplitudes[1])

        self.zotino0.set_dac([self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                             channels=self.coil_channels)

        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_PGC)
        if self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(self.t_PGC_after_loading)

    self.dds_cooling_DP.sw.off()

###############################################################################
# 2. SUBROUTINE FUNCTIONS
###############################################################################

@kernel
def record_chopped_optical_pumping(self):
    """
    D1: since 2026-10-19 Node1's D1 AOM runs from an external RF source
    switched by ttl_D1_pumping_switch (active low: off() = D1 on), not by the
    Urukul channel dds_GRIN2.

    :param self:
    :return:
    """

    n_chop_cycles = int(self.t_pumping/self.t_OP_chop_period + 0.5)
    assert n_chop_cycles >= 1, "t_pumping should be > t_OP_chop_period"

    # todo: use duty cycle ExperimentVariables
    OP_pulse = self.t_OP_chop_period * self.duty_cycle_OP
    FORT_pulse = self.t_OP_chop_period * self.duty_cycle_FORT

    n_depump_chop_cycles = int(self.t_depumping /self.t_OP_chop_period)

    self.core.reset()

    with self.core_dma.record("chopped_optical_pumping"):

        # hardcoded for offsets for now, but this gives decent separation between the FORT and OP
        # pulses
        start = now_mu()
        period_mu = self.core.seconds_to_mu(self.t_OP_chop_period)

        OP_pulse_length_mu = self.core.seconds_to_mu(OP_pulse)
        FORT_pulse_length_mu = self.core.seconds_to_mu(FORT_pulse)
        FORT_on_mu = self.core.seconds_to_mu(0.0)
        OP_on_mu = self.core.seconds_to_mu(self.t_OP_chop_offset)

        if not (self.pumping_light_off or self.D1_off_in_OP_phase):
            for i in range(n_chop_cycles):
                at_mu(start+i*period_mu+FORT_on_mu)
                self.dds_FORT.sw.off()
                delay_mu(OP_pulse_length_mu)
                self.dds_FORT.sw.on()
                at_mu(start+i*period_mu+OP_on_mu)
                self.ttl_D1_pumping_switch.off()  # D1 on
                delay_mu(OP_pulse_length_mu)
                self.ttl_D1_pumping_switch.on()  # D1 off
        else:
            for i in range(n_chop_cycles):
                at_mu(start + i * period_mu + FORT_on_mu)
                self.dds_FORT.sw.off()
                delay_mu(OP_pulse_length_mu)
                self.dds_FORT.sw.on()

        if n_depump_chop_cycles > 0:

            # turn off the pumping repump
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            self.ttl_pumping_repump_switch.on()
            start = now_mu() + self.core.seconds_to_mu(0.5 * us)

            if not (self.pumping_light_off or self.D1_off_in_depump_phase):
                for i in range(n_depump_chop_cycles):
                    at_mu(start + i * period_mu + FORT_on_mu)
                    self.dds_FORT.sw.off()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_FORT.sw.on()
                    at_mu(start + i * period_mu + OP_on_mu)
                    self.ttl_D1_pumping_switch.off()  # D1 on
                    delay_mu(OP_pulse_length_mu)
                    self.ttl_D1_pumping_switch.on()  # D1 off
            else:
                for i in range(n_depump_chop_cycles):
                    at_mu(start + i * period_mu + FORT_on_mu)
                    self.dds_FORT.sw.off()
                    delay_mu(OP_pulse_length_mu)
                    self.dds_FORT.sw.on()

@kernel
def chopped_optical_pumping(self):
    """
    Optical pumping with FORT chopped.

    Turns on:
        Fiber AOMs 5 and 6
        GRIN1and2 DDS
        GRIN1 AOM

    Turns off at the end:
        D1 DP
        Pumping Repumper
        Fiber AOMs 5 and 6
        GRIN1and2 DDS
        GRIN1 AOM
        cooling DP

    :param self:
    :return:
    """

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    if self.t_depumping + self.t_pumping > 1*ms:
        delay(1500 * us)  # we need extra slack

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light
    delay(10 * us)
    self.dds_AOM_A1.sw.off()
    self.dds_AOM_A2.sw.off()
    self.dds_AOM_A3.sw.off()
    self.dds_AOM_A4.sw.off()
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    delay(10 * us)

    if self.which_node == "alice":
        self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
        ### This is split into GRIN1 and GRIN2 in Alice
    else:
        self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(-5.0))
        ### This is just for GRIN1 in Bob
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))
    delay(1 * us)

    ### Tunring on pumping RP:
    self.ttl_pumping_repump_switch.off()
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    delay(1*ms)

    ### so that D1 can pass
    self.dds_GRIN1.sw.on()
    self.ttl_GRIN1_switch.off()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    delay(1 * ms)  # coil relaxation time

    self.core_dma.playback_handle(op_dma_handle)
    delay(self.t_depumping)
    self.dds_GRIN2.sw.off() ### turning off D1 DP
    self.ttl_pumping_repump_switch.on() ### turning off pumping RP

    delay(2*us)
    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    delay(100 * us)

    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    # self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
    delay(1*ms)
    self.dds_GRIN1.sw.off()
    self.ttl_GRIN1_switch.on()
    delay(10*us)
    # self.dds_cooling_DP.sw.off()

@kernel
def CW_optical_pumping_node1(self):
    """
    - Includes turning on or off some AOMs, or setting their powers, setting the coils, etc.
    - D1 is on GRIN1
    - I am bypassing GRIN1 AOM to increase the D1 power at the atom to speed up O.P.
    - 2026-10-19: D1 now runs from an external RF source switched by
      ttl_D1_pumping_switch and goes through the GRIN1 AOM; this function does
      not open GRIN1 yet.
    - Not used since the 2026-10-20 Node1 rewire (Node1's GRIN AOMs are wired
      like Node2's); kept for rollback.
    """

    self.ttl_repump_switch.on()  # turns off the MOT RP AOM
    self.ttl_exc0_switch.on() # turns off the excitation
    self.dds_cooling_DP.sw.off()  # no cooling light
    delay(5 * us)

    ### Turning on fiber AOMs for delivery of the pumping repump
    # self.dds_AOM_A1.set(frequency=self.AOM_A1_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A3.set(frequency=self.AOM_A3_freq,amplitude=dB_to_V(-5.0))
    # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))

    delay(5 * us)
    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    # ### so that D1 can pass
    # self.dds_GRIN1.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))
    # delay(5 * us)
    # self.dds_GRIN1.sw.on()
    # self.ttl_GRIN1_switch.off()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    ### Physical coil-relaxation time
    # delay(0.4 * ms)  # coil relaxation time
    delay(self.t_OP_coil_relax_time)

    ### Optical pumping phase
    ### D1: since 2026-10-19 Node1's D1 AOM runs from an external RF source
    ### switched by ttl_D1_pumping_switch (active low: off() = D1 on), not
    ### by the Urukul channel dds_GRIN2.
    self.ttl_pumping_repump_switch.off()
    self.ttl_D1_pumping_switch.off()  # D1 on

    delay(self.t_pumping)

    self.ttl_D1_pumping_switch.on()  # D1 off
    delay(self.t_depumping)
    self.ttl_pumping_repump_switch.on()

    delay(2*us)

    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()

    delay(1 * us)

    # self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=self.stabilizer_AOM_A1.amplitude)
    # self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=self.stabilizer_AOM_A2.amplitude)
    # self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)
    # self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=self.stabilizer_AOM_A4.amplitude)
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

    delay(5*us)
    # self.dds_GRIN1.sw.off()
    # self.ttl_GRIN1_switch.on()

@kernel
def CW_optical_pumping_node2(self):
    """
    *** previously named as "optical_pumping_GRIN1"
    *** OP with GRIN1
    optical pumping without chopping the FORT - swtching ON/OFF with TTL

    ** OP AOM is driven with external RF source.
    ** ONLY GRIN1 AOM is used to turn on/off the OP.

    Note: the dds channels used for GRIN1 and GRIN2 kept their Node1 names
    (GRIN1and2_dds, dds_D1_pumping_DP) until 2026-10-20; they are now named
    after what they drive here:
    dds_GRIN1   >>>>   GRIN1 dds
    dds_GRIN2   >>>>   GRIN2 dds

    Eunji - for Node2

    """
    self.dds_cooling_DP.sw.off()  # no MOT cooling light
    self.ttl_repump_switch.on()   # no MOT RP AOM
    self.ttl_exc0_switch.on()     # no excitation
    self.ttl_D1_pumping_switch.on()      # no D1
    delay(5 * us)

    ### Turning on fiber AOMs 5 & 6 for delivery of the pumping repump
    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq,amplitude=dB_to_V(-5.0))
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq,amplitude=dB_to_V(-5.0))
    delay(5 * us)

    self.dds_AOM_A5.sw.on()
    self.dds_AOM_A6.sw.on()

    ### set coils for pumping
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
        channels=self.coil_channels)
    ### Physical coil-relaxation time
    # delay(0.395 * ms)  # coil relaxation time
    delay(self.t_OP_coil_relax_time)

    ### so that D1 can pass
    self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
    self.dds_GRIN1.sw.on()  ## GRIN1 RF ON but not yet activated with external switch
    delay(5 * us)

    with parallel:
        self.ttl_pumping_repump_switch.off()  ## pumping repump ON
        self.ttl_D1_pumping_switch.off()  ##turning D1 ON
        self.ttl_GRIN1_switch.off()  ## turning D1 ON

    delay(self.t_pumping)     ## pumping time

    with parallel:
        self.ttl_D1_pumping_switch.on()  ## turning D1 OFF
        self.ttl_GRIN1_switch.on()  ## GRIN OFF

    delay(self.t_depumping)    ## depumping time

    self.ttl_pumping_repump_switch.on()  ## pumping repump OFF
    delay(2 * us)

    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1]) # back to science setpoint
    # delay(1*us)

    self.dds_AOM_A5.sw.off()
    self.dds_AOM_A6.sw.off()
    self.dds_GRIN1.sw.off()  ## GRIN1 RF OFF
    delay(1 * us)

    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
    delay(5 * us)

@kernel
def record_CW_optical_pumping_node1(self):
    """Not used since the 2026-10-20 Node1 rewire: both nodes play
    record_CW_optical_pumping_node2, which keeps Node1's A1-A4 repump. Kept
    for rollback.
    """
    with self.core_dma.record("CW_optical_pumping_node1"):
        ### Turn off unused light
        with parallel:
            self.ttl_repump_switch.on()
            self.ttl_exc0_switch.on()
            self.dds_cooling_DP.sw.off()

        ### Configure pumping-repump delivery AOMs
        ### These must remain sequential if they share the same Urukul SPI bus
        self.dds_AOM_A1.set(frequency=self.AOM_A1_freq,amplitude=dB_to_V(-5.0))
        self.dds_AOM_A2.set(frequency=self.AOM_A2_freq,amplitude=dB_to_V(-5.0))
        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq,amplitude=dB_to_V(-5.0))
        self.dds_AOM_A4.set(frequency=self.AOM_A4_freq,amplitude=dB_to_V(-5.0))
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

        with parallel:
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        ### Program the coil voltages forward in time
        self.zotino0.write_dac(self.coil_channels[0], self.AZ_bottom_volts_OP)
        self.zotino0.write_dac(self.coil_channels[1], -self.AZ_bottom_volts_OP)
        self.zotino0.write_dac(self.coil_channels[2], self.AX_volts_OP)
        self.zotino0.write_dac(self.coil_channels[3], self.AY_volts_OP)
        delay(1.5 * us)
        self.zotino0.load()

        ### Physical coil-relaxation time
        # delay(0.4 * ms)
        delay(self.t_OP_coil_relax_time)

        ### Optical pumping
        ### D1: since 2026-10-19 Node1's D1 AOM runs from an external RF source
        ### switched by ttl_D1_pumping_switch (active low: off() = D1 on), not
        ### by the Urukul channel dds_GRIN2.
        with parallel:
            self.ttl_pumping_repump_switch.off()
            self.ttl_D1_pumping_switch.off()  # D1 on

        delay(self.t_pumping)

        self.ttl_D1_pumping_switch.on()  # D1 off

        delay(self.t_depumping)

        ### Turn off all pumping RF/light paths
        with parallel:
            self.ttl_pumping_repump_switch.on()
            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A1.set(frequency=self.AOM_A1_freq, amplitude=self.stabilizer_AOM_A1.amplitude)
        self.dds_AOM_A2.set(frequency=self.AOM_A2_freq, amplitude=self.stabilizer_AOM_A2.amplitude)
        self.dds_AOM_A3.set(frequency=self.AOM_A3_freq, amplitude=self.stabilizer_AOM_A3.amplitude)
        self.dds_AOM_A4.set(frequency=self.AOM_A4_freq, amplitude=self.stabilizer_AOM_A4.amplitude)
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

@kernel
def record_CW_optical_pumping_node2_without_coil(self):
    with self.core_dma.record("CW_optical_pumping_node2_without_coil"):
        ### Turn off unused light
        with parallel:
            self.dds_cooling_DP.sw.off()
            self.ttl_repump_switch.on()
            self.ttl_exc0_switch.on()
            self.ttl_D1_pumping_switch.on()

        ### Configure pumping-repump delivery AOMs
        ### These must remain sequential if they share the same Urukul SPI bus
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))

        with parallel:
            self.dds_AOM_A5.sw.on()
            self.dds_AOM_A6.sw.on()

        # ### Program the coil voltages forward in time
        # self.zotino0.write_dac(self.coil_channels[0], self.AZ_bottom_volts_OP)
        # self.zotino0.write_dac(self.coil_channels[1], -self.AZ_bottom_volts_OP)
        # self.zotino0.write_dac(self.coil_channels[2], self.AX_volts_OP)
        # self.zotino0.write_dac(self.coil_channels[3], self.AY_volts_OP)
        # delay(1.5 * us)
        # self.zotino0.load()
        #
        # ### Physical coil-relaxation time
        # delay(0.395 * ms)
        # # delay(self.dummy_variable)

        ### Configure GRIN1 DDS
        self.dds_GRIN1.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        self.dds_GRIN1.sw.on()
        delay(2*us)

        ### Optical pumping
        with parallel:
            self.ttl_pumping_repump_switch.off()
            self.ttl_D1_pumping_switch.off()
            self.ttl_GRIN1_switch.off()

        delay(self.t_pumping)

        ### Turn off D1 pumping light
        with parallel:
            self.ttl_D1_pumping_switch.on()
            self.ttl_GRIN1_switch.on()

        delay(self.t_depumping)

        ### Turn off all pumping RF/light paths
        with parallel:
            self.ttl_pumping_repump_switch.on()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
            self.dds_GRIN1.sw.off()

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        delay(1.7 * us)  ### to match the timing with node1

@kernel
def test_shot(self):
    """
    Based on non-chopped atom readout. Can be used as a middle RO betweek RO1 and RO2 to opimize the readout settings.

    Turns on:
        Cooling DP
        MOT RP
        All 6 fiber AOMs

    Turns off at the end:
        Cooling DP
        MOT RP
    """
    self.ttl_repump_switch.off()  ### turn on MOT RP
    self.dds_cooling_DP.sw.on()  ### Turn on cooling
    delay(5 * us)

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
    delay(5 * us)

    with parallel:
        self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_test_shot)
        self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_test_shot)

    self.SPCM0_test_RO = self.ttl_SPCM0_counter.fetch_count()
    self.SPCM1_test_RO = self.ttl_SPCM1_counter.fetch_count()
    self.SPCM0_OtherNode_test_RO = self.ttl_SPCM0_OtherNode_counter.fetch_count()
    self.SPCM1_OtherNode_test_RO = self.ttl_SPCM1_OtherNode_counter.fetch_count()
    self.AllSPCMs_test_RO = self.SPCM0_test_RO + self.SPCM1_test_RO + self.SPCM0_OtherNode_test_RO + self.SPCM1_OtherNode_test_RO

    delay(10 * us)
    self.dds_cooling_DP.sw.off()  ### turn off cooling
    self.ttl_repump_switch.on()  ### turn off MOT RP
    delay(10 * us)


###############################################################################
# 3. EXPERIMENT FUNCTIONS
# These are the experiments we run, and the name of each should end with
# experiment in order to have it show up in GeneralVariableScans
###############################################################################

@kernel
def MOT_loading_experiment(self):
    self.core.reset()

    for measurement in range(self.n_measurements):
        self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs
        load_MOT_and_FORT(self)

@kernel
def single_photon_experiment_3_atom_loading_advance_node2_AllSPCMs(self):
    """

    This is based on load_MOT_and_FORT_until_atom_and_recycle.
    Does not check for atom after each excitation attempt:

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

    AllSPCMs_RO_atom_check_array = [0]

    delay(100 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles / self.atom_check_every_n)
        tStamps_t1 = [0.0] * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        delay(100 * ms)  ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)
            delay(1 * ms)

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################
        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### 1) use ttl_D1_pumping_switch/ttl_excitation to send D1/Exc to GRIN1/GRIN2
        ### 2) use ttl_GRIN1_switch/ttl_GRIN2_switch to send D1/EXC to atoms

        self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5 * us)
        self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 1  ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            delay(10 * ms)

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.
            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                CW_optical_pumping_node2(self)
                delay(10 * us)

                # ############ microwave phase - ONLY USED FOR VERIFYING OP.
                # if self.t_microwave_pulse > 0.0 and self.verify_OP_in_photon_experiment:
                #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                #     self.ttl_microwave_switch.off()  # todo: switching on with external RF switch creates a lag.
                #     delay(1 * us)
                #     t1 = now_mu()
                #     self.dds_microwaves.sw.on()
                #
                #     at_mu(t1 + int(self.t_microwave_pulse / ns))
                #     self.dds_microwaves.sw.off()
                #     self.ttl_microwave_switch.on()
                #
                #     delay(0.1 * ms)
                #     # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                # ############ blow-away phase - push out atoms in F=2 only
                # if self.t_blowaway > 0.0 and self.verify_OP_in_photon_experiment:
                #     self.core_dma.playback_handle(ba_dma_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(5 * us)
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                # # ############## previous code
                # # FORT_offset = 145
                # # FORT_offset = 90   ## as of 2026/03/08 Excelitas SPCMs
                # FORT_offset = 280   ## as of 2026/05/27 same Excelitas SPCMs; drifted.. why!
                #
                # t1 = now_mu()
                # at_mu(t1 + int(self.t_excitation_offset_mu))
                # self.ttl_GRIN2_switch.off()  # turns on excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                # self.ttl_GRIN2_switch.on()  # turns off excitation
                #
                # #todo: Turn FORT OFF earlier than 100ns. fall time is ~100ns
                # # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 145)  # turn off the FORT 145ns earlier than EXC pulse
                #
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # self.dds_FORT.sw.off()  ### turns FORT off
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset + int(self.t_photon_collection_time / ns)) # turning FORT ON 50ns later than gate end
                # self.dds_FORT.sw.on()  ### turns FORT on
                #
                # ## open SPCM gate 50ns earlier
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.gate_start_offset_mu) -50)  # turn the gate on 50ns earlier than EXC pulse
                # with parallel:
                #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                #
                # # ############## newer code 2026-05-28

                ###todo: 2026/05/28 Note - Eunji
                ### To make FORT, EXC, gate to happen at the same time, I need the offset to be:
                ### FORT at t_FORT
                ### EXC  at t_FORT - 280ns; t_EXC_hardware_offset_mu = -280
                ### GATE at t_FORT + 465ns;  t_GATE_hardware_offset_mu  = +465
                ### These variables do not change so I have it hard coded here.
                ### ### CASE - tested with this single photon code
                ### **** To turn EXC pulse on 50ns after FORT is turned off, I need to set
                ###      t_excitation_offset_mu = -280 + 50 = -230
                ### **** To turn GATE pulse on 25ns after FORT is turned off, I need to set
                ###      gate_start_offset_mu = +465 + 25 = 490 which is 50ns before excitation pulse
                ### **** To turn FORT ON after gate pulse is done,
                ###      Turn off at t_FORT t_photon_collection_time + 50(additional delay to account for gate delay)
                ############

                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on 50ns after photon collection ended

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation after t_excitation_offset_mu

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))  ## turns on gate after gate_start_offset_mu
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1


                # at_mu(t1 + 30000)
                tStamps_t1[
                    excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
                delay(0.4 * ms)  ### coils relaxation time

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()

                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()

                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(0.4 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())

                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

                delay(10 * us)
                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)

            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()

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

        delay(5 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle / self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1) * self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_Photon_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_Photon_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_Photon_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_Photon_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(15 * ms)

@kernel
def single_photon_experiment_3_atom_loading_advance_AllSPCM(self):
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

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        delay(1 * us)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        # ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        # ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        # self.GRIN1and2_dds.sw.on()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

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
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)

                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            # delay(10*us)

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)

            # self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            # delay(2 * ms)
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0 events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1 events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                # at_mu(t1 + 30000)
                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
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

                # with parallel:
                #     self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                #     self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                #
                # SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                # SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                # AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                # AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_3_atom_heat_test(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    This is similar to single_photon_experiment_3 that I used for photon generation, but for the purpose of optimizing the cycle.
    Here I am going to see how many excitation + OP cycle I can run (max_excitation_cycles) on the atom before loosing the atom
    without recooling, or keep max_excitation_cycles fixed, and change t_recooling and see how that affects the retention from
    the 2nd readout (RO2).

    Akbar 2026-03-18
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    # record_chopped_optical_pumping(self)
    # delay(100 * ms)
    #
    # if self.verify_OP_in_photon_experiment:
    #     record_chopped_blow_away(self)
    #     ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    #
    #     self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    #     delay(10 * ms)
    #     self.dds_microwaves.sw.on()
    #     delay(100 * ms)
    #
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
        # tStamps_t1 = [0.0] * (self.max_excitation_cycles * self.n_excitation_attempts)
        # SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        # SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 1  ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):

            self.core.break_realtime()

            ### low level pumping sequnce is more time efficient than the prepackaged chopped_optical_pumping function.

            # ############################### copped optical pumping phase - pumps atoms into F=1,m_F=0
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
            #     self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5V for optical pumping
            #     self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
            #     self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
            #     delay(1 * us)
            #
            #     ### Tunring on pumping RP:
            #     self.ttl_pumping_repump_switch.off()
            #     self.dds_AOM_A5.sw.on()
            #     self.dds_AOM_A6.sw.on()
            #
            #     delay(1 * us)
            #
            #     self.ttl_GRIN1_switch.off()  ### turn ON GRIN1 AOM for D1
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
            #     # delay(1 * us)
            #
            #     self.ttl_GRIN1_switch.on()  ### turn OFF GRIN1 AOM
            #     delay(10 * us)

            ### with cw pumping:
            if self.t_pumping > 0.0:
                delay(2 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(2 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(100 * us)

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

                # at_mu(t1 + int(self.gate_start_offset_mu))
                # with parallel:
                #     t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                #     t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                #
                # ### timestamping SPCM0 events
                # while SPCM0_click_counter < max_clicks:
                #     SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                #     if SPCM0_click_time == -1.0:
                #         break
                #     # SPCM0_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                #     #     SPCM0_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)
                #     SPCM0_click_counter += 1
                #
                # ### timestamping SPCM1 events
                # while SPCM1_click_counter < max_clicks:
                #     SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                #     if SPCM1_click_time == -1.0:
                #         break
                #     # SPCM1_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                #     #     SPCM1_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)
                #     SPCM1_click_counter += 1

                # at_mu(t1 + 30000)
                # tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC_optimization, -self.AZ_bottom_volts_PGC_optimization,
                     self.AX_volts_PGC_optimization, self.AY_volts_PGC_optimization],
                    channels=self.coil_channels)

                ampl_cooling_DP_PGC_optimization = self.ampl_cooling_DP_MOT * self.p_cooling_DP_PGC_optimization
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC_optimization, amplitude=ampl_cooling_DP_PGC_optimization)
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

            delay(10 * us)

        delay(1 * ms)

        self.GRIN1and2_dds.sw.off()

        delay(0.1 * ms)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        # delay(0.1 * ms)
        self.core.break_realtime()

        end_measurement(self)

        delay(5 * ms)

        delay(1 * ms)
        # for i in range((excitation_cycle + 1) * self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle+1)

        delay(1 * ms)


    delay(15 * ms)

@kernel
def single_photon_experiment_4_atom_loading_advance(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    This is similar to single_photon_experiment_3_atom_loading_advance but with modified OP to speed up the rate.
    If max_excitation_cycles is large (like 3000), the timestamp arrays get too large and artiq may freeze
    randomly after sometime. So, keep max_excitation_cycles <= 2000 and n_measurements <=100.

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

    AllSPCMs_RO_atom_check_array = [0]

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

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

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

                ### D1: since 2026-10-19 Node1's D1 AOM runs from an external RF
                ### source switched by ttl_D1_pumping_switch (active low: off() =
                ### D1 on), not by the Urukul channel dds_D1_pumping_DP. This
                ### experiment is Node1 only.
                self.ttl_pumping_repump_switch.off()
                self.ttl_D1_pumping_switch.off()  # D1 on

                delay(self.t_pumping)

                self.ttl_D1_pumping_switch.on()  # D1 off
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
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

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
def single_photon_experiment_5_atom_loading_advance_AllSPCM(self):
    """
    ### this is specifically for node1

    This is similar to single_photon_experiment_4_atom_loading_advance but I have reduced the size of the timestamp array
    to avoid freezing artiq. For this:
        1- I am not looking for 2nd clicks in each SPM window. Effectively, this is like max_clicks = 1.
        2- I am not appending the timestamps when a photon is not registered.
        3- I initialize and use a shorter array than max_excitation_cycles because we detect 5% of the times on average.


    IMPORTANT: for this function to work, you need to modify BaseExperiment to be like the following:
        # self.experiment.set_dataset("SPCM0_SinglePhoton_tStamps", [[0.0,0.0]], broadcast=True)
        # self.experiment.set_dataset("SPCM1_SinglePhoton_tStamps", [[0.0,0.0]], broadcast=True)
        self.experiment.set_dataset("SPCM0_SinglePhoton_tStamps", [0.0], broadcast=True)
        self.experiment.set_dataset("SPCM1_SinglePhoton_tStamps", [0.0], broadcast=True)

    """

    self.core.reset()
    delay(1 * ms)

    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

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

        Any_SPCM_click_counter = 0

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * short_timestamps_length
        SPCM0_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_timestamps = [-1.0]  * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0]  * short_timestamps_length

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light

        self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()

        # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

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

                ### D1: since 2026-10-19 Node1's D1 AOM runs from an external RF
                ### source switched by ttl_D1_pumping_switch (active low: off() =
                ### D1 on), not by the Urukul channel dds_D1_pumping_DP. This
                ### experiment is Node1 only.
                self.ttl_pumping_repump_switch.off()
                self.ttl_D1_pumping_switch.off()  # D1 on

                delay(self.t_pumping)

                self.ttl_D1_pumping_switch.on()  # D1 off
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
            self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation
                self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


                if SPCM0_click_time > 0.0 or SPCM1_click_time > 0.0 or  SPCM0_OtherNode_click_time > 0.0 or SPCM1_OtherNode_click_time > 0.0:
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0.0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)

                        if SPCM1_click_time > 0.0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)

                        if SPCM0_OtherNode_click_time > 0.0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)

                        if SPCM1_OtherNode_click_time > 0.0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)

                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

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
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                SPCM0_OtherNode_RO_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
                SPCM1_OtherNode_RO_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

                AllSPCMs_RO_atom_check = int(SPCM0_RO_atom_check + SPCM1_RO_atom_check + SPCM0_OtherNode_RO_atom_check + SPCM1_OtherNode_RO_atom_check)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                    self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        i = 0
        while i < len(tStamps_t1) and tStamps_t1[i] > 0.0:
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])
            i += 1


        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_6_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_5_atom_loading_advance. But works for any nodes.
    """

    self.core.reset()
    delay(1 * ms)

    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    delay(100*ms)

    while self.measurement < self.n_measurements:
        Any_SPCM_click_counter = 0

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * short_timestamps_length
        SPCM0_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_timestamps = [-1.0]  * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0]  * short_timestamps_length

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])

            # ### Changing the bias field to test the effect of Zeeman shift on the photons
            # self.zotino0.set_dac(
            #     [self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave, self.AX_volts_microwave, self.AY_volts_microwave],
            #     channels=self.coil_channels)


            # delay(2 * ms)
            # self.core.break_realtime()

            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


                if SPCM0_click_time > 0.0 or SPCM1_click_time > 0.0 or  SPCM0_OtherNode_click_time > 0.0 or SPCM1_OtherNode_click_time > 0.0:
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0.0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)

                        if SPCM1_click_time > 0.0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)

                        if SPCM0_OtherNode_click_time > 0.0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)

                        if SPCM1_OtherNode_click_time > 0.0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)

                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

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
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()
                else:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()

                delay(self.t_recooling)

                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()

                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
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
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                SPCM0_OtherNode_RO_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
                SPCM1_OtherNode_RO_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

                AllSPCMs_RO_atom_check = int(SPCM0_RO_atom_check + SPCM1_RO_atom_check + SPCM0_OtherNode_RO_atom_check + SPCM1_OtherNode_RO_atom_check)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    break

                delay(10 * us)
                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

            delay(10 * us)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        i = 0
        while i < len(tStamps_t1) and tStamps_t1[i] > 0.0:
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])
            i += 1

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_7_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_6_atom_loading_advance.
    But uses DMA for
    * OP
    * recooling
    """

    self.core.reset()
    delay(1 * ms)

    #### recording DMA
    record_CW_optical_pumping_node2(self)
    record_recooling(self)

    #### getting DMA handles
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
    recooling_dma_handle = self.core_dma.get_handle("recooling")

    delay(1 * ms)
    self.core.break_realtime()
    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    delay(100*ms)

    while self.measurement < self.n_measurements:
        Any_SPCM_click_counter = 0

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * short_timestamps_length
        SPCM0_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_timestamps = [-1.0]  * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0]  * short_timestamps_length

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # delay(1 * ms)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(5 * us)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    delay(10 * us)
                    CW_optical_pumping_node1(self)
                    delay(10 * us)
                else:
                    self.core_dma.playback_handle(CW_OP_node2_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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

                ######### time stamping the photons. Counting to be done in analysis.
                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)


                if SPCM0_click_time > 0.0 or SPCM1_click_time > 0.0 or  SPCM0_OtherNode_click_time > 0.0 or SPCM1_OtherNode_click_time > 0.0:
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0.0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)

                        if SPCM1_click_time > 0.0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)

                        if SPCM0_OtherNode_click_time > 0.0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)

                        if SPCM1_OtherNode_click_time > 0.0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)

                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

                delay(10 * us)

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.core_dma.playback_handle(recooling_dma_handle)

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
                    self.dds_AOM_A5.sw.on()
                    self.dds_AOM_A6.sw.on()

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                SPCM0_OtherNode_RO_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
                SPCM1_OtherNode_RO_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

                AllSPCMs_RO_atom_check = int(SPCM0_RO_atom_check + SPCM1_RO_atom_check + SPCM0_OtherNode_RO_atom_check + SPCM1_OtherNode_RO_atom_check)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    break

                delay(10 * us)
                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

            delay(10 * us)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
        delay(0.1 * ms)

        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
        delay(5 * us)

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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        i = 0
        while i < len(tStamps_t1) and tStamps_t1[i] > 0.0:
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])
            i += 1

        ### Number of compact timestamp rows saved during this measurement
        self.append_to_dataset('n_photon_events', Any_SPCM_click_counter)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def single_photon_experiment_8_atom_loading_advance_AllSPCM(self):
    """
    This is similar to single_photon_experiment_7_atom_loading_advance. with DMA in excitation phase
    But uses DMA for
    * OP
    * recooling
    * EXC
    """

    self.core.reset()
    delay(1 * ms)

    #### recording DMA
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)
    record_recooling(self)
    record_excitation_and_photon_collection(self)

    #### getting DMA handles
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
    recooling_dma_handle = self.core_dma.get_handle("recooling")
    excitation_dma_handle = self.core_dma.get_handle("excitation_and_photon_collection")

    delay(1 * ms)
    self.core.break_realtime()
    AllSPCMs_RO_atom_check_array = [0]

    ### short_timestamps_length is a fraction of the max_excitation_cycles. A photon is registered only in 5% of the
    ### excitation cycles. So no need to keep a large array. For assurance, we can use an array 10% of the max_excitation_cycles.
    short_timestamps_length = int(self.max_excitation_cycles * 0.1)

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    delay(100*ms)

    while self.measurement < self.n_measurements:
        Any_SPCM_click_counter = 0

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * short_timestamps_length
        SPCM0_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_timestamps = [-1.0]  * short_timestamps_length
        SPCM0_OtherNode_timestamps = [-1.0]  * short_timestamps_length
        SPCM1_OtherNode_timestamps = [-1.0]  * short_timestamps_length

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation
        # delay(1 * ms)

        # load_MOT_and_FORT_until_atom(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # delay(1 * ms)

            first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(5 * us)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            self.core.break_realtime()

            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    # CW_optical_pumping_node1(self)
                    self.core_dma.playback_handle(CW_OP_node1_handle)
                else:
                    self.core_dma.playback_handle(CW_OP_node2_handle)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                dma_start_mu = now_mu()
                self.core_dma.playback_handle(excitation_dma_handle)

                t1 = dma_start_mu + self.core.seconds_to_mu(5 * us)
                t_end_SPCM = t1 + int(self.gate_start_offset_mu) + self.core.seconds_to_mu(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM)

                if SPCM0_click_time > 0.0 or SPCM1_click_time > 0.0 or  SPCM0_OtherNode_click_time > 0.0 or SPCM1_OtherNode_click_time > 0.0:
                    if Any_SPCM_click_counter < short_timestamps_length:
                        if SPCM0_click_time > 0.0:
                            SPCM0_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_click_time)

                        if SPCM1_click_time > 0.0:
                            SPCM1_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_click_time)

                        if SPCM0_OtherNode_click_time > 0.0:
                            SPCM0_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)

                        if SPCM1_OtherNode_click_time > 0.0:
                            SPCM1_OtherNode_timestamps[Any_SPCM_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)

                        tStamps_t1[Any_SPCM_click_counter] = self.core.mu_to_seconds(t1)
                        Any_SPCM_click_counter += 1

                delay(5 * us)

            ############################ atom cooling phase with PGC settings
            if self.t_recooling > 0:
                self.core_dma.playback_handle(recooling_dma_handle)

            ############################# readout to see if the atom survived every self.atom_check_every_n
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                # self.zotino0.set_dac(
                #     [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                #     channels=self.coil_channels)
                #
                # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                # delay(0.4*ms)
                ### todo: this is inside recooling block. I am recooling every cycle so do not have set coils here

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

                delay(1 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                SPCM0_OtherNode_RO_atom_check = self.ttl_SPCM0_OtherNode_counter.fetch_count()
                SPCM1_OtherNode_RO_atom_check = self.ttl_SPCM1_OtherNode_counter.fetch_count()

                AllSPCMs_RO_atom_check = int(SPCM0_RO_atom_check + SPCM1_RO_atom_check + SPCM0_OtherNode_RO_atom_check + SPCM1_OtherNode_RO_atom_check)
                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.single_atom_threshold:
                    delay(100 * us)
                    break

                delay(10 * us)
                self.dds_cooling_DP.sw.off()
                self.ttl_repump_switch.on()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                delay(1 * us)
                if not self.PGC_and_RO_with_on_chip_beams:
                    self.dds_AOM_A5.sw.off()
                    self.dds_AOM_A6.sw.off()
                delay(1 * us)
                ##############################

            delay(5 * us)

        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()
        self.dds_D1_pumping_DP.sw.off()
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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        i = 0
        while i < len(tStamps_t1) and tStamps_t1[i] > 0.0:
            self.append_to_dataset('SPCM0_SinglePhoton_reduced_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_reduced_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_reduced_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_reduced_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])
            i += 1

        ### Number of compact timestamp rows saved during this measurement
        self.append_to_dataset('n_photon_events', Any_SPCM_click_counter)

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def atom_photon_parity_2_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Used to get nice parity oscillation before state rotation.

    A simple parity oscillation experiment. To speed up the experiment, I reuse the atom and
    do OP after say 5 excitation attempts and try excitation again. I repeat the loop as long as
    there is an atom measured every atom_check_every_n.
    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection of a photon (either from SPCM0 or 1) proceeds with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- Blow away F=2 and measure retention
    7- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    # AllSPCMs_RO_atom_check_array = [0]

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        # AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        # SPCM0_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        # SPCM1_SinglePhoton_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)

        self.core.break_realtime()
        # delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(1000 * us)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off() ### alowing D1 to pass
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on() ### blocking D1
                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

            for excitation_attempt in range(self.n_excitation_attempts):
                delay(200 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))

                ######### Using the edge_counter (works well):
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
                SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

                delay(20 * us)

                if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    delay(5 * us)
                    self.ttl_microwave_switch.off()
                    delay(self.t_microwave_11_pulse)
                    self.ttl_microwave_switch.on()
                    delay(5 * us)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(200 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    delay(20 * us)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(100*us)
            ############################ atom cooling phase with PGC settings
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
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

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


            ############################# readout every self.atom_check_every_n to see if the atom survived
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(1 * ms)

                self.dds_cooling_DP.sw.on()
                self.ttl_repump_switch.off()
                delay(1 * us)
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
                delay(10 * us)

                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                    self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

                SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

                delay(10*us)

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
                if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

                else:
                    atom_loaded = False

            excitation_cycle +=1


        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1*ms)
        self.GRIN1and2_dds.sw.off()

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        self.measurement -= 1
        end_measurement(self)


        delay(5 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python stops the loop one count earlier.
        # for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
        #     self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        # delay(1 * ms)
        # for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    delay(15 * ms)

@kernel
def atom_photon_parity_3_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Similar to parity_2_experiment, but with MW+RF pi/2 pulse to rotate to x-basis.

    To speed up the experiment, I reuse the atom and
    do OP after say 5 excitation attempts and try excitation again. I repeat the loop as long as
    there is an atom measured every atom_check_every_n.
    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    record_chopped_optical_pumping(self)
    delay(200*ms)

    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(10 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switches.

    # self.dds_MW_RF.set_phase_mode(PHASE_MODE_TRACKING)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    delay(1 * ms)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:
        SPCM0_SinglePhoton = 0
        SPCM1_SinglePhoton = 0

        self.core.break_realtime()
        # delay(100 * ms) ### with n_excitation_attempts = 5, 30ms delay is not enough

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(1*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(1000 * us)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
            delay(2*us)

            self.n_excitation_attempts = 1

            for excitation_attempt in range(self.n_excitation_attempts):
                delay(20 * us)
                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(20 * us)
                t1 = now_mu()

                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))

                ######### Using the edge_counter (works well):
                with parallel:
                    self.ttl_SPCM0_counter.gate_rising(self.t_photon_collection_time)
                    self.ttl_SPCM1_counter.gate_rising(self.t_photon_collection_time)
                SPCM0_SinglePhoton = self.ttl_SPCM0_counter.fetch_count()
                SPCM1_SinglePhoton = self.ttl_SPCM1_counter.fetch_count()

                delay(2 * us) ### 1us is not enough.
                self.ttl_microwave_switch.off()
                delay(self.t_microwave_11_pulse)
                self.ttl_microwave_switch.on()

                delay(2*us)

                if self.t_MW_RF_pulse > 0:
                    self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))
                    delay(5 * us)

                    with parallel:
                        self.ttl_microwave_switch.off()  ### turn on MW
                        self.dds_MW_RF.sw.on()  ### turn on RF

                    delay(self.t_MW_RF_pulse)

                    with parallel:
                        self.ttl_microwave_switch.on()  ### turn off MW
                        self.dds_MW_RF.sw.off()  ### turn off RF

                delay(2 * us)

                if SPCM0_SinglePhoton>0 or SPCM1_SinglePhoton>0:
                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(20 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(20 * us)
                    atom_parity_shot(self)
                    self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton)
                    self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton)
                    self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    delay(20 * us)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(100*us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            ############################ atom cooling phase with PGC settings
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
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

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


            ############################# readout every self.atom_check_every_n to see if the atom survived
            if (excitation_cycle + 1) % self.atom_check_every_n == 0:
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)

                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                delay(1 * ms)

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
                AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

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
                if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    atom_loaded = True

                else:
                    atom_loaded = False

            excitation_cycle +=1


        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1*ms)
        self.GRIN1and2_dds.sw.off()

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        self.measurement -= 1
        end_measurement(self)


        delay(5 * ms)

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python stops the loop one count earlier.
        # for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
        #     self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        # delay(1 * ms)
        # for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
        #     self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
        #     self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
        #     self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)
        delay(1 * ms)

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def atom_photon_parity_4_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Atom-photon parity oscillation in X-basis. Compared to #3, here I am timetagging the photons and starting the MW pulses
    at a fixed time relative to photon detection.

    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment

    """

    self.core.reset()
    delay(1 * ms)

    self.set_dataset("Slack_tracking", [0.0], broadcast=True)
    delay(1*ms)

    record_chopped_optical_pumping(self)
    delay(200*ms)

    record_chopped_blow_away(self)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on() ### turns on the DDS not the switch
    self.ttl_microwave_switch.on() ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    delay(1*ms)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

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

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
                self.ttl_microwave_switch.on()  ### close the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
            self.ttl_microwave_switch.on()  ### close the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()
        delay(10*us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment does not advance after a while!! I am guessing the excitation
            ### dds turns off in the middle of the first iteration, after 20-30 measurements!!
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(1 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(1 * us)

                ### Tunring on pumping RP:
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(10*us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

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

                at_mu(t1 + int(self.gate_start_offset_mu))

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)

                if SPCM0_click_time>0 and SPCM1_click_time<0:
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

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase = 0.0)

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase = 0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)



                    # self.core.break_realtime()
                    # self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    # self.append_to_dataset('SPCM0_SinglePhoton', 1)
                    # self.append_to_dataset('SPCM1_SinglePhoton', 0)
                    # self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    # self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    # delay(20 * us)
                    # self.core.break_realtime()


                    delay(1*ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time<0 and SPCM1_click_time>0:
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

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)



                    # self.core.break_realtime()
                    # self.append_to_dataset('AllSPCMs_parity_RO', self.AllSPCMs_parity_RO)
                    # self.append_to_dataset('SPCM0_SinglePhoton', 0)
                    # self.append_to_dataset('SPCM1_SinglePhoton', 1)
                    # self.append_to_dataset('angle_780_HWP', self.target_780_HWP)
                    # self.append_to_dataset('angle_780_QWP', self.target_780_QWP)
                    # delay(20 * us)
                    # self.core.break_realtime()

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

            if self.measurement == self.n_measurements:
                break

            delay(20*us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            ####################################### atom check
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(1 * ms)

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
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

            delay(1 * ms)

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
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

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

            excitation_cycle +=1


        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        second_shot(self)

        delay(1*ms)
        self.GRIN1and2_dds.sw.off()

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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_5_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Basically the same as parity_4_experiment with minor differences/improvements

    Atom-photon parity oscillation in X-basis. Compared to #3, here I am timetagging the photons and starting the MW pulses
    at a fixed time relative to photon detection.

    self.measurement advances only after single photon detection. Thus, n_measurements = 100, means
    100 detected photons from either SPCM0 or 1.

    1- Load an atom
    2- OP
    3- Excite
    4- gate the SPCMs, conditioned on detection from SPCM0 proceed with microwave mapping
    5- microwave mapping from F=1,mF=1 to F=2,mF=1
    6- apply a pi/2 MW+RF two-photon pulse
    7- Blow away F=2 and measure retention
    8- Change the 780 waveplates in GVS and repeat the experiment


    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_optical_pumping(self)
    delay(100 * ms)

    record_chopped_blow_away(self)

    op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")
    delay(10 * ms)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_ABSOLUTE)
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

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
                self.ttl_microwave_switch.on()  ### close the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
            self.ttl_microwave_switch.on()  ### close the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

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
        self.GRIN1and2_dds.sw.on()

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:

            delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment freezes and does not advance after a while!!
            ### The problem is in get_handle in the function that messes up with RTIO timeline in the loop.
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)
            #     self.GRIN1and2_dds.sw.on() ### turning back on for excitation after chopped_optical_pumping

            if self.t_pumping > 0.0:
                self.ttl_repump_switch.on()  # turns off the MOT RP AOM
                self.ttl_exc0_switch.on()  # turns off the excitation
                self.dds_cooling_DP.sw.off()  # no cooling light
                delay(10 * us)
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(10 * us)

                ### set coils for pumping
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_OP, -self.AZ_bottom_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                    channels=self.coil_channels)
                delay(1 * ms)  # coil relaxation time. 0.4ms was not enough based on oscilloscope.

                self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(5.0))  ### set to 5dBm for optical pumping
                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=dB_to_V(-5.0))
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=dB_to_V(-5.0))
                delay(10 * us)

                ### Tunring on pumping RP:
                self.ttl_pumping_repump_switch.off()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()

                delay(1 * ms)

                self.ttl_GRIN1_switch.off()
                delay(10 * us)

                self.core_dma.playback_handle(op_dma_handle)
                delay(self.t_depumping)
                self.dds_D1_pumping_DP.sw.off()  ### turning off D1 DP
                self.ttl_pumping_repump_switch.on()  ### turning off pumping RP

                delay(2 * us)
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()
                delay(100 * us)

                self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
                self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)
                delay(1 * ms)

                self.ttl_GRIN1_switch.on()
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))

            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(10 * us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

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

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
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

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 3000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 25000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    FORT_ramp2_smoothstep(self, direction="up")
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
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
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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
            delay(10 * us)

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)

            SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

            delay(1 * ms)

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
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
        self.GRIN1and2_dds.sw.off()

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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_6_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

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

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)
                self.ttl_microwave_switch.on()  ### close the switch
                delay(20*us)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)
            self.ttl_microwave_switch.on()  ### close the switch
            delay(20*us)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

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
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

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
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

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
                    self.core_dma.playback_handle(ba_dma_handle)

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
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
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
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

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
                    self.core_dma.playback_handle(ba_dma_handle)

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
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
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
            AllSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)

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
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_6_node2_AllSPCM_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    Trying to reduce the delay after photon detection. Do we need to set all the phases to 0.0?

    """

    self.core.reset()
    delay(100 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_blow_away(self)
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    # ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    # delay(10 * ms)

    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.off()  ### turns off the DDS
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
        self.dds_microwaves.sw.off()

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    SPCM0_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1.0] * self.n_measurements
    angle_780_QWP = [-1.0] * self.n_measurements

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()
        delay(100 * ms)

        self.ttl_exc0_switch.on()  # turns off the excitation
        self.dds_microwaves.sw.off()  ### turns off the DDS not the switch
        self.ttl_microwave_switch.on()  ### close the switch
        delay(10 * ms)

        # load_MOT_and_FORT_until_atom_recycle(self)
        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(1 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(1 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        ### 1) use ttl_D1_pumping_switch/ttl_excitation to send D1/Exc to GRIN1/GRIN2
        ### 2) use ttl_GRIN1_switch/ttl_GRIN2_switch to send D1/EXC to atoms

        # self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        # self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        # self.GRIN1and2_dds.sw.on()      # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        # self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        # delay(10*ms)
        self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
        delay(5*us)
        self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
        delay(5 * us)
        self.GRIN1and2_dds.sw.on()      # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)
        self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
        delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            atom_loaded = True
        else:
            atom_loaded = False

        while atom_loaded:
            delay(1 *ms)
            # delay(20 * ms)

            ############################### optical pumping phase - pumps atoms into F=1,m_F=0
            ### strange that with chopped_optical_pumping function, the experiment freezes and does not advance after a while!!
            ### The problem is in get_handle in the function that messes up with RTIO timeline in the loop.
            # if self.t_pumping > 0.0:
            #     chopped_optical_pumping(self)
            #     delay(1 * ms)
            #     self.GRIN1and2_dds.sw.on() ### turning back on for excitation after chopped_optical_pumping

            if self.t_pumping > 0.0:
                CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
            delay(5*us)
            self.core.break_realtime()

            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(100 * us)

            for excitation_attempt in range(self.n_excitation_attempts):
                # delay(100 * us)
                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()
                # todo: node2- change this to turning on the ttl;
                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)
                # #### node2
                # FORT_offset = 90   ## as of 2026/03/08 Excelitas SPCMs
                #
                # t1 = now_mu()
                # at_mu(t1 + int(self.t_excitation_offset_mu))
                # self.ttl_GRIN2_switch.off()  # turns on excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                # self.ttl_GRIN2_switch.on()  # turns off excitation
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset - 100)  # turn off the FORT 100ns earlier than EXC pulse
                # self.dds_FORT.sw.off()  ### turns FORT off
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + FORT_offset + int(self.t_photon_collection_time / ns))
                # self.dds_FORT.sw.on()  ### turns FORT on
                #
                # at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.gate_start_offset_mu))  # turn the gate on simultaneously with the EXC pulse
                t1 = now_mu()
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
                self.dds_FORT.sw.on()  ### turns FORT on 50ns after photon collection ended

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation after t_excitation_offset_mu

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                at_mu(t1 + int(self.gate_start_offset_mu))  ## turns on gate after gate_start_offset_mu

                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)

                # self.ttl_microwave_switch.off()
                # delay(5*us)
                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()   #### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()  #### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF


                    delay(1*ms)
                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(100 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:

                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()   #### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()  #### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF


                    delay(1*ms)
                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(100 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:

                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.ttl_microwave_switch.off()
                    self.dds_microwaves.sw.on()
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + int(self.t_microwave_11_pulse / ns))
                    # self.ttl_microwave_switch.on()
                    self.dds_microwaves.sw.off()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 5000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds,
                                                amplitude=dB_to_V(self.p_microwaves))

                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 10000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(
                            SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            # self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_microwaves.sw.on()
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(
                            self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            # self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_microwaves.sw.off()
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    self.ttl_microwave_switch.on()  ### turn off MW

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(10 * ms)

                    self.measurement += 1

                    break

                # print("what happens here")
            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            # FORT_ramp2_smoothstep(self, direction="up")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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
            delay(10 * us)

            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)
                self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_second_shot)

            # SPCM0_RO_atom_check = self.ttl_SPCM0_counter.fetch_count()
            # SPCM1_RO_atom_check = self.ttl_SPCM1_counter.fetch_count()
            # BothSPCMs_RO_atom_check = int((SPCM0_RO_atom_check + SPCM1_RO_atom_check) / 2)
            AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                      self.ttl_SPCM1_counter.fetch_count() + \
                                      self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                      self.ttl_SPCM1_OtherNode_counter.fetch_count())


            delay(1 * ms)

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
            if AllSPCMs_RO_atom_check / self.t_SPCM_second_shot > self.single_atom_threshold:
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

        delay(10 * ms)
        self.GRIN1and2_dds.sw.off()

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

    delay(50 * ms)
    self.core.break_realtime()
    for i in range(self.n_measurements):
        delay(10*ms)
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton[i])
        self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])

    delay(50 * ms)

@kernel
def atom_photon_parity_7_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    In parity_6_experiment, we lose most of the atoms in blowaway phase. In this experiment
    I am trying to debug this. This experiment does not care about the photon detection
    to run faster.

    """

    self.core.reset()
    delay(1 * ms)

    move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    delay(10 * ms)
    self.core.reset()

    record_chopped_blow_away(self)

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

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                delay(1 * ms)
                self.ttl_microwave_switch.on()  ### close the switch
                delay(20 * us)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            delay(1 * ms)
            self.ttl_microwave_switch.on()  ### close the switch
            delay(20 * us)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

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

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            # ### Changing the bias field.
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # for excitation_attempt in range(self.n_excitation_attempts):
            slack = now_mu() - self.core.get_rtio_counter_mu()
            if slack < 1e5:
                # self.print_async("slack added in measurement:", self.measurement)
                self.core.break_realtime()

            self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
            delay(5 * us)
            # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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


            delay(self.t_start_MW_mapping_mu*ns)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

            # delay(self.t_start_MW_mapping_mu*ns + 1000*ns)
            # self.ttl_microwave_switch.off()
            # delay(self.t_start_MW_mapping_mu*ns + 1000*ns + self.t_microwave_11_pulse)
            # self.ttl_microwave_switch.on()

            ############################ blow-away phase - push out atoms in F=2 only
            delay(20 * us)
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

            delay(10 * us)
            self.core_dma.playback_handle(ba_dma_handle)

            delay(10 * us)
            atom_parity_shot(self)

            delay(1 * ms)
            AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
            SPCM0_SinglePhoton[self.measurement] = 1.0
            SPCM1_SinglePhoton[self.measurement] = 1.0
            angle_780_HWP[self.measurement] = self.target_780_HWP
            angle_780_QWP[self.measurement] = self.target_780_QWP
            delay(1 * ms)

            self.measurement += 1

            if self.measurement == self.n_measurements:
                break

            delay(20 * us)
            FORT_ramp2_smoothstep(self, direction="up")

            ### stopping the excitation cycle after the atom is lost
            if self.AllSPCMs_parity_RO / self.t_SPCM_second_shot > self.single_atom_threshold:
                delay(100 * us)  ### Needs a delay of about 100us or maybe less
                # atom_loaded = True
                atom_loaded = False

            else:
                atom_loaded = False

            # ################################### atom cooling phase with PGC settings
            # if self.t_recooling > 0:
            #     self.zotino0.set_dac(
            #         [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            #         channels=self.coil_channels)
            #
            #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)
            #     delay(1 * ms)  ### coils relaxation time
            #
            #     self.dds_cooling_DP.sw.on()
            #     self.ttl_repump_switch.off()
            #     delay(1 * us)
            #     self.dds_AOM_A1.sw.on()
            #     self.dds_AOM_A2.sw.on()
            #     self.dds_AOM_A3.sw.on()
            #     self.dds_AOM_A4.sw.on()
            #     if not self.PGC_and_RO_with_on_chip_beams:
            #         self.dds_AOM_A5.sw.on()
            #         self.dds_AOM_A6.sw.on()
            #     else:
            #         self.dds_AOM_A5.sw.off()
            #         self.dds_AOM_A6.sw.off()
            #
            #     delay(self.t_recooling)
            #
            #     self.dds_cooling_DP.sw.off()
            #     self.ttl_repump_switch.on()
            #     delay(1 * us)
            #     self.dds_AOM_A1.sw.off()
            #     self.dds_AOM_A2.sw.off()
            #     self.dds_AOM_A3.sw.off()
            #     self.dds_AOM_A4.sw.off()
            #     self.dds_AOM_A5.sw.off()
            #     self.dds_AOM_A6.sw.off()
            #     delay(1 * us)

            excitation_cycle += 1

        # delay(20 * us)
        # self.ttl_exc0_switch.on()  # block Excitation

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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_8_experiment(self):
    """
    TODO: MUST add _OtherNode SPCMs now that we have a 50:50 BS in the BSA setup.

    this is yet another debugging experiment similar to atom_photon_parity_7_experiment, but based on Rabi_2
    with added excitation. The goal is to see why blowaway removes the atoms from the trap after excitation.
    All atoms should be in f=1 manifold and should survive the blowaway.
    """

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0
    self.SPCM1_RO1 = 0
    self.SPCM1_RO2 = 0

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    record_chopped_blow_away(self)

    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")

    self.core.break_realtime()

    if self.enable_laser_feedback:
        ### todo: set cooling_DP frequency to MOT loading in the stabilizer.
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    # delay(1 * ms)
    # move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
    # move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)
    # delay(10 * ms)
    # self.core.reset()

    self.dds_microwaves.set(frequency=self.f_microwaves_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()
    delay(1 * ms)

    self.ttl_GRIN2_switch.on()  # turns off excitation switch
    self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
    delay(5 * us)
    self.GRIN1and2_dds.sw.on()

    self.measurement = 0
    while self.measurement < self.n_measurements:
        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)

            first_shot(self)

        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        ### with cw pumping:
        if self.t_pumping > 0.0:
            delay (10 * us)
            if self.which_node == "alice":
                CW_optical_pumping_node1(self)
            else:
                CW_optical_pumping_node2(self)
            delay(10*us)

        ############################
        # excitation phase
        ############################

        self.GRIN1and2_dds.sw.on() ### require this here to turn on the channel. It will be off without this.

        self.ttl_GRIN2_switch.on()
        delay(5*us)

        # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
        # delay(1*us)
        # delay(self.t_delay_between_shots)
        # self.ttl_exc0_switch.on()  # turns off the excitation0 AOM

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

        ############################
        # microwave phase
        ############################
        delay(10*us)
        if self.t_microwave_pulse > 0.0:
            # ### Changing the bias field
            # self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                       self.AX_volts_microwave, self.AY_volts_microwave],
            #                      channels=self.coil_channels)
            # delay(1 * ms)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            # FORT_ramp2_smoothstep(self, direction="down")
            delay(1 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(20 * us)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            # FORT_ramp2_smoothstep(self, direction="up")
            delay(10 * us)

        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        end_measurement(self)
        delay(5 * ms)  ### hopefully to avoid underflow.

    delay(1 * ms)
    self.dds_microwaves.sw.off()

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def atom_photon_parity_9_AllSPCM_experiment(self):
    """
    *** compatible with either node1 or node2

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
    SPCM0_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0
    angle_780_QWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                self.ttl_microwave_switch.on()  ### close the switch
                delay(20*us)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            self.ttl_microwave_switch.on()  ### close the switch
            delay(20*us)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                delay(10 * us)
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10 * us)

            # if self.which_node == "alice":
            #     ### Changing the bias field.
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            #     delay(1 * ms)

            #
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            # delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                slack = now_mu() - self.core.get_rtio_counter_mu()
                if slack < 1e5:
                    # self.print_async("slack added in measurement:", self.measurement)
                    self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                # delay(5 * us)
                self.ttl_exc0_switch.off()  # turns on the excitation0 AOM
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
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM0)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM1)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)

                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM0_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

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
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM1_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

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
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM0_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10 * us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break

                if SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu)
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000)
                    self.ttl_microwave_switch.off()
                    at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 1000 + int(self.t_microwave_11_pulse / ns))
                    self.ttl_microwave_switch.on()

                    if self.t_MW_RF_pulse > 0:
                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 6000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_microwaves.set(frequency=self.f_microwaves_m11_dds, amplitude=dB_to_V(self.p_microwaves))

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 15000 + int(self.t_microwave_11_pulse / ns))
                        self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds), phase=0.0)

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns))
                        with parallel:
                            self.ttl_microwave_switch.off()  ### turn on MW
                            self.dds_MW_RF.sw.on()  ### turn on RF

                        at_mu(SPCM1_OtherNode_click_time + self.t_start_MW_mapping_mu + 20000 + int(self.t_microwave_11_pulse / ns) +
                              int(self.t_MW_RF_pulse / ns))

                        with parallel:
                            self.ttl_microwave_switch.on()  ### turn off MW
                            self.dds_MW_RF.sw.off()  ### turn off RF

                    ############################ blow-away phase - push out atoms in F=2 only
                    delay(10*us)
                    # FORT_ramp2_smoothstep(self, direction="up")
                    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    delay(10 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(10 * us)
                    atom_parity_shot(self)

                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 0.0
                    SPCM1_SinglePhoton[self.measurement] = 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP
                    delay(1 * ms)

                    self.measurement += 1

                    break


            delay(20 * us)
            excitation_cycle += 1

        delay(1 * ms)
        second_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        # self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        # self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton[i])
        self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(50 * ms)

@kernel
def atom_photon_parity_10_AllSPCM_experiment(self):
    """
    *** compatible with either node1 or node2
    *** DMA

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

    if self.enable_laser_feedback:
        delay(0.1 * ms)
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_RO setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)
        self.dds_microwaves.sw.on()

    self.core.reset()
    self.dds_microwaves.set_phase_mode(PHASE_MODE_CONTINUOUS)
    self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
    delay(1 * ms)
    self.dds_microwaves.sw.on()  ### turns on the DDS not the switch
    self.ttl_microwave_switch.on()  ### close the switch

    self.dds_MW_RF.set_phase_mode(PHASE_MODE_ABSOLUTE)
    self.dds_MW_RF.set(frequency=self.f_MW_RF_dds, amplitude=dB_to_V(self.p_MW_RF_dds))
    self.dds_MW_RF.sw.off()
    delay(1 * ms)

    #### recording DMA
    record_chopped_blow_away(self)
    record_CW_optical_pumping_node1(self)
    record_CW_optical_pumping_node2(self)
    # record_recooling(self)
    record_excitation_and_photon_collection(self)
    record_parity_mapping(self)

    #### getting DMA handles
    ba_dma_handle = self.core_dma.get_handle("chopped_blow_away")
    CW_OP_node1_handle = self.core_dma.get_handle("CW_optical_pumping_node1")
    CW_OP_node2_handle = self.core_dma.get_handle("CW_optical_pumping_node2")
    # recooling_dma_handle = self.core_dma.get_handle("recooling")
    excitation_dma_handle = self.core_dma.get_handle("excitation_and_photon_collection")
    parity_mapping_dma_handle = self.core_dma.get_handle("parity_mapping")
    delay(10 * ms)

    self.measurement = 0  # advances in end_measurement

    AllSPCMs_parity_RO = [-1] * self.n_measurements
    SPCM0_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_SinglePhoton = [-1.0] * self.n_measurements
    SPCM0_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    SPCM1_OtherNode_SinglePhoton = [-1.0] * self.n_measurements
    angle_780_HWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0
    angle_780_QWP = [-1.0] * self.n_measurements  ## in node2, this needs to be -1.0

    self.core.break_realtime()

    while self.measurement < self.n_measurements:

        self.core.break_realtime()

        self.ttl_exc0_switch.on()  # turns off the excitation
        delay(1 * ms)

        if self.advance_if_atom_in_first_shot:
            ### Keep a measurement only once the FIRST SHOT really saw an atom.
            ### Safe to retry because first_shot writes NO datasets -- every
            ### append happens in end_measurement -- so a rejected attempt
            ### leaves no trace in any histogram or h5.
            atom_confirmed = False
            confirm_tries = 0
            while (not atom_confirmed
                   and confirm_tries < self.max_loading_rounds):
                confirm_tries += 1

                load_until_atom_smooth_FORT_recycle(self)
                # load_MOT_and_FORT_until_atom_recycle(self)

                delay(1 * ms)
                self.ttl_microwave_switch.on()  ### close the switch
                delay(20*us)
                self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

                first_shot(self)

                ### counts/s on BOTH sides: AllSPCMs_RO1 is raw counts over
                ### t_SPCM_first_shot and the threshold is a rate.
                if (self.AllSPCMs_RO1 / self.t_SPCM_first_shot
                        > self.single_atom_threshold):
                    atom_confirmed = True
                else:
                    ### The loader's already-loaded pre-check reads AllSPCMs_RO2.
                    ### The shot above just proved the trap empty, so a stale value
                    ### would make the loader return without loading and this loop
                    ### would spin. Only on a retry: clearing it before attempt 1
                    ### would defeat the recycle path and force a full MOT load
                    ### every measurement.
                    self.AllSPCMs_RO2 = 0
                    ### the rejected attempt's host RPCs advanced the wall clock
                    ### while the timeline cursor stood still
                    self.core.break_realtime()

            ### Loading probability is 1 / mean(this series): the confirm loop
            ### makes n_atom_loaded_per_iteration count loads across attempts
            ### rather than measurements, so this is where that goes.
            self.append_to_dataset('n_first_shot_confirm_tries', confirm_tries)
            self.core.break_realtime()

        else:
            load_until_atom_smooth_FORT_recycle(self)
            # load_MOT_and_FORT_until_atom_recycle(self)

            delay(1 * ms)
            self.ttl_microwave_switch.on()  ### close the switch
            delay(20*us)
            self.dds_microwaves.sw.on()  ### turns on the DDS not the switch

            first_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        delay(10 * us)
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()
        delay(10 * us)

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?

        if self.which_node == "alice":
            # self.ttl_GRIN2_switch.on()  # turns off excitation ## make this global?
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            delay(5 * us)
            self.GRIN1and2_dds.sw.on()  # GRIN1 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON
            delay(5 * us)

        excitation_cycle = 0  ### just for initialization.

        if self.AllSPCMs_RO1 / self.t_SPCM_first_shot > self.single_atom_threshold:
            delay(1 * ms)

            ### with cw pumping:
            ### TODO: Make sure this works fine in this experiment. I have not run this experiment with CW OP yet.
            ### Akbar 2026-03-26
            if self.t_pumping > 0.0:
                # delay(10 * us)
                if self.which_node == "alice":
                    # CW_optical_pumping_node1(self)
                    self.core_dma.playback_handle(CW_OP_node1_handle)
                else:
                    # CW_optical_pumping_node2(self)
                    self.core_dma.playback_handle(CW_OP_node2_handle)
            # t_after_OP_dma = now_mu()

            # if self.which_node == "alice":
            #     ### Changing the bias field.
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            #     delay(1 * ms)

            #
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            # delay(10 * us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            for excitation_attempt in range(self.n_excitation_attempts):
                # self.ttl_exc0_switch.off()  # turns on the excitation0 AOM

                # slack = now_mu() - self.core.get_rtio_counter_mu()
                # if slack < 1e5:
                #     # self.print_async("slack added in measurement:", self.measurement)
                #     self.core.break_realtime()

                self.dds_microwaves.set(frequency=self.f_microwaves_11_dds, amplitude=dB_to_V(self.p_microwaves))
                delay(5 * us)

                dma_start_mu = now_mu()
                self.core_dma.playback_handle(excitation_dma_handle)

                t1 = dma_start_mu + self.core.seconds_to_mu(5 * us)
                t_end_SPCM = t1 + int(self.gate_start_offset_mu) + self.core.seconds_to_mu(self.t_photon_collection_time)

                SPCM0_click_time = self.ttl_SPCM0.timestamp_mu(t_end_SPCM)
                SPCM1_click_time = self.ttl_SPCM1.timestamp_mu(t_end_SPCM)
                SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM)
                SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM)

                detector = -1
                mapping_click_time = 0

                if SPCM0_click_time > 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 0, SPCM0_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time > 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 1, SPCM1_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time > 0 and SPCM1_OtherNode_click_time < 0:
                    detector, mapping_click_time = 2, SPCM0_OtherNode_click_time
                elif SPCM0_click_time < 0 and SPCM1_click_time < 0 and SPCM0_OtherNode_click_time < 0 and SPCM1_OtherNode_click_time > 0:
                    detector, mapping_click_time = 3, SPCM1_OtherNode_click_time

                ###todo: somehow there are some early counts detected. if this can be sorted out, can resolve underflow caused by this?
                ### probably no! early counts happen in single photon experiments but not here. because after this sequence,
                ### the left over counts will be collected in atom_parity_shot???

                if detector >= 0:
                    at_mu(mapping_click_time + self.t_start_MW_mapping_mu)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
                    ##todo: I am setting FORT drop in dma because FORT power barely drifts.

                    # delay(2 * us)
                    self.core_dma.playback_handle(parity_mapping_dma_handle)

                    # delay(5 * us)
                    # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
                    ##todo: I am setting FORT drop in dma because FORT power barely drifts.
                    # t_after_parity_mapping = now_mu()

                    self.core.break_realtime()

                    delay(5 * us)
                    self.core_dma.playback_handle(ba_dma_handle)

                    delay(5 * us)
                    atom_parity_shot(self)

                    #todo: test if this reduces the mapping time & underflow errors
                    ### Drain additional timestamp events from the previous photon gate
                    if SPCM0_click_time > 0: self.ttl_SPCM0.count(t_end_SPCM)
                    if SPCM1_click_time > 0: self.ttl_SPCM1.count(t_end_SPCM)
                    if SPCM0_OtherNode_click_time > 0: self.ttl_SPCM0_OtherNode.count(t_end_SPCM)
                    if SPCM1_OtherNode_click_time > 0: self.ttl_SPCM1_OtherNode.count(t_end_SPCM)

                    self.core.break_realtime()
                    delay(1 * ms)
                    AllSPCMs_parity_RO[self.measurement] = self.AllSPCMs_parity_RO
                    SPCM0_SinglePhoton[self.measurement] = 1.0 if detector == 0 else 0.0
                    SPCM1_SinglePhoton[self.measurement] = 1.0 if detector == 1 else 0.0
                    SPCM0_OtherNode_SinglePhoton[self.measurement] = 1.0 if detector == 2 else 0.0
                    SPCM1_OtherNode_SinglePhoton[self.measurement] = 1.0 if detector == 3 else 0.0
                    angle_780_HWP[self.measurement] = self.target_780_HWP
                    angle_780_QWP[self.measurement] = self.target_780_QWP

                    delay(1 * ms)

                    # t_delay_checker = (t_after_parity_mapping - t_after_OP_dma)
                    # print("t_delay_checker: ", t_delay_checker)

                    self.measurement += 1
                    break

            delay(20 * us)
            excitation_cycle += 1


        delay(5*us)
        ### Restore feedback amplitudes while RF switches are off
        self.dds_AOM_A5.set(frequency=self.AOM_A5_freq, amplitude=self.stabilizer_AOM_A5.amplitude)
        self.dds_AOM_A6.set(frequency=self.AOM_A6_freq, amplitude=self.stabilizer_AOM_A6.amplitude)

        delay(1 * ms)
        second_shot(self)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

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
        self.append_to_dataset('AllSPCMs_parity_RO', AllSPCMs_parity_RO[i])
        # self.append_to_dataset('SPCM0_SinglePhoton_parity', SPCM0_SinglePhoton[i])
        # self.append_to_dataset('SPCM1_SinglePhoton_parity', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_SinglePhoton', SPCM0_SinglePhoton[i])
        self.append_to_dataset('SPCM1_SinglePhoton', SPCM1_SinglePhoton[i])
        self.append_to_dataset('SPCM0_OtherNode_SinglePhoton', SPCM0_OtherNode_SinglePhoton[i])
        self.append_to_dataset('SPCM1_OtherNode_SinglePhoton', SPCM1_OtherNode_SinglePhoton[i])
        self.append_to_dataset('angle_780_HWP', angle_780_HWP[i])
        self.append_to_dataset('angle_780_QWP', angle_780_QWP[i])
    delay(5 * ms)

@kernel
def trap_frequency_experiment(self):
    """
    For spectroscopy of the trap vibrational frequencies with the Rigol D11022Z function generator.

    One way to use this is with GeneralVariableScan and scan_variable1 = f_Rigol_modulation,
    where f_Rigol_modulation must be within the boundaries Rigol_carrier_frequency +/- Rigol_FM_deviation
    specified in ExpermientVariables (which in turn must be set on the Rigol D11022Z).

    :param self: an experiment instance.
    :return:
    """

    set_RigolDG1022Z(frequency=self.Rigol_carrier_frequency,
                     vpp=self.Rigol_V_pp,
                     vdc=self.Rigol_V_DC)

    self.core.reset()

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0

    self.require_D1_lock_to_advance = False # override experiment variable

    # self.set_dataset(self.SPCM0_rate_dataset,
    #                  [0.0],
    #                  broadcast=True)

    self.measurement = 0
    while self.measurement < self.n_measurements:

        #TODO: just set the rigol frequency using pyvisa

        if self.enable_laser_feedback:
            run_feedback_and_record_FORT_MM_power(self)

        load_MOT_and_FORT(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        # set the FORT AOM to the science setting. this is only valid if we have run
        # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT
        self.dds_FORT.set(frequency=self.f_FORT,
                                amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        if not self.no_first_shot:
            first_shot(self)

            if self.t_recooling_after_first_shot > 0:
                recooling_after_first_shot(self)

        ##############################################################################
        # modulate the FORT
        ##############################################################################

        self.FORT_mod_switch.on()  # toggle the modulation to the VCA
        delay(10 * ms)
        self.FORT_mod_switch.off()
        delay(1 * ms)

        second_shot(self)

        end_measurement(self)

    self.dds_FORT.sw.off()

# ###############################################################################
# # 3. DIAGNOSTIC FUNCTIONS
# ###############################################################################

@kernel
def FORT_monitoring_with_Luca_experiment(self):
    """
    A modified version of atom_loading_experiment for monitoring the FORT scattering.
    This is most easily run with tests/MonitorFORTPowerWithLuca.py

    Load a MOT, load a single atom, readout, wait self.t_delay_between_shots, readout again.

    For this experiment, we were still using the MM fiber monitor (i.e. after the polarizer) to
    feed back to the FORT power. We used this experiment to characterize how well the signal
    from an APD monitoring the scattered 852 nm light in the chamber correlated with the signal
    seen by the Luca in an ROI restricted to a section of the parabolic mirror tube.

    The assumed connection of the Sampler card and channel used to monitor the APD is given
    in load_MOT_and_FORT_for_Luca_scattering_measurement.

    ----
    Alternatively, if you are using the APD to feedback to the FORT:
    The MM fiber is monitored with the channel specified in load_MOT_and_FORT_for_Luca_scattering_measurement.

    The purpose of this experiment is to be able to compare the normalized FORT power recorded
    by the dds_FORT feedback channel, the scattering seen by the camera, and the voltage of
    the detector monitoring the MM fiber after the polarizer. Rotating the 852 nm
    motorized waveplates, e.g. with the APT program, we can check whether the scattering in the
    chamber seen with the Luca correlates with polarization. This is a way to vet our APD
    feedback scheme.

    ----
    For analysis, see FORT feedback/monitor_FORT_scattering_and_Raman_light.ipynb

    :param self: an experiment instance.
    :return:
    """

    self.core.reset()

    self.set_dataset("SPCM0_FORT_loading", [0.0], broadcast=True)
    self.set_dataset("SPCM0_FORT_and_MOT", [0.0], broadcast=True)
    self.set_dataset("SPCM0_FORT_science", [0.0], broadcast=True)
    self.set_dataset("APD_FORT_volts_loading", [0.0], broadcast=True)
    self.set_dataset("APD_FORT_volts_science", [0.0], broadcast=True)

    self.SPCM0_RO1 = 0
    self.SPCM0_RO2 = 0

    self.require_D1_lock_to_advance = False # override experiment variable
    self.require_atom_loading_to_advance = False # override experiment variable

    # self.set_dataset(self.SPCM0_rate_dataset,
    #                  [0.0],
    #                  broadcast=True)

    delay(100*ms)
    now = now_mu()
    for i in range(self.warm_up_shots):
        at_mu(now+i*self.core.seconds_to_mu(500*ms))
        self.ttl_Luca_trigger.pulse(5 * ms)
        # now = now_mu()

    self.core.wait_until_mu(now_mu())

    self.measurement = 0
    while self.measurement < self.n_measurements:

        if self.enable_laser_feedback:
            self.laser_stabilizer.run(monitor_only=self.no_feedback)  # this tunes the MOT and FORT AOMs

        load_MOT_and_FORT_for_Luca_scattering_measurement(self)

        delay(0.1*ms)
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
            channels=self.coil_channels)

        # set the FORT AOM to the science setting. this is only valid if we have run
        # feedback to reach the corresponding setpoint first, which in this case, happened in load_MOT_and_FORT
        if not self.no_feedback:
            self.dds_FORT.set(frequency=self.f_FORT,
                                    amplitude=self.stabilizer_FORT.amplitudes[1])

        # set the cooling DP AOM to the readout settings
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO,
                                amplitude=self.ampl_cooling_DP_MOT*self.p_cooling_DP_RO)

        # take the first shot
        if not self.MOT_light_off:
            self.dds_cooling_DP.sw.on()
        with parallel:
            self.ttl_Luca_trigger.pulse(5 * ms)
            t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_first_shot)
        self.SPCM0_RO1 = self.ttl_SPCM0.count(t_gate_end)
        delay(1 * ms)
        self.dds_cooling_DP.sw.off()

        delay(self.t_delay_between_shots)

        # take the second shot
        if not self.MOT_light_off:
            self.dds_cooling_DP.sw.on()
        # with parallel:
            # self.ttl_Luca_trigger.pulse(5 * ms)
        t_gate_end = self.ttl_SPCM0.gate_rising(self.t_SPCM_second_shot)
        self.SPCM0_RO2 = self.ttl_SPCM0.count(t_gate_end)

        delay(1 * ms)

        end_measurement(self)

        # update experiment-specific datasets:
        self.append_to_dataset("SPCM0_FORT_loading", self.SPCM0_FORT_loading)
        self.append_to_dataset("SPCM0_FORT_and_MOT", self.SPCM0_FORT_and_MOT)
        self.append_to_dataset("SPCM0_FORT_science", self.SPCM0_FORT_science)
        self.append_to_dataset("APD_FORT_volts_loading", self.APD_FORT_volts_loading)
        self.append_to_dataset("APD_FORT_volts_science", self.APD_FORT_volts_science)

    self.dds_FORT.sw.off()

@kernel
def Pulse_microwave_smooth(self, MW_freq):
    """
    This uses dds RAM profiles prepared in BaseExperiment to generate smooth MW pulses. Works but blocks
    other dds channels. Needs fix.
    Use it in the experiment like:
    Pulse_microwave_smooth(self, self.f_microwaves_00_dds)

    """
    ### predefine these 3 parameters. Needed but not used.
    MW_step_size = 1
    MW_total_points = 1
    MW_amplitudes_list = [0]
    MW_pulse_length = 0.0

    if MW_freq == self.f_microwaves_00_dds:
        MW_step_size = self.MW_00_step_size
        MW_total_points = self.MW_00_total_points
        MW_amplitudes_list = self.MW_00_amplitudes_list
        MW_pulse_length = self.t_microwave_00_pulse

    if MW_freq == self.f_microwaves_01_dds:
        MW_step_size = self.MW_01_step_size
        MW_total_points = self.MW_01_total_points
        MW_amplitudes_list = self.MW_01_amplitudes_list
        MW_pulse_length = self.t_microwave_01_pulse

    if MW_freq == self.f_microwaves_11_dds:
        MW_step_size = self.MW_11_step_size
        MW_total_points = self.MW_11_total_points
        MW_amplitudes_list = self.MW_11_amplitudes_list
        MW_pulse_length = self.t_microwave_11_pulse

    self.dds_microwaves.set_frequency(MW_freq) ### takes 0.7us
    self.dds_microwaves.set_att(0.0) ### takes 1.6us


    if MW_pulse_length > 0.0:

        self.dds_microwaves.set_cfr1(ram_enable=0)  ### disable RAM mode to write the config
        self.dds_microwaves.cpld.io_update.pulse_mu(8)  ### pulse the ttl to update and implement settings
        ### the two lines above take 0.7us

        ### Configures the RAM playback engine. Takes 1.3us
        self.dds_microwaves.set_profile_ram(
            start=0,
            end=MW_total_points - 1,
            step=MW_step_size,
            profile=0,
            mode=RAM_MODE_RAMPUP,
        )

        self.dds_microwaves.cpld.set_profile(0) ### takes 0.4us
        # self.dds.cpld.io_update.pulse_mu(8)
        self.ttl_microwave_switch.off()  ### takes 0us

        self.dds_microwaves.write_ram(MW_amplitudes_list)  ### write the data onto RAM.
        ### Takes 33us with 60 MW_amp_points. Takes 48us with 90 points. Takes 96us with 180 points.

        ### Enabling RAM playback, not playing yet. takes 0.7us
        self.dds_microwaves.set_cfr1(
            internal_profile=0,
            ram_enable=1,
            ram_destination=RAM_DEST_ASF,
        )

        # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope

        self.dds_microwaves.sw.on()
        self.dds_microwaves.cpld.io_update.pulse_mu(8)  ### This runs the RAM

        delay(MW_pulse_length)  ### keep the delay as ramp time

        ### shutting off. Takes 0.7us
        self.dds_microwaves.set_cfr1(ram_enable=0)
        self.dds_microwaves.cpld.io_update.pulse_mu(8)
        self.dds_microwaves.sw.off()

        self.core.reset()
        self.urukul2_cpld.init()
        self.urukul2_ch0.cpld.init()
        self.urukul2_ch1.cpld.init()
        self.urukul2_ch2.cpld.init()
        self.urukul2_ch3.cpld.init()

        self.ttl_microwave_switch.on()
        # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

@kernel
def track_1_microwave_transition_experiment(self):
    """
    Starting with testing the dataset read.
    1- Read the dataset and set it to f01. Use getattr to avoid "is already define message"
    2- Use that as the starting point to find the new res. f.
    3- update the dataset. We never want to run the subroutine before exiting the kernel. GVS exits the kernel between each scan.
    So, the dataset gets updated with the new res freq which can be used for the next run of the subroutine.

    """
    self.core.reset()
    delay(1*ms)



    # try:
    #     f01 = self.f01_resonance
    # except ValueError:
    #     f01 = self.get_dataset("f_microwaves_01_dds")

    # self.print_async(f01)
    # delay(10*ms)


    # self.default_setpoints = np.array([getattr(self, dataset) for dataset in self.setpoint_datasets])

    # f = getattr(self, self.f_microwaves_01_dds)

    # self.set_dataset("f_microwaves_01_dds_test", 0.0, broadcast=True)
    #
    # self.f_microwaves_01_dds = 333.0 * MHz
    # self.append_to_dataset("f_microwaves_01_dds", 0.0)


# ###############################################################################
# # 4. TWO NODES EXPERIMENT FUNCTIONS WITH TTL SYNCHRONIZATION
# # These are functions that are used for two nodes experiment with TTL
# # symnchronization. Working as of 20260810. Soon will be deprecated after
# # upgrading to master-satellite configuration. - Eunji
# ###############################################################################

@kernel
def two_nodes_synchronization(self):
    """
    1. Initial coarse synchronization(Two-way communication)
    - ttl_node1_output1 and ttl_node2_output1 turned ON
    - IF both nodes receives this signal, move on to finer synchronization
    2. Fine synchronization(One-way communication)
    - ttl_node1_output2 is pulsed @ t_node1_ref_0
    - node2 receives this siganl and time tags @ t_node2_ref_0
    3. Sets each node's time cursor to 150us after this.
    """
    # ############################################################
    # # Starting synchronization
    # # Signal to the other node.
    # ############################################################
    other_node_ready = False
    readout = 0
    # sync_starting_at = now_mu()

    # self.core.break_realtime()

    delay(10 * us)

    if self.which_node == "alice":
        self.ttl_node1_output1.on()
    else:
        self.ttl_node2_output1.on()

    delay(1 * us)

    ############################################################
    # 1. Initial coarse synchronization
    # Wait until the other node is also ready.
    # Check every wait_time
    ############################################################
    while not other_node_ready:
        wait_time = 4*us
        delay(wait_time)
        if self.which_node == "alice":
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node2_input1.sample_get())
            ## this can be replaced with sample_get_nonrt; but it has core.break_realtime() so needs to be careful

            if readout == 1:
                other_node_ready = True
            else:
                other_node_ready = False

        else:
            self.ttl_node1_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node1_input1.sample_get())

            if readout == 1:
                other_node_ready = True
            else:
                other_node_ready = False

    delay(10 * us)

    ############################################################
    # 2. Fine synchronization
    # Synchronizing with time tag
    ############################################################
    #todo: check if these timings can be reduced or needs to be increased
    #      now it takes 221us for the entire syncrhonization to be done;
    ##### Syncing by sending a TTL pulse from node1 to node2
    t_node1_ref_0 = -1
    t_node2_ref_0 = -1

    if self.which_node == "alice":
        delay(50*us)
        # delay(15 * us)
        t_node1_ref_0 = now_mu()
        self.ttl_node1_output2.pulse(1*us)
    else:
        t_gate_end = self.ttl_node1_input2.gate_rising(150*us)
        # t_gate_end = self.ttl_node1_input2.gate_rising(30 * us)
        t_node2_ref_0 = self.ttl_node1_input2.timestamp_mu(t_gate_end)

        if t_node2_ref_0 == -1:
            self.print_async("Node2 did not detect pulse from Node1")

    delay(10 * us)
    t_node1_ref_1 = t_node1_ref_0 + self.core.seconds_to_mu(200 * us) + 200
    t_node2_ref_1 = t_node2_ref_0 + self.core.seconds_to_mu(200 * us)

    if self.which_node == "alice":
        at_mu(t_node1_ref_1)
        self.ttl_node1_output1.off()
        # sync_time_took = int((t_node1_ref_1 - sync_starting_at) / 1e3)
    else:
        at_mu(t_node2_ref_1)
        self.ttl_node2_output1.off()
        # sync_time_took = int((t_node2_ref_1 - sync_starting_at) / 1e3)

    delay(10 * us)

    # sync_time_took = int((t_node1_ref_1 - sync_starting_at) / 1e3)
    # self.print_async("sync_time_took :", sync_time_took, " us")
    # self.append_to_dataset("sync_time_took", sync_time_took)

@kernel
def two_nodes_synchronization2(self) -> TBool:
    t_node1_ref_0 = -1
    t_node2_ref_0 = -1

    poll_period = 1 * ms
    timeout_cycles = 10000     # 10 s if poll_period = 1 ms and timeout_cycles = 10000
    ack_delay_mu = self.core.seconds_to_mu(50 * us)
    common_delay_mu = self.core.seconds_to_mu(2 * ms)

    alice_extra_offset_mu = 200  ### Calibrated this on the scope.

    self.core.break_realtime()
    delay(100 * us)

    self.print_async("******* Syncing the nodes **********")

    if self.which_node == "bob":

        #### Open Bob's input gate first.
        t_gate_start = now_mu() + self.core.seconds_to_mu(100 * us)
        at_mu(t_gate_start)
        t_gate_end = self.ttl_node1_input1.gate_rising(10.01 * s)

        ### Now assert READY while the gate is open.
        ### This is scheduled after gate opening but before calling timestamp_mu().
        at_mu(t_gate_start + self.core.seconds_to_mu(10 * us))
        self.ttl_node2_output1.on()

        ### Now wait for Alice's reference pulse.
        t_node2_ref_0 = self.ttl_node1_input1.timestamp_mu(t_gate_end)

        if t_node2_ref_0 < 0:
            self.core.break_realtime()
            self.ttl_node2_output1.off()
            return False

        ### Drop READY after successful sync.
        self.core.break_realtime()
        self.ttl_node2_output1.off()

        ### Put Bob at the common future reference.
        at_mu(t_node2_ref_0 + common_delay_mu)
        return True

    else:
        ### Alice polls Bob's READY level.
        ready = False
        for i in range(timeout_cycles):
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            if int(self.ttl_node2_input1.sample_get()):
                ready = True
                break
            delay(poll_period)

        if not ready:
            self.core.break_realtime()
            return False

        ### Alice sends the reference pulse at a deterministic future time.
        t_node1_ref_0 = now_mu() + ack_delay_mu
        at_mu(t_node1_ref_0)
        self.ttl_node1_output1.pulse(500 * ns)

        ### Put Alice at the corresponding future reference.
        at_mu(t_node1_ref_0 + common_delay_mu + alice_extra_offset_mu)
        return True

@kernel
def two_nodes_synchronization_trial(self) -> TBool:
    """
    Faisal's code based on two_nodes_synchronization2
    """
    delay(1*ms)

    t_node1_ref_0 = np.int64(-1)
    t_node2_ref_0 = np.int64(-1)

    poll_period = 1 * ms
    timeout_cycles = 10000     # 10 s if poll_period = 1 ms and timeout_cycles = 10000
    ack_delay_mu = self.core.seconds_to_mu(50 * us)
    common_delay_mu = self.core.seconds_to_mu(2 * ms)

    alice_extra_offset_mu = np.int64(200)  ### Calibrated this on the scope.

    self.core.break_realtime()

    delay(100 * us)
    self.print_async("******* Syncing the nodes **********")

    if self.which_node == "bob":

        #### Open Bob's input gate first.
        t_gate_start = now_mu() + self.core.seconds_to_mu(100 * us)
        gate_duration_mu = self.core.seconds_to_mu(10*s)
        t_gate_end = t_gate_start + gate_duration_mu


        at_mu(t_gate_start)
        with parallel:
            self.ttl_node1_input1.gate_rising_mu(gate_duration_mu)
            with sequential:
                delay(10*us)
                self.ttl_node2_output1.on()


        ### Now wait for Alice's reference pulse.
        t_node2_ref_0 = self.ttl_node1_input1.timestamp_mu(t_gate_end)

        if t_node2_ref_0 < 0:
            self.core.break_realtime()
            self.ttl_node2_output1.off()
            return False

        ### Drop READY after successful sync.
        self.core.break_realtime()
        self.ttl_node2_output1.off()

        ### Put Bob at the common future reference.
        at_mu(t_node2_ref_0 + common_delay_mu)
        return True

    else:

        self.ttl_node2_input1.count(now_mu())

        ### Alice polls Bob's READY level.
        ready = False
        for i in range(timeout_cycles):
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            if int(self.ttl_node2_input1.sample_get()):
                ready = True
                break
            delay(poll_period)

        if not ready:
            self.core.break_realtime()
            return False

        ### Alice sends the reference pulse at a deterministic future time.
        t_node1_ref_0 = now_mu() + ack_delay_mu
        at_mu(t_node1_ref_0)
        self.ttl_node1_output1.pulse(500 * ns)

        ### Put Alice at the corresponding future reference.
        at_mu(t_node1_ref_0 + common_delay_mu + alice_extra_offset_mu)
        return True

@kernel
def load_atom_in_both_nodes_individually(self):
    """
    Eunji - this is assuming that both nodes can try atom loading/readout simultaneously without distruption
    - this maybe useful when we integrate a fast shutter most likeley a MEMS switch.
    - but we need to think more carefully how much loss it has with coupling/mode matching/pol dependence
    - if we need to use MEMS switch, it might need to be spliced.

    Load/recycle atoms until both nodes have atoms simultaneously.

    Assumption:
        This node's TTL output is a live atom-status indicator.

            HIGH = this node currently has an atom
            LOW  = this node currently does not have an atom

        The other node's TTL input is used to determine whether the other node
        currently has an atom.

    This function exits only when:
        this node has atom AND the other node has atom.

    Datasets used:
        time_without_atom_other_node:
            Live-updated time since the other node was last seen without atom.

        Atom_loading_time_other_node:
            Appended once when the other node becomes ready.
    """

    atom_loaded = False
    other_node_has_atom = False
    both_nodes_have_atom = False
    readout = 0

    ############################################################
    # Timer for how long the OTHER node has been without atom.
    ############################################################
    t_before_other_atom = now_mu()
    t_after_other_atom = now_mu()
    t_no_other_atom = now_mu()

    time_without_other_atom = 0.0
    self.atom_loading_time_other_node = 0.0
    other_node_loading_time_recorded = False
    other_node_was_without_atom = False

    delay(100 * us)

    ############################################################
    # Initial check for THIS node using latest RO2 information.
    ############################################################
    if self.measurement > 0:
        if self.AllSPCMs_RO2 / self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True
            self.print_async("This node starts with atom.")
        else:
            atom_loaded = False
            self.print_async("This node starts without atom.")
    else:
        atom_loaded = False
        self.print_async("No previous measurement. This node needs loading.")

    ############################################################
    # Update this node TTL according to the actual initial status.
    ############################################################
    if atom_loaded:
        if self.which_node == "alice":
            self.ttl_node1_output1.on()
            self.print_async("Node1 atom TTL is ON.")
        else:
            self.ttl_node2_output1.on()
            self.print_async("Node2 atom TTL is ON.")
    else:
        if self.which_node == "alice":
            self.ttl_node1_output1.off()
            self.print_async("Node1 atom TTL is OFF.")
        else:
            self.ttl_node2_output1.off()
            self.print_async("Node2 atom TTL is OFF.")

    delay(100 * us)

    ############################################################
    # Main loop.
    # Only exits when this node has atom AND other node has atom.
    ############################################################
    while not (atom_loaded and both_nodes_have_atom):

        ########################################################
        # If this node does not have atom, load/reload it.
        ########################################################
        if not atom_loaded:
            self.print_async("Loading atom in this node...")

            load_MOT_and_FORT_until_atom_recycle(self)

            atom_loaded = True
            self.print_async("Atom loaded in this node.")

            if self.which_node == "alice":
                self.ttl_node1_output1.on()
                self.print_async("Node1 atom TTL is ON after loading.")
            else:
                self.ttl_node2_output1.on()
                self.print_async("Node2 atom TTL is ON after loading.")

            delay(100 * us)

        ########################################################
        # Check whether the OTHER node has atom.
        ########################################################
        if self.which_node == "alice":
            self.ttl_node2_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node2_input1.sample_get())

            if readout == 1:
                other_node_has_atom = True
                self.print_async("Node2 has atom.")
            else:
                other_node_has_atom = False
                self.print_async("Node2 does NOT have atom.")

        else:
            self.ttl_node1_input1.sample_input()
            delay(1 * us)
            readout = int(self.ttl_node1_input1.sample_get())

            if readout == 1:
                other_node_has_atom = True
                self.print_async("Node1 has atom.")
            else:
                other_node_has_atom = False
                self.print_async("Node1 does NOT have atom.")

        ########################################################
        # Track how long the OTHER node has been without atom.
        ########################################################
        if other_node_has_atom:
            self.set_dataset("time_without_atom_other_node", 0.0, broadcast=True)

            if other_node_was_without_atom and not other_node_loading_time_recorded:
                t_after_other_atom = now_mu()

                self.atom_loading_time_other_node = self.core.mu_to_seconds(
                    t_after_other_atom - t_before_other_atom
                )

                self.append_to_dataset(
                    "Atom_loading_time_other_node",
                    self.atom_loading_time_other_node
                )

                other_node_loading_time_recorded = True
                other_node_was_without_atom = False

                self.print_async("Other node loading time recorded.")

        else:
            ####################################################
            # If this is the first time seeing the other node
            # without atom, start the timer.
            ####################################################
            if not other_node_was_without_atom:
                t_before_other_atom = now_mu()
                other_node_was_without_atom = True
                other_node_loading_time_recorded = False

            delay(0.1 * ms)

            t_no_other_atom = now_mu()

            time_without_other_atom = self.core.mu_to_seconds(
                t_no_other_atom - t_before_other_atom
            )

            self.set_dataset(
                "time_without_atom_other_node",
                time_without_other_atom,
                broadcast=True
            )

        ########################################################
        # If both nodes have atoms, exit the loop.
        ########################################################
        if atom_loaded and other_node_has_atom:
            both_nodes_have_atom = True
            self.print_async("Both nodes have atoms. Exiting loading loop.")
            break
        else:
            both_nodes_have_atom = False

        ########################################################
        # Other node does not have atom yet.
        # Wait 200 ms, then check whether THIS node still has atom.
        ########################################################
        delay(200 * ms)

        ########################################################
        # Read out this node's current atom using second_shot.
        ########################################################
        second_shot(self)

        ########################################################
        # Check if this node still has atom.
        ########################################################
        if self.AllSPCMs_RO2 / self.t_SPCM_second_shot > self.single_atom_threshold:
            atom_loaded = True
            self.print_async("This node still has atom while waiting.")

            ####################################################
            # Keep this node TTL HIGH.
            ####################################################
            if self.which_node == "alice":
                self.ttl_node1_output1.on()
            else:
                self.ttl_node2_output1.on()

            ####################################################
            # Recool after the readout.
            ####################################################
            recooling_after_first_shot(self)

        else:
            atom_loaded = False
            both_nodes_have_atom = False
            self.print_async("This node lost atom while waiting.")

            ####################################################
            # This node no longer has atom, so set TTL LOW.
            ####################################################
            if self.which_node == "alice":
                self.ttl_node1_output1.off()
                self.print_async("Node1 atom TTL is OFF.")
            else:
                self.ttl_node2_output1.off()
                self.print_async("Node2 atom TTL is OFF.")

            ####################################################
            # Important:
            # Do NOT reset time_without_atom_other_node here.
            # That timer belongs to the other node's atom status,
            # not this node's atom status.
            ####################################################

            delay(100 * us)

    ############################################################
    # Function exits only here, when both nodes have atoms.
    ############################################################
    self.print_async("Confirmed: both nodes have atoms.")
    delay(100 * us)

@kernel
def load_until_atom_in_both_nodes_recycle(self):
    """
    Two-node combined-loading version.

    * based on load_MOT_and_FORT_until_atom_recycle
    * with synchronization added in the while loop
    * Checks for atoms in both nodes using two_atom_threshold_in_loading

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
        delay(1 * ms)

        max_tries = 100  ### Maximum number of attempts before running the feedback
        atom_check_time = self.t_atom_check_time
        try_n = 0
        t_before_atom = now_mu()  ### is used to calculate the loading time of atoms by atom_loading_time = t_after_atom - t_before_atom
        t_after_atom = now_mu()
        time_without_atom = 0.0

        AllSPCMs_atom_check_loaded = 0  ### for initilization
        AllSPCMs_atom_check_not_loaded = 0

        shim_tune_runs = 0

        n_synchronization_in_loading = 0
        while True:
            two_nodes_synchronization(self)  ### Syncing the two nodes.

            n_synchronization_in_loading += 1

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

                # self.print_async("n_synchronization_in_loading: ", n_synchronization_in_loading)
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

        delay(1 * ms)
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
def two_node_alternating_shot(self):
    """
    Alternating two-node readout. -Eunji as of 2026.07.14

    Alice windows: even windows
    Bob windows:   odd windows

    Each node gets:
        n_readout_windows_per_node windows
        1 ms counter gate per window

    Timing:
        Alice window starts at 0 ms
        Bob window starts 2.3 ms later
        Alice next window starts 4.6 ms later...

    ** Note that unlike other readout functions, I only save AllSPCM data rather than
        saving each dataset separately for each SPCM. Note that at some point if you
        want to implement this, alice and bob should have separate dataset.
        for example:
            SPCM0_alternating_RO_alice, SPCM0_alternating_RO_bob,
            SPCM1_alternating_RO_alice, SPCM1_alternating_RO_bob,
            SPCM0_OtherNode_alternating_RO_alice, SPCM0_OtherNode_alternating_RO_bob,
            SPCM1_OtherNode_alternating_RO_alice, SPCM1_OtherNode_alternating_RO_bob,
    """

    self.AllSPCMs_alternating_RO_alice = 0
    self.AllSPCMs_alternating_RO_bob = 0

    ### set the coils to PGC settings
    self.zotino0.set_dac(
        [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
        channels=self.coil_channels)
    delay(1 * ms)  ## coils relaxation time

    # Set cooling DP AOM to readout settings
    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)

    # Set FORT AOM to science settings
    self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])

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

    n_readout_windows_per_node = 10

    t_readout_window = 1.0 * ms

    self.ttl_repump_switch.on()  # turn off MOT RP
    self.dds_cooling_DP.sw.off()

    # 2*n_readout_windows_per_node because even windows are Alice, odd windows are Bob
    for i in range(2 * n_readout_windows_per_node):

        alice_window = (i % 2 == 0)

        # Only the active node turns on readout light.
        if alice_window:
            if self.which_node == "alice":
                self.ttl_repump_switch.off()   # turn on MOT RP
                self.dds_cooling_DP.sw.on()    # turn on cooling/readout=
        else:
            if self.which_node == "bob":
                self.ttl_repump_switch.off()   # turn on MOT RP
                self.dds_cooling_DP.sw.on()    # turn on cooling/readout

        # waiting for RO beams to turn on
        delay(0.1 * ms)

        # Open all counters during this readout window
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(t_readout_window)
            self.ttl_SPCM1_counter.gate_rising(t_readout_window)
            self.ttl_SPCM0_OtherNode_counter.gate_rising(t_readout_window)
            self.ttl_SPCM1_OtherNode_counter.gate_rising(t_readout_window)

        delay(0.1 * ms)
        self.ttl_repump_switch.on()  # turn off MOT RP
        self.dds_cooling_DP.sw.off()

        SPCM0_alternating_RO = self.ttl_SPCM0_counter.fetch_count()
        SPCM1_alternating_RO = self.ttl_SPCM1_counter.fetch_count()
        SPCM0_OtherNode_alternating_RO = self.ttl_SPCM0_OtherNode_counter.fetch_count()
        SPCM1_OtherNode_alternating_RO = self.ttl_SPCM1_OtherNode_counter.fetch_count()

        window_count = (SPCM0_alternating_RO + SPCM1_alternating_RO
                        + SPCM0_OtherNode_alternating_RO + SPCM1_OtherNode_alternating_RO)

        if alice_window:
            self.AllSPCMs_alternating_RO_alice += window_count
        else:
            self.AllSPCMs_alternating_RO_bob += window_count

        delay(0.1*ms)

@kernel
def Two_nodes_atom_loading_experiment(self):
    """
    Simple atom loading experiment in both nodes
    - checking for atoms in both nodes at the same time
    - based on atom_loading_2_experiment

    Sequence as of 2026.07.13
    1. Load atoms in both nodes simultaneously
    2. First shot
    3. FORT drop if > 0
    4. testing individual node shot - Alternating RO
    5. Second shot
    6. end_measurement
    """

    self.core.reset()
    self.require_D1_lock_to_advance = False  # override experiment variable

    self.n_feedback_per_iteration = 2  ### number of times the feedback runs in each iteration. Updates in atom loading subroutines.
    ### Required only for averaging RF powers over iterations in analysis. Starts with 2 because RF is measured at least 2 times
    ### in each iteration.
    self.n_atom_loaded_per_iteration = 0

    if self.enable_laser_feedback:
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        run_feedback_and_record_FORT_MM_power(self)

    ############################################################
    # Turn OFF this node's TTL to be ready for synchronization.
    ############################################################
    if self.which_node == "alice":
        self.ttl_node1_output1.off()
    else:
        self.ttl_node2_output1.off()
    delay(50 * ms)

    self.measurement = 0
    while self.measurement < self.n_measurements:
        delay(10 * ms)

        two_nodes_synchronization(self)

        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        if self.t_FORT_drop > 0:
            self.dds_FORT.sw.off()
            delay(self.t_FORT_drop)
            self.dds_FORT.sw.on()

        delay(self.t_delay_between_shots)

        two_nodes_synchronization(self)

        delay(1*ms)
        two_node_alternating_shot(self)
        delay(1 * ms)

        two_nodes_synchronization(self)

        second_shot(self)

        end_measurement(self)

    self.append_to_dataset('n_feedback_per_iteration', self.n_feedback_per_iteration)
    self.append_to_dataset('n_atom_loaded_per_iteration', self.n_atom_loaded_per_iteration)

@kernel
def Two_node_single_photon_experiment(self):
    """
    **** This is based on single_photon_experiment_3_atom_loading_advance
    **** Note that 5 is the latest version as of 2026-06-19
    **** But version 5 has not yet been tested in node2 so I will be using version 3

    **** This experiment is used to match the timing between the single photon arrival
         time from each node. Here, to do so inidividually, I physically block the excitation
         light of the other node while collecting the photon arrival histogram from one.

    **** It starts by loading atoms in both nodes, and then does all the excitation cycles
         in both nodes. This is to make sure that two nodes are running the same sequence
         and being synchronized to the next step correctly.

    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            # self.core.break_realtime()  ## Commenting this out to avoid timing error

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
            delay(1*ms)
            if self.t_pumping > 0.0:
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10*us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
            # delay(10*us)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(5 * us)

            self.core.break_realtime()

            two_nodes_synchronization(self)

            if self.which_node == "bob":
                delay_mu(int(self.t_delay_in_bob_mu))

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
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                # at_mu(t1 + 30000)
                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

            delay(20 * us)
            self.ttl_exc0_switch.on()  # block Excitation

            ############################ atom cooling phase with PGC settings
            #todo: check Node2 + timing
            self.core.break_realtime()
            two_nodes_synchronization(self)

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

            self.core.break_realtime()
            two_nodes_synchronization(self)

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
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())


                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.two_atom_threshold:
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

            self.core.break_realtime()
            two_nodes_synchronization(self)
            delay(10 * us)

        delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.1 * ms)

        two_nodes_synchronization(self)

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
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def Two_node_single_photon_2_experiment(self):
    """
    **** This is for Hong-Ou_Mandel Experiment
    **** This is based on Two_node_single_photon_experiment

    **** In HOM Measurement, we want to vary the excitation timing of one node to see
         see the HOM dip. However while we do that, we might accidentally keep the
         FORT from the other node on while one node is collecting photons, making
         false count measurements. Thus, here I added some timing changes in
         excitation cycle loop.

    ****
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    self.core.break_realtime()

    two_nodes_synchronization(self)
    delay(10*us)

    #### So that the other node is not turning the FORT on while photon collection in one node;
    t_diff_mu = int(self.t_delay_in_bob_mu - 189)

    if self.which_node == "alice":
        # t_diff = int(self.t_delay_in_bob_mu - 189)
        # t_FORT_OFF_mu = 0

        if t_diff_mu < 0:
            t_gate_mu = t_diff_mu
            t_photon_add_mu = -1 * t_diff_mu
            t_FORT_OFF_mu = t_diff_mu
        else:
            t_gate_mu = 0
            t_photon_add_mu = t_diff_mu
            t_FORT_OFF_mu = 0

        # t_photon = abs(t_diff)
        # self.t_photon_collection_time = self.t_photon_collection_time + int(t_photon_add) * ns

    else:
        # t_gate_mu = 0
        # t_photon = 0
        if t_diff_mu < 0:
            t_gate_mu = 0
            t_photon_add_mu = -1 * t_diff_mu
            t_FORT_OFF_mu = 0
        else:
            t_gate_mu = -1 * t_diff_mu
            t_photon_add_mu = t_diff_mu
            t_FORT_OFF_mu = -1 * t_diff_mu

    # t_diff_mu = 0
    # t_gate_mu = 0
    # t_photon_add_mu = 0
    # t_FORT_OFF_mu = 0

    while self.measurement < self.n_measurements:

        AllSPCMs_RO_atom_check_array = [0] * int(self.max_excitation_cycles/self.atom_check_every_n)
        tStamps_t1 = [0.0]  * (self.max_excitation_cycles * self.n_excitation_attempts)
        SPCM0_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM0_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]
        SPCM1_OtherNode_timestamps = [[-1.0] * max_clicks for _ in range(self.max_excitation_cycles * self.n_excitation_attempts)]

        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        delay(10*us)
        two_nodes_synchronization(self)

        ##todo: add synchronization here also??
        ##      also, add how much time it took for synchronization.
        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        for excitation_cycle in range(self.max_excitation_cycles):
            ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
            delay(10*us)
            if self.t_pumping > 0.0:
                if self.which_node == "alice":
                    CW_optical_pumping_node1(self)
                else:
                    CW_optical_pumping_node2(self)
                delay(10*us)

            ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
            # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
            # delay(10*us)

            self.ttl_exc0_switch.off() # turns on the excitation0 AOM
            delay(5 * us)

            self.core.break_realtime()

            two_nodes_synchronization(self)

            for excitation_attempt in range(self.n_excitation_attempts):

                if self.which_node == "bob":
                    delay_mu(int(self.t_delay_in_bob_mu))

                t0 = now_mu()
                t1 = int(t0 + 150)

                at_mu(t1 - int(t_FORT_OFF_mu))
                self.dds_FORT.sw.off()  ### turns FORT off

                at_mu(t1 + 50 + int(self.t_photon_collection_time / ns) + int(t_photon_add_mu))
                self.dds_FORT.sw.on()  ### turns FORT on

                at_mu(t1 + int(self.t_excitation_offset_mu))
                self.ttl_GRIN2_switch.off()  # turns on excitation

                at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
                self.ttl_GRIN2_switch.on()  # turns off excitation

                ######### time stamping the photons. Counting to be done in analysis.
                SPCM0_click_counter = 0
                SPCM1_click_counter = 0
                SPCM0_OtherNode_click_counter = 0
                SPCM1_OtherNode_click_counter = 0

                at_mu(t1 + int(self.gate_start_offset_mu) + int(t_gate_mu))
                with parallel:
                    t_end_SPCM0 = self.ttl_SPCM0.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM1 = self.ttl_SPCM1.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM0_OtherNode = self.ttl_SPCM0_OtherNode.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)
                    t_end_SPCM1_OtherNode = self.ttl_SPCM1_OtherNode.gate_rising(self.t_photon_collection_time + int(t_photon_add_mu) * ns)

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

                ### timestamping SPCM0_OtherNode events
                while SPCM0_OtherNode_click_counter < max_clicks:
                    SPCM0_OtherNode_click_time = self.ttl_SPCM0_OtherNode.timestamp_mu(t_end_SPCM0_OtherNode)
                    if SPCM0_OtherNode_click_time == -1.0:
                        break
                    SPCM0_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM0_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM0_OtherNode_click_time)
                    SPCM0_OtherNode_click_counter += 1

                ### timestamping SPCM1_OtherNode events
                while SPCM1_OtherNode_click_counter < max_clicks:
                    SPCM1_OtherNode_click_time = self.ttl_SPCM1_OtherNode.timestamp_mu(t_end_SPCM1_OtherNode)
                    if SPCM1_OtherNode_click_time == -1.0:
                        break
                    SPCM1_OtherNode_timestamps[excitation_cycle * self.n_excitation_attempts + excitation_attempt][
                        SPCM1_OtherNode_click_counter] = self.core.mu_to_seconds(SPCM1_OtherNode_click_time)
                    SPCM1_OtherNode_click_counter += 1

                tStamps_t1[excitation_cycle * self.n_excitation_attempts + excitation_attempt] = self.core.mu_to_seconds(t1)
                delay(30 * us)  ### 20us is not enough

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(10*us)

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

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(1*us)

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
                    self.ttl_SPCM0_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)
                    self.ttl_SPCM1_OtherNode_counter.gate_rising(self.t_SPCM_recool_and_shot)

                AllSPCMs_RO_atom_check = int(self.ttl_SPCM0_counter.fetch_count() + \
                                          self.ttl_SPCM1_counter.fetch_count() + \
                                          self.ttl_SPCM0_OtherNode_counter.fetch_count() + \
                                          self.ttl_SPCM1_OtherNode_counter.fetch_count())


                AllSPCMs_RO_atom_check_array[int(excitation_cycle / self.atom_check_every_n)] = AllSPCMs_RO_atom_check

                ### stopping the excitation cycle after the atom is lost
                if AllSPCMs_RO_atom_check / self.t_SPCM_recool_and_shot < self.two_atom_threshold:
                    delay(100 * us)  ### Needs a delay of about 100us or maybe less
                    break

                # self.core.break_realtime()
                two_nodes_synchronization(self)
                delay(1*us)
            delay(10 * us)

        # delay(1 * ms)
        # self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.01 * ms)

        # two_nodes_synchronization(self)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        self.core.break_realtime()

        ### only the elements in range [0:excitation_cycle + 1] contain non-zero values because the loop exits after
        ### the atom is lost. +1 is because python sttops the loop one count earlier.
        for val in AllSPCMs_RO_atom_check_array[0:int(excitation_cycle/self.atom_check_every_n)]:
            self.append_to_dataset('AllSPCMs_RO_atom_check', val)

        delay(1 * ms)
        for i in range((excitation_cycle + 1)* self.n_excitation_attempts):
            self.append_to_dataset('SPCM0_SinglePhoton_tStamps', SPCM0_timestamps[i])
            self.append_to_dataset('SPCM1_SinglePhoton_tStamps', SPCM1_timestamps[i])
            self.append_to_dataset('SPCM0_OtherNode_SinglePhoton_tStamps', SPCM0_OtherNode_timestamps[i])
            self.append_to_dataset('SPCM1_OtherNode_SinglePhoton_tStamps', SPCM1_OtherNode_timestamps[i])
            self.append_to_dataset('reference_tStamps_t1', tStamps_t1[i])

        self.append_to_dataset('n_excitation_cycles', excitation_cycle)

        # delay(1*ms)
        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()

@kernel
def Two_node_single_photon_2_optimization_experiment(self):
    """
    **** This is for Hong-Ou_Mandel Experiment
    **** This is based on Two_node_single_photon_experiment

    **** In HOM Measurement, we want to vary the excitation timing of one node to see
         see the HOM dip. However while we do that, we might accidentally keep the
         FORT from the other node on while one node is collecting photons, making
         false count measurements. Thus, here I added some timing changes in
         excitation cycle loop.

    ****
    """

    self.core.reset()
    delay(1 * ms)

    max_clicks = 2  ### maximum number of clicks that will be time tagged in each gate window.
    ### Have to change SPCM0_SinglePhoton_tStamps in BaseExperiment accordingly.

    AllSPCMs_RO_atom_check_array = [0]

    # record_chopped_optical_pumping(self)
    # delay(100*ms)

    # op_dma_handle = self.core_dma.get_handle("chopped_optical_pumping")

    if self.enable_laser_feedback:
        delay(0.1 * ms)  ### necessary to avoid underflow
        ### set the cooling DP AOM to the MOT settings. Otherwise, DP might be at f_cooling_Ro setting during feedback.
        self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        delay(0.1 * ms)
        run_feedback_and_record_FORT_MM_power(self)

    self.measurement = 0  # advances in end_measurement
    self.core.break_realtime()

    t_diff_mu = 0
    t_gate_mu = 0
    t_photon_add_mu = 0
    t_FORT_OFF_mu = 0

    while self.measurement < self.n_measurements:
        self.core.break_realtime()
        self.ttl_exc0_switch.on()  # turns off the excitation

        ### load atoms in both nodes
        load_until_atom_in_both_nodes_recycle(self)

        delay(1 * ms)

        two_nodes_synchronization(self)

        first_shot(self)
        delay(1 * ms)

        if self.t_recooling_after_first_shot > 0:
            recooling_after_first_shot(self)

        ### first_shot doesn't turn off the fiber AOMs. thus, PR was actually being done with all 6 beams!!!! :(
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        delay(0.1 * ms)
        if not self.PGC_and_RO_with_on_chip_beams:
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

        # two_nodes_synchronization(self)

        ########################################################
        # lower level optical pumping and excitation sequence to optimize for speed
        ########################################################

        ### this will stay on for the entire excition + OP loop, because both the D1 and excitation light use it
        ### use GRIN1 and GRIN2 switches to swith on/off D1 or Exc light
        if self.which_node == "alice":
            self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))
            delay(5 * us)
        else:
            self.GRIN1and2_dds.set(frequency=self.f_GRIN1_D1_pumping, amplitude=dB_to_V(self.p_GRIN1_D1_pumping))
            delay(5 * us)
            self.dds_D1_pumping_DP.set(frequency=self.f_GRIN2_excitation, amplitude=dB_to_V(self.p_GRIN2_excitation))
            self.dds_D1_pumping_DP.sw.on()  # GRIN2 RF ON, external sw not activated yet  # GRIN2 DDS ON

        delay(5*us)

        self.GRIN1and2_dds.sw.on()
        excitation_cycle = 1 ### just for initialization.

        ############################### CW optical pumping phase - pumps atoms into F=1,m_F=0
        delay(10*us)
        if self.t_pumping > 0.0:
            if self.which_node == "alice":
                CW_optical_pumping_node1(self)
            else:
                CW_optical_pumping_node2(self)
            delay(10*us)

        ############################### excitation phase - excite F=1,m=0 -> F'=0,m'=0, detect photon
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=self.stabilizer_excitation.amplitudes[0])
        # self.GRIN1and2_dds.set(frequency=self.f_excitation, amplitude=dB_to_V(self.p_excitation))  ###doing this outside of the loop
        # delay(10*us)

        self.ttl_exc0_switch.off() # turns on the excitation0 AOM
        delay(5 * us)

        self.core.break_realtime()

        two_nodes_synchronization(self)

        for excitation_attempt in range(self.n_excitation_attempts):

            if self.which_node == "bob":
                delay_mu(int(self.t_delay_in_bob_mu))

            t1 = now_mu()

            self.dds_FORT.sw.off()  ### turns FORT off

            at_mu(t1 + 50 + int(self.t_photon_collection_time / ns))
            self.dds_FORT.sw.on()  ### turns FORT on

            at_mu(t1 + int(self.t_excitation_offset_mu))
            self.ttl_GRIN2_switch.off()  # turns on excitation

            at_mu(t1 + int(self.t_excitation_offset_mu) + int(self.t_excitation_pulse / ns))
            self.ttl_GRIN2_switch.on()  # turns off excitation


            self.core.break_realtime()
            two_nodes_synchronization(self)

        delay(20 * us)
        self.ttl_exc0_switch.on()  # block Excitation

        delay(1 * ms)
        self.core.break_realtime()

        self.GRIN1and2_dds.sw.off()

        if self.which_node == "bob":
            self.dds_D1_pumping_DP.sw.off()

        delay(0.1 * ms)

        two_nodes_synchronization(self)

        if self.t_microwave_pulse > 0.0:
            # ### Changing the bias field
            # if self.which_node == "alice":
            #     self.zotino0.set_dac([self.AZ_bottom_volts_microwave, -self.AZ_bottom_volts_microwave,
            #                           self.AX_volts_microwave, self.AY_volts_microwave],
            #                          channels=self.coil_channels)
            # delay(1 * ms)

            # self.zotino0.set_dac([3.5], self.Osc_trig_channel)  ### for triggering oscilloscope
            # delay(0.1 * ms)
            # self.zotino0.set_dac([0.0], self.Osc_trig_channel)

            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.p_FORT_holding * self.stabilizer_FORT.amplitudes[1])
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])

            # FORT_ramp2_smoothstep(self, direction="down")
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[2])
            delay(2 * us)

            self.ttl_microwave_switch.off()
            delay(self.t_microwave_pulse)
            self.ttl_microwave_switch.on()
            delay(2 * us)

            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitudes[1])
            # FORT_ramp2_smoothstep(self, direction="up")
            delay(2 * us)


        if self.t_blowaway > 0.0:
            self.core_dma.playback_handle(ba_dma_handle)

        two_nodes_synchronization(self)

        second_shot(self)

        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(0.1 * ms)

        end_measurement(self)

        self.core.break_realtime()

    # delay(15 * ms)
    self.core.break_realtime()