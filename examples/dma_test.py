"""
A Direct Memory Access example.

For sequences of short pulses ~ 1 us or less, we will quickly run into RTIO underflow errors as the CPU tries to
assign timestamps in realtime but can't keep up. To get around this, we can record a waveform directly to RAM onboard
the core device and then play it back later.

Notes:
    From my tests, the maximum number of pulses that can be saved to RAM is 4096, i.e. 2**12. If you try to record more
    than this, the recording will run without printing any errors or warnings. I found that my loop was truncated
    by analyzing the waveform in gtkwave.

    By printing the cursor value at various places we can verify that the playback of the recording lasts for the
    expected amount of time, that the total readout lasts the correct amount of time (i.e. the recording was played
    back the expected number of times within one readout), and that the loop in run does the readout the expected number
    of times. However, if you analyze the waveform in gtkwave, you will see the DMA sequence played only once! To be
    sure that the waveform was really being played as I expected, I analyzed the signal out of the Urukul (actually,
    I analyzed the diffracted laser power after an AOM driven by the Urukul) on an oscilloscope.

Resources:
https://m-labs.hk/artiq/manual/core_drivers_reference.html?highlight=dma#module-artiq.coredevice.dma
https://github.com/m-labs/artiq/blob/e0ebc1b21dcbb250435c3e36dc4026f0b7d607aa/artiq/examples/kc705_nist_clock/repository/dma_blink.py#L5
"""
# from artiq.experiment import *
# import logging
#
#
# class DMATest(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("core_dma")
#         self.setattr_device("urukul0_cpld")
#         self.setattr_device("urukul0_ch0")
#         self.setattr_device("urukul0_ch1")
#         self.setattr_device("core_dma")
#         self.dds_FORT = self.urukul0_ch0
#         self.dds_cooling_DP = self.urukul0_ch1
#
#         self.t_readout = 10*ms
#         self.t_chop_period = 2*us
#         self.n_chop_cycles = int(self.t_readout/self.t_chop_period)
#         self.n_playbacks = int(self.n_chop_cycles/(2**12)) # round down
#         logging.debug(f"requested readout time {self.t_readout/ms:.2f} ms")
#         logging.debug(f"actual readout time {(self.t_readout*self.n_playbacks*2**12/self.n_chop_cycles)/ms:.2f} ms")
#
#     @kernel
#     def record(self):
#         with self.core_dma.record("chopped_readout"):
#             # FORT pulse period is 2*um so this should last 10 ms
#             for i in range(2**12):
#                 self.dds_FORT.sw.off()
#                 delay(200 * ns)
#                 self.dds_cooling_DP.sw.on()
#                 delay(600 * ns)
#                 self.dds_cooling_DP.sw.off()
#                 delay(200 * ns)
#                 self.dds_FORT.sw.on()
#                 delay(self.t_chop_period/2)
#
#     @rpc(flags={"async"})
#     def print_async(self, x):
#         print(x)
#
#     @kernel
#     def chopped_readout(self):
#         handle = self.core_dma.get_handle("chopped_readout")
#         now = now_mu()
#         for i in range(self.n_playbacks):
#             self.core_dma.playback_handle(handle)
#         after = now_mu() # the cursor advances by 8.19 ms as expected
#         # self.print_async(self.core.mu_to_seconds(after - now))
#         # at_mu(self.core.seconds_to_mu(2**12*self.t_chop_period)+now) #  move the cursor
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.record()
#         self.core.break_realtime()
#         now = now_mu()
#         for i in range(5):
#             delay(10*ms)
#             self.chopped_readout()
#             self.print_async(self.core.mu_to_seconds(now_mu() - now))
#         delay(1*ms)
#         self.dds_cooling_DP.sw.on()
#         self.dds_cooling_DP.sw.off()
#         print("finished DMA test")
#












"""
Experiment 3: Can we gate_rising on DMA?
"""

from artiq.experiment import *

class DMATest(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl7")
        self.setattr_device("ttl13")
        self.setattr_device("ttl0")
        self.setattr_device("ttl0_counter")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.dds_FORT = self.urukul0_ch0

    @kernel
    def Gate_rise_TTL0(self):
        with self.core_dma.record("TTL0_gates"):
            for i in range(20):
                self.ttl0_counter.gate_rising(1 * ms)
                delay(1 * ms)

    @kernel
    def TTL0_pulse_sensitivity(self):
        with self.core_dma.record("TTL0_sensivity_pulses"):
            for i in range(1000):
                self.ttl0._set_sensitivity(1)
                delay(1 * us)
                self.ttl0._set_sensitivity(0)
                delay(1 * us)


    @rpc(flags={"async"})
    def print_async(self, x):
        print(x)


    @kernel
    def run(self):
        self.core.reset()
        self.dds_FORT.sw.on()  ### turns FORT on
        delay(1 * ms)
        self.ttl0._set_sensitivity(0)



        #### The following seems to work as we want; sets the sensitivity using DMA and counts the events afterwards.
        self.TTL0_pulse_sensitivity()
        delay(100 * ms)
        TTL0_gates_handle = self.core_dma.get_handle("TTL0_sensivity_pulses")
        delay(100 * ms)

        self.core_dma.playback_handle(TTL0_gates_handle)
        SPCM0_SinglePhoton = self.ttl0.count(now_mu())





        # #### The following counts only for the first event out of 20.
        # self.Gate_rise_TTL0()
        # delay(100 * ms)
        # TTL0_gates_handle = self.core_dma.get_handle("TTL0_gates")
        # delay(100 * ms)
        #
        # self.core_dma.playback_handle(TTL0_gates_handle)
        # SPCM0_SinglePhoton = self.ttl0_counter.fetch_count()




        # for i in range(10):
        #     delay(100 * us)
        #     self.dds_FORT.sw.on()  ### turns FORT on
        #     delay (100 * us)
        #     self.ttl0_counter.gate_rising(1 * ms)
        #     delay(100*us)
        #     self.dds_FORT.sw.off()  ### turns FORT on

        # self.ttl0_counter.gate_rising(20 * ms)
        # SPCM0_SinglePhoton = self.ttl0_counter.fetch_count()




        delay(100*ms)
        self.dds_FORT.sw.off()  ### turns FORT off
        self.print_async(SPCM0_SinglePhoton)

        print("***********  finished DMA test  ***********")















"""
Experiment 2: Sending two series of ttl pulses in parallel through two channels.
"""

from artiq.experiment import *

class DMATest(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl7")
        self.setattr_device("ttl13")

    @kernel
    def record_pulses_TTL7(self):
        with self.core_dma.record("TTL7_pulses"):
            for i in range(5):
                self.ttl7.pulse(1 * us)
                delay(1 * us)

    @kernel
    def record_pulses_TTL13(self):
        with self.core_dma.record("TTL13_pulses"):
            for i in range(5):
                delay(100 * ns)
                self.ttl13.pulse(1 * us)
                delay(900 * ns)


    @rpc(flags={"async"})
    def print_async(self, x):
        print(x)


    @kernel
    def run(self):
        self.core.reset()

        self.record_pulses_TTL7()
        self.record_pulses_TTL13()
        TTL7_pulse_handle = self.core_dma.get_handle("TTL7_pulses")
        TTL13_pulse_handle = self.core_dma.get_handle("TTL13_pulses")

        delay(100 * us)
        for i in range(2):
            delay(10*us)
            with parallel:
                self.core_dma.playback_handle(TTL7_pulse_handle)
                self.core_dma.playback_handle(TTL13_pulse_handle)

        print("***********  finished DMA test  ***********")











# """
# Experiment 1: The simplest DMA test to run pulses on ttl7.
# """
# from artiq.experiment import *
#
# class DMATest(EnvExperiment):
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("core_dma")
#         self.setattr_device("ttl7")
#
#     @kernel
#     def record_pulses(self):
#         with self.core_dma.record("TTL_pulses"):
#             for i in range(10):
#                 self.ttl7.pulse(1000 * ns)
#                 delay(1000 * ns)
#
#     @rpc(flags={"async"})
#     def print_async(self, x):
#         print(x)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.record_pulses()
#         pulses_handle = self.core_dma.get_handle("TTL_pulses")
#         # self.core.break_realtime()
#         now = now_mu()
#
#         for i in range(5):
#             delay(1000*us)
#             self.core_dma.playback_handle(pulses_handle)
#             self.print_async(self.core.mu_to_seconds(now_mu() - now))
#
#         print("***********  finished DMA test  ***********")
