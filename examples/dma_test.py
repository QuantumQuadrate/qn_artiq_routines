"""
A Direct Memory Access example.

For sequences of short pulses ~ 1 us or less, we will quickly run into RTIO underflow errors as the CPU tries to
assign timestamps in realtime but can't keep up. To get around this, we can record a waveform directly to RAM onboard
the core device and then play it back later.

Notes:
    From my tests, the maximum number of pulses that can be saved to RAM is 4096, i.e. 2**12. If you try to record more
    than this, the recording will run without printing any errors or warnings. I found that my loop was truncated
    by analyzing the waveform in gtkwave.

Resources:
https://m-labs.hk/artiq/manual/core_drivers_reference.html?highlight=dma#module-artiq.coredevice.dma
https://github.com/m-labs/artiq/blob/e0ebc1b21dcbb250435c3e36dc4026f0b7d607aa/artiq/examples/kc705_nist_clock/repository/dma_blink.py#L5
"""
from artiq.experiment import *


class DMATest(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("core_dma")
        self.dds_FORT = self.urukul0_ch0
        self.dds_cooling_DP = self.urukul0_ch1

    @kernel
    def record(self):
        with self.core_dma.record("chopped_readout"):
            # FORT pulse period is 2*um so this should last 10 ms
            for i in range(5000):
                self.dds_FORT.sw.off()
                delay(200 * ns)
                self.dds_cooling_DP.sw.on()
                delay(600 * ns)
                self.dds_cooling_DP.sw.off()
                delay(200 * ns)
                self.dds_FORT.sw.on()
                delay(3 * us)

    @kernel
    def run(self):
        self.core.reset()
        self.record()
        handle = self.core_dma.get_handle("chopped_readout")
        self.core.break_realtime()
        for i in range(5):
            delay(1*ms)
            self.core_dma.playback_handle(handle)
        print("finished DMA test")
