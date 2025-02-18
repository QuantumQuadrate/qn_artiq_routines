"""
LED blinking, as an example for analyzing RTIO slack.
See here https://forum.m-labs.hk/d/188-analyzing-rtio-slack-with-gtkwave

To analyze the rtio slack run "artiq_coreanalyzer -w yourfilename.vcd" in
an Anaconda prompt in the artiq-master directory. Then open the waveform
with "gtkwave yourfilename.vcd" (assumes you added gtkwave to path).


Only uncomment the run function that you want to analyze.
"""

from artiq.experiment import *


class BlinkForever(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("led0")

    # this will underflow around 15 us (according to the vcd analysis).
    # underflow occurs because events are being added to the scheduler slower than they
    # are being dispatched, so at some point we're scheduling a device to be run at a
    # time stamp that already happened.
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     while True:
    #         self.led0.pulse(.1*us)
    #         delay(.1*us)

    # this seems to go on "forever" with no problem. terminate it in the dashboard.
    # now, the events are being submitted more slowly than they are actually being run
    # and the slack remains positive (and actually increases). Around 48 ms the
    # slack wraps around, and this continues to happen, so this will go until we
    # stop it or the hardware gives out.
    @kernel
    def run(self):
        self.core.reset()
        while True:
            self.led0.pulse(1 * us)
            delay(1 * us)

    # # more complicated example
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     for i in range(10):
    #         for i in range(10):
    #             self.led0.pulse(10 * ns)
    #             delay(10 * ns)
    #         # delay(50 * us)
    #         self.core.reset() # in lieu of a delay, just cheat and reset the clock?
    #         for i in range(10):
    #             self.led0.pulse(10 * ns)
    #             delay(10 * ns)
    #     print("blink example done.")

    # todo: add an overflow example
