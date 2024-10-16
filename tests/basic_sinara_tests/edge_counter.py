"""
test the edge counter. adapted from https://forum.m-labs.hk/d/501-delayed-response-of-edge-counter
"""

from artiq.experiment import *

class EdgeCounterTest(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ccb")
        self.setattr_device("scheduler")
        self.setattr_device("ttl0_counter")

        self.setattr_argument("t_accu", NumberValue(
            default=100 * ms,
            unit="ms"))

    def prepare(self):
        self.ccb.issue("create_applet", "counts_0",
                       "${artiq_applet}plot_xy test.counts_0", ["Test"])
        self.set_dataset("test.counts_0", [0], broadcast=True)

    def run(self):
        self.run_kernel()

    @kernel
    def run_kernel(self):
        while True:
            self.core.break_realtime()
            with parallel:
                self.ttl0_counter.gate_rising(self.t_accu)

            counts_0 = self.ttl0_counter.fetch_count()
            if self.cb(counts_0):
                break

    def check_stop(self) -> TBool:
        try:
            self.scheduler.pause()
        except TerminationRequested:
            return True
        return False

    def cb(self, counts_0) -> TBool:
        self.append_to_dataset("test.counts_0", counts_0)

        try:
            self.scheduler.pause()
        except TerminationRequested:
            return True
        return False