from artiq.experiment import *


class TestRTIOClockVsReference(EnvExperiment):
    """TestRTIOClockVsReference

    Definitive RTIO-clock check against the lab frequency reference.

    Feed a reference-locked TTL-level square wave of known frequency (e.g.
    a signal generator locked to the house 10 MHz that also disciplines the
    125 MHz Si5324 reference) into the input of the chosen EdgeCounter
    channel, then run this. The gate duration is timed by the RTIO clock, so

        offset_ppm = (measured_count / (f_ref * gate_s) - 1) * 1e6

    directly compares the RTIO clock to the reference. Locked clock:
    ~0 ppm within +-1 count. Free-running crystal: tens of ppm, unmistakable.

    Defaults: 1 MHz into ttl0_counter for 10 s -> 1e7 counts, 0.1 ppm
    resolution. An EdgeCounter channel is required (it counts in gateware);
    a plain TTL input would overflow its RTIO input FIFO at high rates.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_argument("counter_device", StringValue("ttl0_counter"))
        self.setattr_argument("f_ref_Hz", NumberValue(1e6, ndecimals=0, step=1e6))
        self.setattr_argument("gate_s", NumberValue(10.0, ndecimals=1, step=1.0))

    def prepare(self):
        self.counter = self.get_device(self.counter_device)
        self.gate_time = float(self.gate_s)

    @kernel
    def measure(self) -> TInt32:
        self.core.reset()  # also clears any stale input events
        delay(10 * ms)
        self.counter.gate_rising(self.gate_time)
        return self.counter.fetch_count()

    def run(self):
        print("opening a %.1f s gate on %s (silent until it closes)..."
              % (self.gate_time, self.counter_device))
        count = self.measure()
        expected = float(self.f_ref_Hz) * self.gate_time
        ppm = (count / expected - 1.0) * 1e6
        print("counted %d edges | expected %.0f | offset: %+.3f ppm"
              % (count, expected, ppm))
        if count == 0:
            print("ZERO counts -> no signal reaching %s: check cabling, "
                  "level (0-3.3 V into high-Z), and channel name."
                  % self.counter_device)
        elif abs(ppm) < 0.5:
            print("RTIO clock agrees with the reference -> locked.")
        elif count > expected * 1.5:
            print("Gross excess counts -> ringing / multiple triggers per "
                  "edge; fix the signal (slower edges, Schmitt, shorter "
                  "cable), not the clock.")
        else:
            print("RTIO clock does NOT match the reference at this level.")
