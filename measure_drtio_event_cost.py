"""measure_drtio_event_cost.py

How much more lead time does an RTIO event need on the satellite than on the
master? Answering that turns the Node2 underflows from guesswork into a
number you can put in a delay().

WHY THIS PAIR OF CHANNELS
-------------------------
    spi_urukul1  ->  channel 0x00002b  (destination 0, master   = Node1 crate)
    spi_urukul4  ->  channel 0x01002b  (destination 1, satellite = Node2 crate)

Same peripheral, same slot, same local channel number 0x2b -- the ONLY
difference is which side of the fibre it lives on. Anything that shows up as a
difference between the two is the DRTIO crossing and nothing else.
spi_urukul4 is also the exact channel (65579) that underflowed in
load_until_atom_smooth_FORT_recycle.

WHY THIS IS SAFE TO RUN AT ANY TIME
-----------------------------------
The probe is ``bus.set_config_mu(...)``, which is a single ``rtio_output()``
to the SPI core's config register (spi2.py:169). It loads the transfer
parameters and nothing else: no chip select is asserted (``cs=0``), no SCK is
generated, no data is shifted. Not one pin moves on the Urukul, and every
ARTIQ routine that later uses the bus (``cpld.sta_read``, ``dds.set``, ...)
programs the config register itself before transferring, so the register we
leave behind is overwritten before it is ever used.

Note this deliberately does NOT probe with ``cpld.sta_read()``: that looks
read-only but writes ``self.cfg_reg`` back to the CPLD (urukul.py:231) from
the *device_db defaults*, which would silently drop the RF switches on the
probed Urukul.

Nothing else is bound: no coils, no DDS, no TTL, no core.reset() unless you
tick the box.

WHAT IT MEASURES
----------------
Test A -- submission cost. Issue a short burst of events from comfortable
slack and watch how far the RTIO counter moves. That is how fast the CPU can
push events onto that channel, i.e. the rate at which slack is spent. The
burst is kept below the FIFO depth on purpose: a longer burst would stall on
a full FIFO and measure the programmed delay instead of the submission cost.

Test B -- minimum viable lead (the number that matters). Place the timeline
cursor a known distance ahead of the live RTIO counter, submit one event, and
see whether it is accepted. Walking that distance down until RTIOUnderflow
appears gives the smallest lead the channel tolerates. The gap between the
local and remote answers is the extra delay a satellite event needs.

RESULT, 2026-09-10 (RID 38429)
------------------------------
There is NO measurable DRTIO penalty. Master and satellite came out
indistinguishable on both tests:

    submission cost   master 713 mu/event   satellite 711 mu/event
    minimum lead      master 1000 mu        satellite 1000 mu
                      (both 8/8 at 1000 mu, both 0/8 at 500 mu)

So a satellite event is no more demanding than a master one, and the
"crossing the fibre costs slack" theory is dead. ARTIQ buffers events at the
satellite ahead of time -- the DEST#1 boot log even reports "buffer space is
128" -- so the link latency is absorbed rather than charged to the caller.

That reframes the Node2 underflow entirely. The failure had -585808 mu of
slack, roughly 500x past the 1000 mu floor, which no plausible delay() would
have covered. The slack is spent BEFORE the call, by host round trips --
set_dataset, append_to_dataset, print_async -- which advance wall clock while
the timeline cursor stands still. It is slack EROSION across loop iterations,
not per-event cost, and it is not Node2-specific: Node1 hit the same class of
failure on channel 39 (ttl_urukul0_sw0, destination 0) in the feedback path.

Re-run this after any gateware or link change; it is the control that tells
you whether a new underflow is about DRTIO at all.
"""

import numpy as np

from artiq.experiment import *
from artiq.coredevice.exceptions import RTIOUnderflow


LOCAL_TAG = 0
REMOTE_TAG = 1

_TAG_NAMES = {
    LOCAL_TAG: "spi_urukul1  dest 0  master",
    REMOTE_TAG: "spi_urukul4  dest 1  satellite",
}

#: Lead times to try, in machine units (1 mu = 1 ns), largest first. Walking
#: down means every event is submitted later than the previous one, so the
#: per-channel timestamps stay monotonic and we never trip a sequence error.
_LEAD_LADDER_MU = (
    2000000, 1000000, 500000, 200000, 100000, 50000, 20000,
    10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10,
)


class measure_drtio_event_cost(EnvExperiment):
    """measure_drtio_event_cost"""

    def build(self):
        self.setattr_device("core")
        self.setattr_device("spi_urukul1")  # destination 0, master
        self.setattr_device("spi_urukul4")  # destination 1, satellite

        self.setattr_argument(
            "n_repeats",
            NumberValue(8, ndecimals=0, step=1, min=1, max=100),
            tooltip="Attempts per lead time. A lead only counts as viable if "
                    "every attempt is accepted.",
        )
        self.setattr_argument(
            "n_burst",
            NumberValue(16, ndecimals=0, step=1, min=1, max=32),
            tooltip="Events per burst in Test A. Keep below the RTIO output "
                    "FIFO depth or the burst measures FIFO backpressure.",
        )
        self.setattr_argument(
            "reset_core",
            BooleanValue(False),
            tooltip="Clear stale RTIO errors before measuring. WARNING: this "
                    "resets the RTIO core, which drops TTL output states "
                    "(AOM switches). Leave off unless results look nonsense.",
        )

    def prepare(self):
        self.n_repeats = int(self.n_repeats)
        self.n_burst = int(self.n_burst)
        self.lead_ladder = [np.int64(value) for value in _LEAD_LADDER_MU]
        self.ladder_results = {LOCAL_TAG: [], REMOTE_TAG: []}
        self.burst_results = {}

    # ------------------------------------------------------------------ RPC

    def record_lead(self, tag, lead_mu, n_ok, n_try):
        """Called from the kernel after each rung of the ladder."""
        self.ladder_results[tag].append((int(lead_mu), int(n_ok), int(n_try)))
        verdict = "ok" if n_ok == n_try else ("UNDERFLOW" if n_ok == 0 else "marginal")
        print("  %-30s lead %9d mu  %d/%d  %s"
              % (_TAG_NAMES[tag], int(lead_mu), int(n_ok), int(n_try), verdict))

    # --------------------------------------------------------------- kernels

    @kernel
    def clear_rtio(self):
        self.core.reset()
        delay(10 * ms)

    @kernel
    def burst_cost_mu(self, bus, n: TInt32) -> TInt64:
        """RTIO counter advance while submitting n inert config events."""
        self.core.break_realtime()
        delay(1 * ms)
        t_start = self.core.get_rtio_counter_mu()
        for _ in range(n):
            bus.set_config_mu(0, 8, 6, 0)
        t_stop = self.core.get_rtio_counter_mu()
        return t_stop - t_start

    @kernel
    def scan_leads(self, bus, tag: TInt32):
        """Walk the lead ladder down, counting acceptances at each rung."""
        for lead_mu in self.lead_ladder:
            n_ok = 0
            for _ in range(self.n_repeats):
                self.core.break_realtime()
                try:
                    # Absolute placement against the live counter: this is the
                    # slack the event is submitted with, by construction.
                    at_mu(self.core.get_rtio_counter_mu() + lead_mu)
                    bus.set_config_mu(0, 8, 6, 0)
                    # Force the event to retire here so a deferred error is
                    # attributed to this attempt, and so the next timestamp is
                    # guaranteed to be later than this one.
                    self.core.wait_until_mu(now_mu())
                    n_ok += 1
                except RTIOUnderflow:
                    pass
            self.record_lead(tag, lead_mu, n_ok, self.n_repeats)

    # ------------------------------------------------------------------ host

    def run(self):
        if self.reset_core:
            print("resetting RTIO core (TTL outputs dropped)...")
            self.clear_rtio()

        probes = ((LOCAL_TAG, self.spi_urukul1), (REMOTE_TAG, self.spi_urukul4))

        print("Test A -- submission cost of %d inert config events" % self.n_burst)
        for tag, bus in probes:
            samples = [int(self.burst_cost_mu(bus, self.n_burst))
                       for _ in range(self.n_repeats)]
            best = min(samples)
            self.burst_results[tag] = best
            print("  %-30s %7d mu total  %8.1f mu/event  (best of %d)"
                  % (_TAG_NAMES[tag], best, best / self.n_burst,
                     self.n_repeats))

        print("")
        print("Test B -- minimum viable lead time")
        for tag, bus in probes:
            self.scan_leads(bus, tag)

        self.report()

    def _min_viable_lead(self, tag):
        """Smallest lead on the ladder that every attempt survived."""
        viable = [lead for lead, n_ok, n_try in self.ladder_results[tag]
                  if n_ok == n_try]
        return min(viable) if viable else None

    def report(self):
        local_lead = self._min_viable_lead(LOCAL_TAG)
        remote_lead = self._min_viable_lead(REMOTE_TAG)

        print("")
        print("=" * 68)
        for tag in (LOCAL_TAG, REMOTE_TAG):
            lead = self._min_viable_lead(tag)
            self.set_dataset("drtio_probe_burst_mu_%d" % tag,
                             self.burst_results.get(tag, 0), broadcast=True)
            self.set_dataset("drtio_probe_min_lead_mu_%d" % tag,
                             -1 if lead is None else lead, broadcast=True)
            print("%-30s min viable lead: %s"
                  % (_TAG_NAMES[tag],
                     "none on ladder" if lead is None else "%d mu" % lead))

        if local_lead is None or remote_lead is None:
            print("")
            print("At least one channel failed at every lead time on the "
                  "ladder. Either the link is down or something else is "
                  "wrong -- do not read a delay value out of this run.")
            return

        extra = remote_lead - local_lead
        print("")
        print("Extra lead the satellite needs: %d mu (%.3f ms)"
              % (extra, extra / 1e6))

        if extra <= 0:
            print("The satellite is no more demanding than the master here, "
                  "so the underflow is not a plain per-event DRTIO cost. "
                  "Look instead at how much slack is spent before the call: "
                  "host RPCs and dataset writes in the loading loop burn "
                  "wall-clock time while the cursor stands still.")
        else:
            print("So a satellite event submitted with less than ~%.3f ms of "
                  "slack is at risk. The failing call in "
                  "load_until_atom_smooth_FORT_recycle had 0.1 ms of delay "
                  "in front of it." % (remote_lead / 1e6))
            print("Either give remote-bound sequences at least this much "
                  "lead, or use core.break_realtime() where the preceding "
                  "timing does not matter.")
