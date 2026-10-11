"""measure_urukul_sync_windows.py

Map the real SYNC_IN timing window of every Urukul channel, instead of
trusting whatever number tune_sync_delay happened to converge on.

WHY
---
artiq_sinara_tester reports one delay per channel, found by hill-climbing
from the seed stored in EEPROM. tune_sync_delay only scans +-6 taps around
that seed (ad9910.py: ``search_span = 13``, with a FIXME pointing at
sinara-hw/Urukul#16 about Kasli SYNC_IN jitter), so the answer it gives
depends on where it started. A channel whose true optimum is far from its
stored seed reports "no valid window/delay" even though a window exists, and
a channel with a fuzzy window reports a different value every run.

This sweeps each channel from seeds spanning the whole 0..31 delay line, so
every window in range is reachable from at least one starting point, and
reports what was found from each.

READ THE WINDOW, NOT THE DELAY
------------------------------
The 2026-09-10 sweeps settled how to interpret this, and it is not what the
delay spread suggests. EVERY channel, healthy ones included, returns delays
in about three clusters roughly 12 taps apart -- because a clean sampling
point repeats once per SYSCLK period, and tune_sync_delay simply returns
whichever repeat is nearest its seed. (Measured pitch on this system: 11.8
taps per 1 ns SYSCLK period, i.e. ~85 ps per tap against the AD9910's
~75 ps nominal.) Every repeat catches the same SYNC_IN edge, since SYNC_IN
is 62.5 MHz (16 ns) while the whole delay line spans only ~2.5 ns, so they
are equally valid. A wide delay spread across seeds is therefore EXPECTED
and says nothing about health -- and it also explains the scatter
sinara_tester reports, which is just different seeds landing in different
repeats.

The discriminator is the WINDOW. Window 0 means the delay was clean at zero
validation margin but could not hold any wider one -- no setup/hold margin.
Channels like that still pass sometimes, which is exactly the intermittency
that made the old Node2 cards look flaky for a year.

HOW MANY ZEROS MAKE A FAULT
---------------------------
Count the seeds that return window 0:

    0 zeros     healthy
    1 zero      inconclusive -- re-run before concluding anything
    2+ zeros    MARGINAL

Calibrated on three runs. The faulty channels (old urukul4_ch2 and all four
of the old urukul5) returned 4-5 zeros on RID 38430 and 2-6 on RID 38431 --
every run. Healthy channels returned none on those runs, but on RID 38432,
after the replacement, two healthy channels each returned exactly one, both
from seed 3 at the bottom of the delay line -- one of them urukul0_ch0, a
Node1 card that was never touched. So a single zero occurs on healthy
channels, and a faulty one can dip as low as two. A threshold of two
classifies all three runs without error; the two earlier rules ("half the
seeds", then "any zero") each got one of them wrong.

WHAT THE NUMBERS MEAN
---------------------
delay  - position on the AD9910 SYNC_IN delay line, 0..31 taps, nominally
         ~75 ps each (~85 ps measured on this system), so the line spans
         roughly 2.5 ns. This is TIME, not frequency.
window - validation width in the same tap units: how wide a region around
         the sampling point stayed free of the chip's SMP_ERR flag. Bigger
         is better; it is the margin you have against drift.

SAFETY
------
Initializes the Urukuls exactly the way artiq_sinara_tester does (cpld.init
then per-channel init), so treat it like a tester run: RF switch and
attenuator state is left at init defaults, and you should run your normal
experiment afterwards to restore it. Nothing else is bound.

It deliberately does NOT write the EEPROM. The tester rewrites calibration
as it goes, which is what makes a failing run leave half-updated seeds
behind; this one only measures. Per-channel init uses blind=True so that
init's own tune_sync_delay call cannot abort the sweep on a bad channel --
the whole point is to measure the channels that currently fail.

PLL lock is checked and reported separately, so "the PLL is not locked" can
never be mistaken for "the SYNC window is bad".
"""

import numpy as np

from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_pll_lock


#: Seeds spanning the delay line. tune_sync_delay reaches about +-6 taps, so
#: a spacing of 4 leaves no gap and gives overlapping views of each window.
SEEDS = (3, 7, 11, 15, 19, 23, 27)

#: Zero-margin seeds needed to call a channel MARGINAL. A single zero is
#: within the noise of this measurement -- see HOW MANY ZEROS MAKE A FAULT.
MARGINAL_ZEROS = 2


def _zero_margin(windows):
    return sum(1 for window in windows if window == 0)


class measure_urukul_sync_windows(EnvExperiment):
    """measure_urukul_sync_windows"""

    def build(self):
        self.setattr_device("core")
        self.setattr_argument(
            "cards",
            StringValue("0,1,2,3,4,5"),
            tooltip="Urukul card numbers to sweep. 0-2 are Node1 (master), "
                    "3-5 are Node2 (satellite). Keep the healthy cards in "
                    "for comparison -- they are the control.",
        )

    def prepare(self):
        self.card_numbers = [int(part) for part in self.cards.split(",")
                             if part.strip()]
        self.seeds = [np.int32(seed) for seed in SEEDS]

        self.cplds = []
        self.channels = []
        self.channel_names = []
        for card in self.card_numbers:
            self.cplds.append(self.get_device(f"urukul{card}_cpld"))
            for index in range(4):
                name = f"urukul{card}_ch{index}"
                self.channels.append(self.get_device(name))
                self.channel_names.append(name)

        self.results = [[] for _ in self.channel_names]
        self.pll_locked = [-1] * len(self.channel_names)

    # ------------------------------------------------------------------ RPCs

    def record_lock(self, channel_index, locked):
        self.pll_locked[channel_index] = int(locked)

    def record(self, channel_index, seed, found_delay, found_window):
        self.results[channel_index].append(
            (int(seed), int(found_delay), int(found_window))
        )

    def announce(self, channel_index):
        print("  sweeping {}...".format(self.channel_names[channel_index]))

    # --------------------------------------------------------------- kernel

    @kernel
    def run_sweep(self):
        for card_index in range(len(self.cplds)):
            self.core.break_realtime()
            self.cplds[card_index].init()
            delay(1 * ms)

        for i in range(len(self.channels)):
            self.announce(i)

            # blind=True skips init's own tune_sync_delay, which is exactly
            # the call that raises on the channels we are here to measure.
            self.core.break_realtime()
            self.channels[i].init(True)
            delay(1 * ms)

            status = self.channels[i].cpld.sta_read()
            locked = (urukul_sta_pll_lock(status)
                      >> (self.channels[i].chip_select - 4)) & 1
            delay(1 * ms)
            self.record_lock(i, locked)

            for s in range(len(self.seeds)):
                self.core.break_realtime()
                found_delay = -1
                found_window = -1
                try:
                    found_delay, found_window = \
                        self.channels[i].tune_sync_delay(self.seeds[s])
                except ValueError:
                    found_delay = -1
                    found_window = -1
                self.record(i, self.seeds[s], found_delay, found_window)

    # ----------------------------------------------------------------- host

    def run(self):
        print("Sweeping SYNC_IN delay windows for {} channels, seeds {}."
              .format(len(self.channel_names), list(SEEDS)))
        print("This does NOT write the EEPROM.")
        self.run_sweep()
        self.report()

    @staticmethod
    def _verdict(windows):
        """Classify a channel by its setup/hold margin.

        Judged on window width only. The delay values legitimately differ
        between seeds -- clean sampling points repeat every SYSCLK period
        (~12 taps) and the tuner returns whichever repeat is nearest its
        seed -- so delay spread carries no information about health.
        """
        if not windows:
            return "NO WINDOW from any seed"
        zeros = _zero_margin(windows)
        if zeros == 0:
            return "healthy (margin >= {} at every seed)".format(min(windows))
        if zeros < MARGINAL_ZEROS:
            return ("inconclusive -- zero margin at {}/{} seeds, re-run"
                    .format(zeros, len(SEEDS)))
        return ("MARGINAL -- zero setup/hold margin at {}/{} seeds"
                .format(zeros, len(SEEDS)))

    def report(self):
        header = "{:<16}{:<6}".format("channel", "PLL")
        for seed in SEEDS:
            header += "{:<9}".format("s%d" % seed)
        print("")
        print(header + "verdict")
        print("-" * (len(header) + 40))

        suspects = []
        rerun = []
        for i, name in enumerate(self.channel_names):
            row = "{:<16}{:<6}".format(
                name, "ok" if self.pll_locked[i] == 1 else "LOCK?")
            windows = []
            for _, found_delay, found_window in self.results[i]:
                if found_delay < 0:
                    row += "{:<9}".format("-")
                else:
                    row += "{:<9}".format("%d/%d" % (found_delay, found_window))
                    windows.append(found_window)
            print(row + self._verdict(windows))

            zeros = _zero_margin(windows)
            if not windows or zeros >= MARGINAL_ZEROS:
                suspects.append(name)
            elif zeros:
                rerun.append(name)
            self.set_dataset("urukul_sync_windows_%s" % name, windows,
                             broadcast=True)

        print("")
        print("delay/window are in delay-line taps (~85 ps measured here). "
              "Judge a channel by the WINDOW: it is the setup/hold margin, "
              "and 0 means none at all.")
        print("Differing delays between seeds are EXPECTED -- a clean sampling "
              "point repeats every SYSCLK period (~12 taps), and the tuner "
              "returns whichever repeat is nearest its seed. All of them catch "
              "the same 62.5 MHz SYNC_IN edge, so they are equally valid.")
        print("")
        if suspects:
            print("SUSPECT CHANNELS ({}): {}".format(len(suspects),
                                                     ", ".join(suspects)))
            cards = sorted({name.split("_")[0] for name in suspects})
            print("Cards involved: {}. A card where ALL FOUR channels are "
                  "suspect points at the card; a lone channel points at that "
                  "one AD9910.".format(", ".join(cards)))
        else:
            print("No channel is marginal ({}+ zero-margin seeds). SYNC is "
                  "healthy.".format(MARGINAL_ZEROS))
        if rerun:
            print("Single zero -- re-run to confirm ({}): {}. One zero has been "
                  "seen on known-healthy channels; faulty ones showed 2+ on "
                  "every run.".format(len(rerun), ", ".join(rerun)))
        print("Nothing was written to EEPROM. Run your normal experiment to "
              "restore switch and attenuator state.")
