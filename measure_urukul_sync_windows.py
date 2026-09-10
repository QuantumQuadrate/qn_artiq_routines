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

That is why urukul4 came back as 19/9/15/10 on 2026-09-10: four channels on
one card share the same SYSCLK and the same CPLD-distributed SYNC, so their
optima must agree within about a tap. A ten-tap spread (~750 ps) is not
physical skew, it is an unreliable measurement.

This sweeps each channel from seeds spanning the whole 0..31 delay line, so
every window in range is reachable from at least one starting point, and
reports what was found from each. Read the result as:

  * same delay from most seeds, decent window   -> real, trustworthy optimum
  * scattered delays / tiny windows             -> marginal SYNC sampling
  * nothing found from any seed                 -> no usable window at all

WHAT THE NUMBERS MEAN
---------------------
delay  - position on the AD9910 SYNC_IN delay line, 0..31 taps of ~75 ps
         (so the whole range spans ~2.3 ns). This is TIME, not frequency.
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
    def _verdict(hits, best_window):
        """Classify a channel from the delays found across all seeds."""
        if not hits:
            return "NO WINDOW from any seed"
        spread = max(hits) - min(hits)
        if len(hits) >= len(SEEDS) - 1 and spread <= 2:
            return "consistent (spread {} taps)".format(spread)
        if spread > 4:
            return "SCATTERED (spread {} taps) -- marginal".format(spread)
        if len(hits) <= len(SEEDS) // 2:
            return "only {}/{} seeds found a window".format(len(hits), len(SEEDS))
        return "usable (spread {} taps, window {})".format(spread, best_window)

    def report(self):
        header = "{:<16}{:<6}".format("channel", "PLL")
        for seed in SEEDS:
            header += "{:<9}".format("s%d" % seed)
        print("")
        print(header + "verdict")
        print("-" * (len(header) + 40))

        for i, name in enumerate(self.channel_names):
            row = "{:<16}{:<6}".format(
                name, "ok" if self.pll_locked[i] == 1 else "LOCK?")
            hits = []
            best_window = -1
            for _, found_delay, found_window in self.results[i]:
                if found_delay < 0:
                    row += "{:<9}".format("-")
                else:
                    row += "{:<9}".format("%d/%d" % (found_delay, found_window))
                    hits.append(found_delay)
                    best_window = max(best_window, found_window)
            print(row + self._verdict(hits, best_window))

            self.set_dataset("urukul_sync_hits_%s" % name, hits, broadcast=True)

        print("")
        print("delay/window are in ~75 ps taps. A channel is healthy when most "
              "seeds agree to within a tap or two and the window is not 0.")
        print("Channels on one card share SYSCLK and SYNC, so their optima "
              "should agree; a wide spread within a card means the measurement "
              "is unreliable, not that the channels differ.")
        print("Nothing was written to EEPROM. Run your normal experiment to "
              "restore switch and attenuator state.")
