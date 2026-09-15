# Node2 Urukul SYNC window margin

**Date:** 2026-09-10 **Status:** **Resolved 2026-09-15** — both cards replaced.
**Hardware:** Node2 crate, `urukul4` and `urukul5` (satellite, destination 1).
In standalone naming these were Node2's `urukul1` and `urukul2`.

---

## Resolution (2026-09-15)

`urukul4` and `urukul5` were replaced with spare v1.5 and v1.5.6 cards — the
same CPLD generation and `proto_rev` 8, so no software or device_db change.

* `artiq_sinara_tester` passed **on the first run**: no `no valid
  window/delay`, and no re-running until it creeps through. From the common
  seed 15, all four channels of each new card landed within two taps of each
  other (urukul4: 13/15/15/14, urukul5: 18/18/17/19). The old urukul4
  scattered 19/9/15/10 from the same seed.
* `measure_urukul_sync_windows` (RID 38432): the failure signature is gone.
  Every channel on both new cards returned margin from every seed but one —
  `urukul4_ch1` returned window 0 once, from seed 3 at the bottom of the
  delay line. `urukul0_ch0`, a Node1 card never touched and clean on
  2026-09-10, did exactly the same in this run. The old cards returned
  window 0 on 2–6 of 7 seeds, every run. See "How many zeros make a fault"
  below for why a single zero is noise.

**A tester pass is not proof of margin.** On 2026-09-10 the tester passed
three of the five channels the sweep flagged (old urukul4_ch2 `15 3`,
urukul5_ch0 `20 0`, urukul5_ch1 `20 0`), failed on the fourth, and never
reached the fifth. It runs one hill-climb per channel and prints where it
landed; it never reports the window. Certify cards with the sweep.

**Why v1.5.x spares and not the new v1.6 cards.** v1.6 replaces the CPLD with
an FPGA (iCE40), and the current Urukul gateware is `proto_rev` 9. ARTIQ 7's
driver accepts only 8 and raises `Urukul proto_rev mismatch`; v1.6 support
arrived in ARTIQ 9. Whether a `proto_rev` 8 build exists for the v1.6 FPGA
is unconfirmed — ask M-Labs before putting those cards on ARTIQ 7.

---

## Verdict (2026-09-10, before replacement)

`urukul5` has a **card-level fault**: all four channels have zero SYNC
setup/hold margin. `urukul4_ch2` has the same fault on **one channel only**;
its other three channels are as good as Node1's.

Everything else — all of Node1 (`urukul0/1/2`) and, importantly, `urukul3` in
the *same Node2 crate* — is clean.

**Action: replace `urukul5`.** Replacing `urukul4` is optional and only
buys you its one bad channel.

Because `urukul3` shares the Node2 crate, the same SYNC distribution, the
same EEM chain and the same MMCX clock, and measures perfectly, **the crate,
the cabling and the clock distribution are all exonerated.** That is the
swap-test answer without having to pull the rack apart.

---

## Symptom (about one year, since the 2025 incident)

* `artiq_sinara_tester` fails partway through the Urukul test with
  `ValueError: no valid window/delay`.
* Re-running the tester several times gets progressively further each time
  and eventually passes completely.
* It comes back after every power cycle, so checking it became part of the
  boot routine.
* Only ever these two cards. Node1's three were never affected.

It began after an accidental `GeneralVariableScan` over `t_microwave_pulse`
with a range entered in units of 335 MHz, which ARTIQ read as 335 Ms — about
3877 days. The experiment was terminated, but afterwards working code
(AOMsCoils, GVS) stopped running.

The tester failure landed on the microwave DDS channel, which is
`urukul2_ch3` in standalone naming = **`urukul5_ch3`** today.

---

## What the error actually is

`ValueError: no valid window/delay` is raised at `ad9910.py:1038`, inside
`tune_sync_delay()`. That is the **per-channel SYNC_IN delay calibration** —
phase 2 of the tester, after every CPLD has already initialised.

`artiq_sinara_tester.py:281` runs the test in two distinct phases:

1. **Per card** — `cpld.init()` then the attenuator test, printing
   `"urukulN_cpld: initializing CPLD..."`.
2. **Per channel** — `"Calibrating inter-device synchronization..."` then
   `calibrate_urukul()` on each of the 24 channels.

A failure in phase 1 would be a CPLD, SPI or clock problem. This one is
firmly in phase 2. The neighbouring faults have their own distinct messages,
so they can be told apart at a glance:

| Message | Subsystem |
| --- | --- |
| `no valid window/delay` | SYNC_IN delay calibration ← **this issue** |
| `PLL lock timeout` | AD9910 PLL / reference clock |
| `Urukul AD9910 AUX_DAC mismatch` | SPI communication |

---

## Measurement

Run **`measure_urukul_sync_windows`** (repo top level). It sweeps seven
seeds spanning the whole 0–31 delay line so every window in range is
reachable, and reports what each seed found as `delay/window`. It does
**not** write the EEPROM.

Result, 2026-09-10, RID 38430:

```
channel         PLL   s3     s7     s11    s15    s19    s23    s27
urukul0_ch0     ok    0/1    11/2   11/2   13/1   24/2   24/2   25/1
urukul0_ch1     ok    0/1    10/2   11/2   13/1   22/1   24/2   24/2
urukul0_ch2     ok    1/2    12/2   12/2   13/2   13/2   25/2   25/2
urukul0_ch3     ok    0/2    12/2   11/1   12/2   24/2   24/2   24/2
urukul1_ch0     ok    3/2    4/2    15/2   15/2   15/2   26/2   27/2
urukul1_ch1     ok    4/2    5/2    16/2   16/2   17/2   28/2   28/2
urukul1_ch2     ok    4/2    5/2    16/2   16/2   17/2   28/2   28/2
urukul1_ch3     ok    4/2    5/2    15/2   15/2   15/2   26/2   26/2
urukul2_ch0     ok    3/1    2/2    13/2   14/2   14/2   26/2   27/2
urukul2_ch1     ok    2/1    2/1    13/2   13/2   13/2   24/2   25/2
urukul2_ch2     ok    2/2    2/2    13/2   13/2   14/1   24/1   25/2
urukul2_ch3     ok    2/2    2/2    14/2   14/2   15/2   26/2   26/2
urukul3_ch0     ok    1/2    12/2   12/2   13/2   24/2   24/2   26/2
urukul3_ch1     ok    1/2    12/2   12/2   13/2   13/2   25/2   26/2
urukul3_ch2     ok    1/2    1/2    13/2   14/2   14/2   25/2   27/2
urukul3_ch3     ok    1/2    13/2   12/2   13/2   24/2   24/2   25/2
urukul4_ch0     ok    7/2    7/2    9/2    19/1   20/2   22/2   23/1
urukul4_ch1     ok    8/2    8/2    9/2    9/2    21/2   23/1   22/2
urukul4_ch2     ok    1/0    7/0    11/0   11/2   24/2   24/2   27/0   <-- 4 zeros
urukul4_ch3     ok    8/2    8/2    10/2   10/2   21/2   22/2   22/2
urukul5_ch0     ok    6/2    7/0    11/0   20/2   19/0   23/0   25/0   <-- 5 zeros
urukul5_ch1     ok    3/0    7/2    11/0   16/0   20/2   23/0   21/2   <-- 4 zeros
urukul5_ch2     ok    3/0    7/0    11/0   13/0   19/0   22/2   22/2   <-- 5 zeros
urukul5_ch3     ok    5/0    7/0    11/0   9/2    20/2   23/0   24/0   <-- 5 zeros
```

Twenty channels returned a non-zero window from **every** seed. Four
returned window 0 from most seeds. Nothing in between — the separation is
absolute, and it lands exactly on `urukul4_ch2` plus all of `urukul5`.

PLL lock was confirmed on all 24 channels, so this is not a clock or PLL
problem.

### Reproducibility: count the zeros, don't weight them

A second sweep (RID 38431, same session, no power cycle) gave the same
answer only once the verdict rule was fixed:

| channel | run 38430 | run 38431 |
| --- | --- | --- |
| urukul4_ch1 | 0 zeros | 0 zeros |
| urukul4_ch2 | 4 zeros | 6 zeros |
| urukul4_ch3 | 0 zeros | 0 zeros |
| urukul5_ch0 | 5 zeros | 4 zeros |
| urukul5_ch1 | 4 zeros | 6 zeros |
| urukul5_ch2 | 5 zeros | **2 zeros** |
| urukul5_ch3 | 5 zeros | **3 zeros** |

**How many** seeds hit zero is not reproducible — these channels sit right
on the boundary, so the count swings run to run. `urukul5_ch2` moved from
5/7 to 2/7. **Which** channels produce a zero at all looked perfectly
reproducible here: the same five in both runs, and every healthy channel
produced none in either. A third run disproved that — see below.

The tool originally flagged a channel only when at least half the seeds hit
zero, which split these two runs into different verdicts for identical
hardware. It was then changed to **any zero is suspect**, with the warning
that if a future run disagreed, the rule should be suspected before the
hardware. That is what happened.

### How many zeros make a fault (after RID 38432)

The post-replacement sweep (RID 38432) gave exactly one zero on each of two
healthy channels — `urukul0_ch0`, a Node1 card never touched and clean on
RID 38430, and the new `urukul4_ch1` — both from seed 3, at the bottom of
the delay line. So "any zero" raises false alarms. Tested against all three
runs (38431 was only partly captured):

| rule | 38430 | 38431 | 38432 |
| --- | --- | --- | --- |
| zeros ≥ half the seeds | correct | misses urukul5_ch2, ch3 | correct |
| any zero | correct | correct | flags 2 healthy channels |
| **zeros ≥ 2** | correct | correct | correct |

Rules that ignore zeros near the delay-line edge also fit, but only if the
"edge" is drawn in exactly the right place, so they were rejected.

The tool now reports three tiers: **0 zeros = healthy, 1 = inconclusive
(re-run), 2+ = MARGINAL.** The faulty channels showed 2–6 zeros on every
run; healthy channels have shown at most one. That gap is thin — a faulty
channel dipped to exactly two once — so never condemn or clear a card on a
single run.

---

## How to read the numbers

### They are time, not frequency

**`sync_delay_seed`** is a position on the AD9910's SYNC_IN delay line:
0–31 taps. Nominal tap size is ~75 ps; measured here (see below) it is
**~85 ps**, so the full line spans about 2.3 ns. `15` is not a "correct"
value, only ARTIQ's default midpoint.

**`io_update_delay`** (the `0`/`3` that `sinara_tester` prints) is a
different quantity — the alignment between `IO_UPDATE` and the internal
`SYNC_CLK`, in whole SYSCLK cycles. SYNC_CLK is SYSCLK/4, so it only ever
takes 0–3.

**`window`** is the validation width in the same tap units: how wide a
region around the sampling point stayed free of the chip's `SMP_ERR` flag.
**This is the setup/hold margin, and it is the number that matters.**

### Judge by window. Ignore delay spread.

Every channel — healthy ones included — returns delays in about three
clusters, roughly 12 taps apart. That is **not** noise and **not** a fault.

The chip samples SYNC_IN using SYSCLK. Increasing the delay slides the
SYNC_IN edge relative to the SYSCLK sampling edge; slide it a full SYSCLK
period and you are back to the same relative position. So valid and invalid
zones alternate with a one-SYSCLK-period pitch, and `tune_sync_delay`
returns whichever repeat is nearest its seed. Every repeat catches the same
SYNC_IN edge — SYNC_IN is 62.5 MHz (16 ns) while the whole delay line spans
only 2.3 ns — so they are all equally valid.

This is confirmed rather than assumed. Taking the 16 healthy channels and
measuring the gap between landing clusters:

```
n = 32 gaps, mean 11.8 taps, median 11.7, min 10.5, max 13.0
one SYSCLK period at 1 GHz = 1.000 ns
=> implied tap size = 1000 ps / 11.8 = 85 ps
```

which agrees with the AD9910's ~75 ps typical (it varies with process,
voltage and temperature). The spacing independently measures the delay-line
step, which is what makes the periodic interpretation solid.

The decisive evidence that this is real periodicity and not just "the search
found something near where it started": subtracting the seed from each
result gives swings of ±6 with **adjacent seeds running in opposite
directions**, and different seeds converging on identical values. On
`urukul1_ch0`, seed 7 walked *down* to 4 while seed 11 walked *up* to 15 —
four taps apart, fleeing in opposite directions, because there is a dead
zone between them. Seeds 11, 15 and 19 all landed on exactly 15.

> **Consequence:** the scatter `artiq_sinara_tester` reports (e.g. urukul4 as
> 19 / 9 / 15 / 10) is just different stored seeds landing in different
> repeats. It is **not** evidence of a fault. Only the window is.

---

## Why the ritual works, and why it comes back

`tune_sync_delay` scans only **±6 taps around the stored seed**
(`ad9910.py`, `search_span = 13`, carrying a FIXME that cites
[sinara-hw/Urukul#16](https://github.com/sinara-hw/Urukul/issues/16) about
Kasli SYNC_IN jitter). The seed comes from EEPROM, per channel, addressed in
the device_db as e.g. `"sync_delay_seed": "eeprom_urukul4:64"`.

`artiq_sinara_tester.py:301` writes each channel's freshly tuned result back
to that EEPROM **as it goes**, and the loop **aborts on the first failure**.
So a failing run still refreshes every channel *before* the failure, and
those become the seeds for the next run. Each run therefore reaches one or
two channels further than the last.

That is the "run it several times and it eventually passes" behaviour. It is
deterministic seed propagation, not luck. `artiq_coremgmt log` contributes
nothing — it only reads the device log.

Note the margin arithmetic: window pitch ~11.8 taps, search reach ±6. From
the worst starting point — dead centre between two windows, ~5.9 taps from
either — the search *only just* reaches a valid zone. There is essentially
zero margin in the search itself, which is why a seed that drifts slightly
returns `no valid window/delay` outright, and why a cold boot at a different
temperature can undo a calibration that worked yesterday.

---

## Why it cannot be worked around in software

**It is not cosmetic — experiments hit it too.** `AD9910.init()` itself calls
`tune_sync_delay` at `ad9910.py:509`:

```python
if self.sync_data.sync_delay_seed >= 0 and not blind:
    self.tune_sync_delay(self.sync_data.sync_delay_seed)
```

and `BaseExperiment_master_satellite.py:1335` calls `init()` on every DDS
channel on every run (standalone does the same at `DeviceAliases.py:79`).

**SYNC cannot simply be disabled.** Setting `sync_delay_seed: -1` would skip
the whole path, but the code sets deterministic phase modes in 21 places —
`dds_microwaves` and `dds_MW_RF` with `PHASE_MODE_TRACKING` /
`PHASE_MODE_ABSOLUTE` in the parity and microwave-map sequences. Those derive
phase from the RTIO timestamp, which assumes a fixed relationship to the
DDS's internal divider — exactly what SYNC establishes. And the failing card
carries the microwave DDS.

**Freezing the seeds is not enough either.** Writing known-good values as
literals in the device_db would remove the seed-propagation fragility, but it
cannot manufacture margin that is not there. A channel with window 0 will
still fail intermittently.

---

## Wrong turns (do not re-test these)

1. **"It crosses DRTIO, so it needs more time."** Unrelated, and separately
   disproved — see the RTIO underflow work of the same date. `spi_urukul1`
   (destination 0) and `spi_urukul4` (destination 1) need identical lead
   times. This Urukul problem also predates the master-satellite migration
   entirely; it appeared in standalone.
2. **"MMCX clock distribution, or supply droop along the init sequence."**
   Both plausible from the symptom, both killed by `urukul3` sitting in the
   same crate on the same distribution and measuring perfectly.
3. **"The scatter within urukul4 proves the card is marginal."** Wrong — see
   the periodicity section. Scatter is expected on healthy cards too.
4. **"Then it is only bad EEPROM seeds; the cards are fine."** Also wrong.
   Seeds fully explain the *ritual*, but not why four channels have zero
   margin while twenty have margin at every seed.
5. **Physically swapping cards between crates to separate card from slot.**
   Correct in principle, rejected as too invasive here — it needs the rack
   out and the ribbon cables off the back. `urukul3` acting as an in-crate
   control gave the same answer for free.

---

## Next steps — done 2026-09-15

1. ~~Replace `urukul5`.~~ Done, and `urukul4` was replaced too.
2. ~~Re-run `measure_urukul_sync_windows`.~~ Done (RID 38432): no channel
   marginal under the corrected rule.
3. ~~Run `artiq_sinara_tester` once.~~ Done: it passed on the first run and
   wrote fresh EEPROM calibration.
4. Still to confirm: that the boot-time ritual stays gone over the next few
   power cycles. So far there has been only one boot with the new cards.

---

## Reference

| What | Where |
| --- | --- |
| Sweep tool | `measure_urukul_sync_windows.py` (repo top level) |
| Search span limit | `artiq/coredevice/ad9910.py`, `search_span = 13` |
| Error raised | `artiq/coredevice/ad9910.py:1038` |
| `init()` calls the tuner | `artiq/coredevice/ad9910.py:509` |
| Tester writes EEPROM, aborts on failure | `artiq/frontend/artiq_sinara_tester.py:299-302` |
| Experiments call `init()` per channel | `utilities/BaseExperiment_master_satellite.py:1335` |
| Seed source in device_db | `"sync_delay_seed": "eeprom_urukulN:64"` |
