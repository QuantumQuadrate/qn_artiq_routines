# Chopped readout: motivation, benefits, and requirements

**Started:** 2026-10-12 **Status:** **Idea / partially built** — `record_chopped_RO`
exists in `subroutines/experiment_functions.py`, and `use_chopped_readout` is
`False` on both nodes.
**Scope:** single-node first; the two-node version depends on decisions recorded
in [D1_optical_pumping_merge_ideas.md](D1_optical_pumping_merge_ideas.md).

Sections 1–7 are the physics case and the requirements. Section 8 is what the
repository actually contains today, measured, and section 9 is what blocks it.

---

## 1. Purpose

The purpose of chopped readout is to perform fluorescence readout while reducing
the harmful effects of the optical dipole trap (FORT) during near-resonant
photon scattering.

In a normal continuous readout, the atom remains trapped while readout or
cooling light is applied. This can be problematic because the FORT causes
ground- and excited-state light shifts. These shifts make the readout transition
position-dependent and can break the clean cycling character of the transition.
The result can be reduced scattering efficiency, worse cooling, extra heating,
Raman depumping, and atom loss.

The chopped readout idea is to avoid this by alternating the trap and readout
beams out of phase:

```
FORT/trap:      ON   OFF   ON   OFF   ON   OFF
Readout light:  OFF  ON    OFF  ON    OFF  ON
SPCM gate:      OFF  ON    OFF  ON    OFF  ON
```

During the readout-light-on period the trap is off, so the atom scatters photons
without experiencing the FORT-induced excited-state Stark shift. During the
trap-on period the readout light is off, so the atom is recaptured before it
escapes.

## 2. Literature motivation

This idea is used in Saffman-lab neutral-atom readout work. In Kwon, Ebert,
Walker and Saffman, *Parallel Low-Loss Measurement of Multiple Atomic Qubits*,
the probe beams were chopped out of phase with the trap using a 50% duty cycle
at 1.25 MHz. The paper explicitly states that the probe and ODT lights were
chopped out of phase so that the excited state is not populated while the ODT is
on. Their reported readout achieved 97% uncorrected state detection fidelity,
98.7% corrected fidelity, and <2% atom loss after correction for background
losses.

A similar idea appears in Graham et al., *Midcircuit Measurements on a
Single-Species Neutral Alkali Atom Quantum Processor*. Their initial occupancy
measurement used readout light chopped out of phase with the trapping light with
a 50% duty cycle and a chopping frequency of 1.42 MHz. They state that chopping
the trap prevented atoms from experiencing trap-induced Stark shifts while
illuminated by the readout beams.

## 3. Main benefit

The main benefit is not simply speed. The deeper benefit is:

```
reduced trap-induced Stark shift during fluorescence
-> cleaner readout transition
-> less depumping / heating / loss
-> higher retention for a given photon collection
```

It can enable faster usable readout, but only indirectly. Since the readout
light is on only for part of the time, the average scattering rate is reduced by
the duty cycle. For example, with a 50% duty cycle and the same instantaneous
readout intensity, the average scattering rate is roughly halved.

However, because chopped readout avoids the worst trap-induced light-shift
problem, one may be able to use stronger or more resonant readout light without
losing the atom as quickly. In that practical sense, chopped readout can make
the overall readout more efficient.

## 4. Requirements

### 4.1 Chopping must be fast compared to atomic motion

The FORT-off time must be short enough that the atom does not move
significantly before the trap turns back on. A useful rule is:

```
f_chop >> trap frequency
```

For our system the radial trap frequency is roughly 90–105 kHz, corresponding to
a trap period of about 10 µs. Therefore the FORT-off time should be much shorter
than 10 µs. A practical target is:

```
chopping frequency:            0.5-1.5 MHz
FORT-off / readout-on time:    sub-us to ~1 us
```

This is consistent with the Saffman-lab examples at 1.25 MHz and 1.42 MHz.

### 4.2 Trap and readout light must be out of phase

The trap and readout light should not significantly overlap. If readout photons
are scattered while the FORT is on, the main benefit is lost. The intended
timing is:

```
FORT on  -> readout off
FORT off -> readout on
```

The actual optical pulses should be checked on a photodiode, not only inferred
from TTL or RF commands.

### 4.3 AOM/DDS switching must be fast enough

At MHz chopping each half-cycle is only a few hundred nanoseconds. The optical
rise/fall time must be short compared to this window. If the AOM response is too
slow, the effective duty cycle becomes poorly defined and the trap/readout
overlap may become non-negligible.

What should be checked experimentally:

```
FORT optical pulse
readout/cooling optical pulse
relative phase between them
SPCM gate timing
```

### 4.4 SPCM gate should be synchronized to the readout window

For our SPCM-based implementation the cleanest version is to open the SPCM gate
only during the readout-light-on windows. This is slightly different from
camera-based readout, where the camera integrates over the full exposure. The
advantage is that SPCM counts are collected only during the useful fluorescence
windows:

```
readout light on  -> SPCM gate open
readout light off -> SPCM gate closed
```

This suppresses background counts and prevents counting photons during trap-on
or transient periods.

### 4.5 Duty cycle must be optimized

A 50% duty cycle is a natural starting point and is used in the Saffman-lab
examples. But it is not automatically optimal. Increasing readout duty cycle
gives more photons but can increase heating and loss. Increasing trap duty cycle
improves confinement but reduces photon collection rate. The tradeoff is:

```
more readout duty -> more signal, more heating/loss
more trap duty    -> better retention, fewer photons
```

## 5. Difference from alternating two-node readout

Chopped readout and alternating two-node readout are different schemes.

Chopped readout means:

```
FORT/readout/SPCM gate are rapidly modulated within one node's readout.
```

Alternating two-node readout means:

```
Alice and Bob are read out in separate time windows:
Alice RO -> Bob RO -> Alice RO -> Bob RO -> ...
```

The alternating readout separates the nodes in time so that each node has an
individual readout window. The chopped readout separates trap-on and readout-on
periods within a single readout process. Therefore the names should stay
distinct:

```
chopped_readout              -> fast FORT/readout modulation
two_node_alternating_shot    -> Alice/Bob readout windows alternate
```

## 6. Practical recommendation for our setup

For our Rb single-atom node a reasonable starting point is:

```
chopping frequency:        0.5-1.5 MHz
duty cycle:                ~50%
FORT off / readout on:     ~0.3-1 us
SPCM gate:                 aligned to readout-on window
total readout duration:    a few ms to ~10 ms, depending on count rate
```

The key experimental checks should be:

1. Measure the actual optical FORT and readout pulse shapes on a fast photodiode.
2. Verify that the readout pulse occurs when the FORT is off.
3. Verify that the SPCM gate overlaps only the readout-on window.
4. Compare retention between normal continuous readout and chopped readout.
5. Compare total collected photons and atom survival for different chopping
   frequencies and duty cycles.

## 7. Bottom line

Chopped readout is best understood as a low-loss fluorescence readout method,
not simply a faster readout method. Its purpose is:

```
scatter readout photons while the FORT-induced Stark shift is absent,
then restore the FORT before the atom escapes.
```

If implemented correctly it can reduce heating and atom loss, improve retention,
and potentially allow faster usable fluorescence readout for the same survival
probability.

---

## 8. As-built state in this repository (measured 2026-10-12)

### 8.1 It already exists, and it is already DMA

`record_chopped_RO` at `subroutines/experiment_functions.py:1651`, played back by
`first_shot` (`:1788`) and `second_shot` (`:1875`) behind
`use_chopped_readout`, which is `False` on both nodes today. The trace, per
cycle:

```
at_mu(t0 + 100 ns)                       dds_FORT.sw.off()
at_mu(t0 + 100 ns + chop_RO_FORT_off)    dds_FORT.sw.on()
at_mu(t0 + chop_RO_pulse_offset)         dds_cooling_DP.sw.on()
at_mu(t0 + ... + chop_RO_pulse_length)   dds_cooling_DP.sw.off()
at_mu(t0 + chop_RO_gate_offset)          ttl0/ttl1 _set_sensitivity(1)
at_mu(t0 + ... + chop_RO_gate_length)    ttl0/ttl1 _set_sensitivity(0)
delay(0.5 us)
```

with the Node1 values `chop_RO_FORT_off = chop_RO_pulse_length =
chop_RO_gate_length = 500 ns`, `chop_RO_pulse_offset = 100 ns`,
`chop_RO_gate_offset = 900 ns`, and `range(4000)` cycles.

**DMA is not a choice here, it is the only option.** 8 events per cycle over a
1.9 µs period is 32,000 events in 7.6 ms; submitting those live costs
32,000 × 713 ns = **22.8 ms**, i.e. **300% duty**. (713 ns/event is measured on
this hardware — RID 38429, `tests/measure_drtio_event_cost.py`.) For contrast
the *unchopped* alternating readout runs at 0.66% duty and gains nothing from
DMA; see the DMA analysis in the conversation of 2026-10-12.

### 8.2 The calibration itself says the cooling AOM is slow

The readout light is commanded on at `t0 + 100 ns` and off at `t0 + 600 ns`, but
the SPCM gate opens at `t0 + 900 ns` — **300 ns after the light was commanded
off**. The gate is a direct TTL input with no Urukul delay. The Urukul `sw` path
contributes ~150 ns, so if that were the whole story the light would span
250–750 ns and a 900 ns gate would collect nothing.

That the calibration landed at 900 ns implies roughly **800 ns of delay in the
cooling-light path**, about five times the DDS contribution. The remainder is
AOM acoustic transit plus optics.

**This gives a falsifiable test for the planned optics upgrade** (focus the beam
at the DP AOM with a pair of lenses, per section 4.3): after focusing,
re-calibrate `chop_RO_gate_offset`. It should fall substantially. If it does
not, the delay is somewhere other than the AOM and the lenses bought no speed.

Rise time scales as beam diameter over acoustic velocity, τ ≈ d / v_a. For fused
silica (v_a ≈ 4–6 mm/µs) a 1 mm beam gives ~240 ns and a 100 µm beam ~24 ns.
There is a ceiling: focusing raises divergence δθ ≈ λ/(π w₀), and at w₀ = 50 µm
/ 780 nm that is ~5 mrad, comparable to a typical AOM's Bragg acceptance, so
diffraction efficiency starts falling. **100–150 µm is usually the sweet spot.**
Check the AOM's acoustic material first — TeO₂ slow-shear is ~0.6 mm/µs, seven
times worse than fused silica.

### 8.3 The slow AOM is what caps the duty cycle at 26%, not 50%

This is the synthesis worth keeping. The trace's period is set by the gate:
`chop_RO_gate_offset (900) + chop_RO_gate_length (500) + delay (500) = 1.9 µs`,
so 526 kHz with a 500 ns FORT-off — a **26% duty cycle**, not the ~50% of
section 4.5 and the literature.

Reaching 50% at 1 MHz needs a ~1 µs period, which the 900 ns gate offset alone
forbids. So:

```
slow DP AOM  ->  800 ns path delay  ->  900 ns gate offset
             ->  period >= 1.4 us   ->  duty capped near 26%
```

**The optics upgrade is what unlocks the literature's 50% at 1–1.5 MHz.** It is
not an incremental improvement to an already-working scheme; it is the thing
standing between the current 26% and the regime the papers operate in.

### 8.4 The FORT needs the same treatment — but one TTL can serve both

In the trace the FORT and the readout are **exactly complementary and co-timed**:
FORT off and cooling on are both commanded at `t0 + 100 ns`, both released at
`t0 + 600 ns`. Both go through Urukul `sw` lines, so the ~150 ns DDS delay is
common-mode and their *relative* alignment survives. What does not survive is
the AOM, so the FORT path needs focusing too — harder, because it is
high-power 852 nm and damage threshold limits how tightly it can be focused.

Because the two edges are identical and opposite, **one TTL can drive both
external RF switches, with the FORT's wired inverted.** That matters for the TTL
budget in section 9.

The trade-off to accept consciously: it welds `chop_RO_FORT_off` to
`chop_RO_pulse_length`. They are both 500 ns today so nothing is lost now, but
if the FORT ever needs to stay off slightly longer than the light pulse — to let
the light shift settle before photons are collected — the second line is needed
back.

### 8.5 Two traps in the existing path, to fix before reviving it

1. **It reads only 2 of 4 SPCMs, and averages instead of summing.**
   `AllSPCMs_RO1 = int((SPCM0_RO1 + SPCM1_RO1) / 2)` (`:1797`) against the
   unchopped branch's four-detector **sum** (`:1821`). That is roughly **4×
   smaller numbers**, so `single_atom_threshold` and `two_atom_threshold` do
   **not** transfer between chopped and unchopped readout. Same class of bug as
   the atom-check thresholds: a threshold that looks calibrated but is scaled
   wrongly reads as a physics problem.
2. **It uses raw `self.ttl0` / `self.ttl1` and the private
   `_set_sensitivity()`** — the `TTLInOut` path, not the edge counter. The
   two-node AST guard forbids raw physical names for exactly the reason that
   bites here: in two-node mode `self.ttl0` resolves, silently, to the master.
   A two-node chop trace needs the `_NodeX` counters and all four SPCMs.

Note also that `_set_sensitivity` + `count()` generates **one input event per
photon**, unlike `EdgeCounter.gate_rising`, which emits one event per gate. At
these rates the input FIFO is fine, but the two paths are not interchangeable.

### 8.6 Sizing a two-node chopped alternating readout

Per cycle: 2 FORT (both nodes) + 2 cooling (one node lit at a time) + 4 counter
gates ≈ **12–16 events**. For a 10 ms per-node exposure at a 1.9 µs period that
is ~5,300 cycles → **~70–85k events, roughly 1.5–2 MB of SDRAM**. Workable.

Two constraints on the pattern:

* Several timestamp groups land at exactly **4 events**, which is the whole
  four-per-timestamp budget with no headroom. The chop edges would need
  staggering by a few nanoseconds.
* The trace bakes in every geometry knob at record time. `chop_RO_*`,
  `t_alternating_RO_*` and `t_Node2_rtio_offset_mu` are all dashboard-scannable,
  so the trace must be **re-recorded per scan point** or a scan silently runs
  the first point's timing. And the optics block must stay **outside** the
  trace, or `stabilizer_FORT_NodeX.amplitudes[1]` is frozen at record time —
  the same failure mode as the FORT-commanded-to-0.0 incident.

## 9. What blocks it

**The TTL budget, and it is shared with the D1 plan.** The gateware configures
only `ttl4-7` and `ttl12-15` as outputs, and all eight are aliased on both
nodes. `ttl2` and `ttl10` are unaliased and referenced nowhere in the sequence
code, but they sit on the input-configured banks, so using one as an output
needs a gateware/DIO direction change rather than a device_db alias.

The one output that frees up is **`ttl15`** (`ttl_node1_output2` /
`ttl_node2_output1`), driven only by the standalone two-node handshake at
`experiment_functions.py:14270`. Combined with the D1 note:

| TTL | role | depends on |
|---|---|---|
| ttl12 | `ttl_D1_pumping_switch`, both nodes | retiring the standalone handshake |
| **ttl15** | **chop gate: cooling DP + FORT (one line, FORT inverted)** | the same |
| ttl11 | spare | — |

So **both** this project and the D1 merge hang on one decision: whether two-node
*standalone* synchronization is retired. If it is, two output TTLs free up on
both nodes and both unblock together. If not, neither has a symmetric home, and
the one-TTL-inverted trick of section 8.4 becomes necessary rather than merely
convenient.
