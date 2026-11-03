# Two-atom threshold and Node1's readout fluorescence

**Started:** 2026-10-27 **Status:** **Plan only. Nothing has been changed**: no
dataset and no file. **Reviewed** 2026-10-27 against the h5 files, the code and
a fresh copy of the dataset DB; corrections are folded in and listed in §13.
Deferred by the user (2026-10-27): "probably i can try to lower the loading
FORT power so that the fluctuation is smaller. - might not do it now; make a
plan and keep it in upgrade_notes".

**Note on line numbers.** The same day, the shared two-node settings were being
reorganised (the 13 `TWO_NODE_AGREEING_SCALARS` such as `t_SPCM_first_shot`
becoming plain globals, and the HOM loop's per-node pairs becoming
`*_two_nodes` globals read under their bare names). Where this note says
"agreeing scalar" or cites the max-of-two logic at
`experiment_functions_two_nodes.py:2628-2630`, read it against the code as it
stands when you pick this up.
**Scope:** the joint two-node thresholds (`two_atom_threshold`,
`two_atom_threshold_for_loading`), the per-node thresholds they are built
from, Node1's first-shot fluorescence, and a test of a lower Node1 loading
FORT. Out of scope: HOM photon timing and photons per attempt, which have
their own follow-ups.

Units: every number below is a count per **10 ms** on the **sum of all four
SPCMs**. Runs at 20 ms are scaled to 10 ms. Thresholds are rates in c/s, so
42000 c/s means 420 counts in 10 ms. Each node's light reaches all four
detectors through the BSA beamsplitter, so a detector pair never identifies a
node. Measured numbers cite their RIDs. *(inferred)* marks reasoning that no
run measured directly. Line numbers refer to the working tree of 2026-10-27
(branch `20261020_GRIN_new_wiring`, which has uncommitted edits in
`experiment_functions_two_nodes.py`).

---

## 1. Bottom line

1. **The HOM in-loop atom check fails at its first check because the
   threshold is above the "both atoms" peak.** In each of 39553 and 39554,
   499 of 500 measurements stopped at cycle 49; one passed once (427 and 423
   counts) and stopped at cycle 99 (`n_excitation_cycles`,
   `AllSPCMs_RO_atom_check`). `two_atom_threshold` = 42000
   (420 counts). The joint first-shot peaks are now ~150 (no atom), ~230
   (Node1 only), ~280 (Node2 only) and **~360 (both)** (39553/39554). Before
   2026-10-21 the "both" peak was 459-496 (38705-38724) and 481 (39108), and there 420
   worked.
2. **The "both" peak fell because Node1's single-atom signal halved at RID
   39195 → 39196 (2026-10-21, ~12:00).** That RID boundary is also where
   `p_cooling_DP_MOT_Node1` went from −2 to −5 dBm. The readout and PGC
   amplitudes are *derived* from the MOT amplitude
   (`utilities/BaseExperiment_master_satellite.py:1506-1513`). Net counts went
   140 → 57 (science FORT 0.30). From 2026-10-24 they are 75, after the
   science set point went to 0.20.
3. **"Fluctuating" breaks down into three parts.**
   - The step above is the largest change.
   - The one-atom peak is close to Poisson shot to shot (Fano 1.06-1.11 on
     2026-10-26).
   - There is a slow drift of 2-7 % of the net signal that moves together with
     the background. In relative terms it is the same size as Node2's.

   The same absolute noise now matters twice as much because Node1's net
   signal is half what it was.
4. **The loading FORT is not in the single-node readout path.** PGC and the
   readout both run at the science depth. A lower loading FORT should not undo
   the step *(inferred)*. Its only route into the single-node readout is the
   sudden loading→science drop before PGC, and the near-Poisson Fano leaves
   little shot-to-shot excess for it to remove *(inferred)*. It does enter
   three other places:
   - the loading-check background (FORT-induced counts, roughly 200 of Node1's
     ~310 at 1.16 V *(inferred)*);
   - the light shift on the atom's MOT-light signal;
   - two-node PGC, which runs at the loading depth.

   It is worth testing for those (§6).
5. **The joint threshold can be derived from the two single-node peaks.**
   - "both" = A1 + A2;
   - "brighter single" = A1 + A2 − S_min.

   This predicted 39553/39554 to within 3 %. Today's value is **≈ 32000 c/s**
   (§7). The options for the user's decision are in §8.

## 2. Measured history: single-node first shot

Each fit is a two-Gaussian valley on RO1 ∪ RO2, with medians on either side.
Ranges are per-day medians unless marked "runs". Single runs spread wider:
for example Node1 net 125 (38445), 133 (39062) and 176 (39195), and Node2 net
149 (39212).

| Node, period (RIDs) | t_SPCM_first_shot | science FORT (APD V) | p_cooling_DP_MOT | background | one atom | net S | d′ |
|---|---|---|---|---|---|---|---|
| Node1 09-15 → 10-21 11:00 (38445-39195) | 10 ms | 0.30 | −2 dBm | 80-91 | 217-254 | 133-162 | 6.4-7.3 |
| Node1 10-21 12:00-14:00 (39196-39217) | 10 ms | 0.30 | **−5** | 76-93 (runs) | 137-178 (runs) | **59-87** (runs) | — |
| Node1 10-21 14:00 → 10-24 13:00 (39218-39342) | **20 ms** | 0.30 | −5 | 78-82 | 134-145 | 57-63 | 4.7-5.1 at 20 ms |
| Node1 10-24 02:00 interleaved science scan (39334-39339) | 20 ms | 0.13 / 0.20 / 0.30 | −5 | 48 / 65-66 / 76-82 | 139 / 138-142 / 134-136 | 91 / 73-75 / 54-57 | — |
| Node1 10-24 13:00 → 10-26 (39343-39499; 39353 ran at 0.30) | 10 ms | **0.20** | −5 | 59-62 | 134-145 | 75-83 | 4.9-5.4 |
| Node2 09-16 → 10-27 (38452-39567) | 10 ms | 0.13 | −4 | 74-110 | 206-244 | 126-147 | 6.1-7.6 |

An archive scan of all 362 `GeneralVariableScan` h5 files from 09-15 to 10-27
confirms the `p_cooling_DP_MOT` column: Node1's only change is −2 → −5 at
39196. The other readout knobs did not change in this window on either node:

| Node | p_cooling_DP_RO | f_cooling_DP_RO | p_cooling_DP_PGC | PGC_and_RO_with_on_chip_beams | PGC coils | loading APD set point | FORT_APD_monitor |
|---|---|---|---|---|---|---|---|
| Node1 | 0.5 | 120.37 MHz | 0.7 | True | AZ −0.3, AX 0, AY −0.3 | 1.16 | 1.14-1.18 |
| Node2 | 0.8 | 130 MHz | 0.8 | True | AZ −0.205, AX −0.3, AY −0.9 | 0.38 | 0.374-0.388 |

`FORT_APD_monitor` is the APD reading at the loading set point after
feedback, recorded once per scan point. It shows the loading FORT power was
stable to about ±2 %. `FORT_MM_monitor_Node1` read 0.46-0.49 and has read
0.45 since 39482 (2026-10-26 14:00, after the 852 nm waveplate scans).

The single-node thresholds followed the signal:

- `single_atom_threshold_Node1`: 14000 → 17000 (39043) → 15000 / 12000 /
  11000 (10-21 14:00) → 12000 → 11000 / 10000 (10-24).
- `t_SPCM_first_shot_Node1`: 10 → 20 ms at 39218, back to 10 ms at 39343
  (39340 and 39342 still ran at 20 ms and science 0.30).

### Noise and drift in long runs (100-measurement blocks)

| Node, RIDs | Fano within block | excess block-to-block drift (sd, counts) | drift / net | corr(atom, bg) |
|---|---|---|---|---|
| Node1 10-16, parity_9 (39062-39066) | 1.08-1.15 | 3.0-4.0 | 0.02-0.03 | −0.6 … +0.4 |
| Node1 10-22/23, 20 ms (39256, 39283-39316) | 1.24-1.33 | 7-10 per 20 ms | 0.06-0.08 | 0.76-0.96 |
| Node1 10-25/26, parity_11 (39411-39499) | 1.06-1.11 | 1.7-5.1 | 0.02-0.07 | 0.20-0.73 |
| Node2 10-13/14, parity_11 (38830-38902) | 1.13-1.29 | 4.8-7.8 | 0.03-0.06 | −0.06-0.33 |
| Node2 10-26/27, parity_11 (39536-39548) | 1.14-1.17 | 3.4-7.3 | 0.03-0.06 | 0.46-0.72 |

When the drift is large, the atom peak and the background move together. That
points to readout-light intensity or collection, not atom temperature
*(inferred)*.

On 10-23 (20 ms) the Node1 block medians ran 266-345 raw counts within single
runs (39261-39294), with the background following.

### Loading check (MOT light, loading FORT, 10 ms)

| Node, period | MOT-phase background (median / p99, first tries) | loaded median | threshold_for_loading |
|---|---|---|---|
| Node1 09-16 → 10-13 | 296-340 / 360-403 | 469-528 | 440-450 |
| Node1 10-21 → 10-26 | 289-399 / 340-453 | 454-591 | 400-500 |
| Node1 10-25 21:00 episode (39367-39378) | **375-420 / 437-453** | 461-590 | 450 → 500 |
| Node2 10-26/27 (39536-39548) | 239-250 / 280-299 | 499-507 | 430 |

In the 10-25 21:00 episode, 66 % of 39367's first shots were empty. The
likely cause is that the background reached the threshold and the loader
passed on it *(inferred: the first-try p99 was 445 against 450, and a load
takes many tries)*. The fix was
`single_atom_threshold_for_loading_Node1` 45000 → 50000 (39371).

Node1's p99 also has excursions to 458-497 in 39411, 39412, 39431, 39468 and
39469, against a threshold of 500. Node1's MOT-phase background runs about 70
counts above Node2's and moves more from run to run.

## 3. Measured history: joint (two-node) first shot

The fit is a constrained four-peak model: none, Node1 only, Node2 only, both.
The means are B, B+S1, B+S2 and B+S1+S2; each width is √(F·mean).

The node labels come from **brightness, not from detectors**. The fit returns
a dim and a bright single-atom step, and "Node1" is assigned to the dim one
because the single-node runs have S1 < S2 (75 vs 130). Today the two single
peaks are only ~2.3 σ apart, so no single shot can be assigned to a node; the
population fractions are model estimates *(inferred)*. On 10-08 (S 163/154)
the two singles are degenerate and their split means nothing.

| RIDs | date | B | S dim / bright | single / single / both peaks | F | both fraction | thr (counts) | at thr: both→"lost" / brighter single→"both" |
|---|---|---|---|---|---|---|---|---|
| 38705-38724 Two_nodes_atom_loading | 10-08 | 172-189 | 135-175 | 316-349 / both 459-496 | 1.5-2.1 | 0.56-0.84 | 390-400 | ≤1 % / ≤4 % |
| 39108 Two_nodes_atom_loading | 10-18 | 169 | 139 / 173 | 308 / 342 / 481 | 1.8 | 0.65 | 420 | 2 % / 0 % |
| 39553 HOM, Node1 excited | 10-27 | 154 | 83 / 126 | 236 / 279 / 362 | 1.26 | 0.71 | 420 | **100 %** / 0 % |
| 39554 HOM, Node2 excited | 10-27 | 146 | 79 / 133 | 225 / 279 / 358 | 1.42 | 0.72 | 420 | **100 %** / 0 % |

In 39553 and 39554, 0.2 % and 0 % of RO1 exceeded 420, and 69 % exceeded 320.
Re-fitted independently for this review: 39553 B 154, S 83/126, F 1.26, both
0.705; 39554 B 146, S 79/133, F 1.42, both 0.722; 39108 B 169, S 139/173,
F 1.83, both 0.652; 38724 both 0.837. All agree with the table.

**Sum rule, measured.** The joint peaks are the single-node peaks added up:

- B = bg1 + bg2;
- both = A1 + A2;
- Node2 only = bg1 + A2;
- Node1 only = bg2 + A1.

| Check | single-node inputs | predicted | fitted | agreement |
|---|---|---|---|---|
| 10-27 | Node1 10-26 (bg 60, A 134); Node2 10-27 (bg 92, A 222) | 152 / 226 / 282 / 356 | 146-154 / 225-236 / 279 / 358-362 | within 3 % |
| 10-18 | 10-16 singles (Node1 83/225, Node2 91/228) | 174 / 316 / 311 / 453 | 39108: 169 / 342 / 308 / 481 | within 8 % |

**Joint loading check** (MOT light, 10 ms, `two_atom_threshold_for_loading` =
100000 c/s, i.e. 1000 counts):

| RIDs | 0 atoms | 1 atom | 2 atoms |
|---|---|---|---|
| 39553 | ~585 | ~830 | loaded median 1054 (truncated at 1000) |
| 39554 | ~530 | ~780 | loaded median 1025 (truncated at 1000) |
| 38724 (10-08) | ~560 | ~800-830 | loaded median 1078 |

- The first-try p99 is 890-910, below the 1000 threshold.
- Even so, the RO1 four-peak fit finds only **71-72 %** of loads with both
  atoms on 10-27, against 0.84 on 10-08 (38720, 38724). This is a model
  estimate (see the label caveat above).
- The "2 atoms" column is the median of loaded checks, which are cut off at
  the 1000 threshold. It is biased upward and is not a fitted peak.
- Whether the 28 % are false loads or atoms lost before RO1 is not separated.
  The PGC-at-loading-depth gap in §6.1 is one candidate *(inferred)*.

## 4. Where each threshold sits now

Live values are from `dataset_db.pyon`, copied 2026-10-27 at about 11:55 and
re-copied at about 12:30 for the review. No threshold, set point, exposure or
readout amplitude differed between the two copies; only Node2 fiber-AOM powers
had moved.

| Variable (live) | counts / 10 ms | Measured peaks now | Verdict |
|---|---|---|---|
| `two_atom_threshold` 42000 (file default 39000) | 420 | joint: 150 / 230 / 280 / **360**, both-peak sd ≈ 22 | **2.7 σ above "both"**: 998/1000 first in-loop checks failed; RO2 pre-check passed 0/500 (39553) and 1/500 (39554); a first-shot confirm would almost never confirm |
| `two_atom_threshold_for_loading` 100000 (file default 89000) | 1000 | joint MOT-phase: 530-585 / 780-830 / loaded median 1025-1054 (cut off at 1000) | between 1 and 2 atoms; OK *(the 1- and 2-atom positions are approximate)* |
| `single_atom_threshold_Node1` 10000 (file 14000) | 100 | bg 60 (p95 ~75), atom 134 (p05 ~112), valley ~85-90 | OK |
| `single_atom_threshold_Node2` 17000 | 170 | bg 92 (p95 ~115-120), atom 222 (p05 ~195), valley ~145-150 | OK |
| `single_atom_threshold_for_loading_Node1` 50000 (file 44000) | 500 | MOT bg ~300-330, p99 340-412, excursions to 497; loaded ~520-570 | OK, thin on bad hours |
| `single_atom_threshold_for_loading_Node2` 40000 | 400 | MOT bg ~240, p99 ~290; loaded ~500 | OK. It was 43000 for 39536-39553 and changed to 40000 by 39567 (10-27, 09:00-11:00) |
| `alternating_shot_threshold_Node1` 20000 | 200 | last measured 10-12 (38773-38778), before the Node1 step: bg 129-138, atom 270-276 | **stale** *(inferred)*: the Node1 window atom peak is now ≈ 195-205 |
| `alternating_shot_threshold_Node2` 21000 | 210 | 10-12: bg 139-160, atom 259-297 | probably OK |

`two_atom_threshold` is read in four places, and each has a different failure
mode (`subroutines/experiment_functions_two_nodes.py`):

| Where | Line | Check | Today's effect at 420 |
|---|---|---|---|
| RO2 "already loaded" pre-check of the joint loader | :1100 (alternating copy :1428) | `RO2/t_SPCM_second_shot > thr` | almost never passes, so nearly every measurement reloads: 500 loads for 500 measurements in 39553 and 499 in 39554 (`n_atom_loaded_per_iteration`, `Atom_loading_time`). That ends atom reuse |
| first-shot confirm (`advance_if_atom_in_first_shot`) | :2339-2341 (`Two_nodes_atom_loading_experiment`), :2699-2701 (`Two_nodes_HOM_experiment`) | `RO1/t_SPCM_first_shot > thr` | would almost never confirm (1 of 1000 first shots was above 420), so usually `max_loading_rounds` (20) load+shot rounds per measurement *(inferred)* |
| HOM in-loop check | :2887 (t_atom_check = max of the nodes' `t_SPCM_recool_and_shot`, :2628-2630) | below threshold means lost | ends at cycle 49 (499/500 measurements per run) |
| setting the threshold too *low* | (all three above) | — | the pre-check and the loop accept a missing atom |

## 5. What changed on Node1

**The step.** Between 39195 (10-21 11:00, net 176) and 39196 (12:00, net 87),
nothing changed in the archive except `p_cooling_DP_MOT_Node1` −2.0 → −5.0
dBm and the feedback outputs.

- In the master-satellite stack `ampl_cooling_DP_RO = ampl_cooling_DP_MOT ×
  p_cooling_DP_RO` and the same holds for PGC
  (`BaseExperiment_master_satellite.py:1506-1513`). The standalone stack does
  the same (`utilities/BaseExperiment.py:291`, `:426`). So the MOT change
  lowered the readout and PGC RF by 3 dB.
- The fiber-AOM feedback put the MOT-beam power back. The archived
  `p_AOM_A1_Node1` was −11.9 dBm at 39195, −12.6 at the start of 39196 and
  −10.2 by 39200 (A2: −11.8, −12.3, −9.6).
- The readout runs at half the MOT amplitude. There, in the small-signal
  limit, a double-pass AOM loses up to 6 dB of light for 3 dB of RF. Take
  off the ~2.4 dB that the fiber AOMs (on at their fed-back RF during readout)
  put back, and the estimate is about −3.6 dB, or ×0.44. That is close to the
  ×0.49 measured at 39196 (net 176 → 87) *(inferred mechanism; the
  coincidence in time is verified)*.
- **The mechanism does not explain what came next.** Net kept falling,
  87 → 74 → 61 → 59 (39196, 39213, 39214, 39217, all 10 ms). By then the
  fiber AOMs had been raised (A1 −10.2 dBm from 39200 on). The simple picture
  predicts that raising them gives back readout light, not that more is lost.
  No readout knob changed in the archive (§2 scan). The cause of this second fall is
  **unverified**; the `p_cooling_DP_RO_Node1` scan in §7.3 is the direct test.

**The science FORT.** From the interleaved 10-24 02:00 runs (39334-39339):

- going 0.30 → 0.20 → 0.13 raised the net from 54-57 to 73-75 to 91;
- the background fell by about 175 counts per 10 ms per APD volt. That slope
  is the FORT-induced (fiber Raman) part.

The 0.20 set point (from 39343) is why Node1 now reads 75 rather than 57.

**Per-shot vs slow.** See the drift table in §2. The in-loop atom checks
reproduce the first-shot level on both nodes:

| RID | Node | first shot | in-loop check |
|---|---|---|---|
| 38894 | Node2 | 215 | 213 |
| 39511 | Node2 | 213 | 212 |
| 39221 | Node1 | 132 (20 ms, scaled) | 133 |

So thresholds derived from first shots transfer to the in-loop checks.

## 6. A lower Node1 loading FORT: why it might help, what could go wrong, how to test it

### 6.1 Where the loading depth enters

**Single node:**
- The loading check runs at the loading depth
  (`subroutines/experiment_functions.py:414`, in
  `load_MOT_and_FORT_until_atom`; the other single-node loaders are not
  re-checked here).
- The FORT is lowered to science *before* PGC (`:490-491`). This is one
  `dds_FORT.set`, a sudden step of 5.8× on Node1, not a ramp.
- `first_shot` reads at science (`:1047`).

**Two nodes:**
- The loader turns the FORTs on at loading
  (`experiment_functions_two_nodes.py:1144-1150`).
- It runs PGC (`:1259`) **without lowering them**: nothing in `:1160-1322`
  sets a FORT amplitude. This gap is already documented at `:1359-1362`.
- `first_shot` sets science (`:863`).

**Depth ratios.** Loading/science is 1.16/0.20 = **5.8** on Node1 and
0.38/0.13 = **2.9** on Node2. These are APD volts with per-node calibrations,
so only the ratios compare.

### 6.2 Physics

1. **Background.** FORT-induced counts are linear in power. At ~175 per APD
   volt, 1.16 V gives about **200** of Node1's ~300-330 MOT-phase background
   *(inferred extrapolation)*. Going to 0.6 V would remove about 100. That
   lowers the floor, its shot noise (√310 → √210) and the false-load risk seen
   on 10-25.
2. **Light shift.** The differential light shift scales with FORT intensity.
   At 5.8× the science depth, a trapped atom's MOT-light signal depends on
   where it sits and how hot it is, which gives a broad and variable "loaded"
   signal *(inferred)*.
3. **Two-node PGC.** In two-node runs, PGC at 5.8× depth cools under a large
   light shift *(inferred)*. There is no visible penalty yet: in 39553/39554
   Node1's S was 79-83 against 75 single-node.
4. **Not expected.** The single-node RO1 step and the slow readout drift
   should not change, because PGC and the readout run at science
   *(inferred)*. One caveat: the sudden drop sets which atoms survive and how
   hot they are going into PGC. That could show up as false loads (metric 4
   below), but the near-Poisson Fano (1.06-1.11) leaves little shot-to-shot
   readout excess for it to remove.

### 6.3 What could go wrong

- **Loading rate.** A shallower trap captures less. Node1's median
  `Atom_loading_time` was 0.26-0.30 s on 10-26 and is the cost to watch.
- **Loading thresholds.** The MOT-phase atom signal will change, so
  `single_atom_threshold_for_loading_Node1` and
  `two_atom_threshold_for_loading` need to be re-derived.
- **Feedback.** The FORT feedback works on a smaller APD signal. Check that
  `p_FORT_loading_history_Node1` converges.
- **Ramp ratio.** A different loading/science ratio changes survival through
  the ramp, and could change two-atom loading (blockade).
- **The override trap.** 39370 and 39371 (10-25) put
  `'set_point_FORT_APD_loading_Node1': 0.5` in **Node2's** override dictionary
  of a Node1 run. It was silently not applied: `FORT_APD_monitor_Node1` read
  1.145 and 1.151, and the archive kept 1.16. The idle node's dictionary is
  not evaluated (`BaseExperiment_master_satellite.py:38-46`). That the entry
  was meant as a Node1 test is *inferred* from its name.
- **Untested path.** No run has yet applied a FORT set-point override
  successfully, so it is **unverified** that an override reaches the FORT
  stabilizer's set points (`subroutines/aom_feedback.py:398` reads them from
  the experiment). Metric 1 below is the check. If it fails, set the dataset
  for the test and restore it afterwards.
- **No scanning.** FORT set points cannot be GVS-scanned, so the test needs
  one run per value.

### 6.4 Test (single node, Node1, after the dataset decisions)

**Pre-test (no atoms).** Gate all four SPCMs for 10 ms with only Node1's FORT
on (no 780 nm light), at 1.16, 0.6 and 0.20 V. This measures the
FORT-induced background directly and checks the ~175 /V slope.

**Main test.**
- Experiment: `atom_loading_2_experiment` through
  `GeneralVariableScan_master_satellite_single_node`, Node1, `dummy_variable`
  × 10, `n_measurements` 100.
- Set the loading point with `override_ExperimentVariables_Node1 =
  "{'set_point_FORT_APD_loading': X}"`.
- Interleave against the reference: 1.16 → 0.8 → 1.16 → 0.58 (Node2's ratio)
  → 1.16 → 0.4.
- Keep `single_atom_threshold_for_loading_Node1` fixed during the comparison.

**Metrics per run** (all from the h5):

1. `FORT_APD_monitor_Node1` ≈ X. This confirms the override took effect.
2. MOT-phase background: median and p99 of the first-try entries of
   `AllSPCMs_atom_check_in_loading`. The dataset stores pairs (loaded,
   first-try); use the first-try entries that differ from their pair.
3. Loaded median minus that background, i.e. the atom signal under MOT light.
4. False-load fraction: the share of RO1 below the readout valley (~88).
5. Median `Atom_loading_time`.
6. One-atom RO1 peak: mean, within-block Fano, excess block drift, d′.
7. Retention, RO2 / RO1.
8. Two-atom tail of RO1 (above A + 0.6 S).

**Decision rule** (the user sets the numbers): adopt the lowest X at which (2)
and (4) improve and (5) grows by less than an agreed fraction. Then re-derive
the two loading thresholds. For the joint one, run one short
`Two_nodes_atom_loading_experiment` and read the 0/1/2-atom MOT-phase peaks.

## 7. A joint threshold that survives Node1 drifting

### 7.1 Derive it from the two single-node runs (recommended procedure)

`two_atom_threshold` [c/s] = (A1 + A2 − S_min / 2) / t_exposure

- A_i is each node's one-atom first-shot peak and S_i = A_i − bg_i.
- Take both from the latest single-node runs at the **same** exposure,
  readout amplitude and science set point as the two-node run.
- Today: (134 + 222 − 75/2) / 0.01 ≈ **31800 c/s** (call it 32000).
- Model error rates at 10 ms, with the fitted F of 39553/39554: 2-4 % of
  "both" read as lost, and 2-3 % of "Node2 only" read as both.
- Neither error is expensive at these rates *(inferred from the code)*:
  - A false "both" runs up to 50 one-atom cycles until the next check.
  - A false "lost" ends that measurement's loop early. It usually does not
    cost a reload, because the next measurement's RO2 pre-check (`:1100`)
    tests again and will normally see both atoms. With a fixed
    `n_measurements` it does cost excitation attempts per run.
  - So the midpoint is a reasonable default; there is no strong case for
    biasing it either way.

**Tolerance.** A Node1 drift of ±5 % of net (±4 counts) moves the "both" peak
against a ±38 margin, which is fine. A step like 10-21 (−70) breaks it. So
any change to Node1's MOT, RO or PGC power or to either science set point must
trigger a re-derivation. The same applies to the RO2 pre-check, which uses the
same threshold on RO2.

No new ExperimentVariable is needed. This can be a small Analysis helper that
reads the newest single-node h5 per node and prints the four recommended
thresholds; the user then sets the datasets.

### 7.2 Fit the run itself

Fit the four-peak model of §3 (prototype `jointfit.py`, below) to RO1 and RO2
of the latest two-node run. It validates §7.1 and also reports the true
both-loaded fraction, which the joint loader cannot see.

### 7.3 Widen the margin instead of chasing it

The current joint d′ ("both" against "Node2 only") is 75 / √(1.3 × (360 +
280)) ≈ **2.6**. Two ways to raise it:

- **Restore Node1's readout amplitude.** Raising `p_cooling_DP_RO_Node1` from
  0.5 to ~0.71 undoes the −3 dB, and PGC 0.7 → ~0.99 likewise. The expected S1
  is ~140 or more, giving d′ ≈ 4.3 and a threshold of about 357 counts
  *(inferred: assumes the pre-10-21 response returns)*. Measure first with a
  Node1 single-node scan of `p_cooling_DP_RO_Node1` (0.5 / 0.6 / 0.71 / 0.85)
  in `atom_loading_2`: net, background and retention. Note that
  `single_atom_threshold_Node1` moves with it.
- **Use 20 ms joint shots and checks.** Set `t_SPCM_first_shot`,
  `t_SPCM_second_shot` (the RO2 pre-check divides by it) and
  `t_SPCM_recool_and_shot` to 20 ms on **both** nodes. The first two are
  agreeing scalars in two-node mode and raise if the nodes differ
  (`BaseExperiment_master_satellite.py:1395-1416`, `:1407-1408`). The HOM
  check takes the longer of the two nodes' `t_SPCM_recool_and_shot`
  (`experiment_functions_two_nodes.py:2628-2630`). This gives d′ ×
  1.41 for 10 ms more per check.

### 7.4 A per-node check: possible only by gating the light in time

The four SPCMs cannot tell the nodes apart, but alternating readout windows
can: each window lights one node. It worked on 10-12 in
`Two_nodes_alternating_retention_experiment`:

| RIDs | Node1 window | Node2 window |
|---|---|---|
| 38773, 38774, 38778 (re-fitted) | bg 129-138, atom 270-276 | bg 139-160, atom 259-297 |

That was 10 × 1 ms windows per node (`t_alternating_RO_pad` 0.1 ms), cleanly
bimodal, d′ 4.8-6.1 across RO1/RO2 and both nodes (median ≈ 5.1).

- **Gains.** Each node's threshold is independent of the other node's
  brightness, and the check says *which* atom was lost.
- **Costs.**
  - About 26 ms instead of 10 ms per check. The 26 ms is the formula
    2·n·(window + 3·pad) from the alternating *loader*'s docstring
    (`experiment_functions_two_nodes.py:1389-1391`), applied here to the
    readout by analogy *(inferred; not timed)*.
  - Node1's d′ would be ≈ 3.5-3.7 today *(inferred)*, which is no better than
    a 20 ms joint check.
  - `alternating_shot_threshold_Node1` must be recalibrated.
  - The HOM loop's check would have to be replaced.
- The alternating *loader* is broken (`:1327-1331`); the alternating readout
  is not.

### 7.5 Comparing each check to the measurement's own first shot

Self-referencing cancels slow drift but costs √2 in noise, giving d′ ≈ 1.9 at
today's S1. Not until S1 is restored.

## 8. Decisions needed (the user's)

1. **`two_atom_threshold` now.** ~32000 c/s (§7.1), or keep 42000 until §7.3
   is done. At 42000 the HOM loop stops at cycle 49 and nearly every
   measurement reloads.
2. **Node1's readout and PGC amplitude after the 10-21 MOT change.** Restore
   them (after the `p_cooling_DP_RO_Node1` scan) or keep 0.5 and accept the
   halved S1.
3. **Joint shots and checks at 10 or 20 ms.**
4. **The loading-FORT test (§6.4).** Which set points to try, and the
   acceptable rise in loading time.
5. **The per-node alternating check (§7.4).** Pursue it, or stay with the
   joint check plus a derived threshold.
6. **`two_atom_threshold_for_loading`.** Leave it at 100000 until after the
   loading-FORT test; that test will move the joint MOT-phase peaks.

## 9. What makes it safe, and what makes it dangerous

**Safe.** The thresholds are datasets. Changing them changes no code and no
pulse order. GUI arguments and overrides still win for the run that sets them.

**Dangerous:**

- **Too low** a `two_atom_threshold` lets the RO2 pre-check and the HOM loop
  accept a missing atom. That costs nothing visible in the data except fewer
  coincidences.
- **Too high** stops the loop and forces reloads (today), or spins 20 load
  rounds per measurement if `advance_if_atom_in_first_shot` is turned on.
- **Exposure and units.** Thresholds are c/s against counts / t. Changing an
  exposure without re-deriving the threshold is silent.
- **Coupled amplitudes.** A Node1 MOT-power change silently moves its readout
  and PGC (`BaseExperiment_master_satellite.py:1506-1513`).

## 10. Prerequisites

- The user's dataset decisions (no HOM scans until then).
- For §7.1: one fresh single-node run per node at the final readout settings,
  with at least 500 first shots each.
- For §6.4: Node1 loading stable for the session, and the overrides in
  `override_ExperimentVariables_Node1`.

## 11. Latent bugs noticed on the way, deliberately not fixed

1. **The failing in-loop check count is never saved. This is by design, not
   an off-by-one, but it hurts tuning.**
   - Checks run at cycles 49, 99, … (`:2827`), and each is stored at index
     `int(cycle / n)` (`:2884`).
   - The save slice `[0:int(cycle / n)]` (`:2965-2966`) keeps only the
     checks that passed. The comment at `:2964` says so: "only the passed
     checks: the loop exits on the failed one".
   - The consequence is that a run whose in-loop threshold is wrong records
     nothing to re-tune it from. 39553 and 39554 saved one value each (427
     and 423), the single passed check in each.
   - A side effect: when a measurement runs all 2000 cycles without a break,
     the slice `[0:39]` also drops the 40th check, which passed *(inferred
     from the code; such measurements are rare)*.
   - The single-node version has the same slice
     (`experiment_functions.py:6366`, `:6424-6425`). Its comment describes a
     different intent (`[0:excitation_cycle + 1]`).
2. **Readout and PGC amplitudes are slaved to the MOT amplitude**
   (`BaseExperiment_master_satellite.py:1506-1513`;
   `utilities/BaseExperiment.py:291`, `:426`). This is the likely route of
   the 10-21 step, but the size of the drop is not explained (§5,
   *inferred*).
3. **The two-node loaders run PGC at the loading FORT depth**
   (`experiment_functions_two_nodes.py:1144-1150`, `:1259`; already noted at
   `:1359-1362`). The single-node loaders lower the FORT first
   (`experiment_functions.py:490-491`).
4. **`alternating_shot_threshold_Node1` is stale** since the Node1 step
   *(inferred)*.
5. **Orphaned globals.** Five `*_two_node` timing globals are declared but
   read nowhere outside the tests: `t_atom_check_time_two_node`,
   `t_SPCM_first_shot_two_node`, `t_SPCM_second_shot_two_node`,
   `t_delay_between_shots_two_node` and `t_MOT_dissipation_two_node`
   (`ExperimentVariables_master_satellite_global.py:143-177`). Two-node runs
   use the agreeing per-node scalars instead. `max_atom_check_tries_two_node`
   (`:179-185`) **is** read (`experiment_functions_two_nodes.py:1159`,
   `:1487`).
6. **File defaults differ from the live datasets.** This only matters to an
   initializer.

   | Variable | file default | live |
   |---|---|---|
   | `two_atom_threshold` | 39000 | 42000 |
   | `single_atom_threshold_Node1` | 14000 | 10000 |
   | `single_atom_threshold_for_loading_Node1` | 44000 | 50000 |
   | `set_point_FORT_APD_science_Node1` | 0.30 | 0.20 |
   | `two_atom_threshold_for_loading` | 89000 | 100000 |
   | `alternating_shot_threshold_for_loading_Node1` | 62000 | 50000 |
   | `alternating_shot_threshold_for_loading_Node2` | 60000 | 45000 |

7. **The idle node's override dictionary is ignored by design**, which is a
   trap: it caught 39370/39371 (§6.3).

## 12. How the numbers were made

Read-only, from the results h5 files (the file without the long suffix;
`archive` holds the parameters) for 2026-09-15 to 2026-10-27, and from a copy
of `dataset_db.pyon`.

- **Single node.**
  - Fit a two-Gaussian EM on RO1 ∪ RO2, normalised to 10 ms, to get the
    valley.
  - Take medians and robust sd on either side, using RO1 only for the atom.
  - Drift is the variance of 100-measurement block means minus the
    within-block sampling variance.
- **Joint.** Maximum likelihood on the constrained four-peak model.
- **Loading checks.** Use the (loaded, first-try) pairs in
  `AllSPCMs_atom_check_in_loading`. A first try that equals its loaded value
  is excluded from the background.
- **Scripts** (scratch, not in the repo): `inventory.py`, `analyze2.py`,
  `drift.py`, `daily.py`, `loadcheck.py`, `jointfit.py`, `hist_params.py`.
  Review scripts: `sn.py` (single-node fits), `tn.py` (four-peak fits),
  `knobs.py` (archive scan), `misc.py`, `ovr.py`, `alt.py`, `range.py`.
- **Note:** `daily.py` fills the `p_cooling_DP_MOT` column from the RID
  (< 39196 → −2) rather than from each archive. The `knobs.py` scan confirms
  that this matches the archives.

## 13. Review log (2026-10-27)

Read-only. Nothing in the repo, the dataset DB or the hardware was touched.

**Confirmed (re-computed independently):**
- Single-node first shots: 38445, 39062, 39195 (net 176), 39196 (87), 39213,
  39214, 39217, 39221, 39334-39339, 39411, 39499 (Node1); 38894, 39212,
  39511, 39548, 39567 (Node2). All are within or next to the table ranges.
  The interleaved science scan reproduces 48/65-66/76-82 background and
  91/73-75/54-57 net.
- Four-peak fits of 39553, 39554, 39108 and 38724, the 420 crossings, and
  the joint loading-check numbers.
- 39367: 66 % of first shots empty. 39370/39371: override in the Node2
  dictionary, archive 1.16, monitor 1.145/1.151.
- In-loop vs first shot (38894, 39511, 39221). Fano 1.08 (39499). Node2
  loading (39548): 239/280, loaded 505.
- The live thresholds (§4), the readout knobs (§2), and every code line
  cited, except those corrected below.

**Corrected:**
- 1 of 500 measurements per HOM run passed its first check; 39554 skipped
  one reload.
- §11.1 is deliberate (`:2964`), not an off-by-one.
- `max_atom_check_tries_two_node` is read. Two file-default differences were
  added.
- The proposed-threshold error rates (2-4 % / 2-3 %) and the cost reasoning.
- Alternating d′ is 4.8-6.1, not 5.6, and the 26 ms is borrowed from the
  loader.
- The second Node1 fall (87 → 59) contradicts the simple AOM mechanism.
- The node labels in the four-peak model are by brightness only.
- Minor RID boundaries (39343, 39370 at 420) and line ranges
  (:1144-1150, :1359-1362, :6424-6425, :143-177).

**Still unverified:**
- The drift table (§2) and the per-day ranges were not all recomputed.
- The 175/V slope extrapolated to 1.16 V.
- Whether a FORT set-point override reaches the stabilizer.
- The cause of the 10-21 step's size and of the second fall.
