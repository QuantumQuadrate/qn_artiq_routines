# Atom-atom ZZ and parity: the `_2` rename

**Started:** 2026-11-14 **Status:** **DEPLOYED 2026-11-14 (evening), uncommitted**, together
with the GVS health-check direct re-optimization (gvs_health_check_direct_optimization.md).
Built in a scratch copy (scratchpad `ren2/`, then `hcopt/w2`); live tests and a live offline
compile of all four entry points after the copy.

User, 2026-11-14: "if you have made new version for this, do not duplicate.
make a new one with a suffix: Two_nodes_atom_atom_ZZ_2_experiment", and
"Those duplicated definitions should have the suffix _2_".

Line numbers below are the live file
`subroutines/experiment_functions_two_nodes.py` before the edit (sha256
`5424b566...`, saved 11-12 23:08). Where the post-edit number differs, it is
given as "now N".

---

## 1. Bottom line

1. The file defined three names twice. Python keeps the **last** definition,
   so since the 11-11 19:41 deploy the earlier versions could not be reached:

   | Name | earlier (dead today) | later (live today) |
   | --- | --- | --- |
   | `Two_nodes_atom_atom_ZZ_experiment` | 3913, `@kernel` | 5921, host wrapper -> `_run_atom_atom_ZZ_with_cycle_log` (5358) |
   | `_run_atom_atom_phase_scan_latest_ZZ` | 4647, `@kernel` | 5958, `@kernel` |
   | `Two_nodes_atom_atom_parity_experiment` | 5136, host wrapper -> kernel 4647 | 6467, host wrapper -> kernel 5958 |

2. After the rename the **later** definitions get the `_2` names and the
   earlier ones are live again under the original names:

   | Name | Definition | Runs |
   | --- | --- | --- |
   | `Two_nodes_atom_atom_ZZ_experiment` | 3913 `@kernel` | the 11-10/11-11 ZZ (what RID 40946 ran) |
   | `Two_nodes_atom_atom_ZZ_2_experiment` | now 5925 wrapper -> `_run_atom_atom_ZZ_with_cycle_log` | today's ZZ (what RID 41177 ran) |
   | `Two_nodes_atom_atom_parity_experiment` | 5136 wrapper -> `_run_atom_atom_phase_scan_latest_ZZ` (4647) | the 11-10 four-point parity (what 40871-40926 ran) |
   | `Two_nodes_atom_atom_parity_2_experiment` | now 6471 wrapper -> `_run_atom_atom_phase_scan_latest_ZZ_2` (now 5962) | today's four-point parity (never run on hardware) |

3. **To keep running what runs today, submit the `_2` names.** A clone or
   driver template that names `Two_nodes_atom_atom_ZZ_experiment` or
   `Two_nodes_atom_atom_parity_experiment` silently gets the older sequence
   once the edit is saved (section 5).
4. No sequence changes. The edit renames three `def` lines and one call, adds
   one name to one applet list, and corrects comments and two log labels.

## 2. What changes

`subroutines/experiment_functions_two_nodes.py` (+4 lines):

* 5921 `def Two_nodes_atom_atom_ZZ_experiment` -> `def Two_nodes_atom_atom_ZZ_2_experiment`.
* 5958 `def _run_atom_atom_phase_scan_latest_ZZ` -> `def _run_atom_atom_phase_scan_latest_ZZ_2`.
* 6467 `def Two_nodes_atom_atom_parity_experiment` -> `def Two_nodes_atom_atom_parity_2_experiment`.
* 6557 the call `_run_atom_atom_phase_scan_latest_ZZ(self)` inside the later
  parity wrapper -> `_run_atom_atom_phase_scan_latest_ZZ_2(self)`.
  **Required:** a host wrapper looks the kernel up by name when it runs, so
  without this edit `parity_2` would silently run the OLD 4647 kernel.
* Comments in the 11-11 block header (5242, 5244-5245; now 5242-5249): the
  block's public names are the `_2` names, the earlier definitions keep the
  originals; "appending this later also overrides the same names" is dropped.
* Docstring 5408: "Use the public host wrapper Two_nodes_atom_atom_ZZ_2_experiment".
* Log labels (string constants in `print_async`, no timing effect):
  5503 "Two_nodes_atom_atom_ZZ_2_experiment: Node2 delay D =" and
  6041 "Two_nodes_atom_atom_parity_2_experiment (phase scan): ...". Before,
  6041 printed the same text as 4724, so the log could not tell parity from
  parity_2. Both lines are optional and can be dropped from the edit.

`applets_master_satellite.py`: `FIRST_SHOT_ALIVE_EXPERIMENTS` (241) gains
`"Two_nodes_atom_atom_parity_2_experiment"` next to the original parity name.
Without it, a parity_2 run with `advance_if_atom_in_first_shot` on would close
the "First shot alive atom loading time (s)" applet instead of opening it.
The ZZ names stay out, as today (section 7).

Nothing else keys on these names:

* GVS builds its two-node registry from every module function with
  "experiment" in its name, so the `_2` names appear in the dropdown by
  themselves (6 -> 8 names). `_run_*_ZZ_2` has no "experiment" and stays
  hidden. The default stays `Two_nodes_HOM_experiment`.
* The GVS h5 copy is named from `experiment_name[:-11]`, giving
  `..._TwoNodes_Two_nodes_atom_atom_ZZ_2_scan_over_<var>.h5`;
  `Analysis/h5_data_utilities.py` `ms_two_nodes_files` maps that back to the
  `_2` name.
* Base seeds the datasets the same way for every experiment.

Tests (`tests/test_experiment_functions_two_nodes.py`,
`tests/test_general_variable_scan_master_satellite.py`): the name-pinning
tests now list all four names. New helper `_functions_run_by` follows a host
wrapper into the kernel it calls, so a check on what an entry point writes
covers wrapper and kernel versions alike. Two new tests:

* `test_no_function_is_defined_twice`: no top-level name defined twice.
  This is the guard against the next pasted block redefining a name.
* `test_atom_atom_versions_run_their_own_kernels`: each entry point reaches
  its own kernel; parity writes `atomic_parity_code_version`
  `20261110_latest_ZZ_append_only_v1`, parity_2 writes
  `20261111_independent_maintenance_cycle_log_v2`; only the `_2` entries write
  `atom_atom_cycle_log_version`.

## 3. What the original names run after the rename

From the 11-14 behaviour audit (read-only; offline compiles and h5 checks).

**`Two_nodes_atom_atom_ZZ_experiment` (3913).** What RID 40946 ran
(11-11 18:41-19:02, before the 19:41 block deploy), minus the per-herald
`excitation_cycle_of_row` append that was taken out at 19:41. The same as
ZZ_2 in loading and the first-shot confirm, the MW rearm and setup, the DMA
trace order, the herald (exactly two clicks), the any-click stop, the mapping,
the blow-away, the alternating readout and the measurement counting.
Different:

* Maintenance on a check cycle: PGC coils + 0.4 ms, the check, then the
  recool if its own schedule is due (check -> recool). ZZ_2: PGC prep, the
  scheduled recool, the check, a post-check recool (recool -> check ->
  recool). At check 1 / recool 1 (41177's settings) every cycle differs; at
  the stored check 30 / recool 1 only the check cycles do.
* Recool: the DMA trace `recooling_both_nodes` (today's
  `record_recooling_both_nodes`, 368-410: PGC coils, PGC DDS, 0.4 ms, 0.6 ms
  light) with a 20 us lead. About 0.1 ms per cycle faster than ZZ_2.
* No cycle-log datasets: no `atom_atom_cycle_log_version`,
  `atom_atom_maintenance_order`, `excitation_cycle_index_base`,
  `atom_check_every_n_of_row`, `recool_every_n_OP_of_row`.
  `excitation_cycle_of_row` stays at Base's `[0]` seed (one entry, not one per
  herald).
* No host-side check that the intervals are positive integers.
* Drops a passed check on the final cycle of a pass that runs out at
  `max_excitation_cycles` with no herald (rare).
* Its log line says "Two_nodes_atom_atom_parity_experiment: Node2 delay D ="
  (4037, an old mislabel, not changed).

**`Two_nodes_atom_atom_parity_experiment` (5136 -> 4647).** The 11-10 version
that 40871, 40884-40886, 40904 and 40926 ran (`atomic_parity_code_version`
`20261110_latest_ZZ_append_only_v1`), with one exception:
`_record_atom_atom_recooling_after_check_both_nodes` (4567-4620) was edited
11-12 23:08 and that trace has never run on hardware. The cooling DP and the
repump now open at +5 us with A1/A2 instead of +20 us after A5/A6, so the
beams start staggered (A1/A2 about 15 us, A3/A4 10 us, A5/A6 5 us more light)
and end together. That timestamp carries 4 events per node.

* The host forces `atom_check_every_n = recool_every_n_OP = 1` (restored in
  `finally`; local datasets record 1). GUI and stored intervals are ignored.
  parity_2 uses the configured intervals.
* Per cycle: PGC coils + 0.4 ms, the check, then the after-check recool.
  One recool per cycle. parity_2 at 1/1: prep, recool, check, recool.
* Both write the six `parity_MW_RF_*` series and the phase/basis scalars.
  The original writes none of the cycle-log datasets.

## 4. Data written under the original names

Since the 11-11 19:41 deploy until the rename, the original names ran the
later code. Of the finished two-node GVS runs 11-10..11-14 with an atom-atom
function, only **RID 41177** (ZZ name, 11-14 11:49) falls in that window. It
carries `atom_atom_cycle_log_version = 20261111_v2_unseeded`. No parity run
was made after the deploy. So `atom_atom_cycle_log_version` (ZZ) and
`atomic_parity_code_version` (parity) separate the versions in every file
there is.

## 5. What to submit

* Today's ZZ (as 41177): `experiment_function = Two_nodes_atom_atom_ZZ_2_experiment`.
* Today's four-point parity: `Two_nodes_atom_atom_parity_2_experiment`.
* A clone of 40946 under the original ZZ name gets back what 40946 actually
  did.
* Driver templates that name the original ZZ:
  `scratchpad/pz/drv/driver.py:738` (`zz_expid` asserts
  `experiment_function == "Two_nodes_atom_atom_ZZ_experiment"`, template
  RID 40946) and `scratchpad/zz1111/driver.py:604` (finished 11-11). Switch
  them to `Two_nodes_atom_atom_ZZ_2_experiment` before they are reused, or
  they run the older maintenance. The pz "ZZ2" labels
  (`ZZ2_COMPILE_PASS`, `pz/compile/run_zz2_compile_live.py`) name a
  configuration, not the `_2` function: they compile 40946's arguments and
  would now compile the older kernel.
* Runs already prepared keep the module they imported; anything built after
  the save gets the new mapping. Deploy with an empty pipeline.

## 6. What makes it safe, and what does not

* Safe: no statement inside any sequence changes. The two print labels are
  string constants. All four entry points compile offline in the edited copy
  (section 8).
* Not safe by itself: the **meaning of the original names changes**. That is
  the point of the request, and it is silent for anything that submits by
  name (section 5).
* The revived originals run with **today's** helpers: the
  `record_recooling_both_nodes` PGC lines (uncommented since HEAD fa70fe4) and
  the 11-12 edit of the after-check recool trace. The rename does not cause
  this, but neither version of the original parity's after-check trace has
  been checked on hardware.

## 7. Open decisions (user)

1. `FIRST_SHOT_ALIVE_EXPERIMENTS`: the ZZ entries write
   `Atom_loading_time_first_shot_alive` but were never in the list, so a ZZ
   run closes that applet. Kept as is (pinned in
   `test_every_confirm_loop_records_its_first_shot_alive_time` and
   `test_first_shot_alive_applet_only_with_the_confirm_on`). Adding both ZZ
   names would be a deliberate applet change.
2. Comments outside the 11-11 block left alone (comment-only; some files are
   hashed by drivers):
   `utilities/BaseExperiment_master_satellite.py:1815` ("Two_nodes_atom_atom_ZZ_experiment
   resets it to [] (no seed)" - now ZZ_2 and parity_2 do; the originals keep
   the [0] seed), 1800/1807/1888; `applets_master_satellite.py:185`;
   `subroutines/experiment_functions.py:7028` (atom_photon_parity_12
   docstring; the cycle loop it means is now ZZ_2's);
   `applets/plot_retention_and_loading.py:10`; `plan_claude_detail.md:2853, 3657`.
3. Analysis (open in Jupyter, not edited):
   `Analysis/master_satellite/2026_Two_Nodes_Parity_Analysis.ipynb` cell
   a2a00002 line 10 `experiment_functions` must gain
   `Two_nodes_atom_atom_ZZ_2_experiment` and
   `Two_nodes_atom_atom_parity_2_experiment`, or `_2` runs are never listed
   (`ms_two_nodes_files` matches names exactly). Recommended: add
   `atom_atom_cycle_log_version` and `atomic_parity_code_version` to
   `COMBINE_SETTINGS` (line 301), so a combine across versions is flagged.

## 8. Verification (scratch copy `ren2/s`)

* Unit tests (`C:/Python39`, below-normal priority, one at a time).
  Live: 253 tests, 8 failures + 1 error. Edited copy: 255 tests, 2 failures.
  The 2 that remain fail on the live tree too, for reasons unrelated to the
  rename: `test_no_run_of_switch_events_shares_a_timestamp` (same
  offenders, only two line numbers +4) and
  `test_single_node_registry_uses_current_functions_and_exclusions`
  (identical). Fixed by the rename's test updates: the 7 name-pinning
  tests. Mutation checks: the live source fails the duplicate guard; a
  parity_2 wrapper calling the old kernel fails the versions test; dropping
  parity_2 from `FIRST_SHOT_ALIVE_EXPERIMENTS` fails the confirm-loop test.
* Offline compile (nothing executed; `ren2/s_compile/compile_s.log`), all PASS:
  ZZ_2 and the original ZZ with RID 41177's arguments; parity_2 and the
  original parity with RID 40926's four-point arguments. The loaded module
  was the edited copy (sha256 `86dfd917...`).
