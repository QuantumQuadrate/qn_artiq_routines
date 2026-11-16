# GVS health check: direct re-optimization

**Started:** 2026-11-14 **Status:** **DEPLOYED 2026-11-14 (evening), uncommitted**, together
with the atom-atom `_2` rename (atom_atom_ZZ_parity_2_rename.md). Not yet run on hardware:
pilots P1/P2 (section 6) come first.

**User decisions (2026-11-14 ~18:20):** for the ZZ, re-optimize **both nodes** (order as
coded: Node1, then Node2); **cap: after 3 direct scans of one (node, transition) in one
event without a pass (failed or crashed), the GVS ends cleanly before the point**
(`HEALTH_CHECK_DIRECT_MAX_ATTEMPTS` = 3; local datasets `gvs_health_check_stopped_at_scan_point`
and `gvs_health_check_stop_reason`; status line "STOPPED before scan point N"); check mode
stays unbounded. `health_check_optimize_directly` without `run_health_check` is refused at
prepare. Open decisions 2-4 below are settled by this; 1 and 5-9 stay open.

User (2026-11-14): "F11 drift in node2 is mostly due to the FORT polarization
drift due to SM fiber; AC Stark shift due to that is what changes the
frequency. That is why we need health check. Currently there is no option in
healthcheck to do optimization right away without checking. I think it is
important to implement it. And for long experiments like ZZ, we need to
reoptimize for the frequency 11 every ~5points."

**Scope:** the master-satellite GVS health check
(`GeneralVariableScan_master_satellite_mixin.py`, deployed 2026-10-29,
committed 52c831c), its tests, and both plans (section 8.9 and the decision
log of `plan_claude_detail.md`, item 4 of `plan_claude_concise.md`).
Not `MicrowaveScanOptimizer_master_satellite.py`, the standalone HealthCheck
files, the notebooks or ezz.py. "Was" line numbers are the live mixin of
2026-11-14 (sha256 `3a8d73fb...`, which carries the uncommitted
`_MICROWAVE_CALIBRATION` change at 239-243); "now" numbers are the edited
scratch file.

## 1. The change

* New GUI argument `health_check_optimize_directly` (Bool, default False,
  group "Health Check", right after `run_health_check`; now 490-499) and a
  class attribute of the same name, False (now 318), for hand-built scans and
  the test helpers.
* `_health_check_expid` (was 1376-1408, now 1393-1432), when on, differs from
  the check's expid in exactly three values; the keys are the same:

  | Key | Check | Direct |
  | --- | --- | --- |
  | `run_health_check_and_optimize` | True | False (the optimizer's scan branch) |
  | `n_measurements` | `health_check_n_measurements` | `health_check_optimization_n_measurements`, or the check's when that is <= 0 (the optimizer's own rule, `submit_opt_exp_general` 1081-1083) |
  | `enable_geometric_frequency_scan` | True (unused by a check) | False (one fixed list, the 11-14 rule of the optimizer's child, 1102) |

  With `enable_node_specific_adaptive_scan` True, the optimizer's prepare
  overwrites every kHz value from the node's pi time
  (`_size_frequency_scan_from_pi_time`, optimizer 1179-1201), so the direct
  scan scans the same list, around the same stored centre, as the child a
  failed check submits. A test prepares both and compares the lists.
* `_read_health_check_result` (was 1622, now 1650-1655): a lone failure of
  `chain_rids[0]` counts as "crashed" only in check mode (there it is the
  check's own failure, so its scan must have crashed). In direct mode
  `chain_rids[0]` is the scan, and its failed fit check is "failed": the
  same scan is resubmitted at once at 99. Without this, every direct
  fit-check failure would wait out the 60 s crash backoff.
* Wording only: the two tooltips (now 528-538); the block comment (steps 1 and 3, now
  1206-1217); docstrings of `_health_check_expid`, `_read_health_check_result`
  and `_run_one_health_check`; the log line "Node2 freq11 (direct
  re-optimization), forwarded ..." (now 1691); the status header "MW health
  check (direct re-optimization) every ..." (now 1503-1506); the current line
  "re-optimizing since" instead of "checking since" (now 1704). The status
  heading "Resubmissions this run (did not pass / crashed chain):" (now 1513)
  changed for both modes (was "failed check / crashed chain"); nothing parses
  that text.

## 2. Decisions already taken

* No optimizer change. In scan mode with a `parent_rid` it already publishes
  the parent-keyed result from its fit check (880, 1003, 1006-1046), stores
  the centre (round(f0, -3), persisted, node-suffixed) only on a pass (884),
  and on a failure stores and submits nothing (888).
* Same resubmission policy as a failed check: at once, at 99, without limit;
  a crash (no result of its own: a FitError, an exception in prepare or run,
  a deleted run) after `health_check_crash_backoff_s` (60 s).
* Same cadence (never before scan point 0, `_health_check_due`, now 1434),
  order (nodes, then 00, 01, 11, m10, m11), wait and pause, centre re-read,
  records, status applet and Terminate handling.
* No new ExperimentVariable or dataset. The mode is in the h5's archived
  expid arguments and in the status/log text.

## 3. What makes it safe

* Default False. ARTIQ 7 gives an expid without the key the default
  (`artiq/language/environment.py` 217-222, `ProcessArgumentManager.get`), so
  every saved expid, dashboard clone, `artiq_client` submission and scratch
  driver stays check-first; the classifier change is gated by the flag.
* With the flag False the expid is byte-identical to today's (test 2680, now
  `test_expid_carries_the_forwarded_dictionary_to_its_node_only`, unchanged
  and passing); only the status heading text differs.
* The optimizer is untouched, so no paired deploy. An expid carrying the new
  key sent to the OLD mixin is silently check-first (ARTIQ ignores unknown
  keys): deploy first, then submit, and confirm the mode from the status
  header or the "(direct re-optimization)" log line.

## 4. Tests (all run at below-normal priority, one at a time)

New, in `tests/test_general_variable_scan_master_satellite.py`:
`test_direct_expid_is_the_scan_a_failed_check_submits`,
`test_direct_one_scan_per_transition_in_order_never_before_the_first_point`,
`test_direct_pass_stores_and_rereads_the_center`,
`test_direct_failed_scan_is_resubmitted_at_once_at_99` (the regression for
the classifier), `test_direct_crashed_scan_waits_for_the_backoff`,
`test_direct_terminate_during_the_scan_ends_after_it`; the defaults test
asserts the new argument and class attribute; `_make_scan` takes `direct`.
In `tests/test_microwave_scan_optimizer_master_satellite.py`:
`test_gvs_direct_expid_runs_the_scan_a_failed_check_submits` (every key an
argument; prepare in scan mode with n 100; 16 points at 1.5/5.0, equal to the
child's list) and `test_a_frequency_scan_stores_only_a_passing_fit_check`
(F11 Node2, f0 339.2714 MHz: stored 339271000.0 on a pass, nothing on a
failure, nothing submitted either way).

| Tree | Tests | Failures |
| --- | --- | --- |
| live, 11-14 (baseline) | 253 | 8 failures + 1 error, all pre-existing (two-node registry/structure and single-node registry) |
| `hcopt/wt` (this change) | 261 | the same 9, identical set |
| `hcopt/w2` (ren2's `_2` rename, then this change) | 263 | 2: `test_no_run_of_switch_events_shares_a_timestamp`, `test_single_node_registry_uses_current_functions_and_exclusions` (ren2's expected 2) |

The shared GVS test file composes in either order (same bytes,
`hcopt/order_check.py`). Offline compile (`hcopt/compile`, RID 41177's
arguments, three dummy points, health check on, direct, both nodes, F11,
every point): PASS on `wt` with `Two_nodes_atom_atom_ZZ_experiment` and on
`w2` with `Two_nodes_atom_atom_ZZ_2_experiment`, 6 kernel compiles each. No
kernel code changed; the harness's scheduler is RID 0, so the events take the
"no ARTIQ master (RID 0), health check skipped" branch: build, prepare, the
records and the status text ran with the new code, the submission did not
(the unit tests cover it).

## 5. Latent issues noticed, not fixed here

* Frequency scans have no f0-outside-the-scan guard (Time scans have one,
  optimizer 891-897).
* A failed direct scan stores nothing, so its resubmission scans the same
  window: a line beyond +-1 linewidth (Node2 about +-99 kHz) or an OP fault
  loops until relocked, as with check-first.
* An optimizer frequency scan queued behind a long GVS is prepared when the
  GVS starts; its list (optimizer 324-337) is not rebuilt by the run-start
  reload (556-562).
* The in-run centre is recorded only in `gvs_health_check_center_after_Hz`;
  ezz.py and the notebook settings table show the run-start f11.
* `force_other_node_off` (optimizer 577) keeps the other node's FORT off for
  each optimizer run: 6-14 min per direct scan instead of about 1.2 min per
  passing check. Effect on the FORT fibre polarization untested.
* `health_check_n_measurements` below 1 is still refused in direct mode, even
  when the optimization n is above 0 and the check's n is unused.

## 6. Pilot (the user submits)

* P1: two-node GVS, `dummy_variable` `np.arange(0,3,1)`, every 1, Node2,
  F11, direct. Two events; first hardware test of store -> GVS re-read.
* P2: as P1 with `target_fidelity` 0.99 and `np.arange(0,2,1)`: one "did not
  pass ... resubmitting now", then Terminate during the second scan; the scan
  ends after it, nothing is left queued.

## 7. Using it for a long ZZ

`GeneralVariableScan_master_satellite_two_nodes`, `experiment_function` the
ZZ (today `Two_nodes_atom_atom_ZZ_experiment`; **after the `_2` rename, the
same sequence is `Two_nodes_atom_atom_ZZ_2_experiment`** and the old name
runs the earlier kernel), `dummy_variable` `np.arange(0,P,1)`, main
pipeline, a distinct priority below 99. Health Check: `run_health_check` and
`health_check_optimize_directly` on, `health_check_every_n_ite` 5, node
selection Node2 (or both_nodes), only `Frequency_11_Scan`, target 0.80,
optimization n 100, adaptive 1.0/4.0. Do not override or scan
`f_microwaves_11_dds_Node1/2` (prepare refuses it). Run fresh F11 (and F01)
optimizers before the ZZ: points 0-4 get no event.

## 8. Open decisions for the user

1. An event before point 0 in direct mode? (Recommended: no, as now.)
2. One switch, or per node? (Recommended: one; select Node2 only unless
   Node1 drifts.)
3. both_nodes order: Node1 runs first, so Node2's FORT is off for Node1's
   7-8.5 min scan right before Node2 is re-centred. Node2 only, or reverse?
4. Resubmission limit: unbounded, as today, although a direct failure now
   costs a full 6-14 min scan per attempt.
5. F01 too? F11 maps through f01, so a drifted f01 shows up as repeated F11
   fit-check failures.
6. A frequency range guard in the optimizer (deploy O and M together).
7. In-run f11 in ezz.py and the notebook settings table (notebooks are the
   user's).

README row (`documentation/upgrade_notes/README.md`):

| [gvs_health_check_direct_optimization.md](gvs_health_check_direct_optimization.md) | GVS health check option `health_check_optimize_directly`: re-optimize each ticked transition at every event instead of checking first (F11 every ~5 points of a long ZZ) | Built and tested, not deployed (2026-11-14) |
