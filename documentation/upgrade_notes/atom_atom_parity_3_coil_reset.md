# Atom-atom phase scan with ZZ_3's coil reset: `Two_nodes_atom_atom_parity_3_experiment`

**Started:** 2026-11-17 **Status:** staged and tested in a scratch copy (scratchpad `par3/`); deployed by
`par3/deploy.py` at an empty pipeline (this note is copied to `documentation/upgrade_notes/` by the deploy).

User, 2026-11-17: "run atom-atom parity experiment. If the parity experiment hasn't yet been modified to account for
the coil drift, make another version so that it follows how it is being done in ZZ_3."

`Two_nodes_atom_atom_parity_2_experiment` (the four-point equatorial phase scan) runs ZZ_2's maintenance loop
(`_run_atom_atom_phase_scan_latest_ZZ_2`): it had not been given ZZ_3's coil reset
([atom_atom_ZZ_3_coil_reset.md](atom_atom_ZZ_3_coil_reset.md)).

---

## 1. What changed

`subroutines/experiment_functions_two_nodes.py`: one block appended at the end; nothing above it is edited.

- `_run_atom_atom_phase_scan_latest_ZZ_3`: parity_2's kernel with exactly the edit that made ZZ_3's kernel out of
  ZZ_2's. After EVERY non-herald excitation cycle: 1 ms more at the O.P. field, the PGC coils on both nodes
  (`_atom_atom_prepare_PGC_field_both_nodes`, 0.4 ms settle), then a dark wait so that the coils have sat at PGC for
  `ZZ_3_COIL_RESET_DWELL` (42.9 ms, the constant of the ZZ_3 block) when the next O.P. trace starts; the periodic
  recool, the atom check and the post-check recool are subtracted from the wait and still run last. A unit test checks
  that the cycle-loop change parity_2 -> parity_3 equals ZZ_2 -> ZZ_3 statement for statement.
- `Two_nodes_atom_atom_parity_3_experiment`: parity_2's host wrapper calling the new kernel. Same scan (dummy_variable =
  Node2's RF phase in turns, 0 / 0.25 / 0.5 / 0.75; Node1's at 0; pi/2 = half each node's stored t_MW_RF_pi_pulse),
  same datasets and per-herald phase/cycle metadata. It writes `atom_atom_maintenance_order =
  "coil_reset_every_cycle_then_periodic_recool_then_check_then_recool_if_passed"`, the archive-only
  `atom_atom_coil_reset_PGC_dwell_s`, and `atomic_parity_code_version = "20261117_parity_3_coil_reset_v1"`;
  `atomic_parity_basis` stays `equatorial_phase_scan`, so the parity notebook's section F reads it as a phase scan.

`applets_master_satellite.py`: parity_3 added to `FIRST_SHOT_ALIVE_EXPERIMENTS` (its confirm loop, parity_2's, writes
`Atom_loading_time_first_shot_alive`; a test derives that list from the source).

Tests: parity_3 in the entry-point, registry, kernel and applet tables of `tests/test_experiment_functions_two_nodes.py`
and `tests/test_general_variable_scan_master_satellite.py`, and
`test_parity_3_resets_the_coils_after_every_non_herald_cycle`.

## 2. Verification before the deploy

- Unit tests in the staged copy: Ran 269 (live 268 + the new test); failing = the two the live tree already had
  (`test_no_run_of_switch_events_shares_a_timestamp`, `test_single_node_registry_uses_current_functions_and_exclusions`);
  two more fail only in the scratch copy because the AOM-feedback config path is resolved relative to the repository.
- Offline compile (`par3/compile/run_par3_compile.py`, template RID 41333): parity_3 (4 phases x 5, check 1, recool 1),
  the parity_2 control, and parity_3 at check 30 / recool 2 all PASS; parity_3 writes the order, the 42.9 ms dwell and
  the new code version; pi/2 = 128.2 us (Node1) and 52.25 us (Node2) at the stored pi times.

## 3. Cost and expectations

As ZZ_3: a non-herald cycle takes ~46 ms instead of ~23.7 ms (1 attempt, check 1, recool 1, 20 ms check), so the herald
rate drops to ~x0.56-0.80 of parity_2's. Only the coil history is reset (see the ZZ_3 note, section 3).
