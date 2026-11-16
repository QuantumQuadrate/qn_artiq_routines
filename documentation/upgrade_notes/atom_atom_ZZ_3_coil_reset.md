# Atom-atom ZZ with a parity_11 coil reset: `Two_nodes_atom_atom_ZZ_3_experiment`

**Started:** 2026-11-16 **Status:** **DEPLOYED 2026-11-16 (night), uncommitted**. Built in a scratch copy
(scratchpad `zz3/`), offline-compiled there and again on the live tree after the copy; unit tests on both.

User, 2026-11-16 01:19: "rerun parity13 in node1 with atom_check_n =1, excitationattempt = 1. If the fidelity is
less than 91%, i think we should modify ZZ experiment so that it mimics what parity11 does - to reset the coils";
01:28, for that case: "I deploy and start ZZ".

Trigger: Node1 parity13 at QWP 73 with the ZZ's per-node settings (RID 41305: 1 attempt, check + recool every
cycle, 20 ms shots, HWP 0..70, n 100) gave a = 0.859 +- 0.044, against parity_11 at the same QWP (RID 41301, 10 ms
shots) a = 0.961 +- 0.042.

---

## 1. What changed

`subroutines/experiment_functions_two_nodes.py`: one block appended at the end; nothing above it is edited.

- `ZZ_3_COIL_RESET_DWELL = 42.9e-3` (module constant, s; no new ExperimentVariable).
- `_run_atom_atom_ZZ_3_with_cycle_log`: the kernel of `Two_nodes_atom_atom_ZZ_2_experiment`
  (`_run_atom_atom_ZZ_with_cycle_log`) with one change. After EVERY non-herald excitation cycle:
  1. 1 ms more at the O.P. field (parity_11's delay before its second_shot);
  2. the PGC coils on both nodes (`_atom_atom_prepare_PGC_field_both_nodes`, 0.4 ms settle). ZZ_2 loads them only
     when a recool or an atom check is due;
  3. a dark wait at the PGC field, so that the coils have sat at PGC for `ZZ_3_COIL_RESET_DWELL` when the next O.P.
     trace starts. The periodic recool, the atom check and the post-check recool are subtracted from the wait, so
     they still run last before the O.P. trace, and the check still catches an atom lost during the wait.
- `Two_nodes_atom_atom_ZZ_3_experiment`: the host wrapper, as ZZ_2's. It writes
  `atom_atom_maintenance_order = "coil_reset_every_cycle_then_periodic_recool_then_check_then_recool_if_passed"`
  and the archive-only `atom_atom_coil_reset_PGC_dwell_s` (= the constant).

Unchanged: cycle 0 of each pass, the herald path (attempts, mapping at the O.P. field, blow-away, alternating
readout), f11 (set once per pass), every dataset name and the per-herald cycle metadata.

`tests/test_experiment_functions_two_nodes.py`: ZZ_3 in the registry and entry-point tables, and
`test_zz_3_resets_the_coils_after_every_non_herald_cycle` (ZZ_3's cycle loop is ZZ_2's plus the reset; the helper
still ends with its 0.4 ms settle). `tests/test_general_variable_scan_master_satellite.py`: ZZ_3 in the two-node
registry list.

## 2. Why 42.9 ms

parity_11 runs one excitation cycle per pass, so between two of its O.P. traces the coils sit at the PGC set from
second_shot's coil load to the next pass's O.P. trace: 22.9 ms of programmed delays + t_SPCM_first_shot +
t_SPCM_second_shot (Node1), i.e. 42.9 ms in RID 41301 (10 ms shots); its pass rate gives ~43-46 ms in practice.
ZZ_2 at check 1 / recool 1 with a 20 ms check: ~21.9 ms. With the ZZ's 20 ms shots parity_11 would sit 62.9 ms:
change the constant AND the `assertAlmostEqual(dwell, 42.9e-3)` line in
`test_zz_3_resets_the_coils_after_every_non_herald_cycle`.

## 3. Cost and expectations

- A non-herald cycle takes ~46 ms instead of ~23.7 ms (1 attempt, check 1, recool 1, 20 ms check): the herald rate
  drops to roughly x0.56-0.80 of ZZ_2's; ~22 ms more dark trap time per cycle.
- Meant for atom_check_every_n 1 and recool_every_n_OP 1: sparser schedules also get the full dwell every cycle.
- The evidence for a gain is weak. Cycle 0 of a parity13 pass already has parity_11's coil history: in 41305 it gave
  a = 1.05 +- 0.16 (43 heralds) against 0.853 +- 0.046 for cycles >= 1; in six earlier 20 ms-check runs (11-15)
  0.835 +- 0.068 against 0.834 +- 0.018. Node2's coil-memory model gives ~5 kHz of f11 shift for ZZ_2 at check 1
  and ~1.3-1.5 kHz at 42.9 ms; Node1's coil memory is not measured. Only the coil history is reset: the light
  history between cycles stays ZZ_2's (one 20 ms check plus recooling, against parity_11's two readouts), and the
  20 ms vs 10 ms shot length is not addressed. A ZZ_3 result no better than ZZ_2 rules out the coils only.

## 4. Use

- GVS two-node, experiment function `Two_nodes_atom_atom_ZZ_3_experiment`; arguments as for ZZ_2.
- The dashboard needs "Scan repository" before ZZ_3 shows in the list.
- `Analysis/master_satellite/2026_Two_Nodes_Parity_Analysis.ipynb` lists its files only after
  `Two_nodes_atom_atom_ZZ_3_experiment` is added to its experiment-function list (not edited here; the notebook is
  the user's). Combining ZZ_2 and ZZ_3 files is NOT refused: combine_parity_runs merges them and only prints a
  WARNING naming experiment_function. Combine series-B (ZZ_3) files only with each other.
- The driver's own E_ZZ for series B (p4/ezz_B.py, ezz_pool_B.py) uses the notebook's thresholds, Node1 14000 and
  Node2 17000 c/s; series A's driver numbers used the stored Node2 21000.
