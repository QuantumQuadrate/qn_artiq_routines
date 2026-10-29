# Master-Satellite System: Startup Brief

File `plan_claude_concise.md`; full specification `plan_claude_detail.md`.
State: 2026-10-16, branch `20261006_two_nodes_experiment`; HEAD is the
commit that last updated this file.
"Section N", the rules R1-R26 and the open items O-1 to O-36 all refer to
`plan_claude_detail.md`.

## 1. Purpose and precedence

This brief is for a new AI assistant or lab member. It covers what the
system is, the rules that bite most often, the current status, and where to
read more. It holds no fact the detailed file lacks. If the two disagree,
the detailed file wins and this brief gets corrected (section 1.3).

Order of authority (section 1.2): the user's latest explicit instruction,
then `plan_claude_detail.md`, then this brief, then code comments, then
assumptions (which must be reported before use).

For what actually runs, the code on disk and the hardware evidence (h5 files,
RIDs, logs) are authoritative. When they disagree with the plan, report it
and fix the side that is wrong. Several code comments are stale (O-25).

These two files replace `PLAN.md`, `PLAN_FROM_SCRATCH.md`,
`plan_codex_detail.md` and `plan_codex_consise.md` (section 1.4), which are
deleted and kept only in git history. The code that cited them now cites
these files (section 1.6, O-24).

## 2. Physical system in brief (section 2)

- **The system:** two neutral Rb-87 single-atom quantum-network nodes on
  ARTIQ 7 / Sinara. The two former Kasli-SoC systems now run as ONE DRTIO
  system.
- **Nodes:** Node1 (Alice) is the DRTIO master, at 192.168.1.75 (a v1.1
  board since 2026-09-09). Node2 (Bob) is the satellite, destination 1. It
  runs satman, so it has no network; its console is UART on COM8.
- **Link:** fiber SFP0 to SFP0, 125 MHz RTIO, external 125 MHz clock.
- **Per node:** 2 DIO cards with edge counters, 3 Samplers, 1 Zotino and
  3 AD9910 Urukuls.
- **Urukuls:** Node2's `urukul4`/`urukul5` were replaced on 2026-09-15 to
  fix a SYNC fault. Never use v1.6 cards: they are proto_rev 9, and ARTIQ 7
  needs 8 (section 2.10).
- **Device database:** a single file, `device_db/device_db_master_satellite.py`
  (= `artiq-master/device_db.py`). `channel >> 16` gives the destination.
  Node2's devices are `ttl16`-`31`, `sampler3`-`5`, `zotino1` and
  `urukul3`-`5` (sections 2.4, 2.5).
- **Detectors:** all four SPCMs (H1/V1/H2/V2 = ttl0/1/8/9) are master-local
  Node1 inputs. Each one sees both traps (section 2.6).
- **GRIN aliases** (renamed 2026-10-20 from `GRIN1and2_dds` and
  `dds_D1_pumping_DP`): `dds_GRIN1` drives GRIN1 (D1 pumping) and
  `dds_GRIN2` GRIN2 (excitation) on both nodes since the 2026-10-20 Node1
  GRIN rewire; the code was changed ahead of the cables (section 2.8).
- **USB devices:** all connect to the master host. One `k10cr1_ndsp`
  controller drives the 8 rotators (section 13).

## 3. Rules and working conventions (section 3)

Priorities, in order: working hardware control, easy debugging, minimal-risk
changes, preserving validated sequences, then easy rollback. All of these
come before cleanup (sections 3.1, 3.3).

**Only on explicit instruction:**
- R2: commit or push. The user verifies on hardware first. Never rewrite
  history.
- R3: flash gateware or operate outputs.
- R6: add, rename or move files, imports or modules. Explain why, with
  alternatives, first.
- R7: make a cross-cutting change. Write a note in
  `documentation/upgrade_notes/` instead.
- R20, R24: change applets, build InfluxDB/Grafana, or redesign the
  dashboard.

**Always:**
- **R26, saved means deployed.** The master serves the working tree: it runs
  without `--git`, and every h5 `expid` has `repo_rev` "N/A". So a saved
  edit runs on the next submission, committed or not. Never leave
  half-finished code on disk while the user may be running. "Uncommitted"
  does not mean "not yet run".
- **R1:** keep the standalone stack working, with no dependency on
  master-satellite code. Shared modules get only backward-compatible,
  authorized extensions whose default is standalone behaviour.
- **R4:** check `git status` and preserve unrelated changes. Notebooks may
  be open in Jupyter (section 15.5).
- **R5:** inspect first and cite evidence. Keep facts apart from
  assumptions. Passing tests is not hardware proof. Flag uncertain hardware
  behaviour rather than guessing.
- **R8:** before a bulk edit, pilot one named site, verify it, report, then
  stop (section 3.4).
- **R9:** keep each sequence one linear function in the legacy shape. New
  machinery goes in Base.
- **R10, R11:** never calculate DRTIO channels or offsets. Two-node code
  names each device's node (`self.dds_FORT_Node2`).
- **R12, R13, R14:** avoid new ExperimentVariables. Submitted GUI arguments
  win for that run. Obey the SED limit (section 5).
- **R15, R16:** call `core.break_realtime()` inside shared helpers and after
  every `laser_stabilizer.run()`. Never lazily import `aom_feedback.py`.
- **R17:** a local `set_dataset` after this run broadcast the same name
  deletes it from the master database.
- **R18, R19, R21, R25:** leave historical experiments alone, and do no
  unrelated cleanup. Probe Urukul SPI only with `set_config_mu(...)`, never
  `cpld.sta_read()`. Analysis logic goes in `Analysis/h5_data_utilities.py`.
- **R22, R23 (settled):** one `k10cr1_ndsp` controller; no `best_effort`;
  no two_nodes mode for MicrowaveScanOptimizer. The legacy
  `which_node == "two_nodes"` never returns.

## 4. Execution model (sections 4 to 8)

**Stacks** (sections 4.1-4.3):
- **Standalone:** `standalone/` (17 entry points),
  `utilities/BaseExperiment.py` and `utilities/DeviceAliases.py`, on
  standalone gateware. Kept for rollback, but not re-run on that gateware
  since the migration (O-20).
- **Master-satellite:** the top-level `*_master_satellite*` files,
  `utilities/BaseExperiment_master_satellite.py` ("Base"),
  `utilities/DeviceAliases_master_satellite.py` and
  `subroutines/experiment_functions_two_nodes.py`.
- **Shared:** `subroutines/experiment_functions.py`,
  `subroutines/aom_feedback.py` and `utilities/config/alice|bob/`.

**Modes** (sections 5.6, 6.5, 6.6):
- **single_node (Node1 or Node2):** the selected node gets the old bare
  names (`self.dds_FORT`, `self.f_FORT`), so legacy functions run
  unchanged. `force_other_node_off()` makes the idle node safe once per run
  (AOMsCoils never calls it).
- **two_nodes:** suffixed names only, so a bare device name fails to
  compile. Functions come only from `experiment_functions_two_nodes.py`
  (4 registered). Settings both nodes share are plain globals, or
  `<x>_two_nodes` globals read bare in two_nodes mode (2026-10-27).
- Neither mode is the removed legacy `two_nodes` node mode (section 4.6).

**Explorer** (section 4.5):
- GVS `_single_node`/`_two_nodes` and AOMsCoils `_Node1`/`_Node2`.
- ExperimentVariables `_Node1`/`_Node2`/`_global`.
- The single-node experiments MicrowaveScanOptimizer,
  FORT_Polarization_Optimizer, AtomLoadingOptimizer,
  SamplerMOTCoilAndBeamBalanceTune and MonitorMagnetometer.
- `tests/MonitorSPCMinApplet` (Node1).
- `_mixin` files are private (section 4.4).

**Naming** (sections 5.1, 5.2, 6.9, 7.1):
- The node is a SUFFIX, never a prefix (`dds_FORT_Node2`). Physical names
  (`urukul4_ch3`) never get one.
- Node1 = legacy `alice`, Node2 = `bob`.
- Each node has 213 `*_NodeX` variables; joint quantities are 18
  unsuffixed globals.
- Node calibrations are suffixed in single_node mode. Run state such as
  `measurements_progress` is unsuffixed, and so are two-node results,
  because loading is one joint measurement.

**GVS flow** (sections 8.3, 8.4):
1. Prepare: `configure_execution` -> `base.prepare()` ->
   `prepare_laser_stabilizer()` -> resolve the scans and the global and
   per-node override dictionaries (shared variables in
   `override_ExperimentVariables_global` since 2026-10-29).
2. Run: reload if queued -> apply the overrides, recording them first
   (R17) -> seed the result state -> force the idle node off.
3. Each scan point: `initialize_hardware()` (core reset) -> reset the
   results -> call the function -> save a named h5 copy.

## 5. Essential invariants

- **Single physical resolution** (sections 5.4, 5.5, 5.7): a node plus a
  logical device gives the same physical device in either mode. The only
  translation source is `utilities/config/master_satellite/device_aliases.json`
  (43 explicit entries per node; the "+16" pattern is never computed). DDS
  alias and default-variable resolution happen in one step.
- **SPCM rule** (section 5.8): canonical detectors bind through the Node1
  resolver, never through the selected node's TTL translation. `ttl_SPCM0/1`
  and `ttl_SPCM0/1_OtherNode` (and their `_counter`s) map to H1/V1/H2/V2 the
  same way on both nodes. Two-node code gates all four.
- **ExperimentVariables** (sections 6.9, 7):
  - Persisted values beat code defaults, so editing a default in code
    changes nothing on the master.
  - A new Variable breaks EVERY master-satellite experiment until its
    initializer runs (R12). Prefer, in order: reuse an existing variable;
    derive and broadcast it in Base; a module constant. If a new one is
    unavoidable, hand-seed it first (`artiq_client set-dataset --persist`)
    and check that the compile sweep's ACTION NEEDED block is empty.
  - Both stacks share one `dataset_db.pyon`, so persistent names must be
    globally unique. Names shared with standalone (`n_measurements`,
    `two_atom_threshold`, `two_atom_threshold_for_loading`) take whichever
    initializer wrote last.
  - A standalone run de-persists `n_measurements`; every master-satellite
    run re-persists it.
  - Seed datasets unconditionally in each worker, never "only if absent".
- **Timing** (section 12):
  - At most 8 events may share a timestamp; the 9th is silently DISCARDED,
    leaving only a core-log line. Put `delay(5*us)` after every four
    switching events, never inside gated windows or physics-length
    intervals (section 12.5).
  - Use `with parallel` only for per-node stages of different duration,
    never for same-timestamp pairs (section 12.6).
  - Underflows are slack erosion from host RPCs. DRTIO adds nothing per
    output event (RID 38429; sections 12.3, 12.4).
  - Remote input, SPI, Sampler and TTL-timestamp operations are
    latency-sensitive. Cross-node skew is unmeasured (section 12.2).
  - Hardware init is not passive. Never use GVS as the first check after a
    flash or repair (section 6.8).

## 6. Current status (2026-10-14; section 17)

**Works on hardware:**
- **DRTIO:** link up since 2026-09-09; per-event output cost measured;
  Urukul SYNC clean since 2026-09-15.
- **Single-node (both nodes):** GVS (atom loading, microwave, sp9, parity
  11), MicrowaveScanOptimizer, AtomLoadingOptimizer, SamplerMOT, AOMsCoils,
  FORT optimizer (last run RID 38772, 2026-10-12), MonitorSPCMinApplet.
- **Two-node:** namespace sanity, joint atom loading, alternating retention,
  sequential feedback per crate.
- **Not yet run:** master-satellite MonitorMagnetometer and any two-node
  photon experiment. The TTL skew probe ran on 2026-10-18 (RIDs 39100,
  39102); its cables are gone since 2026-10-19.

**Recent changes:**
- **2026-10-07:** two_nodes in Base and feedback (`e768d50`); two-node atom
  loading; `with parallel`; the SED four-event rule.
- **2026-10-11:** PR #68 merged into `main`.
- **2026-10-12:** alternating retention; one two-node loader with the
  confirm flag (`c1ddbd2`).
- **2026-10-13:** confirm-by-first-shot at 31 single-node sites (section
  9.7); `advance_if_atom_in_first_shot` live True; `max_loading_rounds` 20.
  Parity 12/13 and `single_photon_experiment_10` removed.
- **2026-10-14:** `ef55cae` records applied overrides in the h5, first used
  in RIDs 38952-38957, with a loader fallback for older files; `5b2f618`
  updates the notebooks.

- **2026-10-16 (committed):**
  - `619d3ec`: MicrowaveScanOptimizer override-order and scan-restore fixes,
    with tests (section 11.1). Running since 2026-10-14, but neither fixed
    path has been exercised yet (O-17).
  - `bfca64d`: magnetometer monitoring (O-35). The optimizer creates the
    `Magnetometer_*_NodeX` datasets; a `break_realtime` in both readers
    fixes the RID 38968 underflow, and Node1 monitored health checks now
    complete (RIDs 38976, 38980, 38982). Three applet specs (Zero / OP /
    MOT) are opened by hand only (section 14.2). Includes the user's 2 ms
    coil-settle delays (O-34).
  - `a369392`: the plan consolidation (these two files).
  - Node2's A1 photodiode moved to `sampler4` ch6 because `sampler3` ch7 is
    bad (bob `feedback_channels.json`, section 10.1), and the new
    `tests/MOTBeamPowerCheck_master_satellite.py` (section 16.4).

- **2026-10-19 (uncommitted):** `Analysis/master_satellite/
  2026_Two_Nodes_Atom_Loading.ipynb`, ported from the standalone two-node
  notebook, with two-node loaders appended to `Analysis/h5_data_utilities.py`
  (section 15.2).
- **2026-10-19 (`aab3f8c`):** the handshake cables were removed. Every
  legacy handshake alias now points at its node's dummy TTL (input /
  output: Node1 ttl3 / ttl12, Node2 ttl24 / ttl31), and init no longer
  touches Node1 ttl11 or Node2 ttl25 (section 2.7). device_db is
  unchanged. Node1 ttl15 is now Node1's `ttl_D1_pumping_switch` (also
  `ttl_D1_pumping_switch_Node1`), blocked at init like Node2's ttl12; the Node1
  sequences open it since 2026-10-20. `ttl_D1_pumping` was renamed `ttl_D1_pumping_switch`
  everywhere. AOMsCoils Node2's `run_feedback` now re-applies the combined
  GRIN modes (new test), and the skew probe no longer drives ttl15.
  Node1's D1 AOM now runs from an external source switched by ttl15:
  AOMsCoils Node1 `D1_pumping_switch_ON` opens `ttl_D1_pumping_switch` and keeps
  the Urukul D1 channel off (section 11.7); the Node1 sequences followed on
  2026-10-20 (below).
- **2026-10-19:** AOMsCoils `_Node1` (`51e169c`) and `_Node2` are
  self-contained files written box by box like standalone `AOMsCoils.py`,
  with the same GUI arguments and hardware operations; the AOMsCoils mixin
  is deleted (section 11.7). After laser feedback the stabilized AOMs now
  come back on at the corrected power: `turn_on_AOMs()` no longer sets the
  box DDSes, as in standalone.
- **2026-10-20:** the AOMsCoils D1 box is now `D1_pumping_switch_ON`, in each
  file's TTL-gates group (was `D1_pumping_DP_ON`, among the DDSes); what it
  switches is unchanged (section 11.7).
- **2026-10-20:** two-node feedback writes its dB histories as
  `p_AOM_A1_history_NodeX`, as single-node does, so the `feedback RF NodeX`
  applets also update during two-node runs (before: `p_AOM_A1_NodeX_history`,
  which nothing read; section 10.3).
- **2026-10-20:** the Node1 sequences switch D1 with
  `ttl_D1_pumping_switch` (ttl15) instead of `dds_GRIN2`: the Node1
  pumping functions, single photon 4/5, `measure_GRIN1`, the chopped-pumping
  DMA, the two-node Node1 side, standalone alice init and AOMsCoils, the MOT
  scripts and `tests/A2-SPCMGateTimeVsCounts.py`. `CW_optical_pumping_node1`
  and its DMA version did not open Node1's GRIN1 path (superseded by the
  Node1 GRIN rewire below).
- **2026-10-20:** the functions no longer in use (45 of 103) moved from
  `experiment_functions.py` to `subroutines/experiment_functions_obsolete.py`,
  which nothing imports; the single-node GVS menu has 23 entries (section 9.2).
- **2026-10-20 (branch `20261020_GRIN_new_wiring`):** the aliases are renamed
  `GRIN1and2_dds` -> `dds_GRIN1` and `dds_D1_pumping_DP` -> `dds_GRIN2` in every
  file except the archive `subroutines/experiment_functions_obsolete.py`, alias
  and feedback JSON keys included (values, device_db and variables unchanged), with a guard test that each node's alias lists match its
  `ALIAS_MAP` and every feedback config key. Not yet checked on hardware
  (O-32, section 18.4).
- **2026-10-20 (same branch), Node1 GRIN rewire, code ahead of the cables:**
  `urukul0_ch3` (`dds_GRIN1`) drives Node1's GRIN1 AOM and `urukul2_ch2`
  (`dds_GRIN2`) its GRIN2 AOM, as on Node2. Node1 runs Node2's pumping DMA
  (plus its A1-A4 repump) and excitation arming; the Node1 pumping functions
  are kept but no longer called; AOMsCoils Node1 has combined
  `Node1_GRIN1_AOM_ON` / `Node1_GRIN2_AOM_ON` modes. Node1's GRIN values are
  placeholders to calibrate (sections 2.8, 11.7; D1 upgrade note).
- **2026-10-23 (same branch):** one pumping DMA for every single-node
  experiment: `record_CW_optical_pumping` (trace `CW_optical_pumping`; was
  `record_CW_optical_pumping_node2`), pumping repump on A5/A6 only on both
  nodes (Node1's A1-A4 dropped on 2026-10-21 after a test). Six dead pumping
  helpers moved to `subroutines/experiment_functions_obsolete.py` (section 9.2).

**Uncommitted:** check `git status`; under R26 uncommitted code may already
have run.

## 7. Open issues (section 18; next steps in section 18.5)

- **Physics:**
  - O-1: the alternating loader loads no atoms (disabled).
  - O-2, O-3: cross-node skew and the excitation delay are unmeasured.
  - O-4, O-5: parallel feedback and the two-node FORT polarization optimizer
    are not implemented. The grouping is decided: A1/A2/A5, then A3/A4/A6,
    then FORT (section 10.5).
  - O-6 to O-9: parity 10/11 underflow; the two-node confirm loop defeats
    recycling; two-node PGC runs at loading depth; single-node `first_shot`
    has an A5/A6 asymmetry.
- **Code:**
  - O-10, O-11: unprotected slack sites; AOMsCoils core resets.
  - O-13, O-14: old parity 2-7 compile faults; thin sweep coverage.
  - O-33: `_wait_for_satellite` does not actually wait.
  - O-17, O-34, O-35: the 2026-10-16 commits above (O-17 still
    unexercised).
- **Validation gaps (O-18 to O-23):** no `MonitorMagnetometer_master_satellite`
  run (end-of-measurement magnetometer readings exist, RID 38968),
  no AOMsCoils per-channel walk, no standalone rollback test, no DMA
  revalidation record, and no master-satellite HealthCheck files.
- **Hygiene (O-24 to O-30):** plan citations repointed (O-24); fix stale
  comments,
  delete the orphan `max_loading_rounds_two_node`, and fix the
  `ms_node_of` misclassification risk.
- **Queued, not started:** O-31 (chopped readout), O-32 (D1/GRIN rename).

## 8. How to verify (section 16)

- **Unit tests** (hardware-free): from the repository root, with an
  artiq-free interpreter (`C:/Python39`, the `python` in Git Bash), run
  `python -B -m unittest discover -s tests -p "test_*.py"`. Always use
  discover, never `artiq7`. Expect 161 OK on the 2026-10-20 working tree.
- **Compile sweep:** run `artiq-master/compile_check_master_satellite.py`
  with `artiq7` active after ANY master-satellite edit (about 5 minutes).
  Check the PASS/COMPILED/FAIL counts and that the target kernel is in the
  compiled list. Read the ACTION NEEDED block: a green sweep cannot detect a
  missing persistent dataset (section 16.2).
- **Positive controls:** pair every detector with a control that must fire.
  Before trusting a fix, replay real data through the real code (for example
  `DatasetManager`). Compute throughput from the h5 (sections 9.8, 16.3).
- **Hardware sanity:** run small tests per device before any physics
  (sections 16.4, 16.5).
- **Tools** (section 2.12): read h5 with
  `C:/ProgramData/Anaconda3/envs/artiq7/python.exe -I`. `artiq_client
  submit` needs `-R`. For worker tracebacks, use `artiq_run` on a COPY of
  the dataset database.

## 9. How to work and report (section 20)

- **Before acting:** read this brief and the detail sections for the files
  you will touch, then inspect the code. Find the smallest relevant file set,
  its dependencies and possible regressions. Report what you inspected, the
  facts with evidence, ambiguities, risks and the bounded scope. Give
  alternatives for structural changes (R6); name the pilot site for
  repetitive edits (R8).
- **After editing:**
  - List the changed files and APIs.
  - Report the unit tests, the compile sweep, `git diff --check` and
    `git diff --stat`.
  - Say which shared or standalone files changed and why they stay backward
    compatible (R1), and whether any Variable needs seeding (R12).
  - Say what still needs hardware validation, and how to run it.
  - Update both plan files when the design moves. Never commit
    automatically (R2).
- **Hardware records:** image hashes, routing, cabling, boot logs,
  readiness, the operation tested, RIDs and pass/fail evidence (section 20.3).

## 10. Where to read more in `plan_claude_detail.md`

| Topic | Section |
|---|---|
| Precedence, retired files, citation map | 1.2, 1.4, 1.6 |
| Hardware, channel maps, Urukul SYNC, rollback, ARTIQ notes, reflash | 2.1-2.13 |
| Rules, scope control, multi-site edit method | 3.2-3.4 |
| File layout, Explorer, porting checklist | 4.4, 4.5, 4.7 |
| Mapping JSON, resolver, SPCM binding | 5.4-5.8 |
| Base: modes, idle node, init order, result naming | 6 |
| ExperimentVariables: globals, timing layers, pitfalls | 7.3-7.8 |
| GVS, overrides and their recording; experiment functions, confirm loop | 8, 9 |
| Feedback (single, two-node, parallel design); optimizers, AOMsCoils | 10, 11 |
| DRTIO, skew probe, slack, SED, `with parallel`, DMA | 12 |
| Host devices; applets and dashboard; analysis | 13, 14, 15 |
| Tests, compile sweep, controls; status and history | 16, 17 |
| Open issues O-1 to O-36, next steps | 18 |
| Decisions and superseded claims; reporting contract | 19, 20 |
