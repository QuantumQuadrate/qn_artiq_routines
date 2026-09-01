# Master-Satellite Migration: AI Startup Brief

## Purpose

This repository controls two neutral Rb-87 single-atom quantum-network nodes
using ARTIQ/Sinara. The project is migrating two independent Kasli-SoC systems
to one DRTIO master-satellite system without sacrificing the working
standalone system.

- Node1 = Alice = future DRTIO master at `192.168.1.75`
- Node2 = Bob = future DRTIO satellite, destination 1
- Current standalone addresses are `.75` and `.76`
- The new system has one canonical core and one unified master-side
  `device_db.py`

The motivation is coordinated two-node timing and control through DRTIO while
retaining known-good single-node physics and an easy hardware rollback.

## Highest-level rule

Preserve the standalone architecture and add a separate master-satellite
architecture beside it. Do not gradually turn the standalone files into the
new architecture.

```text
Standalone                         Master-satellite
GeneralVariableScan.py             GeneralVariableScan_master_satellite_mixin.py
ExperimentVariables.py             GeneralVariableScan_master_satellite_single_node.py
BaseExperiment.py                  GeneralVariableScan_master_satellite_two_nodes.py
DeviceAliases.py                   ExperimentVariables_master_satellite_Node1.py
experiment_functions.py            ExperimentVariables_master_satellite_Node2.py
AOMsCoils.py                       ExperimentVariables_master_satellite_global.py
                                   BaseExperiment_master_satellite.py
                                   DeviceAliases_master_satellite.py
                                   experiment_functions_two_nodes.py
                                   AOMsCoils_master_satellite_mixin.py
                                   AOMsCoils_master_satellite_Node1.py
                                   AOMsCoils_master_satellite_Node2.py
```

`_mixin` files hold only private shared implementations and are not Explorer
experiments. Execution mode is fixed per public GVS file; `selected_node` is
the single-node submitted argument. Both standalone and master-satellite GVS
carry optional per-point underflow retry behind `enable_Catch_UnderFlow`
(off by default).

Standalone node selections remain only `alice` and `bob`. The deleted legacy
`which_node == "two_nodes"` mode was not DRTIO and must never return. The
standalone `GeneralVariableScan.py` absorbed
`GeneralVariableScan_CatchUnderflow.py`: per-iteration underflow retry is its
optional `enable_Catch_UnderFlow` argument (off by default).

## Master-satellite execution model

```text
experiment_mode = single_node
    selected_node = Node1 or Node2

experiment_mode = two_nodes
    both nodes active
```

`selected_node` is the authoritative identity. Reused legacy single-node
functions receive `self.which_node = "alice"` or `"bob"` only as a
compatibility presentation.

Node-specific logical attributes and persistent datasets use suffixes:

```text
dds_FORT_Node1, dds_FORT_Node2
f_FORT_Node1, f_FORT_Node2
sampler0_Node1, sampler0_Node2
```

Single-node mode projects the selected authoritative value/device to the old
unsuffixed namespace. The projection is never a second persistent source.

## Essential invariants

1. Physical resolution happens before logical presentation.
2. A logical device on a given node resolves to the same physical device in
   single-node and two-node modes.
3. Physical mappings are exact dictionary lookups from
   `utilities/config/master_satellite/device_aliases.json`; no calculated
   offsets are allowed.
4. Encoded satellite RTIO channels come only from the generated database.
5. There is one canonical core.
6. Node2 hardware is not accessed until destination 1 reports ready.
7. Existing persistent calibration values are never reset by code defaults.
8. Scans and overrides change authoritative in-memory attributes but do not
   persist calibration automatically.
9. Historical independent-Kasli functions remain untouched and are not native
   master-satellite functions.
10. Hardware validation proceeds from passive link tests toward physics.
11. A submitted experiment's own GUI argument values always win over
    persistent datasets for that run.

## Special SPCM rule

All active detectors are master-local Node1 inputs:

```text
SPCM_H1 -> ttl0 / ttl0_counter
SPCM_V1 -> ttl1 / ttl1_counter
SPCM_H2 -> ttl8 / ttl8_counter
SPCM_V2 -> ttl9 / ttl9_counter
```

Even in single Node2 mode, legacy SPCM aliases bind to these same Node1
devices. H1/V1/H2/V2 are detector identities, not Kasli ownership.

## Current software status

Implemented and hardware-free tested:

- explicit physical translation;
- master-satellite resolver and Base;
- fixed SPCM policy;
- Node1/Node2/global ExperimentVariables;
- authoritative variable loading, projection, reload, and cache refresh;
- mode-specific GVS registries and scan resolution;
- optional built-in underflow retry (`enable_Catch_UnderFlow`) in both GVS
  modes;
- run-local GUI-argument precedence over persistent datasets
  (`n_measurements`);
- native namespace sanity function;
- minimal magnetometer compatibility and node-suffixed results;
- deterministic per-node manual AOM/DDS, optical-gate, microwave/RF, and
  `disable_coils` MOT-coil control through
  `AOMsCoils_master_satellite_Node1/Node2.py`;
- single-node laser-feedback compatibility (`run_laser_feedback`) and K10CR1
  waveplate control (one `k10cr1_ndsp` controller, node-suffixed axis names);
- the single-node FORT polarization optimizer with node-suffixed waveplate
  calibration (the two-node variant is documented in the file and waits on
  parallel FORT feedback);
- derived RF amplitudes and the single-node atom-result state in Base;
- `MicrowaveScanOptimizer_master_satellite` (`which_node = Node1/Node2`,
  never two_nodes): suffixed microwave calibrations and `health_check_uw_*`
  fidelities, scan definitions read from the untouched standalone source,
  self-resubmission carries `which_node`; the standalone optimizer and
  HealthCheck files stay untouched, and future master-satellite HealthCheck
  files must submit the new class with `which_node`.

The manual utilities are node-split. Each public experiment binds the two-node
device superset during build, then configures Base as `single_node` for its
fixed node during prepare and loads only that node's variables plus globals.
Only the selected node is initialized and driven: a second explicit OFF/0 V
pass precedes the requested state, and Base initialization runs with
`reset_core=False` so the other node's established outputs are never touched.
Making a node safe means running that node's own experiment with its controls
unchecked. Coils follow the standalone `disable_coils` semantics: unless
disabled, the four coils are driven to the node's persistent MOT calibration
voltages; nothing is entered manually or persisted. Microwave/RF defaults OFF
and requires confirmation. K10CR1 780/852 rotations act on the selected
node's axes through the single `k10cr1_ndsp` controller (node-suffixed
nicknames, bound only when a waveplate action is selected).
`run_laser_feedback` matches the standalone utility (feedback, or monitor
when unticked, only with all six fiber AOMs plus cooling DP on) and persists
through node-suffixed datasets via unchanged aom_feedback code. Node2 also
carries the standalone combined `Node2_GRIN1/GRIN2_AOM_ON` modes, applied
last and authoritative over the individual GRIN controls. Old
independent-Kasli networking and Rigol are excluded.

ARTIQ repository examination supplies `None` for submitted arguments and may
run before any new persistent datasets exist. Every public master-satellite
`build()` must therefore succeed with both conditions. GVS and AOMsCoils bind
the suffixed Node1/Node2 device superset in build, but defer strict mode
validation, authoritative dataset loading, selected-node projection, and DDS
preparation to prepare. Runtime missing datasets still fail before hardware.

Explorer labels come from the first class-docstring line and must match the
recognizable filename. Imported `EnvExperiment` subclasses must use private
module aliases; otherwise ARTIQ can rediscover them as duplicate experiments.
Shared implementations live in `_mixin` files that define no public
experiment class.

The software is ready for controlled gateware integration testing, not an
immediate physical upgrade.

## Current hardware gate

Both physical nodes still run standalone gateware. Before flashing:

- make verified rollback images of both SD cards;
- preserve the actual active standalone device databases and configuration;
- confirm carrier/peripheral revisions;
- identify the exact ARTIQ-Zynq revision;
- build/obtain matching master and satellite `boot.bin` files;
- generate and verify the matching unified database;
- establish the master SFP port, satellite upstream port, cable/transceivers,
  and routing table;
- archive all images, JSONs, build logs, revisions, and checksums together.

Do not run the current GVS as the first post-flash check. Its Base
initialization changes DDS, TTL, Sampler, and Zotino state.

## Safe next sequence

1. Secure rollback artifacts.
2. Prepare an exact gateware/database release bundle.
3. Boot and inspect master/satellite logs.
4. Run a passive destination-1 readiness kernel.
5. Verify namespace resolution without initializing outputs.
6. Test one harmless local and one harmless remote RTIO operation.
7. Use `AOMsCoils_master_satellite_Node1/Node2.py` to verify one AOM and coil
   at a time on the corresponding node, beginning with all outputs OFF/0 V.
8. Run `test_magnetometer_experiment` on Node1, then Node2.
9. Add atom-loading/feedback compatibility.
10. Validate the master-satellite microwave optimizer on hardware.
11. Validate DMA and timing-critical single-photon work last.

## How a new AI should work

Read `plan_codex_detail.md` before implementation. Inspect first, state scope,
and distinguish proven facts from assumptions. Do not modify hardware, commit,
push, or broaden scope without explicit authorization. After changes, run
tests, inspect the diff, report hardware validation still needed, and confirm
that standalone files were untouched.
