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
GeneralVariableScan.py             GeneralVariableScan_master_satellite.py
ExperimentVariables.py             ExperimentVariables_Node1.py
BaseExperiment.py                  ExperimentVariables_Node2.py
DeviceAliases.py                   ExperimentVariables_master_satellite.py
experiment_functions.py            BaseExperiment_master_satellite.py
                                    DeviceAliases_master_satellite.py
                                    experiment_functions_master_satellite.py
AOMsCoils.py                        AOMsCoils_master_satellite.py
```

Standalone node selections remain only `alice` and `bob`. The deleted legacy
`which_node == "two_nodes"` mode was not DRTIO and must never return.

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
- master-satellite CatchError GVS;
- native namespace sanity function;
- minimal magnetometer compatibility and node-suffixed results;
- deterministic two-node manual AOM/DDS, optical-gate, microwave/RF, and
  independent coil control through `AOMsCoils_master_satellite.py`.

The manual utility always binds Base in `two_nodes` mode so an unselected node
can be actively made safe. Its selector is `which_node = node1 / node2 /
two_nodes`. Selecting one node respects only that node's checkboxes; all
controlled DDS outputs on the other node are switched OFF and all four of its
coils are explicitly written to 0 V. Coil voltages are run-local GUI values
seeded from node-specific MOT calibration and are never persisted. Microwave/
RF defaults OFF and requires confirmation. Feedback, old independent-Kasli
networking, K10CR1, and Rigol are excluded.

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
7. Use `AOMsCoils_master_satellite.py` to verify one Node1/Node2 AOM and coil
   at a time, beginning with all outputs OFF/0 V.
8. Run `test_magnetometer_experiment` on Node1, then Node2.
9. Add atom-loading/feedback compatibility.
10. Create a separate master-satellite microwave optimizer.
11. Validate DMA and timing-critical single-photon work last.

## How a new AI should work

Read `plan_codex_detail.md` before implementation. Inspect first, state scope,
and distinguish proven facts from assumptions. Do not modify hardware, commit,
push, or broaden scope without explicit authorization. After changes, run
tests, inspect the diff, report hardware validation still needed, and confirm
that standalone files were untouched.
