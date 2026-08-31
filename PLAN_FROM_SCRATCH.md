# Master-Satellite Migration Plan

## Objective

Migrate two independent ARTIQ Kasli-SoCs to one DRTIO master-satellite system
while preserving the validated standalone path and a reliable rollback path.

- Node1 = Alice = future DRTIO master at `192.168.1.75`
- Node2 = Bob = future DRTIO satellite, destination 1
- Current standalone addresses: Node1 `.75`, Node2 `.76`
- One canonical master core and one unified master-side `device_db.py`

Correctness, hardware safety, and rollback take priority over cleanup.

## Operating rules

1. Read this file completely before acting.
2. Inspect first; do not infer mappings or hardware state.
3. Keep the existing standalone stack unchanged unless explicitly authorized.
4. Do not alter historical experiments merely because they mention two nodes.
5. Do not flash hardware or operate outputs without explicit authorization.
6. Do not commit or push unless explicitly requested.
7. Preserve unrelated worktree changes.
8. Use exact generated RTIO channels; never calculate satellite offsets.
9. Separate repository evidence from assumptions requiring hardware tests.
10. Run hardware-free tests after each software phase.

## Architectures that must coexist

Preserve the standalone stack:

```text
GeneralVariableScan.py
ExperimentVariables.py
utilities/BaseExperiment.py
utilities/DeviceAliases.py
subroutines/experiment_functions.py
subroutines/aom_feedback.py
utilities/config/alice/
utilities/config/bob/
```

Standalone selections remain only `alice` and `bob`. The removed legacy
`which_node == "two_nodes"` mode must not return.

Use a separate master-satellite stack:

```text
GeneralVariableScan_master_satellite.py
GeneralVariableScan_CatchError_master_satellite.py
AOMsCoils_master_satellite.py
ExperimentVariables_Node1.py
ExperimentVariables_Node2.py
ExperimentVariables_master_satellite.py
utilities/BaseExperiment_master_satellite.py
utilities/DeviceAliases_master_satellite.py
subroutines/experiment_functions_master_satellite.py
utilities/config/master_satellite/device_aliases.json
```

Execution modes are `single_node` with `selected_node = Node1/Node2`, and
`two_nodes`. `selected_node` is authoritative. Only reused single-node
functions receive the compatibility presentation `Node1 -> alice` and
`Node2 -> bob` through `self.which_node`.

## Hardware reference

Treat these as software mapping references until exact matching images and
physical hardware have been verified:

```text
master_satellite_test_files/node1_master_kasli_soc_20260715.json
master_satellite_test_files/node2_satellite_kasli_soc_20260715.json
master_satellite_test_files/device_db_two_nodes_20260715.py
```

Both nodes have two DIO cards with edge counters, three Samplers, one Zotino,
and three AD9910 Urukuls. The unified database assigns Node1 first:

```text
Node1: ttl0-15, sampler0-2, zotino0, urukul0-2
Node2: ttl16-31, sampler3-5, zotino1, urukul3-5
```

Sequential Python names do not imply sequential encoded RTIO channels.

## Physical resolution invariant

Resolve physical identity before choosing the experiment-facing name. For a
given node and logical device, single-node and two-node modes must select the
same unified physical device.

```text
Bob dds_FORT -> standalone urukul1_ch3 -> unified urukul4_ch3

single Node2: self.dds_FORT       -> urukul4_ch3
two nodes:    self.dds_FORT_Node2 -> urukul4_ch3
```

`utilities/config/master_satellite/device_aliases.json` is the single source
of physical translations. Every mapping is an exact dictionary entry; never
derive Urukul, TTL, destination, or RTIO offsets.

## Fixed SPCM policy

Active photon detectors are always the master-local Node1 inputs:

```text
SPCM_H1 -> ttl0 / ttl0_counter
SPCM_V1 -> ttl1 / ttl1_counter
SPCM_H2 -> ttl8 / ttl8_counter
SPCM_V2 -> ttl9 / ttl9_counter
```

In single-node compatibility mode, for either selected node:

```text
ttl_SPCM0           -> SPCM_H1
ttl_SPCM1           -> SPCM_V1
ttl_SPCM0_OtherNode -> SPCM_H2
ttl_SPCM1_OtherNode -> SPCM_V2
```

Counter aliases follow the same policy. Do not route these aliases through
the selected node's generic TTL translation. Node2's old SPCM-looking inputs
remain physical devices but are not active detectors.

## Experiment-variable ownership

Node-calibratable variables are independent persistent datasets named
`<name>_Node1` and `<name>_Node2`.

```text
ExperimentVariables_Node1.py             owns *_Node1
ExperimentVariables_Node2.py             owns *_Node2
ExperimentVariables_master_satellite.py  owns globals only
```

Approved globals include unsuffixed `n_measurements`, `t_delay_in_bob_mu`, and
`parallel_AOM_feedback`. Execution selections are GVS arguments, not datasets.

Code defaults seed missing datasets only. Running an ExperimentVariables file
must not overwrite an existing persistent calibration. Normal experiments
must fail clearly when a required authoritative dataset is missing.

Single-node compatibility is an in-memory projection, never a second dataset:

```text
f_FORT_Node2 dataset -> self.f_FORT_Node2 -> self.f_FORT
```

## Responsibility boundary

`DeviceAliasesMasterSatellite` owns exact resolution, device binding, DDS
two-stage resolution, caller-selected presentation names, and per-node DDS
bookkeeping. It must not reset core, wait for DRTIO, initialize hardware,
choose modes, implement SPCM policy, or handle USB devices.

`BaseExperimentMasterSatellite` owns mode-dependent variable loading,
compatibility projection, canonical core lifecycle, satellite readiness,
SPCM policy, static wiring metadata, CPLD/DDS/TTL/Sampler/Zotino initialization,
safe output state, result state, scan-target resolution, variable reload, and
dependent-cache refresh.

The readiness-polling argument is `satellite_ready_poll_interval`; it is not
the physical DRTIO communication latency.

## GVS rules

The master-satellite GVS must:

- use mode-specific callable registries instead of `eval()` for functions;
- reuse locally defined legacy functions only in `single_node` mode;
- exclude historical independent-Kasli two-node functions;
- use only `experiment_functions_master_satellite.py` in `two_nodes` mode;
- delegate scan/override target resolution to Base;
- mutate authoritative attributes, then refresh projections and caches;
- keep scan points and overrides non-persistent;
- reload queued persistent values without rebuilding devices;
- never call `build()` or `prepare()` per scan point.

The CatchError variant inherits this behavior and retries only a failed scan
point on `RTIOUnderflow`. It must not hide DRTIO, SPI, mapping, or other errors.

## Manual hardware control

`AOMsCoils_master_satellite.py` is separate from GVS and from the untouched
standalone utility. It always binds Base in `two_nodes` mode so an unselected
node can be actively forced safe. Its own selector is:

```text
which_node = node1 | node2 | two_nodes
```

For `node1`, all controlled Node2 DDS outputs are forced OFF and all Node2
coils are written to 0 V; `node2` is symmetric; `two_nodes` respects both
nodes independently. It controls FORT, cooling DP, AOM A1-A6, D1 pumping,
GRIN1and2, active-low optical gates, microwave/RF, and four independent coils
per node. DDS values come only from authoritative node datasets. Node2 D1 uses
the GRIN2 D1 frequency/power variables.

Coil voltages are run-local GUI values seeded from node MOT calibration and
never persisted. Explicit channels are Node1 `[0, 1, 13, 14]` and Node2
`[0, 1, 2, 3]`. OFF always writes 0 V. Microwave/RF defaults OFF and requires
confirmation. Feedback, old networking TTLs, K10CR1, Rigol, and other USB
devices are excluded.

Safe startup switches DDS outputs OFF before attenuation/profile programming,
zeroes both Zotinos, performs another explicit controlled-output OFF/0-V pass,
then applies the hard-gated requested state. Commit `078dc0e` implemented this
with 51 passing combined hardware-free tests.

## Hardware transition gate

The physical nodes still run standalone gateware. Before any flash:

1. Make and verify byte-for-byte backups of both current SD cards, or obtain
   the exact standalone `boot.bin` for each node.
2. Preserve the experiment computer's actual standalone device databases,
   configuration, startup/idle kernels, ARTIQ version, and boot logs.
3. Confirm Kasli-SoC and peripheral hardware revisions.
4. Establish the exact ARTIQ/ARTIQ-Zynq revision for both new builds.
5. Build or obtain matching master and satellite `boot.bin` files.
6. Generate the unified database from those exact JSON inputs with Node2 at
   destination 1.
7. Determine the master SFP port, satellite upstream port, transceivers,
   cable, and routing table.
8. Archive JSONs, images, database, routing data, build logs, revisions, and
   checksums as one release bundle.
9. Verify the rollback media where practical.

The ARTIQ-Zynq `boot` key is write-only. Do not assume the current boot image
can be recovered later through `artiq_coremgmt`.

## First hardware validation

Do not begin with the current GVS or a physics experiment:

1. Node1 master boot/network sanity.
2. UART/core-log inspection on both nodes.
3. A bounded passive kernel checking only
   `core.get_rtio_destination_status(1)`.
4. Host-side unified namespace resolution with no hardware initialization.
5. One harmless Node1 local RTIO operation.
6. One harmless/read-only Node2 remote RTIO operation.
7. Manual control with all outputs OFF/0 V, then one AOM/coil at a time on
   Node1 followed by Node2.
8. `test_magnetometer_experiment`, single Node1, `n_measurements = 1`.
9. The same magnetometer test for single Node2.

Current Base initialization is not passive: it resets core, changes TTL
directions/states, initializes CPLDs and DDSs, sets DDS attenuation/frequency/
amplitude/switch state, changes Sampler gains, and can write Zotino channels to
zero. Confirm all initialization states are electrically safe before GVS use.

## Deferred validation

Handle these only after basic DRTIO/device access succeeds:

- AD9910 synchronization recalibration;
- satellite-involving DMA traces;
- remote SPI latency and SED-lane capacity;
- AOM/FORT feedback and parallel feedback;
- microwave optimizer compatibility;
- K10CR1 and other host-controlled devices;
- photon timestamps and herald/mapping timing.

Do not claim the standalone `MicrowaveScanOptimizer.py` is compatible. A
separate master-satellite optimizer is required because the old one uses the
standalone Base, unsuffixed calibration targets, legacy feedback/result state,
repeated `prepare()`, and unvalidated DMA traces.

## Required reporting

Before editing, report inspected files, established facts, ambiguity, risk,
and exact proposed scope. After editing, report changed files, APIs, tests,
hardware validation still needed, `git diff --stat`, and confirmation that the
standalone stack was untouched.

Never commit or push automatically.

## Recommended implementation order

1. Read-only repository/hardware audit.
2. Approve explicit physical translation.
3. Implement and test resolver.
4. Design, implement, and test Base lifecycle.
5. Classify and implement node/global ExperimentVariables.
6. Implement authoritative projection and refresh APIs.
7. Audit function discovery and implement GVS registries.
8. Add a native namespace sanity function.
9. Add legacy compatibility narrowly, one target experiment at a time.
10. Prepare rollback artifacts and an exact gateware release bundle.
11. Perform passive DRTIO validation.
12. Validate manual Node1 then Node2 AOM/coil access.
13. Validate Node1 then Node2 magnetometer access.
14. Add atom-loading/feedback compatibility.
15. Add a separate master-satellite microwave optimizer.
16. Validate single-photon and timing-critical operation last.
