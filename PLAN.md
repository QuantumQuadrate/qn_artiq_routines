# AGENTS.md

# Master-Satellite Migration Context

## 1. Project Overview

This repository controls two neutral Rb-87 single-atom quantum-network nodes using ARTIQ/Sinara.

Physical node naming:

- Node1 = Alice
- Node2 = Bob

The current goal is to migrate the two previously independent Kasli systems into a true ARTIQ master-satellite architecture while preserving the currently working standalone single-node system.

The most important constraint is:

> Preserve the existing standalone experiment stack and add the master-satellite stack alongside it.

This migration should prioritize:

1. working hardware control,
2. easy debugging,
3. minimal-risk changes,
4. preserving validated experiment sequences,
5. easy rollback to the old standalone hardware/software path,

over large-scale architectural cleanup.

Do not redesign unrelated parts of the repository simply because a cleaner abstraction is possible.

---

## 2. Development Environment

Migration development is being performed on a separate Git branch:

`20260826_master_satellite_codex`

Most code inspection, restructuring, and static development may be done on a non-experiment computer.

Actual ARTIQ hardware validation must be performed on the experiment computer.

Therefore:

- do not assume ARTIQ hardware is available during local development,
- do not assume successful Python/static checks imply hardware correctness,
- explicitly identify changes that require experiment-system validation,
- do not commit, push, merge, or modify Git history unless explicitly requested.

---

## 3. Two Hardware Architectures Must Coexist

The repository must support two distinct hardware architectures.

### 3.1 Existing standalone architecture

This is the currently validated path.

It uses the existing files:

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

The valid standalone node selections are:

```text
alice
bob
```

Mapping:

```text
alice -> Node1
bob   -> Node2
```

The existing standalone stack must remain usable with the old standalone gateware and old standalone `device_db.py`.

Do not modify this path unless a change is explicitly required and proven safe.

The migration goal is NOT to turn all existing files into generic master-satellite files.

### 3.2 New master-satellite architecture

The new master-satellite path should be implemented alongside the standalone stack.

Preferred new files:

```text
GeneralVariableScan_master_satellite.py
utilities/BaseExperiment_master_satellite.py
utilities/DeviceAliases_master_satellite.py
ExperimentVariables_Node1.py
ExperimentVariables_Node2.py
ExperimentVariables_master_satellite.py
subroutines/experiment_functions_master_satellite.py
utilities/config/master_satellite/device_aliases.json
```

The master-satellite path uses:

- one master Kasli,
- one satellite Kasli,
- one canonical ARTIQ core,
- one unified master-side `device_db.py`,
- one shared physical mapping from the old standalone namespaces to the unified device namespace.

The standalone stack and master-satellite stack must remain clearly separated.

---

## 4. Rollback / Compatibility Principle

The old standalone stack should remain directly usable if the experiment hardware is returned to the old standalone gateware/device database.

Conceptually:

```text
OLD STANDALONE HARDWARE
    old gateware
    old standalone device_db
            |
            v
    GeneralVariableScan.py
            |
            v
    existing BaseExperiment / DeviceAliases / experiment_functions
```

The new path is:

```text
MASTER-SATELLITE HARDWARE
    master boot image
    satellite boot image
    unified device_db
            |
            v
    GeneralVariableScan_master_satellite.py
            |
            +----------------------+
            |                      |
            v                      v
       single_node              two_nodes
       Node1/Node2              Node1+Node2
```

Do not create unnecessary dependencies from the old standalone path onto the new master-satellite path.

---

## 5. Legacy `which_node == "two_nodes"` Mode

The old node-selection value:

```text
which_node == "two_nodes"
```

is obsolete and should remain removed.

It was not a true master-satellite implementation.

It must not be reintroduced.

Only:

```text
alice
bob
```

are valid for the old standalone `BaseExperiment.py` path.

The new master-satellite architecture may use an experiment-mode value such as:

```text
single_node
two_nodes
```

inside `GeneralVariableScan_master_satellite.py`, but this is a completely different concept from the removed legacy:

```text
which_node == "two_nodes"
```

Do not confuse the two.

---

## 6. Node Naming Convention

For NEW master-satellite logical names, use the node name as a suffix.

Preferred pattern:

```text
<existing_logical_name>_Node1
<existing_logical_name>_Node2
```

Examples:

```text
dds_FORT_Node1
dds_FORT_Node2

dds_cooling_DP_Node1
dds_cooling_DP_Node2

ttl_repump_switch_Node1
ttl_repump_switch_Node2

sampler0_Node1
sampler0_Node2

zotino0_Node1
zotino0_Node2

f_FORT_Node1
f_FORT_Node2

p_FORT_loading_Node1
p_FORT_loading_Node2

FORT_monitor_Node1
FORT_monitor_Node2
```

Do NOT prefer the opposite style:

```text
dds_Node1_FORT
Node1_f_FORT
```

Do not rename validated legacy standalone code merely to enforce this convention.

The suffix convention primarily applies to:

- new master-satellite logical experiment attributes,
- master-satellite node-specific datasets,
- new USB/host device nicknames,
- new master-satellite infrastructure.

---

## 7. One Unified Master-Side `device_db.py`

The master-satellite system uses one unified master-side `device_db.py`.

There is not a separate runtime `device_db.py` for Node1 and Node2.

The unified device database is assumed to expose both local master devices and remote satellite devices.

Important rules:

1. Keep one unified master-side `device_db.py`.
2. Keep one canonical ARTIQ `core`.
3. Use the RTIO channel/address data generated for the master-satellite system.
4. Do not manually derive satellite RTIO channels.
5. Do not invent a separate destination argument for normal device access.
6. Sequential Python device names are acceptable, but their RTIO channel integers must not be assumed to be flat sequential integers.

Example:

```text
urukul0_ch0
urukul3_ch0
```

may be sequential-looking Python names, but do not infer RTIO addressing from those names.

The generated gateware/database information is the source of truth for encoded RTIO channels and DRTIO destination information.

---

## 8. Authoritative Master-Satellite Hardware Reference

A previously generated master-satellite hardware configuration is stored at:

`master_satellite_test_files/`

Reference files:

```text
master_satellite_test_files/device_db_two_nodes_20260715.py
master_satellite_test_files/node1_master_kasli_soc_20260715.json
master_satellite_test_files/node2_satellite_kasli_soc_20260715.json
```

For current software development, assume this hardware configuration and unified device database are correct unless explicitly told otherwise.

Current assumptions:

- Node1 = master
- Node2 = satellite
- Node2 uses DRTIO destination 1
- edge counters are present on both nodes

Do not modify these reference files unless explicitly requested.

Hardware validity will be verified later on the experiment system.

---

## 9. Authoritative ARTIQ Reference

This repository currently targets ARTIQ 7.

The authoritative ARTIQ reference for this migration is:

`documentation/ARTIQ manual 7.82.pdf`

For questions involving:

- DRTIO,
- device database behavior,
- RTIO channel mapping,
- DMA,
- core-device behavior,
- satellite readiness,
- SED lanes,
- AD9910 synchronization,

consult this manual instead of guessing from general knowledge.

Relevant areas include:

- ARTIQ Real-Time I/O concepts,
- DMA,
- the device database,
- DRTIO,
- core-device drivers,
- scheduling / SED behavior.

---

## 10. Master-Satellite ARTIQ Constraints

The following conclusions are authoritative for this migration:

1. Keep one unified master-side `device_db.py`.
2. Keep one canonical `core`.
3. Let generated gateware/database information determine encoded RTIO channels.
4. Do not invent a separate destination argument for normal device access.
5. Sequential Python device names are acceptable, but do not assume flat sequential RTIO channel integers.
6. Apply `_Node1` / `_Node2` primarily to logical experiment attributes and globally shared datasets.
7. Explicitly wait for the satellite destination to become available before initializing satellite hardware.
8. Treat remote input and remote SPI operations as latency-sensitive.
9. Revalidate all DMA sequences that include satellite hardware.
10. Recalibrate AD9910 synchronization after final gateware and clocking are finalized.
11. Inspect SED lane configuration before implementing broad parallel RTIO operations.
12. Preserve edge counters on both nodes.

These points should constrain the software architecture.

---

## 11. Physical Device Names vs Logical Experiment Names

Keep a clear distinction between:

1. low-level physical names in the unified `device_db.py`,
2. experiment-facing logical names.

For example:

```text
unified device_db physical names:
    urukul0_ch0
    urukul3_ch0
```

may be presented to experiments as:

```text
dds_FORT_Node1
dds_FORT_Node2
```

Do not assume that every low-level `device_db.py` entry must itself be renamed with `_Node1` / `_Node2`.

Prefer keeping the generated low-level database close to the generated hardware description and applying human-readable node naming in the logical mapping layer.

---

## 12. Existing Alice/Bob `device_aliases.json` Files

The existing files:

```text
utilities/config/alice/device_aliases.json
utilities/config/bob/device_aliases.json
```

remain the source of truth for existing logical DDS aliases within each standalone node.

Conceptually:

```text
logical DDS alias
    ->
standalone node physical DDS name
```

Example:

```text
dds_FORT
    ->
urukul1_ch3
```

depending on the node configuration.

Do not duplicate these logical DDS definitions in the new master-satellite mapping.

---

## 13. New Master-Satellite `device_aliases.json`

Create/use:

```text
utilities/config/master_satellite/device_aliases.json
```

This file is the single source of truth for translating each node's old standalone physical ARTIQ namespace into the unified master-satellite `device_db.py` namespace.

Conceptually:

```text
standalone Node1 physical name
    ->
unified master-satellite physical name

standalone Node2 physical name
    ->
unified master-satellite physical name
```

The file may contain all node-specific ARTIQ physical-device translations required by the master-satellite software, including:

- DDS / Urukul channels,
- Urukul CPLDs,
- Samplers,
- Zotino,
- TTL devices,
- edge counters,
- other node-specific ARTIQ devices required by BaseExperiment logic.

Do NOT create a second duplicated `device_map.json`.

Do NOT encode translation rules in Python using arithmetic such as:

```text
urukul index + 3
ttl index + 16
RTIO channel + offset
```

The mapping must be explicit data derived from the unified reference `device_db.py`.

Do NOT include the canonical `core` as a per-node mapping.

USB/host-controlled devices are not required to live in this ARTIQ physical-device mapping.

---

## 14. Physical Mapping Is Shared by Single-Node and Two-Node Master-Satellite Modes

### Single physical resolution invariant

For a given node and logical device, physical-device resolution MUST be identical
between master-satellite single-node mode and true two-node mode.

The execution mode must affect only the experiment-facing attribute name, not
the physical device selected.

For example, if Node1 FORT resolves to `urukul0_ch1` and Node2 FORT resolves to
`urukul4_ch3`, then:

two_nodes mode:
    self.dds_FORT_Node1 -> urukul0_ch1
    self.dds_FORT_Node2 -> urukul4_ch3

single_node mode, Node1 selected:
    self.dds_FORT -> urukul0_ch1

single_node mode, Node2 selected:
    self.dds_FORT -> urukul4_ch3

Do NOT maintain a separate physical mapping for single-node compatibility mode.

The implementation should first resolve the node-specific physical device from
the shared master-satellite mapping and only then choose how that device is
presented to the experiment namespace.

There must be exactly one master-satellite physical translation source.

Do NOT maintain separate physical translation tables for:

- master-satellite single-node mode,
- master-satellite two-node mode.

Both modes consume:

```text
utilities/config/master_satellite/device_aliases.json
```

The difference is only the logical presentation.

### 14.1 Master-satellite single-node compatibility

When running an existing single-node experiment under the new master-satellite gateware, select either Node1 or Node2 and expose the selected physical node using the old unsuffixed logical namespace.

Examples:

```text
Node1 physical FORT DDS -> self.dds_FORT
Node2 physical FORT DDS -> self.dds_FORT
```

Likewise:

```text
selected node sampler0 -> self.sampler0
selected node zotino0  -> self.zotino0
selected node TTL       -> existing unsuffixed logical name
```

This is required so existing single-node experiment functions can continue to run on the new master-satellite hardware without broad rewrites.

Even a Node2-only experiment uses the one canonical master `core` and accesses Node2 hardware through DRTIO.

### 14.2 True master-satellite two-node mode

When running a true two-node master-satellite experiment, expose both nodes simultaneously with suffixed logical names.

Examples:

```text
dds_FORT_Node1
dds_FORT_Node2

sampler0_Node1
sampler0_Node2

zotino0_Node1
zotino0_Node2
```

Thus:

> Physical resolution is shared; logical presentation depends on execution mode.

---

## 15. Master-Satellite Device Resolver

Keep the existing:

`utilities/DeviceAliases.py`

as the standalone implementation.

Create a new:

`utilities/DeviceAliases_master_satellite.py`

for the master-satellite hardware architecture.

Do not force the existing standalone `DeviceAliases.py` to become the master-satellite implementation unless a very small compatibility change is explicitly required.

The new resolver should:

1. consume the existing Alice/Bob logical DDS configuration,
2. consume `utilities/config/master_satellite/device_aliases.json`,
3. resolve old standalone physical names to unified physical device names,
4. support a legacy unsuffixed presentation for master-satellite single-node mode,
5. support `_Node1` / `_Node2` presentation for true two-node mode,
6. maintain node-specific DDS bookkeeping without collisions,
7. resolve node-specific DDS default variables without collisions,
8. fail clearly if a required mapping is missing,
9. avoid resetting the canonical core independently for each node,
10. avoid assuming fixed Urukul/TTL/RTIO offsets.

The exact class API should be chosen after inspecting the current implementation.

---

## 16. DDS Default Variable Names

The existing Alice/Bob DDS configs refer to unsuffixed experiment-variable names such as:

```text
f_FORT
p_FORT_loading
f_cooling_DP_MOT
p_cooling_DP_MOT
```

In master-satellite mode, node-specific versions must remain distinguishable.

Examples:

```text
f_FORT_Node1
f_FORT_Node2

p_FORT_loading_Node1
p_FORT_loading_Node2
```

The new master-satellite resolver/API must account for this.

Do not assume DDS alias resolution can be separated entirely from node-specific default-variable resolution.

---

## 17. BaseExperiment Architecture

Keep:

`utilities/BaseExperiment.py`

for the old standalone hardware architecture.

Create:

`utilities/BaseExperiment_master_satellite.py`

for the new master-satellite hardware architecture.

Do NOT make the old `BaseExperiment.py` the primary master-satellite implementation.

The new master-satellite BaseExperiment should support two execution modes:

```text
single_node
two_nodes
```

### 17.1 Master-satellite single-node mode

Select:

```text
Node1
or
Node2
```

and expose the selected node through the old unsuffixed logical namespace so existing experiment functions can be reused where possible.

### 17.2 Master-satellite two-node mode

Expose Node1 and Node2 simultaneously through suffixed logical names.

### 17.3 Responsibilities

`BaseExperiment_master_satellite.py` should own or coordinate:

- the one canonical core,
- satellite-readiness checks/waits,
- unified physical device registration,
- selected execution mode,
- logical namespace presentation,
- non-aliased device registration,
- hardware initialization ordering,
- core reset exactly once where appropriate,
- CPLD initialization,
- Sampler/Zotino/TTL exposure,
- integration with the master-satellite device resolver.

Do not make `DeviceAliases_master_satellite.py` independently reset the canonical core once per node.

---

## 18. ExperimentVariables Architecture

Keep:

`ExperimentVariables.py`

for the old standalone stack.

For the master-satellite stack, preferred files are:

```text
ExperimentVariables_Node1.py
ExperimentVariables_Node2.py
ExperimentVariables_master_satellite.py
```

### 18.1 Node-specific files

Node1 and Node2 files contain node-specific:

- hardware settings,
- calibrations,
- node-specific experimental parameters.

Persistent dataset names must be globally unique.

Examples:

```text
f_FORT_Node1
f_FORT_Node2

p_FORT_loading_Node1
p_FORT_loading_Node2

A_X_Node1
A_X_Node2
```

Do not allow Node1 and Node2 settings to overwrite each other in one master dataset database.

### 18.2 `ExperimentVariables_master_satellite.py`

This file should contain parameters specific to:

- master-satellite operation,
- simultaneous two-node behavior,
- cross-node timing,
- synchronization,
- heralding,
- global master-satellite settings,
- feedback execution mode.

Do not duplicate Node1/Node2 calibration values here.

Include a boolean such as:

```python
parallel_AOM_feedback = True
```

where:

```text
False -> sequential Node1 / Node2 feedback
True  -> concurrent Node1 / Node2 feedback
```

---

## 19. GeneralVariableScan Architecture

Keep the existing:

`GeneralVariableScan.py`

unchanged or as close to unchanged as possible for the old standalone hardware architecture.

Create:

`GeneralVariableScan_master_satellite.py`

for the new master-satellite hardware architecture.

Do NOT extend the old `GeneralVariableScan.py` into a large hardware-architecture switch.

This separation is intentional so the old standalone system remains directly usable.

### 19.1 Master-satellite experiment mode

`GeneralVariableScan_master_satellite.py` should support:

```text
experiment_mode = single_node
experiment_mode = two_nodes
```

### 19.2 Single-node mode

When:

```text
experiment_mode = single_node
```

allow:

```text
which_node = Node1
which_node = Node2
```

The selected node should be exposed through the existing unsuffixed logical namespace.

This allows existing single-node physics experiment functions to be reused on the new master-satellite hardware where compatible.

### 19.3 Two-node mode

When:

```text
experiment_mode = two_nodes
```

both nodes are available simultaneously through `_Node1` / `_Node2` logical names.

Only new master-satellite experiment functions should depend on this true two-node namespace.

### 19.4 GUI organization

The master-satellite GVS may maintain separate GUI state for single-node and two-node experiment configuration if useful, including:

```text
experiment_function
scan_variable1_name
scan_sequence1
scan_variable2_name
scan_sequence2
override_experiment_variables
```

Do not make the old standalone GVS depend on these new fields.

---

## 20. Existing `experiment_functions.py`

All functions currently in:

`subroutines/experiment_functions.py`

were written for the previous standalone/non-master-satellite architecture.

This includes historical functions that coordinated two independent Kasli systems.

Important:

> Historical "two-node" experiment functions are NOT master-satellite two-node experiment functions.

Therefore:

- do NOT migrate historical two-node functions automatically,
- do NOT copy the large existing `experiment_functions.py`,
- do NOT delete historical two-node functions merely because the architecture changed,
- do NOT rewrite validated single-node physics sequences unnecessarily.

Existing single-node experiment functions should be reused under master-satellite single-node compatibility mode where practical.

A future/new file:

`subroutines/experiment_functions_master_satellite.py`

should contain only NEW functions intentionally written for the true master-satellite architecture.

Writing new atom-atom/two-node physics functions is not the first infrastructure priority.

---

## 21. AOM Feedback Architecture

Keep the existing standalone feedback path for the old standalone architecture.

The master-satellite path must support independent Node1 and Node2 stabilizers, conceptually:

```text
laser_stabilizer_Node1
laser_stabilizer_Node2
```

Node-specific feedback calibration/configuration should continue to reuse:

```text
utilities/config/alice/feedback_channels.json
utilities/config/bob/feedback_channels.json
```

Do not duplicate node calibration information unnecessarily.

### 21.1 Parallel feedback

Parallel AOM feedback is a desired/required master-satellite feature because Node1 and Node2 optical beams are independent.

Support:

```text
parallel_AOM_feedback = False
```

for sequential feedback and:

```text
parallel_AOM_feedback = True
```

for concurrent feedback.

### 21.2 ARTIQ timing caution

Do NOT assume two existing feedback kernels can simply be wrapped in:

```python
with parallel:
    ...
```

Before implementing broad parallel feedback, inspect:

- kernel boundaries,
- RTIO timing,
- `break_realtime()`,
- Sampler acquisition,
- DDS/SPI updates,
- delays,
- host/kernel interactions,
- timeline cursor behavior,
- DRTIO latency,
- shared resources,
- SED lane configuration.

The final implementation must respect ARTIQ scheduling behavior.

---

## 22. Satellite Readiness

The master-satellite initialization path must explicitly ensure the satellite destination is available before attempting to initialize or access remote hardware.

Do not assume the satellite is ready simply because the master core is reachable.

Initialization order should be designed deliberately.

Remote:

- SPI operations,
- input operations,
- Sampler operations,
- TTL input/timestamp operations,

must be treated as latency-sensitive where applicable.

---

## 23. DMA

Existing DMA sequences were validated under the standalone architecture.

Do not assume that a DMA trace involving satellite devices behaves identically after migration.

Revalidate any DMA sequence that includes Node2 satellite hardware.

Preserve existing standalone DMA implementations unless a master-satellite-specific change is required.

---

## 24. AD9910 Synchronization

AD9910 timing/synchronization must be recalibrated after the final master-satellite gateware and clocking configuration are fixed.

Do not treat old standalone timing/alignment calibration values as guaranteed final master-satellite values.

Do not optimize this prematurely before hardware/clock configuration is stable.

---

## 25. Edge Counters / SPCM Inputs

Preserve edge-counter support on both nodes.

Do not accidentally replace gateware edge-counter channels with plain TTL input channels during migration.

SPCM / edge-counter functionality must be included in hardware sanity tests.

### 25.1 Canonical Master-Satellite SPCM Mapping

For the NEW master-satellite architecture, the four canonical SPCM detector names use the exact four Node1/master-side TTL input and edge-counter channels already used by the existing Node1-side SPCM aliases.

These are not new physical channels. They are new fixed detector/optical-identity names for the existing Node1/master-local channels:

```text
SPCM_H1 -> Node1 ttl0  / ttl0_counter
SPCM_V1 -> Node1 ttl1  / ttl1_counter
SPCM_H2 -> Node1 ttl8  / ttl8_counter
SPCM_V2 -> Node1 ttl9  / ttl9_counter
```

Here:

```text
SPCM_H1 = existing Node1 SPCM0
SPCM_V1 = existing Node1 SPCM1
SPCM_H2 = existing Node1 SPCM0_OtherNode
SPCM_V2 = existing Node1 SPCM1_OtherNode
```

`H1`, `V1`, `H2`, and `V2` describe detector/optical identity, not Kasli ownership.

The old Node2/satellite-side SPCM TTL inputs and edge counters are not used for active photon detection in the master-satellite architecture. They may remain available as ordinary Node2 physical devices for hardware registration or testing, but they must not back the canonical SPCM aliases.

For master-satellite single-node compatibility mode, install the fixed legacy aliases:

```text
SPCM0           -> SPCM_H1
SPCM1           -> SPCM_V1
SPCM0_OtherNode -> SPCM_H2
SPCM1_OtherNode -> SPCM_V2
```

Apply the corresponding fixed mapping to both the TTL input attributes and their edge-counter attributes.

This alias mapping must be identical whether Node1 or Node2 is selected. Do not reorder or reinterpret the SPCM aliases for Node2. In particular, Node2 single-node compatibility mode must continue to use the four Node1/master-local physical detector channels listed above.

For new true two-node master-satellite experiments, expose and use:

```text
SPCM_H1
SPCM_V1
SPCM_H2
SPCM_V2
```

directly.

The single physical resolution invariant applies: execution mode and selected node may change the experiment-facing alias presentation, but they must not change which physical TTL input or edge counter represents a canonical detector.

---

## 26. USB / Host-Controlled Hardware

ALL Node1 and Node2 USB / host-controlled hardware will be connected to the MASTER-SIDE HOST.

This decision is final for this migration.

Do NOT design a distributed USB architecture.

This includes, but is not limited to:

- K10CR1 rotation stages,
- NDSP/controller devices,
- other host-controlled instruments.

Use globally unique host-device nicknames with the node suffix convention.

Examples:

```text
780_QWP_Node1
780_HWP_Node1
852_QWP_Node1
852_HWP_Node1

780_QWP_Node2
780_HWP_Node2
852_QWP_Node2
852_HWP_Node2
```

Reuse the existing multi-rotor K10CR1 NDSP infrastructure if practical.

USB/host devices do not need to be forced into the ARTIQ `master_satellite/device_aliases.json` physical mapping if they are controlled through a separate host/controller mechanism.

---

## 27. Hardware Sanity Tests

Before attempting a real master-satellite physics experiment, independently verify at minimum:

```text
canonical core
satellite destination readiness

Node1 DDS
Node2 DDS

Node1 Urukul CPLDs
Node2 Urukul CPLDs

Node1 TTL
Node2 TTL

Node1 edge counters / SPCM inputs
Node2 edge counters / SPCM inputs

Node1 Zotino
Node2 Zotino

Node1 Sampler
Node2 Sampler

Node1 USB devices
Node2 USB devices
```

Tests should be intentionally small.

Do not debug the infrastructure first using a large physics experiment.

Failures should be isolatable to:

- one node,
- one physical device,
- one mapping entry,
- one DRTIO destination,
- one controller,
- one timing assumption.

---

## 28. Migration Order

Use approximately this order unless actual repository dependencies require adjustment.

### Phase 0 — Preserve old standalone system

1. Keep existing standalone files intact.
2. Keep `GeneralVariableScan.py`.
3. Keep `BaseExperiment.py`.
4. Keep `DeviceAliases.py`.
5. Keep existing `ExperimentVariables.py`.
6. Keep existing `experiment_functions.py`.
7. Keep the obsolete legacy `which_node == "two_nodes"` mode removed.

### Phase 1 — Establish master-satellite physical mapping

8. Treat the reference master/satellite JSON and unified device DB as correct for software development.
9. Create/verify `utilities/config/master_satellite/device_aliases.json`.
10. Explicitly map Node1 standalone physical devices to unified physical names.
11. Explicitly map Node2 standalone physical devices to unified physical names.
12. Include required node-specific ARTIQ devices.
13. Do not encode arithmetic offset rules.

### Phase 2 — Master-satellite resolver

14. Create `utilities/DeviceAliases_master_satellite.py`.
15. Support legacy unsuffixed single-node presentation.
16. Support suffixed true two-node presentation.
17. Handle node-specific DDS defaults/bookkeeping.
18. Fail clearly on missing physical mappings.

### Phase 3 — Master-satellite BaseExperiment

19. Create `utilities/BaseExperiment_master_satellite.py`.
20. Register one canonical core.
21. Wait explicitly for satellite readiness.
22. Register/expose Node1 and Node2 hardware through the unified device database.
23. Support `single_node` and `two_nodes` execution modes.
24. Reset/init shared hardware in a deliberate order.
25. Preserve edge counters.

### Phase 4 — Experiment variables

26. Create `ExperimentVariables_Node1.py`.
27. Create `ExperimentVariables_Node2.py`.
28. Create `ExperimentVariables_master_satellite.py`.
29. Use globally unique Node1/Node2 persistent dataset names.
30. Add sequential/parallel feedback selection.

### Phase 5 — Master-satellite GeneralVariableScan

31. Create `GeneralVariableScan_master_satellite.py`.
32. Do not replace `GeneralVariableScan.py`.
33. Support master-satellite `single_node` mode with Node1/Node2 selection.
34. Support master-satellite `two_nodes` mode.
35. Keep old standalone GVS independent.

### Phase 6 — Feedback

36. Reuse Alice/Bob feedback calibration/configuration.
37. Create independent Node1/Node2 stabilizer contexts.
38. Verify sequential feedback first.
39. Inspect SED lane constraints.
40. Implement/verify parallel feedback.

### Phase 7 — Compatibility tests

41. Under old standalone gateware/device_db, verify old GVS still runs Node1.
42. Under old standalone gateware/device_db, verify old GVS still runs Node2.
43. Under master-satellite gateware/device_db, verify master-satellite GVS single-node Node1.
44. Under master-satellite gateware/device_db, verify master-satellite GVS single-node Node2.
45. Verify selected-node compatibility with existing single-node experiment functions.

### Phase 8 — Master-satellite timing validation

46. Revalidate satellite DMA sequences.
47. Revalidate remote input/SPI timing.
48. Recalibrate AD9910 synchronization if required.
49. Verify broad parallel scheduling only after SED lane inspection.

### Phase 9 — New master-satellite physics functions

50. Create `subroutines/experiment_functions_master_satellite.py`.
51. Write only NEW master-satellite-specific experiment functions.
52. Start with minimal simultaneous two-node tests.
53. Proceed toward atom-atom experiments.

---

## 29. Applets / Dashboard / Grafana

Applet restructuring is LAST PRIORITY.

Do not modify, consolidate, duplicate, or redesign applets during the current infrastructure migration unless an applet directly blocks testing.

Possible future work may include:

- reducing dashboard applet count,
- combining diagnostics,
- moving long-term telemetry to InfluxDB/Grafana.

InfluxDB/Grafana integration is NOT part of the current critical migration path.

---

## 30. Scope Control

Do NOT:

- redesign the entire repository,
- turn all existing standalone files into generic master-satellite files,
- rewrite `BaseExperiment.py` from scratch,
- rewrite `GeneralVariableScan.py` for master-satellite operation,
- migrate every historical experiment,
- migrate historical two-node experiments automatically,
- copy the entire `experiment_functions.py`,
- rewrite analysis code unnecessarily,
- remove old standalone compatibility,
- optimize unrelated physics sequences,
- build a highly generalized framework beyond current needs,
- modify applet/dashboard architecture now,
- implement Grafana/InfluxDB now,
- guess DRTIO destination/channel mapping,
- guess RTIO offsets,
- guess Urukul/TTL index offsets,
- assume local and remote timing are identical,
- assume successful static checks imply hardware correctness,
- duplicate master-satellite physical mapping tables.

Prefer:

- separate standalone and master-satellite execution stacks,
- one source of truth for master-satellite physical translation,
- explicit mappings,
- small changes,
- easy-to-review diffs,
- easy rollback,
- easy hardware debugging,
- preservation of validated physics sequences,
- small hardware tests before large experiments.

---

## 31. Codex Working Rules

Before modifying code:

1. Read this `AGENTS.md` completely.
2. Inspect the actual current repository implementation.
3. Consult the ARTIQ 7.82 manual for ARTIQ-specific behavior when relevant.
4. Use `master_satellite_test_files/` as the assumed hardware reference.
5. Do not assume this document perfectly describes every current code detail.
6. Identify the smallest relevant set of files.
7. Explain significant architectural changes before implementing them.
8. Identify dependencies and possible regressions.
9. Prefer the smallest safe implementation.
10. Do not modify unrelated files.
11. Do not commit or push unless explicitly requested.

For architectural changes, first report:

1. current behavior,
2. relevant code paths,
3. proposed behavior,
4. files that need modification,
5. compatibility impact,
6. hardware validation required.

Do not blindly normalize legacy naming.

Do not make old standalone files depend on the new master-satellite stack unless explicitly required.

Do not silently migrate historical experiment logic.

When uncertain about hardware semantics, explicitly flag the uncertainty instead of guessing.

---

## 32. Current Immediate Priority

The immediate priority is the master-satellite physical/logical mapping and resolver layer.

Do not start by rewriting physics experiment functions.

Current focus:

1. inspect the unified reference `device_db`,
2. finalize `utilities/config/master_satellite/device_aliases.json`,
3. ensure it is the single physical translation source for Node1 and Node2,
4. design `DeviceAliases_master_satellite.py`,
5. ensure the resolver supports both:
   - unsuffixed master-satellite single-node compatibility,
   - `_Node1` / `_Node2` true two-node presentation,
6. then build `BaseExperiment_master_satellite.py`.

Do not modify the reference gateware/device-db files.

Do not modify applets.

Do not migrate historical two-node experiment functions.

Do not commit or push unless explicitly requested.

---

## 33. Implemented Manual Master-Satellite Hardware Control

`AOMsCoils_master_satellite.py` is the dedicated post-gateware manual-control
utility. Standalone `AOMsCoils.py` remains untouched. The new utility is
separate from GVS and always uses Base in `two_nodes` mode because an
unselected node must remain accessible to be actively forced safe.

Its selector is deliberately:

```text
which_node = node1 | node2 | two_nodes
```

This is a hard desired-state gate:

- `node1`: honor Node1 controls; force controlled Node2 DDS outputs OFF and
  every Node2 coil to 0 V;
- `node2`: honor Node2 controls; force controlled Node1 DDS outputs OFF and
  every Node1 coil to 0 V;
- `two_nodes`: honor both independently.

It controls FORT, cooling DP, AOM A1-A6, D1 pumping DP, GRIN1and2, relevant
active-low optical TTL gates, microwave DDS, MW RF, and four independent coils
per node. DDS frequency/power comes through existing resolver bindings from
authoritative node-suffixed ExperimentVariables. Node2 D1 specifically uses
`f_GRIN2_D1_pumping_Node2` and `p_GRIN2_D1_pumping_Node2`.

Each coil has an independent checkbox and run-local voltage seeded from MOT
calibration. These values are never persisted, and OFF always writes 0 V.

```text
Node1: AZ bottom 0, AZ top 1, AX 13, AY 14
Node2: AZ bottom 0, AZ top 1, AX 2,  AY 3
```

Microwave/RF defaults OFF and requires explicit confirmation. Laser feedback,
old independent-Kasli networking, K10CR1, Rigol, and other USB functions are
excluded from the first version.

Safe initialization switches each DDS OFF immediately after initialization
and before attenuation/profile programming, zeroes both Zotinos, explicitly
forces controlled hardware OFF/0 V again, and only then applies requested
state. Base initialization is still not passive because it initializes all
CPLDs/DDS channels, TTL directions/states, Samplers/gains, and Zotinos.

Commit `078dc0e` implements the utility, tests, and DDS safety ordering. The
combined suite passed 51 hardware-free tests. Physical validation is pending.

---

## 34. Updated Immediate Priority

1. Secure verified standalone rollback artifacts.
2. Build/archive matching master and satellite images plus unified database.
3. Perform passive boot, log, and destination-1 checks.
4. Verify namespace and harmless local/remote RTIO access.
5. Start the manual utility with every output OFF and every coil at 0 V.
6. Validate one Node1 AOM and coil at a time, then Node2.
7. Confirm TTL polarity, coil identity/polarity, DRTIO readiness, and remote
   SPI slack.
8. Validate magnetometer Node1 then Node2.
9. Add feedback and atom-loading compatibility only afterward.

AD9910 synchronization, DMA, timing-critical SPCM/herald/mapping behavior,
parallel feedback, K10CR1, and microwave optimizer migration remain deferred.
