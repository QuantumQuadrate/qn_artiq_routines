# Master-Satellite Migration: Detailed Authoritative Specification

## 1. Purpose and authority

This is the implementation specification for migrating two independent ARTIQ
Kasli-SoCs to a true DRTIO master-satellite system. It reorganizes the
decisions accumulated during development into one handoff.

When instructions conflict, use this order:

1. the user's latest explicit instruction;
2. this detailed specification;
3. `PLAN.md` for historical rationale;
4. repository behavior and tests;
5. assumptions, which must be reported before use.

## 2. Physical context

```text
Node1 = Alice = current standalone .75 = future DRTIO master .75
Node2 = Bob   = current standalone .76 = future satellite destination 1
```

Both physical systems still run standalone gateware and have the same
peripheral arrangement. The target has one master Ethernet/core endpoint, one
canonical ARTIQ core, and one unified master-side `device_db.py`.

The purpose is coordinated two-node timing/control while retaining known-good
single-node physics and a dependable rollback path.

## 3. Non-negotiable rules

1. Preserve unrelated worktree changes.
2. Inspect and report evidence before editing.
3. Do not flash hardware without phase-specific authorization.
4. Do not commit or push without explicit instruction.
5. Do not modify applets unless requested.
6. Avoid unrelated cleanup, renaming, and formatting.
7. Leave historical/deprecated experiments untouched merely because they
   mention two nodes or raw TTLs.
8. Never calculate DRTIO channels, destination bits, Urukul offsets, or TTL
   offsets.
9. Treat generated gateware and its database as one matched release.
10. Identify everything requiring experiment-computer validation.

## 4. Separate software architectures

### 4.1 Validated standalone stack

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
`which_node == "two_nodes"` mode was independent-Kasli coordination, not
DRTIO, and must not return.

### 4.2 New master-satellite stack

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

Do not create dependencies from the old stack onto the new one.

## 5. Execution identity and naming

The new GVS owns:

```text
experiment_mode = single_node
selected_node = Node1 or Node2

experiment_mode = two_nodes
```

These are execution arguments, not calibration datasets. `selected_node` is
authoritative. Reused single-node functions alone receive:

```text
Node1 -> self.which_node = "alice"
Node2 -> self.which_node = "bob"
```

New logical names and node-specific datasets use suffixes:

```text
dds_FORT_Node1
dds_FORT_Node2
sampler0_Node1
zotino0_Node2
f_FORT_Node1
p_FORT_loading_Node2
```

## 6. Hardware reference

Use these as software mapping references until matching images and physical
hardware are verified:

```text
master_satellite_test_files/node1_master_kasli_soc_20260715.json
master_satellite_test_files/node2_satellite_kasli_soc_20260715.json
master_satellite_test_files/device_db_two_nodes_20260715.py
```

Do not modify them.

Each JSON describes two DIO cards with edge counters, three Samplers, one
Zotino, and three synchronized AD9910 Urukuls at 125 MHz RTIO frequency with
Urukul `clk_sel = 2`.

Node1 declares `drtio_master`, base `master`, carrier revision `v1.0`,
and address .75. Node2 declares `drtio_satellite`, base `satellite`, and
carrier revision `v1.1`. Verify revisions physically.

Destination 1 is not a satellite JSON field. It comes from combined database
generation and routing configuration.

## 7. Unified low-level namespace

| Subsystem | Node1/master | Node2/satellite |
|---|---|---|
| DIO TTLs | `ttl0`-`ttl15` | `ttl16`-`ttl31` |
| Counters bank 1 | `ttl0_counter`-`ttl3_counter` | `ttl16_counter`-`ttl19_counter` |
| Counters bank 2 | `ttl8_counter`-`ttl11_counter` | `ttl24_counter`-`ttl27_counter` |
| Samplers | `sampler0`-`sampler2` | `sampler3`-`sampler5` |
| Zotino | `zotino0` | `zotino1` |
| CPLDs | `urukul0_cpld`-`urukul2_cpld` | `urukul3_cpld`-`urukul5_cpld` |
| DDS | `urukul0_ch0`-`urukul2_ch3` | `urukul3_ch0`-`urukul5_ch3` |

Node2 direct channels encode destination 1. Python name sequencing does not
prove encoded channel sequencing.

Generated SPI/CNV/LDAC/CLR/switch/sync/IO-update/EEPROM devices stay database
implementation details unless active code directly references them. Search
before excluding them from a required translation.

The unified reference omits standalone `led0`, `led1`, and `k10cr1_ndsp`.
K10CR1 is host-controlled and must not be forced into RTIO mapping.

## 8. Physical translation

`utilities/config/master_satellite/device_aliases.json` is the single source
translating old standalone physical namespaces to unified physical names,
separated by Node1 and Node2.

It includes required node-specific RTIO devices: DDS channels, CPLDs, Samplers,
Zotino, TTLs, and edge counters. It excludes the canonical core, logical DDS
aliases, USB devices without a reason, and all arithmetic rules.

Alice/Bob JSON files retain:

```text
logical DDS alias -> standalone DDS + default variable names
```

The shared master-satellite file adds:

```text
standalone physical name -> unified physical name
```

## 9. Single physical resolution invariant

Resolve physical identity first, then presentation. A node/logical-device pair
must resolve identically in both modes.

```text
Bob dds_FORT -> standalone urukul1_ch3 -> unified urukul4_ch3

single Node2: self.dds_FORT       -> urukul4_ch3
two nodes:    self.dds_FORT_Node2 -> urukul4_ch3
```

Never maintain a separate single-node physical map.

## 10. Canonical SPCM exception

Active detectors are always master-local Node1 devices:

| Detector | TTL | Counter |
|---|---|---|
| `SPCM_H1` | `ttl0` | `ttl0_counter` |
| `SPCM_V1` | `ttl1` | `ttl1_counter` |
| `SPCM_H2` | `ttl8` | `ttl8_counter` |
| `SPCM_V2` | `ttl9` | `ttl9_counter` |

H1/V1/H2/V2 are detector identities, not Kasli ownership. For either selected
node, single-node compatibility is fixed:

```text
ttl_SPCM0                   -> SPCM_H1
ttl_SPCM1                   -> SPCM_V1
ttl_SPCM0_OtherNode         -> SPCM_H2
ttl_SPCM1_OtherNode         -> SPCM_V2
```

Counter aliases follow the same order. Never route them through Node2 mappings
such as `ttl0 -> ttl16`. Node2's old SPCM-looking inputs physically remain
but are not active detectors.

Do not rewrite deprecated/test functions that use raw TTLs. Document the
distinction in the new Base binding layer.

## 11. DeviceAliasesMasterSatellite

The resolver owns:

- loading/validating the exact physical mapping;
- exact dictionary lookup;
- generic binding under a caller-supplied presentation name;
- logical DDS -> standalone DDS -> unified DDS resolution;
- caller-selected DDS presentation;
- node-specific DDS default lookup;
- independent per-node DDS device/frequency/power bookkeeping;
- clear errors for missing aliases, mappings, devices, or defaults.

It does not reset core, wait for DRTIO, initialize hardware, choose modes,
implement SPCM policy, or handle USB devices.

## 12. BaseExperimentMasterSatellite

### Build

Build validates mode/node, registers one core/cache/DMA, loads variables,
creates single-node projections, creates resolver contexts, binds node devices
and SPCM aliases, and installs fixed wiring metadata.

Active contexts:

```text
single Node1: Node1 ordinary devices + Node1 shared SPCM map
single Node2: Node2 ordinary devices + Node1 shared SPCM map
two nodes:    Node1 ordinary + Node2 ordinary + master SPCM map
```

### Prepare

Prepare is one-shot and owns DDS binding/preparation plus initial caches. Never
call it per scan point.

### Repeatable APIs

`resolve_experiment_variable_target(name)` resolves globals to themselves,
maps single-mode legacy names to the selected suffix, accepts selected
authoritative names, rejects other-node names, requires explicit suffixes in
two-node mode, and rejects unknown/ambiguous names clearly.

`refresh_compatibility_variables()` copies authoritative selected-node
attributes to legacy in-memory names without dataset I/O.

`refresh_variable_dependent_state()` refreshes actual cached derived state,
initially DDS frequency/power caches, without prepare, rebinding, hardware
initialization, core reset, or dataset I/O.

`reload_experiment_variables()` rereads relevant authoritative and global
datasets, refreshes projections/caches, and fails if required data disappeared.
It never creates/writes datasets or rebuilds hardware.

Result initialization/reset is separate from variable-dependent refresh.

### Hardware initialization order

1. one Base-owned core reset;
2. destination-1 readiness wait if Node2 is active;
3. canonical SPCM input configuration;
4. active Node1 then Node2 CPLD initialization;
5. active Node1 then Node2 DDS initialization;
6. TTL directions and safe states;
7. Sampler initialization/gains;
8. Zotino initialization/safe state;
9. final slack restoration.

The argument is `satellite_ready_poll_interval`: startup polling, not DRTIO
communication latency.

## 13. Static wiring metadata

Single Node1:

```text
coil_names = ["AZ bottom", "AZ top", "AX", "AY"]
AZ_bottom_Zotino_channel = 0
AZ_top_Zotino_channel = 1
AX_Zotino_channel = 13
AY_Zotino_channel = 14
coil_channels = [0, 1, 13, 14]
UV_trig_channel = [8]
Osc_trig_channel = [10]
FORT_MM_sampler_ch = 7
GRIN1_sampler_ch = 4
Magnetometer_X_ch = 1
Magnetometer_Y_ch = 2
Magnetometer_Z_ch = 3
```

Single Node2 differs only in:

```text
AX_Zotino_channel = 2
AY_Zotino_channel = 3
coil_channels = [0, 1, 2, 3]
```

These are fixed wiring facts, not ExperimentVariables or scan-dependent state.

## 14. ExperimentVariables

Ownership:

```text
ExperimentVariables_Node1.py             owns *_Node1
ExperimentVariables_Node2.py             owns *_Node2
ExperimentVariables_master_satellite.py  owns globals only
```

The global file must not initialize node files. Fresh database workflow is
Node1 once, Node2 once, then global once. Normal experiments only read them and
fail clearly when required datasets are absent.

Existing persistent calibration always wins. Defaults seed missing datasets
only. Preserve asymmetric/historical node variables, and create both suffixes
for conceptually node-calibratable quantities.

Confirmed values and rules:

```text
t_coil_relaxation_time_OP_Node1 = 1 ms seed
t_coil_relaxation_time_OP_Node2 = 0.4 ms seed
n_measurements = 100 global
t_delay_in_bob_mu = 189 global, exact name
parallel_AOM_feedback = True global
```

Keep historical D1 variables but do not invent `set_point_D1_DP`. Active
feedback is AOM1-AOM6 plus FORT. Preserve independent loading/science/holding
FORT setpoints for both nodes.

Bob `dds_D1_pumping_DP` uses `f_GRIN2_D1_pumping_Node2` and
`p_GRIN2_D1_pumping_Node2`.

Do not migrate runtime/results such as `fit_parameter_*`,
`feedbackchannels`, `parent_rid`, or `child_rid`.

## 15. GeneralVariableScanMasterSatellite

Single-node discovery accepts locally defined functions satisfying
`inspect.isfunction`, matching module ownership, and containing
`"experiment"` in the name.

Exclude exactly:

```text
Two_nodes_atom_loading_experiment
Two_node_single_photon_experiment
Two_node_single_photon_2_experiment
Two_node_single_photon_2_optimization_experiment
```

Do not build a large allowlist. True two-node discovery uses only the new
master-satellite function module. Function selection uses registry lookup, not
`eval()`. Expression evaluation for scan sequences/overrides may remain.

Run-wide override flow:

1. resolve through Base;
2. mutate authoritative attribute;
3. refresh compatibility;
4. refresh dependent state.

Scan-point flow:

1. mutate authoritative target(s);
2. refresh compatibility;
3. refresh dependent state;
4. initialize hardware;
5. reset transient result state;
6. invoke selected callable.

Scans/overrides are non-persistent. Use `reload_experiment_variables()` for
queued freshness. Never rebuild or prepare per point.

`n_measurements` loads from the global dataset and may be a current-run GUI
override without changing persistent calibration.

The CatchError subclass retries only the failed point on `RTIOUnderflow`.
It uses bounded retries and host backoff, and must propagate DRTIO, SPI,
resolution, and other non-underflow failures.

## 16. Experiment functions

Do not copy the large legacy module. Reuse compatible functions only in
single-node mode. New DRTIO/two-node functions live only in the new module.
Historical independent-Kasli functions remain untouched.

`master_satellite_namespace_sanity_experiment` is host-side and
non-destructive by itself. Through current GVS it is not passive, because Base
hardware initialization happens first.

## 17. Magnetometer compatibility

The first reused target is `test_magnetometer_experiment`, only after passive
DRTIO validation.

Physical results use node-suffixed transient names following standalone
spelling, for example:

```text
Magnetometer_MOT_X_Node1
Magnetometer_MOT_X_Node2
```

Legacy hard-coded unsuffixed magnetometer append names are redirected only in
the new stack. Do not create duplicate unsuffixed result datasets.
`measurements_progress` stays unsuffixed as run-global state.

Result initialization/reset must not persist calibration, reload variables,
reset core, or rebind devices. Preserve existing `break_realtime()` calls and
final MOT coil behavior. Add no speculative post-reset readiness wait or slack
until hardware demonstrates a failure.

## 18. Feedback scope

Active feedback is AOM1-AOM6 and FORT with loading/science/holding setpoints.
D1 feedback is inactive; retain history without expanding it.

Full standalone feedback compatibility is not implemented in the new Base:
`laser_stabilizer`, feedback histories, amplitude arrays, and related derived
state remain future work. Do not implement parallel feedback before
single-node feedback works and SED-lane constraints are checked.

## 19. Microwave optimizer

The standalone `MicrowaveScanOptimizer.py` is incompatible and remains
untouched. A separate master-satellite optimizer is required because it:

- constructs standalone Base;
- treats legacy `which_node` as authoritative;
- scans/persists unsuffixed node calibration;
- repeats `prepare()`;
- uses old queued rebuild behavior;
- expects full standalone results and feedback objects;
- needs derived cooling/FORT amplitudes;
- uses unvalidated DMA traces;
- has internal resets without destination checks;
- resubmits jobs without master-satellite mode/node state.

Port it only after atom loading, feedback, derived amplitudes, result state, and
DMA work in master-satellite single-node mode.

## 20. DRTIO/timing constraints

- Remote devices are accessed normally through unified database names.
- Do not invent destination arguments for ordinary access.
- Wait for destination 1 before Node2 access.
- Remote input/SPI operations are latency-sensitive.
- RTIO underflow is a scheduling-slack issue, not readiness polling.
- Revalidate every DMA trace involving satellite hardware.
- Inspect SED lanes before broad parallel RTIO.
- Recalibrate AD9910 synchronization after final gateware/clocking.
- Preserve both nodes' edge counters.
- Defer herald/mapping, photon timestamp, and physics-load timing tests.

## 21. Host-controlled devices

All Node1/Node2 USB and host-controlled devices will connect to the master-side
host. K10CR1/NDSP, Rigol, and similar devices do not belong in per-node RTIO
translation without a concrete reason. Eventually use suffixed global
nicknames and selected-node legacy projection where needed. This work is
deferred and must not block basic DRTIO/magnetometer validation.

## 22. Gateware and rollback gate

Before any flash:

1. capture and verify full SD-card images or exact standalone `boot.bin` for
   both nodes;
2. preserve actual active standalone databases/configuration, startup/idle
   kernels, controllers, datasets, ARTIQ version, and logs;
3. record carrier/peripheral revisions and physical port wiring;
4. identify exact ARTIQ and ARTIQ-Zynq build revisions;
5. build/obtain matched master-runtime and satellite-satman images;
6. generate unified database from those exact JSONs with destination 1;
7. determine master downstream SFP port and routing table;
8. verify transceivers/cable for the 125 MHz/2.5 Gbit/s link;
9. archive images, JSONs, DB, route, logs, revisions, and hashes together;
10. verify rollback media.

The repository has no deployable 20260715 image/build directory, route table,
or proof of the master SFP cage. ARTIQ-Zynq's `boot` key is write-only, so
rollback cannot depend on reading it later.

Safe staging concept after approval: stop services; stage satellite image on
.76 and master image on .75 while both are reachable; do not reboot Node2
prematurely; power down; connect verified DRTIO; boot with UART monitoring,
operationally satellite then master; verify .75 and destination 1 before the
dashboard.

## 23. First hardware-validation sequence

1. master boot/network sanity;
2. both UART/core logs;
3. bounded passive `core.get_rtio_destination_status(1)` kernel;
4. host-only namespace resolution without hardware initialization;
5. one harmless Node1 local operation;
6. one harmless/read-only Node2 remote operation;
7. magnetometer single Node1, `n_measurements = 1`;
8. magnetometer single Node2, `n_measurements = 1`.

On failure report exact mode/node, exception, device/operation, reset stage,
destination readiness, resolved alias, and whether slack is relevant before
proposing a targeted fix.

## 24. Initialization hazards

Current Base initialization can reset core, change SPCM/TTL directions and
levels, initialize CPLDs/DDSs, set DDS attenuation/frequency/amplitude/switch
state, change Sampler gains, and write Zotino channels to zero. Thus current
GVS is not a passive first test. Confirm every output state is electrically
safe first.

## 25. Completed status

### Manual master-satellite hardware control

`AOMsCoils_master_satellite.py` is separate from GVS and from the untouched
standalone `AOMsCoils.py`. Build registers the suffixed two-node device
superset without configuring execution or loading datasets. Prepare strictly
configures Base in `two_nodes` mode so both physical nodes remain accessible.
Its GUI selector is deliberately `which_node = node1 | node2 | two_nodes` and
is a hard desired-state gate:

| Selection | Node1 | Node2 |
|---|---|---|
| `node1` | honor controls | force DDS OFF and coils 0 V |
| `node2` | force DDS OFF and coils 0 V | honor controls |
| `two_nodes` | honor independently | honor independently |

Unselected-node checkboxes can never preserve or enable output.

Controlled DDS aliases and exact unified devices are:

| Logical control | Node1 | Node2 |
|---|---|---|
| FORT | `urukul0_ch0` | `urukul4_ch3` |
| cooling DP | `urukul0_ch1` | `urukul4_ch2` |
| AOM A1 | `urukul1_ch2` | `urukul3_ch0` |
| AOM A2 | `urukul1_ch0` | `urukul3_ch1` |
| AOM A3 | `urukul1_ch1` | `urukul3_ch2` |
| AOM A4 | `urukul2_ch0` | `urukul3_ch3` |
| AOM A5 | `urukul2_ch1` | `urukul4_ch0` |
| AOM A6 | `urukul1_ch3` | `urukul4_ch1` |
| D1 pumping DP | `urukul2_ch2` | `urukul5_ch0` |
| GRIN1and2 | `urukul0_ch3` | `urukul5_ch2` |
| microwaves | `urukul2_ch3` | `urukul5_ch3` |
| MW RF | `urukul0_ch2` | `urukul5_ch1` |

Frequency and power come through resolver bindings to authoritative
`*_Node1`/`*_Node2` attributes. Node2 D1 must use
`f_GRIN2_D1_pumping_Node2` and `p_GRIN2_D1_pumping_Node2`. Existing
active-low polarity is retained for repump, pumping-repump, excitation0,
GRIN1, GRIN2, microwave, and Node2's dedicated D1 gate. Node1 D1 shares the
GRIN1 optical gate.

Each AZ-bottom, AZ-top, AX, and AY coil has an independent checkbox and
run-local voltage argument. Explorer metadata uses a safe explicit 0 V default
because persistent MOT datasets may not exist during repository examination.
Values are limited to +/-10 V and never persisted. Disabled coils are
explicitly 0 V.

| Coil | Node1 Zotino | Node2 Zotino |
|---|---:|---:|
| AZ bottom | 0 | 0 |
| AZ top | 1 | 1 |
| AX | 13 | 2 |
| AY | 14 | 3 |

Microwave/RF defaults OFF, obeys hard gating, and requires confirmation. The
first version excludes laser feedback, old independent-Kasli networking,
K10CR1, Rigol, and other USB devices.

Safe startup now switches each requested-off DDS OFF immediately after
`init()` and before attenuation/profile programming. Base still initializes
CPLDs, all DDS channels, TTL directions/states, Samplers/gains, and all Zotino
outputs, so it is not passive. The utility performs a second explicit OFF/0-V
pass before applying requested state.

Commit `078dc0e` implements the utility, focused tests, and DDS ordering.
Kernel compilation, active-low polarity, coil identity/polarity, DRTIO
readiness, remote SPI slack, and AD9910 synchronization still require hardware
testing.

### Repository-examination lifecycle

ARTIQ repository examination supplies `None` for argument values and must work
with a completely empty master-satellite dataset namespace. This is required
so Explorer can discover the ExperimentVariables initializers themselves.

Public GVS and AOMsCoils build behavior is therefore limited to argument
metadata and binding the suffixed Node1/Node2 physical-device superset. It does
not select a runtime mode, load node/global calibration datasets, publish
single-node compatibility, or prepare DDS defaults.

At prepare:

1. submitted `experiment_mode`/`selected_node` or AOMsCoils `which_node` is
   strictly validated;
2. Base configures `single_node Node1`, `single_node Node2`, or `two_nodes`;
3. only the required node variables plus globals are loaded for GVS;
4. AOMsCoils loads both node namespaces plus globals because it must force the
   unselected node safe;
5. selected-node compatibility and DDS bindings are published;
6. missing datasets fail clearly before hardware initialization.

For deferred single-node configuration, both nodes remain registered but the
unselected node's CPLD, Sampler, Zotino, TTL, and safety-state lists are cleared
from the hardware lifecycle. Thus superset registration does not broaden
single-node hardware access.

GVS metadata uses fixed code defaults where persistent values cannot safely be
read during examination: `n_measurements = 100` and the normal single-node scan
presentation. The function selector exposes the union of the two registries
because Explorer metadata cannot dynamically change after mode selection;
prepare enforces mode membership before execution.

### Explorer naming and public-class imports

ARTIQ uses the first line of each public `EnvExperiment` class docstring as the
Explorer label. The expected labels are:

```text
ExperimentVariables_Node1
ExperimentVariables_Node2
ExperimentVariables_master_satellite
GeneralVariableScan_master_satellite
GeneralVariableScan_CatchError_master_satellite
AOMsCoils_master_satellite
```

Do not publicly import one dashboard `EnvExperiment` class into another
repository file. CatchError previously imported `GeneralVariableScanMasterSatellite`
under its public name, so ARTIQ rediscovered the base class and generated
`GeneralVariableScan_master_satellite1`. Importing it as
`_GeneralVariableScanMasterSatellite` preserves inheritance while hiding only
the accidental duplicate. The CatchError subclass remains intentional and
public. The namespace sanity function remains a GVS callable, not an Explorer
experiment.

Commits `a3eb7fe` and `f06f9be` implement these lifecycle/discovery fixes. The
combined hardware-free suite passed 56 tests, including every public build
with `None` arguments and empty dataset storage, plus separate strict runtime
missing-dataset failures.

Implemented and hardware-free tested on branch
`20260826_master_satellite_codex` through commit `f06f9be`:

- legacy selector removal;
- explicit physical mapping;
- resolver and Base lifecycle;
- destination readiness and fixed SPCM policy;
- Node1/Node2/global variables;
- authoritative loading/projection/reload/target/cache APIs;
- mode-specific registries and native namespace sanity;
- GVS and CatchError variant;
- selected-node magnetometer wiring and suffixed results;
- deterministic manual AOM/DDS, TTL optical-gate, microwave/RF, and
  independent coil control;
- DDS switch-OFF-before-programming safety ordering;
- examination-safe deferred Base configuration with empty dataset storage;
- stable Explorer class names and no duplicate public GVS import;
- hardware-free tests.

This is readiness for controlled gateware integration, not proof that physical
flashing or physics is safe.

## 26. Remaining order

1. secure rollback artifacts;
2. prepare matched gateware/database release bundle;
3. pass passive DRTIO validation;
4. validate manual outputs with all OFF first, then one AOM/coil at a time;
5. validate magnetometer Node1 then Node2;
6. add derived cooling/FORT amplitudes and result state;
7. add active AOM/FORT feedback compatibility;
8. validate normal atom loading;
9. create master-satellite microwave optimizer;
10. recalibrate AD9910;
11. revalidate DMA/remote SPI timing;
12. add native two-node physics;
13. validate timing-critical photon/herald/mapping work;
14. update applets after result names stabilize.

## 27. Reporting contract

Before editing report files, evidence, ambiguity, risk, and bounded scope.
After editing list files/APIs, run syntax/tests, run `git diff --check`, show
`git diff --stat`, identify hardware validation, and confirm standalone files
were untouched. Never commit or push automatically.

For hardware work record image hashes, route/configuration, cabling, boot logs,
destination readiness, tested operation, and pass/fail evidence. Local Python
tests never substitute for experiment-system validation.
