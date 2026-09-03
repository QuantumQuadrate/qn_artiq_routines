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
standalone/GeneralVariableScan.py        (all 16 standalone Explorer
standalone/ExperimentVariables.py         entry points live here)
standalone/AOMsCoils.py
standalone/MicrowaveScanOptimizer.py
utilities/BaseExperiment.py
utilities/DeviceAliases.py
subroutines/experiment_functions.py
subroutines/aom_feedback.py
utilities/config/alice/
utilities/config/bob/
```

Standalone selections remain only `alice` and `bob`. The removed legacy
`which_node == "two_nodes"` mode was independent-Kasli coordination, not
DRTIO, and must not return. `GeneralVariableScan.py` absorbed
`GeneralVariableScan_CatchUnderflow.py`: per-iteration underflow retry is its
optional `enable_Catch_UnderFlow` argument (off by default; the retry values
in the "Catch Underflow" GUI group take effect only when it is enabled).

All standalone Explorer experiment entry points live in the `standalone/`
folder; ARTIQ scans the repository recursively, so they appear in Explorer
under a `standalone` group while the master-satellite stack keeps the
repository top level (the primary namespace after the migration). Shared
infrastructure (`utilities/`, `subroutines/`, `fitting/`, configs) and the
already-foldered standalone consumers (`tests/`, `MOT_experiments/`) stay in
place. Import resolution is cwd-based
(`sys.path.append(cwd + "\\repository\\qn_artiq_routines")`), so the move
changed no `utilities.*`/`subroutines.*` imports; the files that import
`ExperimentVariables` or `MicrowaveScanOptimizer` as top-level modules
(`utilities/BaseExperiment.py`, `MOT_experiments/Coils-SPCMCounts.py`, the
two HealthCheck files) additionally append the `standalone` folder to
`sys.path`. Hardcoded scheduler expid `file` paths in the standalone
optimizer and HealthCheck files carry the `qn_artiq_routines\\standalone\\`
prefix, and the master-satellite optimizer reads the standalone scan
definitions from `standalone/MicrowaveScanOptimizer.py`.

### 4.2 New master-satellite stack

```text
GeneralVariableScan_master_satellite_mixin.py
GeneralVariableScan_master_satellite_single_node.py
GeneralVariableScan_master_satellite_two_nodes.py
AOMsCoils_master_satellite_mixin.py
AOMsCoils_master_satellite_Node1.py
AOMsCoils_master_satellite_Node2.py
FORT_Polarization_Optimizer_master_satellite.py
MicrowaveScanOptimizer_master_satellite.py
ExperimentVariables_master_satellite_Node1.py
ExperimentVariables_master_satellite_Node2.py
ExperimentVariables_master_satellite_global.py
utilities/BaseExperiment_master_satellite.py
utilities/DeviceAliases_master_satellite.py
subroutines/experiment_functions_two_nodes.py
utilities/config/master_satellite/device_aliases.json
```

The `_mixin` files hold only private shared implementations, define no public
EnvExperiment class, and never appear in ARTIQ Explorer. Every public
experiment is one class in one file whose name matches its Explorer label.

Do not create dependencies from the old stack onto the new one.

## 5. Execution identity and naming

Execution mode is fixed per public experiment file:

```text
single_node  GeneralVariableScan_master_satellite_single_node
             selected_node = Node1 or Node2 (submitted argument)

two_nodes    GeneralVariableScan_master_satellite_two_nodes
```

The mode is a class constant; `selected_node` is an execution argument, not a
calibration dataset, and is authoritative. Reused single-node functions alone
receive:

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

Build validates mode/node when configured at construction, registers one
core/DMA/scheduler, loads variables, creates single-node projections, creates
resolver contexts, binds node devices and SPCM aliases, and installs fixed
wiring metadata. When constructed without a mode for repository examination,
build binds the suffixed superset and defers validation, loading, and
projection to `configure_execution()`/`prepare()` (see section 25).

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
ExperimentVariables_master_satellite_Node1.py    owns *_Node1
ExperimentVariables_master_satellite_Node2.py    owns *_Node2
ExperimentVariables_master_satellite_global.py   owns globals only
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

The shared implementation lives in
`GeneralVariableScan_master_satellite_mixin.py`
(`_GeneralVariableScanMasterSatelliteMixin`). The two public experiments fix
their mode as a class constant and expose only that mode's function registry:
the single-node file discovers reusable legacy functions from
`subroutines/experiment_functions.py`; the two-node file uses only
`subroutines/experiment_functions_two_nodes.py`.

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

Submitted GUI values always win for the run: `n_measurements` shares its name
with the persistent global dataset, so prepare captures the submitted value
before `configure_execution()` loads datasets and re-asserts it afterward,
and the queued-run reload restores it again. Persistent calibration is never
changed.

Underflow retry is built into the mixin behind the optional
`enable_Catch_UnderFlow` argument (off by default), mirroring the standalone
`GeneralVariableScan.py`. When enabled, one scan point is wrapped in bounded
retries with host backoff using the "Catch Underflow" GUI values. Only
`RTIOUnderflow` is caught, and a retried point re-runs the full scan-point
flow, including hardware initialization and its core reset. DRTIO, SPI,
resolution, and all other non-underflow failures always propagate; with the
argument off, the point executes with no exception handling at all.

## 16. Experiment functions

Do not copy the large legacy module. Reuse compatible functions only in
single-node mode. New DRTIO/two-node functions live only in
`subroutines/experiment_functions_two_nodes.py`. Historical independent-Kasli
functions remain untouched.

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

Single-node feedback compatibility is implemented without modifying
`subroutines/aom_feedback.py`: every dataset access in AOMPowerStabilizer
goes through the experiment object, so
`BaseExperimentMasterSatellite.prepare_laser_stabilizer()` installs a name
map (from the selected node's `feedback_channels.json`) and the
`_DatasetRedirectMixin` routes those reads and writes to the node-suffixed
authoritative datasets (`p_AOM_A1_NodeX`, `*_history`, `feedbackchannels`).
The stabilizer itself is constructed from the projected legacy namespace
(`fast_feedback_dds_list`, samplers, set points, `dds_defaults`,
`which_node`), so shared feedback code stays byte-identical to main.

GVS-level feedback integration and hardware validation remain outstanding.
Do not implement parallel feedback before single-node feedback works on
hardware and SED-lane constraints are checked; Node2 feedback exercises
remote Sampler SPI over DRTIO and is validated only after Node1.

## 19. Microwave optimizer

The standalone `MicrowaveScanOptimizer.py` remains untouched and keeps
serving the standalone systems together with the standalone HealthCheck GVS
files, whose expids point at the standalone class.

`MicrowaveScanOptimizer_master_satellite.py` is the port
(`which_node = Node1 | Node2`; two_nodes is deliberately unsupported and
never planned):

- `scan_options`/`scan_dict` are extracted from the standalone source as
  pure literals (no import of the standalone module), so the scan
  definitions cannot drift;
- Base provides the derived RF amplitudes, the traced single-node
  atom-result datasets/buffers, `print_async`, and `write_results` that the
  reused experiment functions need;
- microwave calibrations and the `health_check_uw_*` fidelities persist
  node-suffixed through the dataset redirect
  (`MICROWAVE_CALIBRATION_DATASETS`); `fit_parameter_*` and `parent_rid`
  stay unsuffixed run state;
- scans and overrides mutate the suffixed authoritative attributes and
  refresh the projection instead of re-running the one-shot Base prepare;
  queued freshness uses `reload_experiment_variables()`;
- health-check failures resubmit THIS class with `which_node` carried
  forward. Future master-satellite HealthCheck files must submit the same
  expid shape; they do not exist yet.

Hardware validation is pending: Node1 first; Node2 scans exercise DMA
playback and remote Sampler SPI over DRTIO (section 20). The port pulls
sections 26.6/26.9 ahead of the atom-loading validation order by explicit
decision.

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
host. Host devices do not belong in per-node RTIO translation.

K10CR1 is implemented: one `k10cr1_ndsp` controller entry in the unified
database owns all eight rotators through one USB hub, with node-suffixed
axis nicknames (`780_HWP_Node1`, ..., `852_QWP_Node2`) and the combined
serial list. The launcher (`ndsp/k10cr1/launcher_multi_rotor.py`) reads
`sn_list`/`nickname_list` from that entry, so the launcher, driver, and
`subroutines/k10cr1_functions.py` are all unchanged. In single-node mode
Base publishes `k10cr1_ndsp` as a lazy node-aware proxy
(`_NodeAxisK10CR1Proxy`): reused legacy code that addresses axes by their
bare standalone names ('852_HWP', ...) is resolved to the selected node's
nicknames, already-suffixed names pass through unchanged, and the real
sipyco client - which connects the moment the device is created - is
fetched only on the first actual call, so experiments that never rotate
never connect. Do not use `best_effort` for the rotators (silent `None`
returns would drop motion commands). After editing a deployed database, run `artiq_client
scan-devices` (no gateware impact). Rigol and similar devices remain
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

Manual control is node-split: `AOMsCoils_master_satellite_mixin.py` holds the
private shared implementation, and the public experiments are
`AOMsCoils_master_satellite_Node1` and `AOMsCoils_master_satellite_Node2`,
separate from GVS and from the untouched standalone `AOMsCoils.py`. Build
registers the suffixed two-node device superset without configuring execution
or loading datasets. Prepare configures Base in `single_node` mode for the
class's fixed node and loads only that node's variables plus globals.

| Experiment | Selected node | Other node |
|---|---|---|
| `AOMsCoils_master_satellite_Node1` | explicit OFF/0 V pass, then honor controls | never initialized or written |
| `AOMsCoils_master_satellite_Node2` | explicit OFF/0 V pass, then honor controls | never initialized or written |

The per-node experiment never assigns the other node an output state, and it
calls Base initialization with `reset_core=False` so outputs already
established on the other node through the shared RTIO fabric are preserved.
Making a node safe therefore means running that node's own experiment with its
controls unchecked; there is no cross-node force-safe gate.

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

Node2 additionally carries the standalone combined `Node2_GRIN1_AOM_ON` /
`Node2_GRIN2_AOM_ON` modes, applied after the individual controls exactly as
in the standalone ordering, so on Node2 they are authoritative: their OFF
branches override `GRIN1and2_ON`, the GRIN gate controls, and the D1 DDS
(`D1_pumping_DP_ON` on Node2 therefore effectively drives only
`ttl_D1_pumping`, matching the original). The GRIN2 mode re-programs
`dds_D1_pumping_DP` to `f/p_GRIN2_excitation` before switching it on.

Coil control matches the standalone utility: a single `disable_coils`
checkbox (default False) drives the selected node's four coils to the
persistent `*_volts_MOT` calibration values, or to 0 V when set. Voltages are
read from the authoritative node datasets at prepare, never entered in the
GUI, and never persisted. Repository examination is unaffected because the
argument is a plain boolean.

| Coil | Node1 Zotino | Node2 Zotino |
|---|---:|---:|
| AZ bottom | 0 | 0 |
| AZ top | 1 | 1 |
| AX | 13 | 2 |
| AY | 14 | 3 |

Microwave/RF defaults OFF, obeys hard gating, and requires confirmation.

K10CR1 780/852 waveplate rotations are ported from the standalone utility
(section 21 has the controller design): thirteen booleans in the 780/852 GUI
groups drive this node's axes only, via node-suffixed axis names. The 852
target moves are gated on the 852 booleans, fixing the standalone copy-paste
bug that gated them on the 780 booleans.

`run_laser_feedback` mirrors the standalone utility: only when every fiber
AOM (A1-A6) and the cooling DP are requested on does the utility run the
laser stabilizer (or `monitor()` when the box is unticked), then re-apply the
requested DDS and optical-gate states. Feedback results persist to the
node-suffixed datasets through the redirection layer (section 18). Old
independent-Kasli networking, Rigol, and other USB devices remain excluded.

Safe startup now switches each requested-off DDS OFF immediately after
`init()` and before attenuation/profile programming. Base still initializes
the selected node's CPLDs, DDS channels, TTL directions/states,
Samplers/gains, and Zotino outputs (the unselected node's groups are removed
from the lifecycle), so it is not passive. The utility performs a second
explicit OFF/0-V pass before applying requested state.

Commit `078dc0e` implemented the original combined utility; commit `a7d8dff`
split it into the per-node experiments and added the Base `reset_core` option.
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

1. submitted `selected_node` (single-node GVS) is strictly validated; the GVS
   execution mode and the AOMsCoils node are fixed per class;
2. Base configures `single_node Node1`, `single_node Node2`, or `two_nodes`;
3. only the required node variables plus globals are loaded for GVS;
4. AOMsCoils loads only its own node's namespace plus globals; the other node
   is out of scope;
5. selected-node compatibility and DDS bindings are published;
6. missing datasets fail clearly before hardware initialization.

For deferred single-node configuration, both nodes remain registered but the
unselected node's CPLD, Sampler, Zotino, TTL, and safety-state lists are cleared
from the hardware lifecycle. Thus superset registration does not broaden
single-node hardware access.

GVS metadata uses fixed code defaults where persistent values cannot safely be
read during examination: `n_measurements = 100` and the normal scan
presentation. Each public GVS file has a fixed mode, so its function selector
exposes only that mode's registry.

### Explorer naming and public-class imports

ARTIQ uses the first line of each public `EnvExperiment` class docstring as the
Explorer label. The expected labels are:

```text
ExperimentVariables_master_satellite_Node1
ExperimentVariables_master_satellite_Node2
ExperimentVariables_master_satellite_global
GeneralVariableScan_master_satellite_single_node
GeneralVariableScan_master_satellite_two_nodes
AOMsCoils_master_satellite_Node1
AOMsCoils_master_satellite_Node2
```

Do not publicly import one dashboard `EnvExperiment` class into another
repository file; a public import makes ARTIQ rediscover the class as a
duplicate experiment (historically `GeneralVariableScan_master_satellite1`).
Shared implementations therefore live in the `_mixin` modules as private
classes, and every public experiment composes them in its own file. The
namespace sanity function remains a GVS callable, not an Explorer experiment.

Commits `a3eb7fe` and `f06f9be` implemented these lifecycle/discovery fixes;
`025c027` split GVS by execution mode, and `a7d8dff` split manual control by
node. The combined hardware-free suite passes 58 tests, including every public
build with `None` arguments and empty dataset storage, strict runtime
missing-dataset failures, submitted-GUI `n_measurements` precedence, and
optional underflow-retry behavior.

Implemented and hardware-free tested on branch
`20260826_master_satellite_codex`:

- legacy selector removal;
- explicit physical mapping;
- resolver and Base lifecycle;
- destination readiness and fixed SPCM policy;
- Node1/Node2/global variables;
- authoritative loading/projection/reload/target/cache APIs;
- mode-specific registries and native namespace sanity;
- mode-split GVS files with built-in optional underflow retry
  (`enable_Catch_UnderFlow`);
- run-local GUI-argument precedence over persistent datasets
  (`n_measurements`);
- selected-node magnetometer wiring and suffixed results;
- deterministic per-node manual AOM/DDS, TTL optical-gate, microwave/RF, and
  `disable_coils` MOT-coil control;
- single-node laser-feedback compatibility (`run_laser_feedback`) through
  unchanged aom_feedback code with node-suffixed persistence;
- K10CR1 waveplate control through one `k10cr1_ndsp` controller with
  node-suffixed axis names;
- the single-node FORT polarization optimizer
  (`FORT_Polarization_Optimizer_master_satellite`, `which_node = node1 |
  node2`), reusing the unchanged legacy record/rotation helpers through the
  dataset redirect (`POLARIZATION_RESULT_DATASETS`) so monitors and the
  `best_852*` calibrations store node-suffixed; its two-node design is
  documented in the file and waits on parallel FORT feedback (section 18);
- derived RF amplitudes and the traced single-node atom-result
  datasets/buffers in Base, plus `print_async`/`write_results`;
- `MicrowaveScanOptimizer_master_satellite` (`which_node = Node1 | Node2`)
  with node-suffixed microwave calibrations and `health_check_uw_*`
  fidelities and a self-resubmission expid carrying `which_node`
  (section 19);
- `_mixin` naming for shared non-experiment modules;
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
6. validate on hardware the derived amplitudes and single-node result state
   now provided by Base;
7. validate active AOM/FORT feedback on hardware and extend it to GVS runs;
8. validate normal atom loading;
9. validate `MicrowaveScanOptimizer_master_satellite` on hardware (Node1
   first; Node2 exercises DMA and remote SPI over DRTIO);
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
