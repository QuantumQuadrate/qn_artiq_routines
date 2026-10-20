# Master-Satellite System: Detailed Specification, Status and History

File `plan_claude_detail.md`; companion brief `plan_claude_concise.md`.
State described: 2026-10-16. Branch `20261006_two_nodes_experiment`; HEAD is
the commit that last updated this file (unpushed; origin is
behind). `main` holds everything through `19dcee9` via PR #68 (`94d0ef6`,
merged 2026-10-11). Committed on 2026-10-16 at the user's request: the
MicrowaveScanOptimizer fixes (`619d3ec`; section 11.1, O-17), the magnetometer
monitoring work together with the user's delay edits in
`subroutines/experiment_functions.py` (`bfca64d`; O-34, O-35), this plan
consolidation (`a369392`), and Node2's A1 photodiode move to `sampler4` ch6
with the MOT beam power check (sections 10.1, 16.4). The analysis notebooks
stay uncommitted until the user commits
them. Run `git status` rather than trusting this paragraph. The ARTIQ master
runs experiments straight from the working tree, so uncommitted code is not
undeployed code (rule R26).

---

## 1. Purpose, authority and precedence

### 1.1 What this file is

The repository controls two neutral Rb-87 single-atom quantum-network nodes
with ARTIQ 7 / Sinara. The two formerly independent Kasli-SoC systems now run
as one DRTIO master-satellite system, while the original standalone software
stack stays usable for rollback. This file is the single authoritative
specification, status record and history of that system. It is written for
two readers: a new AI assistant starting work, and a lab member who needs to
know what the system does and why.

### 1.2 Precedence

For instructions and design intent, in order:

1. the user's latest explicit instruction;
2. this file;
3. `plan_claude_concise.md`;
4. code comments and docstrings;
5. assumptions, which must be reported before they are used.

For facts about what currently runs, the code at HEAD and hardware evidence
(h5 files, RIDs, logs) describe reality. When this file and the code
disagree, do not silently follow either: report the discrepancy, establish
which side is wrong, and fix that side (and this file, if the design moved).
Several code comments are known to be stale (section 18, item O-25).

### 1.3 How the two files relate

The concise file is the startup brief: what the system is, the rules that bite
most often, a one-screen status, and pointers into this file by section
number. It must never contain a fact that this file lacks. If they conflict,
this file wins and the concise file is corrected.

### 1.4 Retired files

These four files are replaced by the pair above and were deleted on
2026-10-14 (committed with this file on 2026-10-16); they live
only in git history, last at `8963734` (`PLAN.md`), `fd10950`
(`PLAN_FROM_SCRATCH.md`), `0e54c84` (`plan_codex_detail.md`) and `1a1bb6c`
(`plan_codex_consise.md`).
Every still-valid item from them has been carried over; superseded claims are
recorded in section 19.

| Retired file | What it was |
|---|---|
| `PLAN.md` | Titled "AGENTS.md / Master-Satellite Migration Context"; the original Codex agent brief (1457 lines, 35 sections, last changed 2026-08-30). Rules, architecture, a 53-step phased migration order. |
| `PLAN_FROM_SCRATCH.md` | A compressed re-derivation of the migration plan (341 lines, 2026-08-30; one paragraph added 2026-10-07). Added IP addresses, the 9-step flash gate, the 9-step first-validation sequence and a 16-step order. |
| `plan_codex_detail.md` | The authoritative detailed spec until now (908 lines, 27 sections, last changed 2026-10-12). |
| `plan_codex_consise.md` | The authoritative startup brief until now (230 lines, last changed 2026-09-03; the filename typo "consise" was never fixed). |

### 1.5 Keeping this file current

- Update this file in the same change that moves the design. Put the date and
  commit next to new facts.
- Never delete a superseded statement silently: move it to section 19 with its
  date and reason.
- Keep section numbers stable, because code cites them. Add subsections at the
  end of a section rather than renumbering.
- Cite code by file plus function, class or constant name, not by line
  number (line numbers drift; one citation into the retired detail file
  drifted by 31 lines in five days).
- Hardware claims need evidence: a RID, a log, or a commit body.
- Mark uncommitted work as uncommitted and hardware-unverified work as such,
  and keep the two apart: the master runs the working tree, so uncommitted
  code may already have run (R26).
- After editing, check that `plan_claude_concise.md` still agrees.

### 1.6 Citation map for code that cites the retired files

Repointed on 2026-10-14 (committed with this file; O-24), comment and docstring text
only. Kept so a reader holding an old citation can find its new home:

| Citing site | Old citation | New target |
|---|---|---|
| `FORT_Polarization_Optimizer_master_satellite.py` module docstring ("FUTURE REFERENCE - two_nodes design"), TWO sites | "see plan_codex_detail.md section 18" (header of the design block) and a bare "once section 18 allows it" (its feedback bullet; no filename, so it would now silently resolve to this file's section 18, Open issues) | both to `plan_claude_detail.md` section 10 (Feedback scope), especially 10.5 and 10.8 |
| `utilities/BaseExperiment_master_satellite.py` `_install_wiring_metadata` docstring | "plan_codex_detail.md already fixes that policy for measurements_progress" | section 6.9 (Result state and dataset naming) |
| `tests/measure_drtio_ttl_skew.py` module docstring | `plan_codex_detail.md:569` and `PLAN.md:945-952` | section 12.2 (Remote input and TTL-timestamp latency). The `:569` citation was already stale: that line had drifted to the microwave-optimizer text. |
| `standalone/README.md` (last paragraph) | `plan_codex_detail.md` / `plan_codex_consise.md` | `plan_claude_detail.md` / `plan_claude_concise.md` |
| memory `gui-arguments-win-over-datasets.md` | "invariant 11 in plan_codex_consise.md" | rule R13 (section 3.2) |
| memory `master-satellite-mixin-conventions.md` | names the `plan_codex_*` files authoritative | section 1.2 of this file |

### 1.7 Conventions used here

- "mu" is machine units: 1 ns at `ref_period` 1 ns. The RTIO coarse clock is
  125 MHz (8 ns).
- ExperimentVariables store times in seconds; the `unit` kwarg only scales the
  GUI display.
- `NodeX` means `Node1` or `Node2`. "Idle node" means the node a single-node
  run does not select.
- Dates are absolute. Commit hashes are short hashes on this repository.

---

## 2. Physical system and hardware reference

### 2.1 Nodes and link

| | Node1 | Node2 |
|---|---|---|
| Name | Alice (`which_node = "alice"` in legacy code) | Bob (`"bob"`) |
| DRTIO role | master | satellite, destination 1 |
| Address | 192.168.1.75 (`core_addr`) | .76 was its standalone address; silent by design now (satman has no network stack) |
| Board | Kasli-SoC v1.1 since 2026-09-09 (MAC fc-0f-e7-0f-69-13) | Kasli-SoC, v1.1 per its JSON |
| Firmware | `boot.bin` sha 2840a2f2, ident `7+ea380abc;drtio_master` | `boot.bin` sha c1412ddc |
| Console | network plus UART | UART only: FTDI FT4232H (VID_0403 PID_6011), third port (COM8), 115200 8N1 |

- Link: fiber SFP0 (master) to SFP0 (satellite), 125 MHz RTIO, 2.5 Gbit/s.
  The master log on 2026-09-09 showed `[LINK#0] ... remote replied`,
  `[DEST#1] destination is up`, `buffer space is 128`.
- Clock: `rtio_clock = ext0_synth0_125to125` (external lab 125 MHz).
- Host ARTIQ: 7.8143.75d75cc (conda env `artiq7`). Gateware ident
  `7+ea380abc`.
- The original Node1 board (v1.0) had a dead GTX/MGTREFCLK path and is the
  RMA item (section 17.2).

### 2.2 Peripherals (identical per node)

Two DIO cards with edge counters (each card: low bank input, high bank
output, `edge_counter` true), three Samplers, one Zotino, three AD9910 Urukuls
(`synchronization` true, `clk_sel = 2`). RTIO frequency 125 MHz.

### 2.3 Hardware reference files

`master_satellite_test_files/`:

```text
node1_master_kasli_soc_20260715.json       variant drtio_master, base master, hw_rev v1.0
node2_satellite_kasli_soc_20260715.json    variant drtio_satellite, base satellite, hw_rev v1.1
device_db_two_nodes_20260715.py            generated unified database
kasli-soc-standalone_node1_with_edgecounters.json   differs from the master JSON only in base/variant
                                                    (and it has no core_addr)
```

- Do not modify these files. Both images were built from these JSONs and run
  on hardware, and the deployed database is byte-identical to
  `device_db_two_nodes_20260715.py` (verified with `cmp`, commit `1382a49`).
- Documented exception: a host-side `k10cr1_ndsp` controller entry was
  appended on 2026-08-31 on explicit request (section 13).
- `hw_rev` in the JSON does not affect `kasli_soc` gateware: the unchanged
  "v1.0" master image runs on the v1.1 replacement board. The standalone JSON
  differing only in base/variant proved `hw_rev` was not involved in the
  bring-up clock failure.
- Destination 1 is not a field in the satellite JSON; it comes from the
  combined database generation and routing configuration.
- The folder README names a non-existent `device_db_nodes12_20260715.py`; the
  real file is `device_db_two_nodes_20260715.py`.

### 2.4 Deployed device database

- `artiq-master/device_db.py` is a symlink to
  `repository/qn_artiq_routines/device_db/device_db_master_satellite.py`
  (tracked since `1382a49`, 2026-09-10). 166 entries (163 local, 3
  controllers; counted 2026-10-14); 57 RTIO channels at
  destination 0 and 57 at destination 1; `core_addr` 192.168.1.75;
  `ref_period` 1 ns; target `cortexa9`.
- `artiq-master/device_db_standalone.py` is kept for rollback (symlink to
  `device_db_node1_with_edgecounters_calibrated.py`).
- Channel encoding: `channel >> 16` is the destination. Node1 channels are
  0x000000-0x000038, Node2 channels 0x010000-0x010038. Examples:
  `ttl0` 0x000000, `ttl0_counter` 0x000008, `ttl16` 0x010000,
  `ttl16_counter` 0x010008, `spi_urukul1` 0x00002b, `spi_urukul4` 0x01002b
  (= 65579), `spi_zotino1` 65569, `ttl_urukul4_sw1` 0x01002f. Use this to
  decode channel numbers in RTIO error messages.
- `led0`, `led1` are absent from the unified database. `ccb` is a worker
  virtual device, not a database entry.
- Generated SPI/CNV/LDAC/CLR/switch/sync/IO-update/EEPROM devices are
  database implementation details unless active code references them
  directly. Search before excluding one from a required translation. Urukul
  RF switches are reached through the channel objects (`.sw`).
- After editing a deployed database, run `artiq_client scan-devices` (no
  gateware impact).

### 2.5 Unified low-level namespace

| Subsystem | Node1 / master | Node2 / satellite |
|---|---|---|
| DIO TTLs | `ttl0`-`ttl15` | `ttl16`-`ttl31` |
| Counters bank 1 | `ttl0_counter`-`ttl3_counter` | `ttl16_counter`-`ttl19_counter` |
| Counters bank 2 | `ttl8_counter`-`ttl11_counter` | `ttl24_counter`-`ttl27_counter` |
| Samplers | `sampler0`-`sampler2` | `sampler3`-`sampler5` |
| Zotino | `zotino0` | `zotino1` |
| CPLDs | `urukul0_cpld`-`urukul2_cpld` | `urukul3_cpld`-`urukul5_cpld` |
| DDS | `urukul0_ch0`-`urukul2_ch3` | `urukul3_ch0`-`urukul5_ch3` |

Sequential Python names do not imply sequential encoded RTIO channels.

### 2.6 Canonical SPCM detectors

All four active photon detectors are master-local Node1 inputs, whichever node
an experiment runs:

| Detector | TTL | Counter | Old standalone Node1 name |
|---|---|---|---|
| `SPCM_H1` | `ttl0` | `ttl0_counter` | SPCM0 |
| `SPCM_V1` | `ttl1` | `ttl1_counter` | SPCM1 |
| `SPCM_H2` | `ttl8` | `ttl8_counter` | SPCM0_OtherNode |
| `SPCM_V2` | `ttl9` | `ttl9_counter` | SPCM1_OtherNode |

H1/V1/H2/V2 are detector (optical) identities, not Kasli ownership. Node2's
old SPCM-looking inputs (`ttl16`/`17`/`24`/`25` and counters) remain physical
devices but are not active detectors. Each detector sees light from both
traps, which is why joint thresholds are global (section 7.4). Binding rules
are in section 5.8.

### 2.7 Per-node TTL roles

From `BaseExperimentMasterSatellite.NODE_TTL_ALIASES` (standalone-local
indices; Node2's resolve to `ttl16`-`ttl31`):

| Alias | Node1 | Node2 |
|---|---|---|
| `ttl_microwave_switch` | ttl4 | ttl4 |
| `ttl_repump_switch` | ttl5 | ttl5 |
| `ttl_exc0_switch` | ttl6 | ttl6 |
| `ttl_pumping_repump_switch` | ttl7 | ttl7 |
| `ttl_GRIN2_switch` | ttl13 | ttl13 |
| `ttl_GRIN1_switch` | ttl14 | ttl14 |
| `ttl_D1_pumping_switch` | ttl15 (since 2026-10-19; blocked at init, no sequence opens it yet) | ttl12 (real gate) |
| dummy input / output (no alias name of their own) | ttl3 / ttl12 | ttl8 / ttl15 (= ttl24 / ttl31) |
| legacy handshake | `ttl_node2_input1` ttl3, `ttl_node1_output1` and `_output2` ttl12 (dummies) | `ttl_node1_input1` and `_input2` ttl8, `ttl_node2_output1` ttl15 (dummies) |

- Inputs: Node1 ttl3, its dummy (the group must not be empty or the kernel
  does not compile); Node2 ttl0-3, ttl8.
- Safe state ON (active-low optical gates blocked): Node1 ttl4-7, 13, 14,
  15; Node2 ttl4-7, 12, 13, 14. Safe state OFF: Node1 ttl12; Node2 ttl15.
- The gateware configures only ttl4-7 and ttl12-15 as outputs on each node;
  ttl2/ttl10 sit on input banks.
- The handshake cables were removed on 2026-10-19 (user). They ran master
  ttl12 -> satellite ttl25, master ttl15 -> satellite ttl26 and satellite
  ttl31 -> master ttl11, plus BNC T loopbacks ttl12 -> ttl10 and ttl31 ->
  ttl27 (section 12.2). Each legacy handshake alias, and each single-node
  `COMPATIBILITY_DUMMY_TTL_ALIASES` entry, now points at its own node's
  dummy, which has nothing connected.
  Before this, single-node Node2's `ttl_node1_output1/2` were on ttl5, its
  repump switch. Node1 ttl10, ttl11 and Node2 ttl9, ttl10 (ttl25, ttl26)
  are now referenced only by the bind-everything `TTLS` and `EDGE_COUNTERS`
  tuples. device_db is unchanged: renaming a key there makes `build()` fail
  for every experiment.
- Node1 ttl15 became Node1's `ttl_D1_pumping_switch` later on 2026-10-19 (user);
  it used to be the dummy ttl3, set in `COMPATIBILITY_DUMMY_TTL_ALIASES`
  and so present only in single-node mode. It now exists in both modes,
  as `ttl_D1_pumping_switch` and `ttl_D1_pumping_switch_Node1`. Nothing writes it yet:
  Node1's optical pumping gates D1 with `dds_D1_pumping_DP`, only the
  `*_node2` pumping functions use the alias and Node1 runs dispatch to the
  `*_node1` ones; AOMsCoils Node1 opens it when D1 is requested (section
  11.7). Init blocks it (safe state ON), like Node2's ttl12, at the user's
  request that master-satellite init turn off everything standalone does.
  ttl15 does gate Node1's D1 (external source since 2026-10-19), so every
  Node1 sequence that pumps gets NO D1 light until it opens ttl15: the runs
  complete and the data look normal, but the atoms are not pumped.
  Standalone (`utilities/BaseExperiment.py`, the user's edits) matches:
  alice has D1 on ttl15 blocked at init, output dummies on ttl12, input
  dummies on ttl3; bob's dummies are its ttl8 (input) and ttl15 (output),
  node-local names for ttl24/ttl31. The `dummy_ttl_input`/`_output`
  aliases added earlier that day were removed again: nothing uses them.
- `ttl_D1_pumping` was renamed `ttl_D1_pumping_switch` on 2026-10-19 (user),
  in both stacks and every file that used it (60 occurrences in 12 files,
  including `ttl_D1_pumping_Node1/_Node2` -> `ttl_D1_pumping_switch_Node1/
  _Node2`). `subroutines/experiment_functions.py` and
  `subroutines/cost_functions.py` carry the Windows read-only attribute;
  it was cleared for that edit and set again, origin unknown.

### 2.8 DDS channel map

Logical alias (alice/bob `ALIAS_MAP`) to unified device, all 24 verified
against both mapping layers:

| Logical control | Alias | Node1 | Node2 |
|---|---|---|---|
| FORT | `dds_FORT` | `urukul0_ch0` | `urukul4_ch3` |
| cooling DP | `dds_cooling_DP` | `urukul0_ch1` | `urukul4_ch2` |
| AOM A1 | `dds_AOM_A1` | `urukul1_ch2` | `urukul3_ch0` |
| AOM A2 | `dds_AOM_A2` | `urukul1_ch0` | `urukul3_ch1` |
| AOM A3 | `dds_AOM_A3` | `urukul1_ch1` | `urukul3_ch2` |
| AOM A4 | `dds_AOM_A4` | `urukul2_ch0` | `urukul3_ch3` |
| AOM A5 | `dds_AOM_A5` | `urukul2_ch1` | `urukul4_ch0` |
| AOM A6 | `dds_AOM_A6` | `urukul1_ch3` | `urukul4_ch1` |
| D1 pumping DP | `dds_D1_pumping_DP` | `urukul2_ch2` | `urukul5_ch0` |
| GRIN1and2 | `GRIN1and2_dds` | `urukul0_ch3` | `urukul5_ch2` |
| microwaves | `dds_microwaves` | `urukul2_ch3` | `urukul5_ch3` |
| MW RF | `dds_MW_RF` | `urukul0_ch2` | `urukul5_ch1` |

The same alias names mean physically different things on the two nodes
(recorded in `documentation/upgrade_notes/D1_optical_pumping_merge_ideas.md`):

- `GRIN1and2_dds`: Node1 = one DDS feeding both GRINs (excitation, defaults
  `f/p_excitation`); Node2 = GRIN1 only (D1 pumping, defaults
  `f/p_GRIN1_D1_pumping`).
- `dds_D1_pumping_DP`: Node1 = a real D1 DP AOM at 373.075 MHz; Node2 =
  GRIN2 (defaults `f/p_GRIN2_D1_pumping`; the Node2 combined GRIN2 mode
  reprograms it to `f/p_GRIN2_excitation`).
- Node1 D1 shares the GRIN1 optical gate; Node2 has a dedicated D1 gate
  (`ttl_D1_pumping_switch`, ttl12).

### 2.9 Fixed wiring metadata

`BaseExperimentMasterSatellite.SINGLE_NODE_WIRING_METADATA` (fixed wiring
facts, not ExperimentVariables or scan state):

```text
                          Node1             Node2
coil_names                [AZ bottom, AZ top, AX, AY]
AZ_bottom_Zotino_channel  0                 0
AZ_top_Zotino_channel     1                 1
AX_Zotino_channel         13                2
AY_Zotino_channel         14                3
coil_channels             [0, 1, 13, 14]    [0, 1, 2, 3]
UV_trig_channel           [8]               [8]
Osc_trig_channel          [10]              [10]
FORT_MM_sampler_ch        7                 7
GRIN1_sampler_ch          4                 4
Magnetometer_X/Y/Z_ch     1 / 2 / 3         1 / 2 / 3
```

Published bare in single_node mode and `_NodeX`-suffixed in two_nodes mode
(section 6.9). The magnetometer is read from the node's `sampler2` (Node1
physical `sampler2`, spi 30 / cnv 32; Node2 physical `sampler5`, spi 65566 /
cnv 65568). The FORT MM monitor is `sampler1_NodeX` ch 7; the FORT APD is
`sampler1_Node1` / `sampler0_Node2` ch 6 (a legacy alice/bob asymmetry).

### 2.10 Urukul cards and the SYNC fix

- 2026-09-15: Node2's `urukul4` and `urukul5` were replaced with spare v1.5 /
  v1.5.6 cards (proto_rev 8). `artiq_sinara_tester` passed first time and the
  SYNC-window sweep (RID 38432) shows no marginal channel. New seeds: urukul4
  13/15/15/14, urukul5 18/18/17/19. History and measurement rules: section
  12.8 and `documentation/debug_notes/2026-09-10_urukul_sync_window_margin.md`.
- Naming: `urukul4` = Node2's old standalone `urukul1`; `urukul5` = old
  `urukul2`. The Node2 microwave DDS is `urukul5_ch3`.
- Do NOT use the new v1.6 Urukul cards: they carry an iCE40 with proto_rev 9
  gateware, and ARTIQ 7 accepts only proto_rev 8 ("Urukul proto_rev
  mismatch"). v1.6 support arrived in ARTIQ 9. Whether a proto_rev 8 build
  for v1.6 exists is unknown (ask M-Labs).
- A SYNC fault cannot be fixed in software: `AD9910.init()` calls
  `tune_sync_delay`, and Base calls `init()` on every DDS channel in every
  run. SYNC cannot be disabled (`sync_delay_seed = -1`) because
  `experiment_functions.py` sets `PHASE_MODE_TRACKING`/`ABSOLUTE` on
  `dds_microwaves` and `dds_MW_RF`: 11 live `set_phase_mode` calls in 9
  functions on 2026-10-14 (the CORPSE map and parity 4-11 including
  6_node2; parity 3's call is commented out). The earlier figure of "21
  sites" is stale. Freezing seeds cannot create margin.

### 2.11 Rollback artifacts and boot facts

- Master image and boot log: `artiq-master/firmware_master_satellite/`
  (`boot.bin` sha 2840a2f2, `boot_log_after_power_cycle.txt`).
- Standalone rollback image (shared by both nodes for recovery):
  `standalone-kasli_node1_with_edgecountersboot_20241011.bin` in OneDrive
  `.../ARTIQ gatewares/.../old_boot_files/preston/`. Swapped-role images
  (Node2 master / Node1 satellite) are in OneDrive
  `.../Master Satellite/node2_1_master_satellite/`.
- `artiq-master/device_db_standalone.py` and
  `node1_artiq_conda_packages_2024-06-25.txt` are kept.
- The SD card needs only `boot.bin` and `CONFIG.TXT` (`ip=192.168.1.75`,
  `rtio_clock=ext0_synth0_125to125` on Node1). A copied `boot.bin` may be
  read-only; clear `IsReadOnly`.
- The ARTIQ-Zynq `boot` config key is write-only: rollback must never depend
  on reading the current image back.
- `artiq_coremgmt config write -f boot <boot.bin>` works remotely and applies
  at the next power cycle; a `TimeoutError` may still have completed.
  `artiq_coremgmt reboot` is NOT implemented on this Kasli-SoC build
  (UnrecognizedPacket): power-cycle physically.
- Not recorded anywhere: byte-for-byte SD-card backups and a single release
  bundle with checksums (the retired gate asked for both).
- The archive location of the satellite `boot.bin` (sha c1412ddc) is not
  recorded: `artiq-master/firmware_master_satellite/` holds only the master
  image and its boot log. Locate it and archive it next to the master image.

### 2.12 Operating notes for ARTIQ on this machine

- `artiq_client submit` needs `-R` here; without it the RID appears, goes to
  "preparing" and vanishes with no h5 and no error (`artiq_client ls` printing
  `C:\` and `D:\` is the tell). Working form, options before the file:
  `artiq_client submit -R -c <Class> qn_artiq_routines/<File>.py
  'experiment_function="..."' 'create_applets=False' 'n_measurements=20'`.
  ARTIQ 7 subcommands: submit, delete, set-dataset, del-dataset, show,
  scan-devices, scan-repository, ls. There is no `log`; `show log` streams
  and blocks.
- Worker tracebacks are not reachable from the CLI. Use
  `artiq_run --dataset-db <COPY of dataset_db.pyon> -c <Class> <file> <args>`
  (drives real hardware; cannot take the core while `artiq_master` runs).
- `artiq_ctlmgr` must run alongside the master for moninj and controllers.
  Its `aqctl_corelog` sends PullLog on connect, which replays and drains the
  core log (WARNING and above, no device timestamps), so old boot errors look
  current; `artiq_coremgmt log` is then empty. Use `aqctl_corelog -v` for
  INFO.
- From a non-conda shell, prepend
  `C:\ProgramData\Anaconda3\envs\artiq7\Library\bin` to PATH (provides
  `ld.lld`).
- The master binds every port on `::1` and `127.0.0.1`; Python 3.10's
  Proactor permanently closes a listener after an accept error. Symptom:
  dashboard `ConnectionRefusedError` WinError 1225. Fix: restart
  `artiq_master` when the schedule is empty, or `artiq_dashboard -s
  127.0.0.1` (its state file is keyed by address, so the layout starts fresh).
- `artiq_client show schedule` streams; for one read use sipyco
  `Client('::1', 3251, None)`, `select_rpc_target('master_schedule')`,
  `get_status()`.
- A kernel hang (RID 38433, in a `cpld.init` loop) could not be stopped from
  the dashboard; killing the worker process (`artiq.master.worker_impl`)
  sufficed, no reset or power cycle. Root cause unknown; capture the ctlmgr
  ERROR lines if it recurs.
- Boot-race link errors: when both boards power-cycle together the master
  logs `[LINK#0] ping failed`, then one `error(s) found (0x03)` after link
  init (bit0 unknown packet type, bit1 truncated packet; sticky flags read
  once the link is up). Benign only if it does not recur. Worry if 0x03
  repeats after boot, the link re-pings mid-run, or bit2 (timeout getting
  remote buffer space) appears.
- For read-only h5 inspection use
  `C:/ProgramData/Anaconda3/envs/artiq7/python.exe -I` (has h5py).
- The master serves the repository WORKING TREE: `artiq_master` runs with no
  `--git` option (process command line checked 2026-10-14), and every h5
  `expid` carries `repo_rev` "N/A". Each submission loads the files as they
  are on disk at that moment, committed or not (rule R26). Example: the
  MicrowaveScanOptimizer fix saved at 19:44 on 2026-10-14 ran in RIDs
  38942-38951 and 38961-38965 two days before it was committed (`619d3ec`).

### 2.13 Reflash and rollback procedure (for any future flash)

The 2026-09 gate is done, but the same steps apply to any future flash (for
example a rollback to standalone, or the return of the RMA board). Only with
explicit authorization (R3):

1. Capture and verify full SD-card images, or the exact `boot.bin` currently
   on each node (the `boot` key cannot be read back).
2. Preserve the active device databases and `CONFIG.TXT`, startup/idle
   kernels, controllers, `dataset_db.pyon`, ARTIQ version and logs.
3. Record carrier/peripheral revisions and port wiring: the master
   downstream SFP port, the satellite upstream port (today SFP0 to SFP0) and
   the routing table. Verify the SFP transceivers and fiber for the 125 MHz /
   2.5 Gbit/s link.
4. Pin the exact ARTIQ and ARTIQ-Zynq revisions; build or obtain matched
   master-runtime and satellite-satman images; regenerate the unified
   database from the exact JSONs (Node2 at destination 1); never hand-edit
   channel numbers.
5. Archive images, JSONs, database, routing, build logs, revisions and hashes
   together. Verify that the rollback media boot, where practical.
6. Staging: stop services. Node1 can be staged remotely (`artiq_coremgmt
   config write -f boot <boot.bin>`, applied at the next power cycle). Node2
   cannot: under satman it has no network stack, and ARTIQ 7
   `artiq_coremgmt` has no option to address a satellite (only
   `--device-db` and `-D/--device`), so write Node2's `boot.bin` to its SD
   card physically. (In 2026-09 remote staging worked for both only because
   Node2 was still standalone at .76.) Do not reboot the satellite
   prematurely; power down; connect the verified link; boot with UART
   monitoring (satellite, then master); verify .75 and destination 1 before
   starting the dashboard.
7. Validate passively first: master boot/network, both UART logs, a bounded
   kernel that only calls `core.get_rtio_destination_status(1)`, host-only
   namespace resolution, one harmless local and one harmless remote RTIO
   operation (`artiq-master/drtio_passive_check.py` does the last steps with
   no outputs), then AOMsCoils per node with every control off, then one AOM
   and coil at a time. Do not use GVS as the first check (section 6.8).
8. After any gateware, clock or card change, recalibrate AD9910 sync first:
   run `artiq_sinara_tester` (it writes the EEPROM sync seeds), then the
   Urukul SYNC sweep (section 12.8) and `tests/measure_drtio_event_cost.py`.
   Re-run the event-cost probe after any link change.

---

## 3. Non-negotiable rules and working conventions

### 3.1 Priorities

In order: (1) working hardware control, (2) easy debugging, (3) minimal-risk
changes, (4) preserving validated experiment sequences, (5) easy rollback to
the standalone hardware and software. All of these come before architectural
cleanup. Do not redesign unrelated parts of the repository because a cleaner
abstraction is possible.

### 3.2 Rules

Rule IDs are stable; cite them as "rule R7".

- **R1. The standalone stack stays functional.** Rollback to standalone
  gateware must keep working with `standalone/` plus the shared modules.
  Never create a dependency from the standalone stack onto master-satellite
  code. Shared modules (section 4.3) may receive only backward-compatible
  extensions whose default reproduces standalone behaviour, and only when
  the user authorizes it. Do not gradually turn standalone files into the new
  architecture.
- **R2. Commit and push only on explicit instruction.** The user verifies
  changes on hardware first. Never rewrite history.
- **R3. No flashing and no operating outputs without explicit,
  phase-specific authorization.**
- **R4. Preserve unrelated worktree changes.** Check `git status` before
  editing; notebooks may be modified on disk by a live Jupyter session.
- **R5. Inspect first, report evidence, separate facts from assumptions.**
  Passing Python or static checks is not hardware proof. Flag everything that
  needs experiment-computer validation. When hardware semantics are
  uncertain, say so instead of guessing.
- **R6. Explain before structural edits.** Before adding, renaming or moving
  files, imports or modules, give a one-paragraph rationale plus alternatives
  with trade-offs and let the user choose. File layout has operational
  meaning (Explorer discovery). Plain content edits inside an agreed change
  need no preface.
- **R7. Document large changes first.** When reconnaissance uncovers a
  cross-cutting change (hundreds of sites, a device-config migration), do not
  edit: write a queued note in `documentation/upgrade_notes/` (named by
  topic, no date prefix; record edit sites, decisions taken, what makes it
  safe or dangerous, prerequisites, latent bugs left alone) and add a row to
  that folder's README index. Post-mortems go to `documentation/debug_notes/`
  as `YYYY-MM-DD_topic.md`, keeping the wrong turns.
- **R8. Pilot one site before bulk edits.** For repetitive multi-site
  changes, change one named site, verify it (parse, AST diff, unit tests,
  compile sweep, hardware if possible), report, and stop to ask before the
  next. Method: section 3.4.
- **R9. Keep experiment sequences linear.** An experiment sequence stays one
  top-to-bottom function in the legacy `experiment_functions.py` shape; do not
  decompose it into mini-helpers. New machinery goes in the binding layer
  (Base), not between the reader and the pulse order. A linearity break is
  allowed only when unavoidable (for example per-node durations) and carries
  a comment. (User, 2026-10-07: "I want to follow the previous experiment
  format and just make it compatible.")
- **R10. Never calculate DRTIO channels, destination bits, or Urukul/TTL/RTIO
  offsets.** Mappings are exact dictionary entries derived from the generated
  database. Treat gateware and its database as one matched release. Never
  invent a `destination` argument for ordinary device access.
- **R11. Two-node code names every device's node explicitly**
  (`self.dds_FORT_Node2`, never a bare `self.dds_FORT`). Bare device names do
  not exist in two_nodes mode and fail to compile (section 6.6).
- **R12. Avoid new ExperimentVariables.** A newly declared Variable has no
  persistent dataset until its initializer runs, and until then EVERY
  master-satellite experiment fails in prepare (including AOMsCoils). The
  user does not want to run initializers ("Whatever you did breaks even
  AOMsCoils. revert it. i don't like to run the test code anyways.",
  2026-10-08). Preference order: (1) reuse an existing variable; (2) derive in
  Base and publish with `set_dataset(broadcast=True)`, recomputed per scan
  point; (3) a module-level constant; (4) a new Variable only after telling
  the user the cost. If one is unavoidable, hand-seed it on the master first
  (section 7.8). Removing a declaration is safe.
- **R13. Submitted GUI arguments win over persistent datasets for that run.**
  When prepare/run loads datasets over attribute names that are also GUI
  arguments, capture the submitted value first and re-assert it after every
  load, including queued-run reloads. Never persist GUI values back, except in
  the ExperimentVariables initializers (and the `n_measurements` broadcast
  of section 6.9).
- **R14. SED: at most 8 events may share a timestamp; the 9th is silently
  DISCARDED.** House rule: put `delay(5*us)` after every four timeline events
  (`dds.sw.on/off`, `ttl.on/off`, `dds.set`, `zotino.set_dac`, `set_att`), not
  after each one. Section 12.5.
- **R15. Restore slack in shared helpers, not at call sites.** Host RPCs in
  loops erode slack. Put `core.break_realtime()` inside the shared helper,
  outside any `which_node` guard, where the preceding timing is irrelevant,
  and after every `laser_stabilizer.run()`. Section 12.4.
- **R16. Never lazily import `subroutines/aom_feedback.py`** (or any module
  that captures `os.getcwd()` at import) inside prepare/run: the worker
  chdirs into `results/<date>/<hh>/` first, so `feedback_channels.json` is
  not found. Base imports `AOMPowerStabilizer` at module level; keep it there.
- **R17. Local write after broadcast deletes the dataset.** In ARTIQ 7
  `DatasetManager.set`, `set_dataset(..., broadcast=False)` on a key this run
  already broadcast deletes the key from the master database
  (`dataset_db.pyon`). Before any local write, check nothing earlier in the
  run broadcast that name.
- **R18. Leave historical experiments alone.** Do not migrate, copy, delete
  or rewrite historical independent-Kasli two-node functions or deprecated
  raw-TTL functions merely because they mention two nodes or raw TTLs.
- **R19. No unrelated cleanup, renaming or reformatting.**
- **R20. Applets change only on request.** (Many applet changes have been
  requested since 2026-09-10; section 14.)
- **R21. Probe Urukul SPI safely.** Use `bus.set_config_mu(flags=0,
  length=8, div=6, cs=0)` (one rtio_output, no pin moves). Never call
  `cpld.sta_read()` in a probe: it writes `cfg_reg` back from device_db
  defaults and drops that Urukul's RF switches. `core.reset()` drops TTL
  output states; diagnostics should use `break_realtime()`.
- **R22. Settled decisions not to reopen:** one `k10cr1_ndsp` controller for
  all eight rotators (never per-node controller processes); no `best_effort`
  on the rotators; keep the separate standalone
  `GeneralVariableScan_HealthCheck.py` / `_HealthCheck_CatchError.py` pair;
  `MicrowaveScanOptimizer_master_satellite` never gets a two_nodes mode.
- **R23. The legacy `which_node == "two_nodes"` mode never returns**
  (section 4.6).
- **R24. Do not build InfluxDB/Grafana or redesign the dashboard as part of
  this work** unless asked.
- **R25. Analysis code lives in `Analysis/h5_data_utilities.py`, not in
  notebooks**, and notebook edits respect the live Jupyter session
  (section 15.5).
- **R26. Saved means deployed.** The master serves the working tree (no
  `--git`; h5 `expid` `repo_rev` = "N/A"; section 2.12). Any saved edit
  under `repository/qn_artiq_routines` is what the next submitted experiment
  runs, committed or not, including edits an assistant makes while the user
  is submitting. Do not leave a half-finished or untested edit on disk while
  the user may be running experiments; finish it, or tell the user before
  saving. Not committing does not keep a change off the hardware, so
  "uncommitted" never means "not yet run": check the h5 files for runs that
  started after the file's modification time. (Recorded 2026-10-14.)

### 3.3 Scope control

Do not: redesign the whole repository; make standalone files generic;
rewrite `utilities/BaseExperiment.py` or the standalone GVS for
master-satellite; migrate every historical experiment or copy all of
`experiment_functions.py`; rewrite analysis code unnecessarily; drop
standalone compatibility; optimize unrelated physics; build an over-general
framework; guess DRTIO destinations, channel mappings or offsets; assume
local and remote timing are identical; duplicate physical mapping tables.

Prefer: separate stacks, one physical translation source, explicit mappings,
small reviewable diffs, easy rollback, easy hardware debugging, preserved
validated physics, small isolated hardware tests before large experiments.
Each failure should isolate to one node, device, mapping entry, destination,
controller or timing assumption.

### 3.4 Method for mechanical multi-site edits

Proven on the confirm-by-first-shot rollout (2026-10-13):

- Never hand-retype and never regex-sweep. Use an explicit hard-coded site
  list of `(line, expected_function, expected_first_line)` triples and a
  generator script.
- Proof is a REVERSE transform: `unwrap(new) == old` byte for byte. Show
  first that the proof fails on deliberate corruptions.
- Per-site assertions: indentation, owning function (from `ast`), enclosing
  measurement loop, no `break`/`continue`/`return` in the span. Global
  assertions: `ast.dump` identical for non-target functions; no line between
  a decorator and its `def` starts with `#`. (ARTIQ's
  `compiler/embedding.py` runs `re.search('def', line)` on the line after the
  decorators: a comment there gives `AttributeError: NoneType.start`, and a
  comment containing "def" silently computes a wrong column. This crashed an
  earlier all-at-once attempt on 2026-10-12, which was fully reverted.)
- Bytes: `subroutines/experiment_functions.py` is CRLF with no final newline
  (do not add one). Read `rb`, split/join on `b'\r\n'`, apply sites
  bottom-up, run on a copy first. `core.autocrlf=true` with no
  `.gitattributes` means `git show HEAD:` returns LF, so take the baseline
  from a worktree byte copy. `ExperimentVariables_master_satellite_global.py`
  is CRLF and `experiment_functions_two_nodes.py` is LF: no tree-wide
  `sed -i`. Interior blank lines stay empty.
- Add the compile-check cases before the edits (section 16.2).

---

## 4. Architecture

### 4.1 Two stacks side by side

Standalone (Explorer group `standalone`, 17 entry points since 2026-09-03):

```text
standalone/GeneralVariableScan.py      standalone/ExperimentVariables.py
standalone/AOMsCoils.py                standalone/MicrowaveScanOptimizer.py
standalone/GeneralVariableScan_HealthCheck.py (+ _CatchError)   ... (17 total)
utilities/BaseExperiment.py            utilities/DeviceAliases.py
MOT_experiments/SamplerMOTCoilAndBeamBalanceTune.py (restored to 70dd16c)
```

Master-satellite (repository top level, the primary namespace):

```text
GeneralVariableScan_master_satellite_mixin.py         (private)
GeneralVariableScan_master_satellite_single_node.py
GeneralVariableScan_master_satellite_two_nodes.py
AOMsCoils_master_satellite_Node1.py                   (self-contained)
AOMsCoils_master_satellite_Node2.py                   (self-contained)
ExperimentVariables_master_satellite_Node1.py
ExperimentVariables_master_satellite_Node2.py
ExperimentVariables_master_satellite_global.py
MicrowaveScanOptimizer_master_satellite.py
FORT_Polarization_Optimizer_master_satellite.py
AtomLoadingOptimizer_load_until_atom_master_satellite.py
SamplerMOTCoilAndBeamBalanceTune_master_satellite.py
MonitorMagnetometer_master_satellite.py
applets_master_satellite.py                           (CCB applet table)
utilities/BaseExperiment_master_satellite.py
utilities/DeviceAliases_master_satellite.py
utilities/config/master_satellite/device_aliases.json
subroutines/experiment_functions_two_nodes.py
device_db/device_db_master_satellite.py               (deployed database)
tests/MonitorSPCMinApplet.py                          (converted in place)
tests/measure_drtio_event_cost.py, measure_drtio_ttl_skew.py,
tests/measure_urukul_sync_windows.py, bringup_test_rtio_clock.py
```

Outside the repository: `artiq-master/compile_check_master_satellite.py`
(section 16.2) and `artiq-master/drtio_passive_check.py`.

### 4.2 Rollback principle

```text
Standalone chain:
  standalone gateware + device_db_standalone.py
    -> standalone/GeneralVariableScan.py
    -> utilities/BaseExperiment.py, DeviceAliases.py, experiment_functions.py

Master-satellite chain:
  master boot.bin + satellite boot.bin + device_db_master_satellite.py
    -> GeneralVariableScan_master_satellite_single_node (Node1 or Node2)
       or GeneralVariableScan_master_satellite_two_nodes (both)
```

No standalone file, `utilities/BaseExperiment.py` or `utilities/DeviceAliases.py`
imports master-satellite code (verified by grep, 2026-10-14). Rollback was
exercised only during diagnosis (Node1 reflashed with the standalone image on
2026-09-08 locked its RTIO PLL); the standalone GVS has not been re-run on
standalone gateware since the migration (section 18, O-20).
`standalone/README.md` says the standalone stack freezes once the migration
is validated.

### 4.3 Shared modules and their sanctioned extensions

These modules serve both stacks. Each change below is backward compatible
(default = standalone behaviour) and was user-authorized:

- `subroutines/aom_feedback.py`: `AOMPowerStabilizer(..., node_suffix='',
  config_node=None)` (`e768d50`, 2026-10-07; on `main` via PR #68); design
  comment for parallel grouping (`12c5314`).
- `subroutines/experiment_functions.py`: confirm-by-first-shot wrapper at 31
  sites, flag default False (2026-10-13); `t_OP_coil_relax_time` (`3257b1c`);
  slack restores in shared helpers (section 12.4).
- `standalone/ExperimentVariables.py` declares `advance_if_atom_in_first_shot`
  (False) and `max_loading_rounds` (20); `utilities/BaseExperiment.py` seeds
  `n_first_shot_confirm_tries` as `[0]` (`2fb2fbf`) so standalone runs the
  same functions.
- `utilities/config/alice|bob/` (`device_aliases.json`,
  `feedback_channels.json`) are read by the master-satellite resolver and
  stabilizers.
- `subroutines/k10cr1_functions.py`, `fitting/`: unchanged, reused.

### 4.4 File layout conventions

- `_mixin` files hold only private `_...Mixin` classes that do not inherit
  `EnvExperiment`, so Explorer never lists them.
- Every public experiment is one class in one file, named for its mode or
  node (`_single_node.py`/`_two_nodes.py`, `_Node1.py`/`_Node2.py`). The first
  line of the class docstring is the Explorer label and equals the filename.
- Never publicly import one `EnvExperiment` class into another repository
  file: ARTIQ rediscovers it as a duplicate (historically
  `GeneralVariableScan_master_satellite1`). Import shared code from mixins or
  under a private alias.
- ARTIQ workers put only the experiment file's own directory on `sys.path`
  (cwd = `artiq-master`). Every entry point outside the repository root needs
  the cwd-based block (`sys.path.append(cwd +
  "\\repository\\qn_artiq_routines")`, plus `standalone` where needed) before
  its first repo-local import. Seven files broke with "No module named
  utilities" when `standalone/` was created and gained the block (`331b669`).
  `utilities/BaseExperiment.py`, `MOT_experiments/Coils-SPCMCounts.py` and
  both HealthCheck files append `standalone/` because they import
  `ExperimentVariables`/`MicrowaveScanOptimizer` as top-level modules.
- Hardcoded scheduler expid `file` paths in the standalone optimizer and
  HealthCheck files carry the `qn_artiq_routines\\standalone\\` prefix.
- Shared infrastructure (`utilities/`, `subroutines/`, `fitting/`, configs)
  and already-foldered consumers (`tests/`, `MOT_experiments/`) stay in
  place. Master-satellite diagnostics live in `tests/`, import nothing
  repo-local, and moved there on 2026-09-15 (`bb5d687`; Explorer paths changed,
  rescan needed).

### 4.5 Explorer entries (2026-10-14)

```text
ExperimentVariables_master_satellite_Node1     (class ExperimentVariablesMasterSatelliteNode1)
ExperimentVariables_master_satellite_Node2
ExperimentVariables_master_satellite_global
GeneralVariableScan_master_satellite_single_node
GeneralVariableScan_master_satellite_two_nodes
AOMsCoils_master_satellite_Node1
AOMsCoils_master_satellite_Node2
MicrowaveScanOptimizer_master_satellite
FORT_Polarization_Optimizer_master_satellite
AtomLoadingOptimizer_load_until_atom_master_satellite
SamplerMOTCoilAndBeamBalanceTune_master_satellite
MonitorMagnetometer_master_satellite
tests/MonitorSPCMinApplet                       (fixed single_node Node1)
```

The namespace sanity function is a GVS callable, not an Explorer experiment.

### 4.6 The removed legacy two_nodes node mode

The old standalone node selection `which_node == "two_nodes"` coordinated two
independent Kaslis through a TTL handshake; it was never DRTIO. It was
removed with its `config/two_nodes` JSONs (`e86046f`, 2026-08-26) and must
never return. `utilities/BaseExperiment.py` accepts only `alice` and `bob`.
The master-satellite experiment modes `single_node`/`two_nodes` are a
different concept. Historical "two-node" functions in `experiment_functions.py`
are NOT master-satellite two-node functions.

### 4.7 Porting checklist (standalone experiment to master-satellite)

Applied to SamplerMOTCoilAndBeamBalanceTune, MonitorMagnetometer and
AtomLoadingOptimizer:

1. Replace `utilities.BaseExperiment` with `BaseExperimentMasterSatellite` and
   add `_DatasetRedirectMixin` to the class bases.
2. Add a node dropdown `EnumerationValue(('Node1', 'Node2'))`; never two
   booleans (they allow neither/both).
3. Construct Base unconfigured in `build()`.
4. In `prepare()`, in order: `base.configure_execution('single_node', node)`;
   `self.which_node = base.NODE_LEGACY_NAMES[node]` (alice/bob, needed because
   `aom_feedback` builds `utilities/config/<which_node>/`); `base.prepare()`;
   `base.prepare_laser_stabilizer()` (raises unless prepare ran first).
5. Delete `base.set_datasets_from_gui_args()` (absent in the master-satellite
   Base; ARTIQ archives arguments under `expid`).
6. Drop `core.reset()` at the top of `run()`: `base.initialize_hardware()`
   owns the reset and the satellite wait. Call `base.force_other_node_off()`
   once per run afterwards.
7. Check every GUI argument name against the per-node and global
   ExperimentVariables; a collision needs the capture-and-re-assert pattern
   (R13).
8. Experiment functions that read `self.stabilizer_*` or host scalars such as
   `AllSPCMs_RO2` need both `prepare_laser_stabilizer()` and
   `initialize_result_state()`, or they fail to COMPILE.
9. Explicitly suffix any best-parameter writes the redirect does not cover
   (use `resolve_experiment_variable_target`).
10. Add a compile-check case and unit tests.

---

## 5. Naming, namespaces and device resolution

### 5.1 Node naming and suffixes

New logical names, node-specific datasets, host nicknames and new
infrastructure use the node as a SUFFIX: `dds_FORT_Node1`,
`dds_cooling_DP_Node2`, `ttl_repump_switch_Node1`, `sampler0_Node2`,
`zotino0_Node1`, `f_FORT_Node2`, `p_FORT_loading_Node1`,
`FORT_MM_monitor_Node2`, `780_QWP_Node1`. Never the prefix style
(`dds_Node1_FORT`). Do not rename validated standalone code to enforce this.
`BaseExperimentMasterSatellite.NODE_LEGACY_NAMES` maps Node1 to `alice` and
Node2 to `bob`.

### 5.2 Physical versus logical names

Low-level physical names (`urukul4_ch3`, `ttl16`) stay exactly as generated in
the unified database and are never renamed with node suffixes. Node naming
is applied in the logical mapping layer. Keep the generated database close to
the generated hardware description.

### 5.3 Logical configuration (per node, shared with standalone)

`utilities/config/alice/device_aliases.json` and
`utilities/config/bob/device_aliases.json` remain the source of truth for
logical DDS aliases: `ALIAS_MAP` (logical alias to standalone physical DDS,
e.g. Bob `dds_FORT -> urukul1_ch3`, Alice `dds_FORT -> urukul0_ch0`) and
`DDS_DEFAULTS` (default frequency/power variable names). Do not duplicate
these in the master-satellite mapping. `feedback_channels.json` in the same
folders holds per-node feedback calibration.

### 5.4 Physical translation

`utilities/config/master_satellite/device_aliases.json` is the single source
translating each node's standalone physical name to the unified physical name
(cc399dc, 2026-08-27). Top-level keys `Node1` and `Node2`, 43 explicit entries
each: 12 DDS channels, 3 CPLDs, 3 Samplers, 1 Zotino, 16 TTLs, 8 edge
counters. Node1 maps to itself; Node2 maps `urukul0/1/2 -> urukul3/4/5`,
`sampler0/1/2 -> sampler3/4/5`, `zotino0 -> zotino1`, `ttlN -> ttl(N+16)`,
`ttl0-3_counter -> ttl16-19_counter`, `ttl8-11_counter -> ttl24-27_counter`.
Each mapping is an explicit entry; the "+16" pattern is a description, not a
rule the code may compute (R10). Excluded: the canonical core (not per node),
logical aliases, host/USB devices, `led0/led1`, any arithmetic. There is no
second `device_map.json`.

### 5.5 Single physical resolution invariant

For a given node and logical device, the physical device is identical in
single-node and two-node mode; only the experiment-facing name changes.
Resolve physical identity first, then presentation. There is exactly one
physical translation source and no per-mode table.

```text
Bob dds_FORT -> standalone urukul1_ch3 -> unified urukul4_ch3
  single Node2: self.dds_FORT        -> urukul4_ch3
  two nodes:    self.dds_FORT_Node2  -> urukul4_ch3
Alice dds_FORT -> standalone urukul0_ch0 -> unified urukul0_ch0
  single Node1: self.dds_FORT        -> urukul0_ch0
  two nodes:    self.dds_FORT_Node1  -> urukul0_ch0
```

(`PLAN.md` gave `urukul0_ch1` for Node1 FORT; that is Node1's cooling DP.)

### 5.6 Presentation by mode

- **single_node**: the selected node is exposed under the old unsuffixed
  names (`self.dds_FORT`, `self.sampler0`, `self.zotino0`, TTL aliases,
  variables such as `self.f_FORT`) so existing single-node functions run
  without rewrites. Even a Node2-only run uses the one master core and reaches
  Node2 through DRTIO. The projection is in memory only, never a second
  dataset: dataset `f_FORT_Node2 -> self.f_FORT_Node2 -> self.f_FORT`.
- **two_nodes**: both nodes are exposed with suffixes only. There are no bare
  device names (section 6.6). Principle: physical resolution is shared;
  logical presentation depends on execution mode.

### 5.7 Resolver: `DeviceAliasesMasterSatellite`

`utilities/DeviceAliases_master_satellite.py` (landed `1f83f11`, 2026-08-28).
One instance per node. Owns: loading and validating the exact mapping; exact
dictionary lookup (`resolve_physical_name`); binding under a caller-supplied
presentation name (`bind_physical_device`); two-stage DDS resolution
(`resolve_dds_physical_name`: logical -> standalone -> unified);
node-specific DDS defaults (`resolve_dds_default_attributes`); `bind_dds`
with independent per-node `dds_list`/frequency/power bookkeeping; clear
`DeviceResolutionError`s for missing aliases, mappings, devices or defaults.
It does NOT reset the core, wait for DRTIO, initialize hardware, choose
modes, implement SPCM policy, or handle USB devices. Standalone
`utilities/DeviceAliases.py` stays the standalone implementation (note: it
silently skips aliases missing from its list; see section 18, O-22).

DDS default variables: the alice/bob configs name unsuffixed variables
(`f_FORT`, `p_FORT_loading`, `f_cooling_DP_MOT`); the resolver and
`Base._default_attribute_name(node, standalone_name)` turn them into the
per-node authoritative names (`f_FORT_Node2`). Bob's `dds_D1_pumping_DP`
defaults resolve to `f_GRIN2_D1_pumping_Node2`/`p_GRIN2_D1_pumping_Node2`.

DDS alias resolution and default-variable resolution are ONE step and must
not be split (retired PLAN.md section 16): the same node context that maps
`dds_FORT -> urukul4_ch3` must map `f_FORT -> f_FORT_Node2`. A resolver that
mapped the device but looked up the default unsuffixed would let Node2
silently use Node1's value, or a stale bare standalone dataset from the shared
database (section 7.3).

### 5.8 SPCM binding exception

Canonical detectors (section 2.6) are bound through the Node1 resolver
(`shared_spcm_resolver`) in every mode, never through the selected node's
generic TTL translation (never `ttl0 -> ttl16`). Single-node compatibility
installs fixed legacy aliases on both the TTL and counter attributes,
identically for Node1 and Node2 and never reordered:

```text
ttl_SPCM0           -> SPCM_H1     ttl_SPCM0_counter           -> SPCM_H1 counter
ttl_SPCM1           -> SPCM_V1     ttl_SPCM1_counter           -> SPCM_V1 counter
ttl_SPCM0_OtherNode -> SPCM_H2     ttl_SPCM0_OtherNode_counter -> SPCM_H2 counter
ttl_SPCM1_OtherNode -> SPCM_V2     ttl_SPCM1_OtherNode_counter -> SPCM_V2 counter
```

Implemented by `Base.SPCM_DEVICES`, `LEGACY_SPCM_ALIASES`,
`_bind_shared_spcms` and `_publish_shared_spcm_compatibility` (which restores
Node2's raw `ttl0/1/8/9` presentation afterwards). The same bare SPCM aliases
are also published in two_nodes mode. Two-node code gates all four canonical
detectors (`SPCM_H1/V1/H2/V2`). Deprecated functions that use raw TTLs are
not rewritten; the distinction is documented in `_bind_shared_spcms`.

---

## 6. BaseExperimentMasterSatellite

`utilities/BaseExperiment_master_satellite.py` (about 2300 lines). Class
constants: `VALID_MODES` (`single_node`, `two_nodes`), `NODE_LEGACY_NAMES`,
`SATELLITE_DESTINATION = 1`, `SPCM_DEVICES`, `LEGACY_SPCM_ALIASES`,
`NODE_TTL_ALIASES`, `EDGE_COUNTERS`, `SINGLE_NODE_WIRING_METADATA`,
`TWO_NODE_AGREEING_SCALARS`, `MAGNETOMETER_RESULT_DATASETS`,
`MICROWAVE_CALIBRATION_DATASETS`, `POLARIZATION_RESULT_DATASETS`.
The standalone `utilities/BaseExperiment.py` is unchanged in role.

### 6.1 Responsibility boundary

Base owns: the canonical core lifecycle (one core, `core_dma`, scheduler);
satellite readiness; execution mode and presentation; authoritative variable
loading and compatibility projection; SPCM policy; wiring metadata; CPLD,
DDS, TTL, Sampler and Zotino initialization and safe states; idle-node
force-off; laser-stabilizer construction; derived RF amplitudes; result
state; scan/override target resolution; variable reload and dependent-state
refresh; the K10CR1 proxy; the dataset redirect. Applet helpers live in
`applets_master_satellite.py`, not in Base.

### 6.2 Construction and build

`BaseExperimentMasterSatellite(experiment, experiment_mode=None,
which_node=None, satellite_wait_attempts=100,
satellite_ready_poll_interval=100*ms)`. Every public experiment constructs it
UNCONFIGURED in `build()`. Build then binds core/core_dma/scheduler,
`print_async` (async RPC) and the `write_results` wrapper, creates a resolver
per node, binds the suffixed Node1+Node2 physical superset and the shared
SPCMs. It does not load datasets, choose a presentation or prepare DDS
defaults. This keeps repository examination safe (section 6.11).

### 6.3 `configure_execution(mode, which_node)`

Called once, in the experiment's `prepare()`, before `base.prepare()`.
Validates mode/node strictly; loads authoritative variables for the active
nodes plus globals (`_load_experiment_variables`, which only reads and raises
via `_raise_missing_datasets` naming the owning initializer file); builds the
single-node `compatibility_variable_map`. In single_node it deactivates the
idle node's init groups (section 6.5), projects the selected node's devices,
TTL aliases and legacy SPCM aliases to bare names, and publishes
`named_devices` and the `k10cr1_ndsp` proxy; in two_nodes it publishes
`named_devices_Node1/_Node2`. Wiring metadata is installed either way.

### 6.4 `prepare()` (one shot)

Binds the 12 DDS aliases per active node (bare or `_NodeX`), builds
`all_dds_channels`, computes derived `ampl_*` per node
(`_compute_derived_amplitudes`) and the initial caches. In two_nodes it also
runs `_publish_broadcast_aliases` (section 6.6). Calling `prepare()` twice
raises. Never call it per scan point.

### 6.5 Single-node mode and the idle node

- `_deactivate_node_hardware_groups` removes the idle node's CPLD, Sampler,
  Zotino, TTL and safety-state lists from `initialize_hardware` but keeps them
  as `_other_node_*`; `_bind_other_node_dds_switches` binds the idle node's 12
  Urukul `sw` TTLs (through `_default_attribute_name`, not
  `_presentation_name`).
- `force_other_node_off()` (kernel, `21f6392`, 2026-09-17): initializes the
  idle node's three CPLDs, switches its 12 RF switches off, blocks its
  active-low optical gates and writes Zotino channels 0-15 (of 32) to 0 V
  (`for channel in range(16)`; every wired coil and trigger channel is below
  16, so this is harmless today), waiting
  for destination 1 first when the idle node is Node2. The CPLD init is
  essential: a Urukul RF switch is the CPLD `cfg` bit OR-ed with the `sw` TTL,
  the register is not carried between experiments and the driver mirror
  reseeds to 0, so `sw.off()` alone can leave a beam on. `CPLD.init()` does
  not run the AD9910 SYNC search. Call it ONCE per run after
  `initialize_hardware()` (core reset can only drive TTLs low; it cannot
  change a held Zotino output).
- Callers: single-node GVS, `MicrowaveScanOptimizer_master_satellite`,
  `FORT_Polarization_Optimizer_master_satellite`,
  `AtomLoadingOptimizer_load_until_atom_master_satellite`,
  `SamplerMOTCoilAndBeamBalanceTune_master_satellite` (`0ff183e`),
  `MonitorMagnetometer_master_satellite`. AOMsCoils never calls it, by
  design. Two-node runs do not call it.
- Consequence: single-node runs DO touch the idle node once per run. Before
  2026-09-17, single_node mode only rebound seven lists to `[]` and the idle
  crate stayed hot.

### 6.6 Two-node mode presentation

- No bare device aliases: every device is `_Node1`/`_Node2`, and a bare name
  fails to compile. The earlier broadcast-wrapper design (one bare object
  driving both crates) was removed in `e768d50` (2026-10-07) after
  `laser_stabilizer_Node2` drove Node1's repump, cooling DP and fiber AOMs
  through it. Despite its name, `_publish_broadcast_aliases` now publishes
  only the legacy SPCM aliases, the agreeing scalars and `which_node =
  'alice'` (the last only so `end_measurement`'s dead magnetometer branch
  type-checks).
- `TWO_NODE_AGREEING_SCALARS` are published bare only while Node1 and Node2
  agree; `_publish_agreeing_scalars` raises `RuntimeError` on divergence. It
  runs at prepare and in every `refresh_variable_dependent_state` (after
  overrides and each scan point), so asymmetric overrides raise and scans of
  an agreeing scalar take effect. The 13 scalars:
  `monitor_magnetometer_in_end_measurement`, `monitors_for_atom_loading`,
  `use_chopped_readout`, `no_first_shot`, `require_D1_lock_to_advance`,
  `enable_laser_feedback`, `Luca_trigger_for_feedback_verification`,
  `t_SPCM_first_shot`, `t_SPCM_second_shot`, `t_atom_check_time`,
  `t_delay_between_shots`, `t_FORT_loading`, `t_MOT_dissipation`. The timing
  ones are joint because one hardware event consumes each (one gate of all
  four SPCMs; one interval on the shared timeline).
- Deliberately per node and addressed explicitly:
  `PGC_and_RO_with_on_chip_beams`, `do_PGC_after_loading`, `t_FORT_drop`
  (removed from the agreeing list in `e768d50`), `t_PGC_after_loading`,
  `single_atom_threshold`, every `f_*`, `p_*` and coil voltage.
- Per-node `ampl_*` and wiring metadata are suffixed.

### 6.7 Repeatable APIs

- `resolve_experiment_variable_target(name, node=None)`: globals resolve to
  themselves; with an owning `node` (from a per-node override dictionary), a
  bare or suffixed name resolves to that node in either mode (`6a89075`,
  2026-10-06); in single_node, bare legacy names map to the selected suffix,
  selected suffixed names are accepted and other-node names rejected; in
  two_nodes without a node, suffixes are required ("Ambiguous ... use
  X_Node1 or X_Node2"); unknown names are rejected clearly.
- `refresh_compatibility_variables()`: single_node only; copies authoritative
  suffixed attributes to legacy names; no dataset I/O.
- `refresh_variable_dependent_state()`: refreshes DDS frequency/power caches
  and per-node `ampl_*`; in two_nodes re-publishes the agreeing scalars and
  WRITES the derived broadcast dataset `t_alternating_RO_exposure =
  n_alternating_RO_windows_per_node * t_alternating_RO_window`. No prepare,
  rebinding, hardware init or core reset. (The retired spec said "no dataset
  I/O"; that is no longer true in two_nodes.)
- `reload_experiment_variables()`: rereads the active nodes' and global
  datasets, refreshes projections and caches, fails if required data
  disappeared; never creates calibration datasets or rebuilds hardware.
- `result_name_tag()`: `Node1`/`Node2` in single_node, `TwoNodes` otherwise.
  It reads `which_node` because experiments spell their selector differently
  (MicrowaveScanOptimizer `Node1`/`Node2`, FORT optimizer `node1`/`node2`).
- `current_iteration_dataset_names()`, `resolve_result_dataset_name(name)`
  (section 6.9).
- Result initialization/reset is separate from variable-dependent refresh.

### 6.8 Hardware initialization

`initialize_hardware(turn_off_dds_channels, turn_off_zotinos,
initialize_spcms, reset_core=True)` (kernel). Order:

1. one Base-owned `core.reset()` (skipped with `reset_core=False`);
2. `_wait_for_satellite()` if Node2 is active: calls
   `core.get_rtio_destination_status(1)` up to 100 times, separated by
   `delay(satellite_ready_poll_interval)` (100 ms), then raises
   `RuntimeError("DRTIO destination 1 did not become available.")`; then
   `break_realtime`. Caution: `delay()` only advances the timeline cursor and
   does not block the CPU, and the status call returns at once, so the 100
   polls take microseconds of wall time, not about 10 s. If destination 1 is
   down the raise is effectively immediate; if it comes up on poll k > 1 the
   cursor is left about k x 100 ms ahead (`break_realtime` only moves it
   forward). Untested on hardware because destination 1 has always been up
   (O-33). `tests/measure_drtio_ttl_skew.py` copies the same loop shape;
3. canonical SPCM input configuration;
4. CPLD init, Node1 then Node2;
5. DDS init, Node1 then Node2 (`_initialize_dds_group`: `init()`, then
   `sw.off()` when outputs are to be off, then `set_att(0.0)` and `set()`, so
   a DDS is OFF before attenuation/profile programming);
6. TTL directions and safe states;
7. Samplers: PGIA gain code 0 (x1) on all 8 channels of both nodes. The
   standalone `set_gain_mu(ch, 8)` idiom also produced x1 (it effectively sent
   0x0000), so every feedback setpoint is calibrated at x1. Do NOT "fix" the 8
   to 3 (x1000);
8. Zotino init and safe state;
9. final slack restoration (`break_realtime`).

GVS calls it at every scan point, so "one reset" means one per call.
`satellite_ready_poll_interval` is startup polling, not DRTIO communication
latency. A reachable master core does not mean the satellite is ready.

Hazard: initialization is not passive. It resets the core, changes SPCM/TTL
directions and levels, initializes CPLDs and DDSs, sets DDS
attenuation/frequency/amplitude/switch state, sets Sampler gains and zeroes
Zotino channels; `force_other_node_off` adds writes to the idle node. Never
use GVS as a first post-flash or post-repair check; confirm every output
state is electrically safe first.

### 6.9 Result state and dataset naming

- `initialize_result_state()` (renamed from
  `initialize_single_node_result_state` in `e768d50`; runs in both modes)
  seeds every host scalar and broadcast series the reused functions need.
  It mirrors the standalone seeding surface (about 117 datasets and 36 host
  scalars); `6329a1a` (2026-09-10) added the 52 datasets that were missing,
  with exact standalone shapes because `[0]`, `[0.0]`, `[[0.0, 0.0]]` and
  `[]` differ. It also seeds `n_first_shot_confirm_tries` as `[0]`, seeds
  `FORT_MM_monitor_NodeX`/`FORT_APD_monitor_NodeX` per node, omits
  `Rigol_modulation_volts` (Rigol out of scope), and in two_nodes broadcasts
  `t_alternating_RO_exposure`.
- `n_measurements` exception: `initialize_result_state` writes the run's
  submitted `n_measurements` with `broadcast=True, persist=True` (`82f5f6f`,
  2026-10-05), its only persisted write. Reasons: applets read only broadcast
  datasets (the retention and "Atom loading time" applets went blank), and
  any standalone run on this master de-persists the shared global (section
  7.8), so every master-satellite run self-heals it.
- Seed datasets unconditionally in each worker, in `prepare()`/run-state
  init, not in kernels: `get_dataset` falls through to the master database,
  but `append_to_dataset` resolves only through this worker's local and
  broadcast dicts, so a "create only if absent" guard is satisfied by a
  leftover master dataset and the append still raises "Cannot mutate
  nonexistent dataset". If that error appears, re-diff the whole seeding
  surface rather than adding one name.
- `reset_result_state_for_scan_point()` (via `initialize_result_datasets()`)
  runs per point in both modes: `measurements_progress = 0.0` (broadcast,
  non-persist) and every `*_current_iteration` dataset cleared from
  `current_iteration_dataset_names()` (before `c1ddbd2` these histograms
  showed accumulated totals); magnetometer results only in single_node. It
  never persists calibration, reloads variables, resets the core or rebinds
  devices. Preserve existing `break_realtime()` calls and the final MOT-coil
  behaviour.
- `resolve_result_dataset_name(name)`: in single_node, suffixes the
  node-specific result and calibration names (magnetometer, polarization,
  microwave calibration, coil calibration, feedback map); identity in
  two_nodes, so both nodes write the same joint result names.
- Run-global names policy. Retired spec, verbatim (the `_install_wiring_metadata`
  docstring cites this): "Legacy hard-coded unsuffixed magnetometer append
  names are redirected only in the new stack. Do not create duplicate
  unsuffixed result datasets. `measurements_progress` stays unsuffixed as
  run-global state." Extension (current code): `measurements_progress` and
  the `*_rate_dataset` / `scan_*` name strings (`SPCM*_counts_per_s`,
  `AllSPCMs_counts_per_s`, `scan_variables`, `scan_sequence1/2`) are
  unsuffixed in BOTH modes, because they name run state, not one node's
  hardware.
- Magnetometer results use node-suffixed transient names in standalone
  spelling (`Magnetometer_MOT_X_Node1`; `MAGNETOMETER_RESULT_DATASETS`, 9
  names), single_node only. They are CREATED only by
  `reset_result_state_for_scan_point`, and `append_to_dataset` cannot create a
  dataset, so every experiment whose sequence can reach `measure_Magnetometer`
  / `measure_Magnetometer_Node2_X_Y_Z` must call it before its first
  experiment call: GVS does, and MicrowaveScanOptimizer does since 2026-10-14
  (`bfca64d`, O-35). The same method also broadcasts
  `magnetometer_axis_labels = ["X", "Y", "Z"]` (single_node, non-persistent),
  the legend of the magnetometer applets (section 14.2).
- Microwave calibrations and `health_check_uw_*` fidelities persist
  node-suffixed (`MICROWAVE_CALIBRATION_DATASETS`, 15 names:
  `f_microwaves_00/01/11/m10/m11_dds`, `t_microwave_00/01/11/m10_pulse`,
  `t_MW_RF_pulse`, `health_check_uw_freq00/01/11/m10/m11`). `fit_parameter_*`
  and `parent_rid` stay unsuffixed run state.
- In two_nodes the result datasets (RO1/RO2, retention, loading) are
  unsuffixed because loading is one joint measurement; FORT_MM/APD monitors
  are suffixed.
- `write_results` saves a named copy per scan point:
  `<RID>-<experiment>_<Node1|Node2|TwoNodes>_<function minus _experiment>_scan_over_<v1>[_and_<v2>].h5`
  (`7ea1476`, 2026-09-17). This guards against ARTIQ h5 corruption during
  worker cleanup, which does happen: the ARTIQ h5 files for RIDs 38816 and
  38826 are unreadable, and 38816 survives only through its named copy.

### 6.10 `_DatasetRedirectMixin` and the K10CR1 proxy

- `_DatasetRedirectMixin` (mixed into every public experiment that uses Base) overrides
  `set_dataset`, `get_dataset` and `append_to_dataset` to pass names through
  `base.resolve_result_dataset_name` and the feedback name map.
  `append_to_dataset` is an async RPC (`@rpc(flags={'async'})`), so the
  master-satellite stack does LESS blocking host work than standalone, whose
  appends were blocking.
- `_NodeAxisK10CR1Proxy`: section 13.

### 6.11 Repository examination lifecycle

ARTIQ repository examination passes `None` for every argument and must work
with an empty master-satellite dataset namespace, because Explorer must
discover the ExperimentVariables initializers that create it. Therefore every
public `build()`: registers argument metadata, binds the suffixed superset
through an unconfigured Base, and does NOT choose a mode, load datasets,
publish compatibility names or prepare DDS defaults. GVS metadata uses code
defaults (`n_measurements = 100`; scan default `f_FORT` single /
`f_FORT_Node1` two). At prepare: strict mode/node validation, load only the
required node(s) plus globals, publish compatibility and DDS bindings, and
fail clearly on missing datasets before any hardware access. Base never
creates calibration datasets. Commits `a3eb7fe`, `f06f9be` (2026-08-30).
Test stubs install an artiq stub only when artiq is genuinely absent (marker
`_qn_hardware_free_stub`), so the shared repository-scan worker is never
poisoned (`c18a6f2`, `3d0ec82`).

---

## 7. ExperimentVariables

### 7.1 Files and ownership

| File | Class | Owns |
|---|---|---|
| `ExperimentVariables_master_satellite_Node1.py` | `ExperimentVariablesMasterSatelliteNode1` | `NODE1_VARIABLES` (213, all `*_Node1`) |
| `ExperimentVariables_master_satellite_Node2.py` | `ExperimentVariablesMasterSatelliteNode2` | `NODE2_VARIABLES` (213, same base names, `*_Node2`) |
| `ExperimentVariables_master_satellite_global.py` | `ExperimentVariablesMasterSatelliteGlobal` | `MASTER_SATELLITE_VARIABLES` (18 globals) |

Renamed from `ExperimentVariables_Node1/_Node2/_master_satellite.py` on
2026-08-31 (`83ae355`); "global" because no node owns it and every mode loads
it. Persistent dataset names did not change. Standalone keeps
`standalone/ExperimentVariables.py`.

### 7.2 Initializer behaviour

`build()` offers every variable as a GUI argument whose default is the
current dataset (code default only when the dataset is absent). `run()` writes
all of them with `broadcast=True, persist=True`, so submitting unchanged
values preserves calibration; per-node `run()` also writes derived
`detuning_*_units_Gamma_NodeX`, `microwave_frequency_GHz_NodeX`,
`f_dds_clock_resonance_NodeX`. The global file never initializes node files.
Running an initializer is non-destructive (absent values are added, existing
ones kept unless edited in the GUI).

Fresh database: run Node1 once, Node2 once, then global once. Normal
experiments only read and fail with "Missing required master-satellite
persistent datasets ... Run the indicated ExperimentVariables
initializer(s): <file>: <names>".

### 7.3 Ownership rules

- Existing persistent calibration always wins; code defaults seed only
  missing datasets.
- One `dataset_db.pyon` serves BOTH stacks, so every persistent name must be
  globally unique (retired PLAN.md 18.1/18.2): node quantities carry
  `_Node1`/`_Node2`, so Node1, Node2, the globals and the standalone stack
  never overwrite each other. The global file never duplicates Node1/Node2
  calibration values. Before declaring any new global (R12 permitting), grep
  `standalone/ExperimentVariables.py` and `dataset_db.pyon`: a shared name
  with a different default makes the live value depend on which initializer
  ran last (today `n_measurements`, `two_atom_threshold`,
  `two_atom_threshold_for_loading`; history: `t_delay_in_bob_mu`, section
  7.5). Bare standalone datasets live in the same database with different
  values (2026-10-14: `single_atom_threshold` 60000 against `_Node1` 14000 /
  `_Node2` 17000; `aom_feedback_iterations` 200; `t_OP_coil_relax_time`
  0.0004 against 0.0015 per node), which is why a bare name can resolve to a
  stale standalone value (section 14.2).
- Create both suffixes for anything conceptually node-calibratable, and
  preserve asymmetric or historical node variables (`Rigol_*` remain declared
  per node, unused).
- Genuinely JOINT two-node quantities are GLOBAL and unsuffixed (one
  comparison over four shared SPCMs, one gate on one timeline), e.g.
  `two_atom_threshold` (global since `5a0bc76`, 2026-09-17; previously the
  bare name did not even exist in two_nodes), `advance_if_atom_in_first_shot`,
  `max_loading_rounds`, the alternating geometry, `dummy_variable` (global
  since `e768d50`; per-node made a scan of `dummy_variable` ambiguous).
- Execution selections (`selected_node`, mode) are arguments, not datasets.
- Never migrate runtime/result state into ExperimentVariables
  (`fit_parameter_*`, `feedbackchannels`, `parent_rid`, `child_rid`).
  `feedbackchannels` is redirected at runtime to `feedbackchannels_NodeX`.
- Feedback variables: keep historical D1 variables but do not invent
  `set_point_D1_DP` (Bob's `feedback_channels.json` still names it; D1 is not
  in the fast list). Master-satellite has `set_point_D1_SP_Node1` (0.2) and
  `_Node2` (0.175). Active feedback is `fast_feedback_dds_list_NodeX` = A1-A6
  plus FORT, with independent loading/science/holding FORT setpoints per node.

### 7.4 Global variables (2026-10-14)

| Name | Declared | Live | Role / status |
|---|---|---|---|
| `n_measurements` | 100 | 100 | Measurements per point; persisted by every run (section 6.9). Standalone declares 50. |
| `t_Node2_rtio_offset_mu` | 0 | 0 | ELECTRICAL cross-node offset (section 7.5); unmeasured. |
| `t_Node2_excitation_delay_mu` | 189 | 189 | OPTICAL photon-coincidence delay (section 7.5); unmeasured. Read by no code yet: a deliberate placeholder for the deferred photon-coincidence work (the unported single-photon code still reads the standalone `t_delay_in_bob_mu`; the skew probe only names it in its docstring). |
| `parallel_AOM_feedback` | True | True | Unread today (only tests name it), but it is the selector the retired PLAN.md designated for the future switch: False = sequential Node1/Node2 feedback, True = concurrent. Keep it declared (section 10.5). |
| `dummy_variable` | 0.0 | | Scan placeholder. |
| `two_atom_threshold` | 39000 | 42000 | Joint two-atom readout threshold (counts/s). Standalone declares 20000. |
| `two_atom_threshold_for_loading` | 89000 | 100000 | Joint loading threshold. Standalone declares 84000. |
| `t_atom_check_time_two_node` | 0.01 s | | UNREAD (joint timing moved to agreeing scalars). |
| `t_SPCM_first_shot_two_node` | 0.01 s | | UNREAD. |
| `t_SPCM_second_shot_two_node` | 0.01 s | | UNREAD. |
| `t_delay_between_shots_two_node` | 0.0 s | | UNREAD. |
| `t_MOT_dissipation_two_node` | 0.05 s | | UNREAD. |
| `max_atom_check_tries_two_node` | 100 | | Two-node loader bound. |
| `max_loading_rounds` | 20 | 20 | Bounds the two-node loader's outer retry and every confirm-by-first-shot loop. Renamed from `max_loading_rounds_two_node` (`2fb2fbf`, 2026-10-13); hand-seeded before the rename landed. |
| `advance_if_atom_in_first_shot` | False | **True** | Confirm-by-first-shot flag (section 9.7). Live True, so the confirm loop is the default behaviour. |
| `n_alternating_RO_windows_per_node` | 10 | | Alternating readout geometry. |
| `t_alternating_RO_window` | 0.001 s | | Alternating readout geometry. |
| `t_alternating_RO_pad` | 0.0001 s | | Alternating readout geometry. |

`t_alternating_RO_exposure` is derived and broadcast by Base, never declared
(the model for R12).

### 7.5 Cross-node timing: two layers

Until 2026-10-06 the single variable `t_delay_in_bob_mu` conflated two
things, and it also collided with `standalone/ExperimentVariables.py` (same
name, default 0), so the live value was whichever initializer ran last (it
was 0). The rename (`e768d50`) split it:

- `t_Node2_rtio_offset_mu` (ELECTRICAL, layer 1): the offset applied to
  Node2-crate events so both crates' output pins move together for one
  scheduled timestamp. Sign: positive when Node2's pin moves LATER; the
  sequence schedules Node2 events at `at_mu(t - offset)`. Measured by
  `tests/measure_drtio_ttl_skew.py` (never run; section 12.2). Copy the
  probe's T2 median in WITH ITS SIGN UNCHANGED: negating it doubles the skew,
  and the only symptom is readout light leaking between alternating windows.
  Consumed by the alternating readout. Replaces the concept behind
  `alice_extra_offset_mu`, a scope-calibrated constant of the retired TTL
  handshake that has no meaning under DRTIO.
- `t_Node2_excitation_delay_mu` (OPTICAL, layer 2): the fine delay that brings
  the two emitted photons into coincidence at the SPCM (fiber and free-space
  path lengths, AOM turn-on). Not measurable by TTL loopback; needs a photon
  coincidence (HOM dip / cross-correlation) measurement, part of the deferred
  herald work.

`standalone/ExperimentVariables.py` still declares `t_delay_in_bob_mu`
(live 0); it is standalone state, not an orphan.

### 7.6 Notable per-node variables

- OP coil relaxation: every optical-pumping sequence delays on
  `t_OP_coil_relax_time` (added `3257b1c`, 2026-09-03; read by
  `CW_optical_pumping_node1/2` and `record_CW_optical_pumping_node1/2`); the
  two-node recorder uses the larger node value. Code seed 0.4 ms (the
  previously hardcoded value; Node2 had 0.395 ms). Live
  `t_OP_coil_relax_time_Node1/_Node2` = 1.5 ms each (checked 2026-10-14), so
  the 1.5 ms tuned on the dead `t_coil_relaxation_time_OP_NodeX` is now in the
  live variable. Standalone live `t_OP_coil_relax_time` = 0.4 ms.
- `t_coil_relaxation_time_OP_Node1/_Node2` were seeded but never read and
  were removed on 2026-10-12 (`0e54c84`). Do not re-seed them. They are no
  longer in `dataset_db.pyon`.
- Alternating thresholds: `alternating_shot_threshold_NodeX` 20000/21000
  (live same); `alternating_shot_threshold_for_loading_NodeX` declared
  62000/60000 (`b6b4de1`) but live 50000/45000, because persisted values win
  and the master was never updated. Only the disabled alternating loader
  reads the loading thresholds. Rationale (`b6b4de1`, 2026-10-09): each
  loading threshold is that node's no-atom floor (50100 / 45700 c/s) plus
  half an atom; the first values, 51000/45000 (`19c289f`), passed 23% / 72%
  of tries on pure background. Halving the joint threshold lands on the
  per-node floor, because each window collects the full common background
  (section 9.5).
- `single_atom_threshold_Node1/_Node2`: live 14000 / 17000 counts/s.
- `aom_feedback_iterations_NodeX`: declared 200 (Node1) / 100 (Node2), live
  100 / 100. It is a cap broken by `in_tol`. `aom_feedback_averages_NodeX`
  = 20. Standalone DECLARES 15 iterations / 10 averages, but the live
  standalone datasets in the shared database are 200 / 20, and persisted
  values win, so a standalone run on this master uses 200 / 20.
- Node2 GRIN amplitudes: `p_GRIN1_D1_pumping_Node2` -8 and
  `p_GRIN2_excitation_Node2` -6 are live; the others are -150 (killed).
  `f_GRIN2_excitation` is 150.60 MHz on Node1 vs 150.50 MHz on Node2.

### 7.7 Unread declarations and orphans

- Seven globals are read by no code: the five `*_two_node` timings above,
  `parallel_AOM_feedback` and `t_Node2_excitation_delay_mu`. The five timings
  could be removed (removing a declaration is safe; adding one is not). Keep
  the other two: `t_Node2_excitation_delay_mu` is a deliberate placeholder,
  and `parallel_AOM_feedback` is the designated selector for the future
  sequential/concurrent feedback switch, so deleting it now and re-adding it
  later would cost an R12 failure window.
- Orphan persistent dataset still in `dataset_db.pyon` (2026-10-14):
  `max_loading_rounds_two_node` = 20. Remove with `artiq_client del-dataset
  max_loading_rounds_two_node` when the user agrees. The per-node
  `two_atom_threshold*_NodeX` and `t_coil_relaxation_time_OP_NodeX` orphans
  recorded in commit bodies are no longer present.
- `tune_852_waveplates_to_target_in_experiment_NodeX` is NOT an orphan.
  It has been declared since `d44d528` (2026-08-30), except for twenty
  minutes on 2026-10-05: `a465351` (15:52) removed it and `a2845e7` (16:12)
  re-added it. It is declared today (line 156 of both node initializers), live
  False on both nodes, and read through the single-node projection by the
  `atom_loading_2` branch (`if self.tune_852_waveplates_to_target_in_experiment:`
  in `experiment_functions.py`). The compile checker's comment calling it
  "NOT a declared master-satellite variable" is stale (O-25).

### 7.8 Persistence pitfalls

- **New Variable breaks everything** until its initializer runs (R12). It
  happened with `19c289f` (run Node1/Node2 initializers), `c1ddbd2` (run the
  global initializer) and `3a90e2d` (broke AOMsCoils; reverted in `ae93327`).
- **Hand-seeding** avoids the failure window: seed on the master BEFORE the
  code lands, `artiq_client set-dataset --persist <name> <pyon value>`
  (used 2026-10-13: `max_loading_rounds 20`). Confirm with a compile sweep
  whose ACTION NEEDED block is empty (section 16.2). Standalone needs no
  seeding (`utilities/BaseExperiment.py` sets an attribute for every dataset
  in the shared database). Per-run series such as
  `n_first_shot_confirm_tries` are seeded in code.
- **Standalone runs de-persist `n_measurements`**: `utilities/BaseExperiment.py`
  broadcasts it without `persist`, which clears the global's persist flag; the
  next master flush drops it and every master-satellite experiment then fails
  in prepare at the next master restart (happened 2026-10-05 after
  `tests/A6-MonitorMagnetometer.py`). It also blanks the "retention and
  loading AllSPCMs" and "Atom loading time (s)" applets
  (`applets/plot_retention_and_loading.py` swallows the KeyError). Fixed by
  the persisted write in section 6.9; manual recovery: one master-satellite
  run, or resubmit the global initializer.
- **Local write after broadcast deletes** (R17). Hit by the GVS override
  record (ordered first, `ef55cae`) and by `MicrowaveScanOptimizer`
  (ordered first, `619d3ec`; section 11.1).
- **Name collisions with standalone**: several globals share names with
  `standalone/ExperimentVariables.py` but have different defaults
  (`n_measurements` 100 vs 50; `two_atom_threshold` 39000 vs 20000;
  `two_atom_threshold_for_loading` 89000 vs 84000). Whichever initializer runs
  last sets the live value, exactly as with `t_delay_in_bob_mu`.
- Persisted values win over code defaults: changing a default in code changes
  nothing on the master until the dataset is updated.

---

## 8. GeneralVariableScan (master-satellite)

### 8.1 Files and GUI

Shared implementation `GeneralVariableScan_master_satellite_mixin.py`
(`_GeneralVariableScanMasterSatelliteMixin(_DatasetRedirectMixin)`); public
files fix `EXPERIMENT_MODE` as a class constant (split `025c027`,
2026-08-31). `_single_node` adds `selected_node` (Node1|Node2); `_two_nodes`
has no node argument. GUI: `n_measurements`, `scan_variable1_name` /
`scan_sequence1`, `scan_variable2_name` / `scan_sequence2`,
`override_ExperimentVariables_Node1`, `override_ExperimentVariables_Node2`,
`control_experiment`, `experiment_function` (enumeration from the mode's
registry), `create_applets` (default True), and the "Catch Underflow" group.
Do not grow the standalone GVS into an architecture switch; the stacks stay
separate.

### 8.2 Function registries

- `build_experiment_function_registry(module, excluded)`: locally defined
  functions (`inspect.isfunction`, module ownership) whose name contains
  "experiment". No large allowlist.
- Single-node: `subroutines/experiment_functions.py` minus
  `HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS` = exactly
  `Two_nodes_atom_loading_experiment`, `Two_node_single_photon_experiment`,
  `Two_node_single_photon_2_experiment`,
  `Two_node_single_photon_2_optimization_experiment`. 51 functions at HEAD.
- Two-node: only `subroutines/experiment_functions_two_nodes.py`: 3 functions
  (`master_satellite_namespace_sanity_experiment`,
  `Two_nodes_atom_loading_experiment` (the native one; same name as an
  excluded historical function), `Two_nodes_alternating_retention_experiment`),
  pinned by `test_registry_contains_exactly_the_intended_experiments`.
- Function selection is registry lookup, never `eval()`. `eval` is used only
  for scan-sequence and override expressions.

### 8.3 `prepare()`

1. Capture the submitted `n_measurements`.
2. `base.configure_execution(EXPERIMENT_MODE, selected_node or None)`;
   re-assert `n_measurements` (R13).
3. Single only: `which_node = alice|bob` for legacy functions.
4. `base.prepare()`; then `base.prepare_laser_stabilizer()` (both modes;
   two_nodes delegates to `prepare_laser_stabilizers_both_nodes`).
5. Resolve scan variables through `resolve_experiment_variable_target` (no
   owning node, so a bare per-node name is ambiguous in two_nodes), evaluate
   scan sequences, merge overrides (section 8.5).
6. Set `needs_experiment_variable_reload` when an earlier RID is queued
   (`_has_earlier_queued_experiment`).

### 8.4 `run()`

`_initialize_run_state()`:

1. `reload_experiment_variables()` if queued behind others, then re-assert
   the submitted `n_measurements`;
2. `_apply_run_wide_overrides()`: resolve, `setattr`, refresh compatibility,
   refresh dependent state, and record each override (section 8.5). MUST stay
   first among dataset writes (R17), pinned by
   `test_override_records_precede_every_broadcast`;
3. `initialize_result_datasets()`, `initialize_result_state()`;
4. broadcast `scan_variables`, `scan_sequence1`, `scan_sequence2` (axis
   labels for the retention applet);
5. `_create_node_applets()` (section 14).

Then `force_other_node_off()` once (single_node only), then per scan point
`_run_scan_point` -> `_execute_scan_point`:

1. `_apply_scan_point`: mutate authoritative target(s), refresh
   compatibility and dependent state;
2. `initialize_hardware()` (core reset every point);
3. `reset_result_state_for_scan_point()`;
4. set the `iteration` dataset;
5. call the selected experiment function;
6. `write_results` (named copy, section 6.9).

Scans and overrides are non-persistent. Never call `build()` or `prepare()`
per point.

### 8.5 Overrides

- One dictionary per node, `override_ExperimentVariables_Node1` and
  `_Node2`; the shared `override_ExperimentVariables` field was removed
  (`5c1f508`, 2026-10-01). Only the active nodes' dictionaries are read (the
  selected node in single_node, both in two_nodes); the idle node's text is
  never evaluated, since its attributes are not even loaded.
- Names may be bare or suffixed and resolve against the owning field's node
  (`6a89075`): `{'t_pumping': ...}` in the `_Node1` field means
  `t_pumping_Node1`, even in two_nodes. Globals go in the running node's
  dictionary (name a global in one dictionary only). A target named in two
  dictionaries raises, listing both sources. Empty means no overrides.
- Merging is shared by GVS and `MicrowaveScanOptimizer` through
  `merge_per_node_override_dictionaries` in the Base module; each caller
  evaluates in its own `globals()` with `{'self': experiment}`.
  `submit_opt_exp_general` forwards both fields; old expids fall back to
  empty dictionaries.
- Recording (`ef55cae`, 2026-10-14): each applied, non-scanned override is
  also written as a LOCAL dataset, `set_dataset(target, value,
  broadcast=False, persist=False)`, under the name the run read (suffixed per
  node, bare global), so the h5 `datasets` group holds the value the run
  actually used, as standalone did. Motivation: RID 38903, a Z-parity run
  with `t_MW_RF_pulse` overridden to 0, was analysed as X parity from the
  archived 69.6 us. Verified by replaying RID 38903 through the real
  `DatasetManager`. Tests: `test_override_records_precede_every_broadcast`,
  `test_overrides_are_recorded_locally_like_standalone`. First hardware use:
  GVS single_node Node2 RIDs 38952, 38954 and 38957 (2026-10-14 20:10-20:28,
  `atom_loading_2_experiment` scanning `target_780_HWP`). Each override in
  `override_ExperimentVariables_Node2` is in the h5 `datasets` group under its
  suffixed name, e.g. `p_GRIN2_excitation_Node2` = -12.0 (archive -6.0),
  `t_MW_RF_pulse_Node2` = 0.0 (archive 69.6 us), `t_FORT_drop_Node2` = 10 us
  (20 us in 38957), `target_780_QWP_Node2`, `t_OP_coil_relax_time_Node2`,
  and the global `advance_if_atom_in_first_shot` bare. The Node1 dictionary
  was not evaluated, as designed (no `_Node1` keys). `n_measurements` = 100
  is still persistent in `dataset_db.pyon` afterwards, so "records first,
  broadcasts after" held on hardware. RID 38904 (before the commit) has none
  of these datasets. Older files: the analysis loader re-applies overrides
  (section 15.2).

### 8.6 GUI precedence for `n_measurements`

`n_measurements` is both a GUI argument and a persistent global. Prepare
captures the submitted value before `configure_execution()` loads datasets
and re-asserts it afterwards; the queued-run reload restores it again; the
run then persists the submitted value (section 6.9). Other calibration is
never changed by a run.

### 8.7 Underflow retry

Built into the mixin behind `enable_Catch_UnderFlow` (default off), mirroring
the standalone GVS (`9904f9b`, 2026-08-31; the separate
`GeneralVariableScan_CatchError_master_satellite_*` files and
`_CatchUnderflowRetryMixin` were deleted the same day). GUI group "Catch
Underflow": `underflow_max_retries` 10, `underflow_backoff_ms` 200,
`skip_only_that_iteration_if_exhausted` True. Only `RTIOUnderflow` is caught;
a retried point re-runs the full scan-point flow including the core reset.
DRTIO, SPI, resolution and all other errors propagate. With the flag off
there is no exception handling.

### 8.8 Two-node specifics

Scan variable names must be suffixed or global. Results are unsuffixed
(joint). Feedback runs per node via `laser_stabilizer_Node1/_Node2`. The
result tag is `TwoNodes`. No `force_other_node_off`.

---

## 9. Experiment functions

### 9.1 Policy

- `subroutines/experiment_functions.py` (17364 lines, 103 top-level
  functions) was written for standalone, including historical functions that
  coordinated two independent Kaslis. Reuse its functions in single_node
  mode; do not copy the module; do not migrate historical two-node functions
  automatically; do not rewrite validated sequences. The four historical
  functions are unchanged since `e86046f` (diff empty, checked 2026-10-14).
- New DRTIO two-node functions live only in
  `subroutines/experiment_functions_two_nodes.py` (renamed from
  `experiment_functions_master_satellite.py` in `025c027`). It imports nothing
  from `experiment_functions.py` and deliberately COPIES `first_shot`,
  `second_shot`, the loader and `end_measurement`: a fix to a single-node
  original does not propagate.
- Shared edits to `experiment_functions.py` are allowed under R1/R8 (slack
  restores, confirm-by-first-shot).

### 9.2 Single-node registry (51 functions, 2026-10-14)

Includes 11 `microwave_*` experiments; `single_photon_experiment_3` (three
variants) and 4-9; `atom_photon_parity_2` to `_11` (with `6_node2`,
`9/10/11_AllSPCM`); `atom_loading(_2/_for_*/_optimizer)_experiment`;
`beam_balancing_*`; `blowaway_fidelity`; `test_*`; `tune_coil_driver`. Count
history: 50, then 51 with `atom_photon_parity_11_AllSPCM` (2026-09-03), 54 on
2026-10-13 (parity 12/13 and single_photon_experiment_10), back to 51
(`2a3592c`).

### 9.3 Two-node module conventions

- Sequences keep the legacy linear shape (R9) with every device naming its
  node (R11). Per-node VALUES (e.g. `f_FORT` 245 vs 240 MHz, coil voltages),
  DURATIONS (`t_PGC_after_loading` 0.6 vs 1.0 ms) and BRANCHES are written out
  per node.
- Stages whose durations differ per node use `with parallel` (section 12.6);
  same-timestamp pairs do not.
- Switching events are grouped in fours with `delay(5*us)` (R14), guarded by
  `test_no_run_of_switch_events_shares_a_timestamp`.
- Tests in `tests/test_experiment_functions_two_nodes.py` pin that
  value-carrying calls name a node, that no raw physical names appear, and the
  SED limits.
- Non-entry helpers: `record_FORT_powers_both_nodes`,
  `record_CW_optical_pumping_both_nodes`, `run_feedback_and_record_FORT_MM_power`
  (sequential per-node feedback), `first_shot`, `second_shot`,
  `recooling_after_first_shot`, `load_until_atom_in_both_nodes_recycle`,
  `load_until_atom_in_both_nodes_alternating_shot_recycle` (disabled),
  `two_node_alternating_shot`, `two_node_alternating_windows`,
  `alternating_first_shot`, `alternating_second_shot`, `end_measurement`.

### 9.4 `master_satellite_namespace_sanity_experiment`

Host-side `hasattr` checks only, non-destructive by itself; through GVS it is
NOT passive because `initialize_hardware` runs first. Ran on hardware
2026-10-06/07.

### 9.5 Two-node atom loading

`Two_nodes_atom_loading_experiment` (ported `af8b343`, 2026-10-07) is
section 4 of `experiment_functions.py` in the same linear shape, with the TTL
handshake deleted. Three hardware bugs were found and fixed in the port:

1. SED discard: 18 switch events on one timestamp; events past 8 were dropped
   (section 12.5).
2. RO1 ran on six beams and RO2 on four: `first_shot` had no else branch for
   A5/A6 while `second_shot` turned them off (with
   `PGC_and_RO_with_on_chip_beams` True), making retention meaningless. Fixed
   in the two-node copy only; still present in the single-node original
   (section 18, O-9).
3. Node1's FORT was commanded to 0 at every readout:
   `FeedbackChannel.amplitudes` is `np.zeros` with only `[0]` seeded, and the
   interleaved feedback ran Node1 at index 2 and Node2 at index 1, so
   `stabilizer_FORT_Node1.amplitudes[1]` stayed 0. Each node now runs both
   indices; every stabilizer run gets its own `break_realtime`.

Hardware: RIDs 38670-38684 (2026-10-07); RID 38710 (2026-10-08, a
`Two_nodes_alternating_shot_experiment` run using this joint loader): median
`Atom_loading_time` 0.65 s, 90th percentile 3.0 s, mean 1.24 s over about
1790 loads, so about 64 tries at 10.1 ms each (re-read from the h5 on
2026-10-14; commit `19c289f` quoted 0.85 s / 3.3 s / 84 tries). Count rates
(2026-10-08, commit `b6b4de1`): joint no-atom 56200 c/s; joint both atoms
106400 c/s (about 25100 c/s per atom); common background (FORT scatter plus
dark, not gated by the cooling DP) about 39900 c/s; half the MOT fluorescence
about 16300 c/s; the two per-node windows summed with no atom 96100 c/s (the
background is counted twice); per-node no-atom floor about 48050 c/s (half of
96100). Since `c1ddbd2`
(2026-10-12) it carries the confirm-by-first-shot flag (section 9.7) and
replaces `Two_nodes_atom_loading_2_experiment`. Known gap: both two-node
loaders run PGC at the LOADING FORT depth, whereas all four single-node
loaders drop the FORT to `stabilizer_FORT.amplitudes[1]` first (O-8).

### 9.6 Alternating (node-resolved) readout

- `Two_nodes_alternating_retention_experiment` (`c1ddbd2`, 2026-10-12)
  replaces both joint shots with `alternating_first_shot` /
  `alternating_second_shot`: interleaved per-node readout windows on one
  timeline. Retention(NodeX) = fraction above threshold in RO2_NodeX among
  those above threshold in RO1_NodeX, against
  `alternating_shot_threshold_NodeX` (20000/21000). Datasets
  `AllSPCMs_alternating_RO1/RO2_NodeX`; the "alternating" infix is
  load-bearing, because `AllSPCMs_RO1_NodeX` would be stripped by the analysis
  loader into the joint `AllSPCMs_RO1` and collide (that spelling is reserved
  for a future node-resolving version of the ordinary shots). Hardware: RIDs
  38773, 38774, 38778 (2026-10-12; 38773: 1000 measurements, mean RO1 counts
  Node1 256 / Node2 266, RO2 213 / 232).
- `single_atom_threshold_NodeX` does NOT transfer to alternating windows: the
  other node's FORT is on and scatters into the shared detectors (a pedestal
  the single-node calibration never saw).
- Gate closing: `gate_rising_mu` emits a close event as well as an open
  event. Closing four counters and switching the light off at the same
  instant is 6 events on one timestamp, so the light-off goes one pad later
  (4 + 2 events) and the counted window is fully lit (`30d7be1`). Function
  entry validates `t_alternating_RO_pad`, `t_alternating_RO_window` and
  `t_Node2_rtio_offset_mu` (a pad smaller than the offset turns an absolute
  `at_mu` into a backwards jump; a pad of 0 collapses the windows).
- `f9a587a` (2026-10-08) added the per-node FORT drop inside `with parallel`
  before the alternating shot.
- The per-node alternating atom check in the loader
  (`load_until_atom_in_both_nodes_alternating_shot_recycle`) reports success
  but loads nothing (RID 38718 vs RID 38710: readout saw 0/100 on Node2 and
  1/100 on Node1 while the loader reported 100/100). Raising the loading
  thresholds to 62000/60000 did not help. Untested suspicion: the window grid
  toggles each node's cooling DP and repump every window and nothing
  re-asserts the MOT AOMs, FORTs or coil DACs afterwards. It is disabled: on
  no compiled path, not type-checked, `TWO_NODE_ALTERNATING_EXPERIMENTS = ()`
  in `applets_master_satellite.py`. If it worked: one try =
  2*n*(window + 3*pad) = 26 ms, median near 2 s, about 50% loading duty.
  `19c289f` proposed a joint pre-check as the cheapest speed-up, on the
  argument that a per-node pass implies a joint pass (51000 + 45000 = 96000 >
  89000). That holds for the declared values (62000 + 60000 = 122000 >
  89000) but NOT for the live ones (50000 + 45000 = 95000 < live
  `two_atom_threshold_for_loading` 100000); recheck with the values in force
  before relying on it.
- `Two_nodes_alternating_shot_experiment` was deleted in `c1ddbd2`; docstrings
  still name it (O-25).

### 9.7 Confirm-by-first-shot advance

Purpose: retry loading until the first shot really sees an atom, so a
measurement does not run on an empty trap; loading-plus-survival probability
= 1 / mean(`n_first_shot_confirm_tries`).

- Flag: global `advance_if_atom_in_first_shot` (declared False, live True);
  bound: global `max_loading_rounds` (20). On exhaustion the measurement is
  kept (simply advances). The try count of each kept measurement is appended
  to `n_first_shot_confirm_tries`.
- Single-node form (pilot `0edd7ce` on `microwave_Rabi_2_experiment`, then
  `f36a225` (10 microwave), `34bff1e` (9 single-photon), `b770aa0` (parity 2,
  3, 6_node2, 8), `5e10655` (parity 4, 5, 6, 7, 9, 10, 11), all 2026-10-13; 31
  sites): `if self.advance_if_atom_in_first_shot:` a while loop around the
  experiment's own load + `first_shot` span, accepting when
  `AllSPCMs_RO1 / t_SPCM_first_shot > single_atom_threshold` (counts/s on
  both sides); else the original span byte for byte. On a REJECTED attempt
  only, it sets `self.AllSPCMs_RO2 = 0` and calls `core.break_realtime()`;
  after the loop it appends the count and calls `break_realtime`. Never reset
  RO2 before attempt 1: the loaders' recycle pre-check reads it, and clearing
  it forces a full MOT load every measurement. Retrying is safe because
  `first_shot` writes no datasets (all appends happen in `end_measurement`).
  Sites that arm the microwaves in the span (`ttl_microwave_switch.on();
  dds_microwaves.sw.on()`) duplicate the arming into both branches.
- Two-node form (`Two_nodes_atom_loading_experiment`): accepts on
  `two_atom_threshold`; it resets `self.AllSPCMs_RO2 = 0` at the top of EVERY
  attempt, so with the flag on the recycle shortcut is never taken (O-7).
- Generated, not hand-edited, with the proofs of section 3.4.
- Hardware evidence (h5 scan 2026-10-14; seed excluded):

| Experiment | Node | RIDs / period | Measurements | Mean tries | Max |
|---|---|---|---|---|---|
| `atom_photon_parity_11_AllSPCM` | Node2 | 15 runs, RIDs 38802-38904 (10-12..14) | 350,416 | 1.109 | 8 |
| `MicrowaveScanOptimizer` (Rabi_2 / map01_map11) | both | all 42 runs with data, 2026-10-12..14 through RID 38947 | 33,239 | 1.080 | 6 |
| `MicrowaveScanOptimizer`, 2026-10-14 only | both | 34 runs, RIDs 38898-38947 | 26,015 | 1.035 | 5 |
| `single_photon_experiment_9` | Node2 | 8 runs, RID 38803 and 38887-38895 | 3,420 | 1.261 | 5 |
| `Two_nodes_atom_loading_experiment` (flag True) | both | RID 38785 (10-12) | 83 | 1.60 | 5 |
| `Two_nodes_atom_loading_2_experiment` (removed) | both | RIDs 38766, 38770 (10-11) | 20,000 | 2.78 | 20 (1.1% at cap) |

  RID 38785 is the only two-node run with the flag, and there were no
  two-node runs at all from RID 38786 through 2026-10-14, so O-7 is
  unmeasured.
  Single-node means of 1.07-1.26 mean 80-93% of first shots confirm on the
  first attempt; about 10% are rejected on Node2 parity 11 at
  `single_atom_threshold_Node2` = 17000. The recycle path is intact (parity
  11 Node2 about 5.6 measurements per atom load, as before the change, RID
  38812). Archived `advance_if_atom_in_first_shot` = True from RID 38784 on;
  `max_loading_rounds` archived from RID 38836 on. Unwrapped functions
  (`atom_loading_2_experiment`, RIDs 38840-38849) show seed-only series, as
  expected; GVS `microwave_Rabi_2_and_EXC` RIDs 38856-38863 are also
  seed-only, cause not verified (O-12).
- Out of scope by decision: the 7 loading/optimizer/beam-balancing helpers (a
  retry loop would destroy their loading-probability measurement).

### 9.8 Loaders, recycling and throughput

- Loaders' already-loaded pre-check reads `AllSPCMs_RO2`; anything that zeroes
  it before the loader defeats atom reuse.
- Measure throughput from the h5 before blaming gateware:
  `meas_per_load = len(AllSPCMs_RO1) / (len(AllSPCMs_atom_check_in_loading)/2)`
  (two atom-check entries per load); seconds per measurement = median diff of
  `atom_loading_wall_clock` (mu / 1e9) divided by meas_per_load.
- Measured 2026-10-13: standalone `parity_9_AllSPCM` 28-29 meas/load at 0.029
  s/meas (RIDs 38187, 38216), falling to 5.1 two days later on the same
  gateware (RID 38247); master-satellite parity 11 Node2 5.6 (RID 38812);
  parity 12 0.79-0.81 (RIDs 38816, 38828) because it reset RO2 before attempt
  1. Load-to-load time is shorter on master-satellite Node2 (0.125 s vs
  0.84-1.04 s). Conclusion: the DRTIO switch cost nothing.
- Parity 11 and `microwave_Rabi_2_and_EXC` loaders use
  `load_MOT_and_FORT_until_atom_recycle` (`43455f7`, `1cb0cc1`).
- Parity 13 showed padding outside the core can be trimmed (28.11 ms to 7.39
  ms per measurement by inlining `end_measurement` and moving 24 appends and
  one blocking round trip to the end-of-point flush); not carried forward.

### 9.9 Other sequence facts

- `atom_photon_parity_11_AllSPCM_experiment` (`4e636da`, 2026-09-03) is
  parity_10 plus complete SPCM-edge draining from sp9: with no valid
  single-detector event parity_10 never drained, so stale FIFO edges
  reappeared as false early counts and underflowed the `at_mu()` mapping
  jump.
- sp9 Node2 optical pumping uses `record_CW_optical_pumping_node2` with the
  coil setting (`e79b0fd`, 2026-10-05); its duplicate OP coil set is
  commented out because the recorded trace already sets those DACs.
- `record_CW_optical_pumping_both_nodes` records ONE DMA trace for both
  nodes. `playback_handle` ends with `delay_mu(advance_mu)`, so `with parallel`
  around two playbacks overlaps only the submission and the traces still run
  in sequence; recording both nodes inside one parallel block per step is the
  only way to move the coil ramps together. Where `t_OP_coil_relax_time`,
  `t_pumping` or `t_depumping` differ, the larger value is used for both.
- Loader re-assertions `dds_microwaves.sw.on(); dds_FORT.sw.on()` after
  feedback ("# bug -- microwave dds and FORT are off after AOM feedback") are
  an old symptom patch for an underflow on channel 56 now fixed by the
  unconditional `break_realtime` in `run_feedback_and_record_FORT_MM_power`.

### 9.10 Removed functions and files (for orientation)

Legacy two_nodes node mode and `config/two_nodes` (`e86046f`, 2026-08-26);
`GeneralVariableScan_CatchError_master_satellite_*` and
`_CatchUnderflowRetryMixin` (`9904f9b`); standalone
`GeneralVariableScan_CatchUnderflow.py` (`1c428e4`); the combined
`AOMsCoils_master_satellite.py` (split per node in `a7d8dff`, then renamed to
`AOMsCoils_master_satellite_mixin.py` in `02b5820`: converted, not deleted;
the mixin itself was deleted 2026-10-19 when both node files became
self-contained, section 11.7); `h5_data_utilities_node1.py`
(`2a5796a`); `ms_filter_by_node` (`0b59304`); the shared override field
(`5c1f508`); the `3a90e2d` alternating diagnostic (reverted `ae93327`);
`t_coil_relaxation_time_OP_NodeX` (`0e54c84`);
`Two_nodes_alternating_shot_experiment` and
`Two_nodes_atom_loading_2_experiment` (`c1ddbd2`);
`single_photon_experiment_10` with `_excitation_attempt_of`,
`excitation_attempt_period_mu`, `record_excitation_group`, and
`atom_photon_parity_12/13` (`2a3592c`, "did not work well in practice"); the
MonitorSPCMinApplet node selector (`50f9049`); `max_loading_rounds_two_node`
(`2fb2fbf`).

---

## 10. Feedback scope

(This section carries the content of the retired detail spec's section 18,
cited by `FORT_Polarization_Optimizer_master_satellite.py`.)

### 10.1 What is fed back

Active feedback is AOM1-AOM6 (fiber AOMs A1-A6) and the FORT with
loading/science/holding setpoints (`set_point_FORT_APD_loading/science/holding`
in `feedback_channels.json`). D1 feedback is inactive: keep its history, do
not expand it, do not invent `set_point_D1_DP`. Every `feedback_channels.json`
entry is `"series": true`.

Node2's AOM A1 photodiode is read on logical `sampler1` ch6 = physical
`sampler4` ch6 since 2026-10-16 13:24 (bob `feedback_channels.json`, a
`"sampler1"` group placed AFTER `"sampler0"`; committed 2026-10-16). Its old
channel, logical `sampler0` ch7 = physical `sampler3` ch7,
is bad (the user's finding): `MOTBeamPowerCheck_master_satellite` RIDs 39018,
39019 and 39021 read a flat 0.0215 V there at every power against a 0.06 set
point, and feedback (RID 39020) could not hold it. On the new channel feedback
holds `MOT1_monitor_Node2` at about 1.00 of the retuned
`set_point_PD1_AOM_A1_Node2` = 0.042, at -16.6 to -18.0 dB (RIDs 39027-39032).
Node2's stabilizer now samples two cards (`sampler3`, `sampler4`), as Node1's
always has; see O-36 for the one limitation that brings.

### 10.2 Single-node path

`Base.prepare_laser_stabilizer()` builds `laser_stabilizer` from the projected
legacy namespace (`fast_feedback_dds_list`, samplers, set points,
`dds_defaults`, `which_node`) with `node_suffix=''`, installs a name map from
the selected node's `utilities/config/<alice|bob>/feedback_channels.json`,
and `_DatasetRedirectMixin` routes reads and writes to the node-suffixed
authoritative datasets (`p_AOM_A1_NodeX`, `*_history`, `feedbackchannels_NodeX`).
It publishes `self.stabilizer_<channel>`.

### 10.3 Two-node path (per crate, sequential)

Since `e768d50` (2026-10-07), `prepare_laser_stabilizer()` in two_nodes
delegates to `prepare_laser_stabilizers_both_nodes()`, which builds one
`AOMPowerStabilizer` per crate (`laser_stabilizer_Node1/_Node2`,
`stabilizer_*_NodeX`, `feedbackchannels_NodeX`) using `node_suffix` /
`config_node`. The suffix is applied to the config directory, sampler/DDS
lookups, set points and every dataset name, so each node writes the same
datasets as a single-node run. Ambient devices (repumps, cooling DP, six MOT
AOMs) are bound through the suffix at construction because `getattr` is
unavailable in a kernel. `run_feedback_and_record_FORT_MM_power` in the
two-node module runs `laser_stabilizer_Node1.run()` then
`laser_stabilizer_Node2.run()`. With suffix `''` standalone is unchanged.
(The retired spec's "aom_feedback.py stays byte-identical to main" and "two-node
mode must not call either" are superseded; section 19.)

The dB history is built as `power_dataset + '_history' + suffix`, so a two-node
run writes `p_AOM_A1_history_Node2`, the name a single-node run writes and the
`feedback RF NodeX` applets plot (2026-10-20). Before, it came out
`p_AOM_A1_Node2_history` (every two-node h5 up to RID 39108), which no applet or
notebook reads, so those applets did not move during two-node runs. Test:
`test_two_node_feedback_writes_the_histories_the_feedback_rf_applets_plot`.

### 10.4 Compile requirement

Kernels read `self.stabilizer_<channel>` and host scalars such as
`AllSPCMs_RO2`, `*_in_health_check`, `*_test_RO`, `atom_loading_time`,
`measurement`; the ARTIQ compiler types them from the host object, so a
missing attribute fails to COMPILE, not to run, and the hardware-free suite
cannot catch it. Every GVS run therefore calls `prepare_laser_stabilizer()` in
`prepare()` and `initialize_result_state()` in `_initialize_run_state()`, in
both modes (`6f4001b`, 2026-09-09; extended to two_nodes in `e768d50`).

### 10.5 Parallel feedback: not implemented

- Feedback is strictly sequential, per node and across nodes. The global
  `parallel_AOM_feedback` (True) is read by no code today, but the retired
  PLAN.md (18.2, 21.1) created it as THE selector for this feature (False =
  sequential Node1/Node2 feedback, True = concurrent); keep it declared until
  the `feedback_group` design decides its fate (section 7.7). Premise:
  parallel feedback was required because Node1's and Node2's beams are
  optically independent, so cross-node concurrency (the two crates' series
  loops overlapping) is allowed even with every channel `series`; the
  per-node grouping below is the further step. The dormant
  `parallel_channels` machinery and a commented-out parallel block (disabled
  2025-06-03) are in `subroutines/aom_feedback.py`; `parallel_channels` has
  always been empty because every channel is `series: true`.
- Gate: do not implement parallel feedback until sequential feedback is
  validated on hardware on both nodes (now true) and the timing analysis
  below is done. Node2 feedback exercises remote Sampler SPI over DRTIO.
- Decided grouping (lab owner, 2026-10-08; comment in `aom_feedback.py`,
  `12c5314`). The constraint is OPTICAL, not electrical: no beam may disturb
  another's photodiode. Group 1 = A1, A2, A5; group 2 = A3, A4, A6; group 3 =
  FORT alone, repeated per setpoint. Each group converges before the next;
  both nodes' members of a group run together. Implementation path: replace
  the boolean `series` with `"feedback_group": 1|2|3` (three ORDERED groups);
  within a group switch on the group's DDSs, `measure()` ONCE and feed back
  every member from the same array. Since `measure()` reads every card in
  `sampler_list` regardless of channel count, this gives about 57% fewer
  measures. Alternatives considered: cross-node sampler rewiring (about 33%,
  hardware plus a merged loop); moving Node1's FORT PD from sampler1 ch6 to a
  free sampler0 ch1/2 (about 33% for one cable; Node2 already has all 7 on
  one card).
- Measure before optimizing: time one `laser_stabilizer_NodeX.run()` against
  7 x iterations x 10 ms to infer the real iteration count. Per average:
  150 us + sample + 100 us = 250 us of settling; 20 averages give 5 ms per
  card per measure. If the loops do not converge early, the setpoints are the
  real problem.
- Timing caution (retired PLAN.md, still binding): do NOT assume two feedback
  kernels can simply be wrapped in `with parallel`. Inspect kernel
  boundaries, RTIO timing, `break_realtime()`, Sampler acquisition, DDS/SPI
  updates, delays, host/kernel interactions, timeline-cursor behaviour, DRTIO
  latency, shared resources and SED lanes (section 12.5).

### 10.6 FORT feedback details

- `FORT_MM_monitor_NodeX` and `FORT_APD_monitor_NodeX` are broadcast-only
  transients; their applets start empty after a master restart until feedback
  republishes them.
- `prepare_laser_stabilizer` logs WARNING lines for `MOT1..6_monitor_NodeX`
  and `FORT_monitor_NodeX` (pre-existing; created on first write).
- The two-node FORT index bug of section 9.5 (item 3).

### 10.7 Slack and feedback

A `laser_stabilizer.run()` returns about 5.8-6 ms behind the wall clock: it
advances the timeline only by coded delays while sampler reads and host RPCs
take real time, and `aom_feedback.py` manages no slack. Master-satellite
settings are heavier than the standalone DECLARED defaults (averages 20 vs
10; iteration cap 100-200 vs 15), but the live standalone datasets are 200
iterations / 20 averages (section 7.6), so on this master the two stacks are
in fact comparable. Rule: `break_realtime` AFTER every `run()`; re-arming BEFORE each run
only moved the underflow (from -4.9 ms inside a warm-up loop to -5.8 ms at
`dds_FORT.sw.off()`). Sites still unprotected: section 18, O-10.

### 10.8 What waits on parallel feedback

`FORT_Polarization_Optimizer_master_satellite` two_nodes design (section
11.2) waits for two-node PARALLEL FORT feedback. Two-node SEQUENTIAL
feedback exists; until parallel exists, a two-node polarization run would
have to feed back sequentially or run at the persistent default FORT powers.

---

## 11. Optimizers and other public experiments

### 11.1 `MicrowaveScanOptimizer_master_satellite`

Committed behaviour (`1277477`, 2026-09-01, with later fixes):

- `which_node = Node1 | Node2` (`VALID_NODES`); two_nodes is deliberately
  unsupported and never planned (R22). Standalone
  `standalone/MicrowaveScanOptimizer.py` stays untouched (last change: the
  2026-09-03 move) and keeps serving the standalone HealthCheck files, whose
  expids point at the standalone class.
- `scan_options`/`scan_dict` are extracted from the standalone source with
  `ast.parse` + `ast.literal_eval` (`_load_standalone_scan_definitions`), with
  no import, so scan definitions cannot drift. `fit_model_dict` is defined
  locally (the module docstring wrongly says it is shared).
- Base supplies the derived RF amplitudes, the result state, `print_async`
  and `write_results`. Calibrations persist node-suffixed
  (`MICROWAVE_CALIBRATION_DATASETS`); `fit_parameter_*` and `parent_rid`
  (persisted, unsuffixed) are run state.
- Scans and overrides mutate suffixed authoritative attributes and refresh
  the projection; queued freshness uses `reload_experiment_variables()`.
  Per-node override dictionaries via `merge_per_node_override_dictionaries`
  (active node only). The optimizer records overrides with `set_dataset(var,
  val)` in `initialize_datasets`.
- Health-check failure resubmits THIS class (`default_expid`: file
  `qn_artiq_routines/MicrowaveScanOptimizer_master_satellite.py`,
  `which_node`, both override dictionaries, `parent_rid`). Future
  master-satellite HealthCheck files must submit the same expid shape; none
  exist yet.
- Applets via `create_applets_for` (so it gets the retirement pass) plus the
  shared set.
- Fixes: `cb7a0ab` (2026-08-31) one `experiment_name` instead of
  `experiment_name_alice/_bob`. `bb72047` (2026-09-02): adaptive-refine center
  shifts used `freq_scan_min_step_size_kHz` as a raw number (about 10 Hz)
  instead of a kHz offset; retention windows made seed-aware
  `[1+i*n : 1+(i+1)*n]` in both optimizers and both standalone HealthCheck
  files (previously every window held one foreign shot). `7ea1476`
  (2026-09-17): centring decided by wing ASYMMETRY, tested before centring,
  depth no longer gating the skew cases, tolerance 0.15 (about 2 sigma at
  n = 100); RID 38471 went from "well centered" to right-skewed (+50 kHz)
  while RID 38470 stays centred; probe indices 2 and 3 sit at
  ±x*(1 - 1/shrink_factor) (±0.6x at 2.5), not ±x/2; the
  `retention_array[0:4]` unpack is guarded (RID 38469 produced 10 windows for
  13 points).
- Hardware: used on both nodes since 2026-09-16. On 2026-10-14 through RID
  38951 there are 44 RIDs, 37 of them with measurements (Node2 38885-38901
  and 38951, Node1 38907-38950; e.g. parent RIDs 38949 Node1, 38901 Node2),
  with the confirm loop; Node1 runs continued with RIDs 38961-38965.

Written 2026-10-14, committed 2026-10-16 as `619d3ec` (the two fixes below)
and `bfca64d` (the magnetometer seeding). Because the master runs the working
tree (R26), it ran on hardware before it was committed: the file was saved at
19:44 and RIDs 38942-38951
(Node1 38942-38950, Node2 38951; 19:54-20:04) and 38961-38965 (Node1,
20:45-20:48) started after that and executed it (38951 records its override
`t_pumping_Node2` = 0.1 ms, and `n_measurements` stayed persisted). Neither
fixed failure path has been exercised: no run overrode `n_measurements`, and
no Frequency_11/m10/m11 scan has completed (38963-38965 are Frequency_11
runs that stopped after 57, 0 and 1 measurements, before any scan point
finished). Whether the fit check now prints the pre-scan calibration must
be confirmed from a completed run; through 2026-10-16 no Frequency_11/m10/m11
scan has completed (the later optimizer runs are Frequency_00/01). The
changes:

- `initialize_datasets` records the applied overrides FIRST, before
  `base.initialize_result_state()` broadcasts and persists `n_measurements`.
  Without this, `{'n_measurements': ...}` in an override dictionary deleted
  `n_measurements` from the master database (reproduced through the real
  `DatasetManager` on 2026-10-14).
- `value_before_scan = getattr(self, scan_variable1)` is captured before the
  scan and restored with `_apply_scan_point` afterwards. Scan points are
  written into the scanned attribute, which for the 11/m10/m11 scans IS the
  calibration the fit check reports as "original": RID 38934 (a Frequency_11
  scan of a 339.213 MHz calibration) printed "Original f0 value: 339.073
  MHz", the refined centre.
- New tests `test_fit_check_sees_the_calibration_from_before_the_scan` and
  `test_no_local_write_follows_a_broadcast_of_the_same_name` (14 to 16 tests).
  Unit tests (152) and the compile sweep pass with these changes. Treat as
  pending until the user verifies the two paths above and commits.
- Magnetometer seeding (added about 21:04, O-35): `reset_datasets` now begins
  with `self.base.reset_result_state_for_scan_point()`, as GVS does, so the
  node's `Magnetometer_*_NodeX` series exist before every experiment call --
  health checks included, since every call is preceded by `reset_datasets`
  -- and are cleared per scan point. The optimizer never created them before:
  RID 38965 (Node1 Frequency_11 health check, 20:48) failed with "Cannot
  mutate nonexistent dataset 'Magnetometer_Zero_X_Node1'" after the
  `end_measurement` it ran called `measure_Magnetometer` (O-34). Test
  `test_magnetometer_results_exist_for_every_experiment_call` (health-check
  and scan paths, both nodes); removing the call reproduces the KeyError.

### 11.2 `FORT_Polarization_Optimizer_master_satellite`

- `614549e` (2026-08-31). `which_node = node1 | node2` (lowercase,
  `NODE_BY_SELECTION`). Scans the selected node's 852 HWP/QWP through the
  single `k10cr1_ndsp` controller (suffixed axis names), reads powers with the
  unchanged legacy helpers, and stores monitors and calibrations node-suffixed
  through `POLARIZATION_RESULT_DATASETS` (`FORT_MM_monitor_NodeX`,
  `best_852HWP_to_max_NodeX`, `best_852QWP_to_max_NodeX`,
  `best_852_power_NodeX`, `best_852_power_ref_NodeX`).
- `tolerance_deg`, `full_range`, `sample_pts` are GUI-wins values re-asserted
  after `base.prepare()`.
- Hardware (9 runs): Node1 RIDs 38539 (2026-10-01), 38583-38585 (10-05);
  Node2 RIDs 38611, 38612 (10-05), 38737 (10-09), 38759 (10-11), 38772
  (10-12, the last run).
- `369f805` (2026-10-05): `break_realtime` after each
  `stabilizer_FORT.run(setpoint_index=0)` (underflow on channel 36), and a
  per-node "waveplate trajectory" applet in a K10CR1 group.
- Future two_nodes design (documented in the module docstring, deliberately
  not implemented; waits on section 10.8): `which_node` gains `two_nodes`;
  both nodes march through one grid by point index; per point all four
  waveplates are set first, then each node records; k10cr1 moves are blocking
  RPCs on one controller (true concurrency needs an async move in
  `ndsp/k10cr1/k10cr1_driver.py` plus a server restart); the legacy record
  helpers hard-code `self.sampler0/1` and `which_node`, so parameterized
  kernels are needed; each node normalizes MM by APD /
  `set_point_FORT_APD_loading_NodeX`, tracks and re-centres on its own best;
  the per-row early stop is per-node state; finish per node and persist the
  four `best_852*_NodeX` values.

### 11.3 `AtomLoadingOptimizer_load_until_atom_master_satellite`

`73a11e3` (2026-10-02), M-LOOP, `which_node` Node1|Node2. Port details:
`which_node` replaced by the legacy alice/bob presentation (because
`aom_feedback` builds `utilities/config/<which_node>/` and no
`config/Node1` exists); best-parameter writes for `set_point_PD*_AOM_A*` are
explicitly suffixed through `resolve_experiment_variable_target` (only coil
volts are covered by the redirect); `t_SPCM_exposure` and `n_measurements` GUI
values captured and re-asserted; `break_realtime` after every stabilizer run;
`run_feedback_and_record_FORT_MM_power` records the FORT powers.
`compare_prev_and_optimized_values` (default on) re-measures both and keeps
the optimized values whenever they beat the previous values; `target_cost` is
reported, not required, since 2026-10-16 at the user's request (it used to be
a second bar, which discarded results like RID 38552's -2793 against a -8000
target). `analyze()` tolerates M-LOOP 3.3.4's 0-d `in_costs`
IndexError (RID 38545 had been marked failed for that alone). Hardware RIDs
38540-38562 (2026-10-01..02). M-LOOP aborts before `run()` in the compile
harness, so `set_experiment_variables_to_best_params` was compiled
separately.

### 11.4 `SamplerMOTCoilAndBeamBalanceTune_master_satellite`

Converted in place 2026-09-15 (`a5e685b`), split out 2026-10-02 (`ecc83ca`)
with a renamed class so ARTIQ does not find two classes of one name;
`MOT_experiments/SamplerMOTCoilAndBeamBalanceTune.py` restored to `70dd16c`
for standalone. One node per run (`selected_node`). Re-arm after each
warm-up feedback run (fixed channel 39, -6.38 ms); the tuning loop's
`delay(100*ms)` is real settling and was left alone; the finishing block is
guarded before the Zotino write. `0ff183e` added `force_other_node_off()`.
Hardware RIDs 38568-38572, 38595-38597.

### 11.5 `MonitorMagnetometer_master_satellite`

`d41b9d6` (2026-10-05). Two independent selections: `selected_node` (the node
that runs; the other is driven dark) and `magnetometer_node` (whose
magnetometer is read; choosing the idle node measures cross-talk). Reads the
node's `sampler2` channels 1/2/3 (section 2.9); `run()` inits the chosen
sampler explicitly. `V_to_mG` defaults to 1.0 (raw volts; the notebook
applies 625). Standalone units disagreed: alice stored mG (x350; RID 34618
about -1651), bob stored volts (RID 38593, -4.64 V). Datasets are unsuffixed
`Magnetometer_X/Y/Z`. No hardware run of it exists in results through
2026-10-14; the only magnetometer run since the migration was standalone
`tests/A6-MonitorMagnetometer.py` (RID 38593).

### 11.6 `tests/MonitorSPCMinApplet.py`

Converted in place 2026-09-16 (`0fea421`): reads the canonical SPCMs, writes
`AllSPCMs_rate_dataset` and the per-SPCM rates (unsuffixed), seeds its
datasets in `prepare()`. Since `50f9049` (2026-10-08) it has no node argument
and is always single_node Node1, because the SPCMs are master-local; choosing
Node2 had made it wait for the satellite and initialize the wrong crate.
Hardware RIDs 38698, 38699, 38702. RID 38457 (pre-port) carries both suffixes
and is refused by the analysis loader.

### 11.7 Manual control: `AOMsCoils_master_satellite_Node1/_Node2`

- One public file per node (`a7d8dff`, 2026-08-31), separate from GVS and
  from standalone `AOMsCoils.py`. No two_nodes mode and no node selector;
  `NODE` is a class constant.
- 2026-10-19 (user's decision): both per-node files are self-contained and
  independent, written box by box like standalone `AOMsCoils.py`; they share
  no code and never import each other, so one node's hardware change cannot
  reach the other. The private `AOMsCoils_master_satellite_mixin.py`
  (`_AOMsCoilsMasterSatelliteMixin`) they used to share is deleted. Both
  conversions keep the mixin versions' GUI arguments (names, groups, order)
  and hardware operations: op logs, delays included, compared over 557
  (Node1) and 696 (Node2) argument sets, each through `apply_manual_state()`
  and `run()`; planted one-line mistakes are caught. One timing difference:
  after laser feedback the mixin paused 1 ms between re-applying the DDS
  group and the gate group, and `turn_on_AOMs()` does not. Both take their
  coil channels from Base's wiring table (`coil_channels`). Each file:
  `everything_off()`, `turn_on_AOMs()` (one block per box),
  `set_microwaves()`, `set_coils()`, `run_feedback()`,
  `k10cr1_operations()`. Tests: `test_node_files_are_self_contained_experiments`;
  `test_every_k10cr1_box_moves_only_its_own_nodes_axis` ties the
  hand-written waveplate axes of both files back to Base; the coil,
  microwave-interlock and feedback-gate tests run on both nodes.
  Hardware: not yet run in this form. The last AOMsCoils runs (Node1 RIDs
  39127-39134, 2026-10-19 17:47-18:01) used the mixin.
  The rest of this section describes the current behaviour.
- 2026-10-19 (user): after laser feedback the stabilized AOMs come back on
  at the corrected power, as in standalone. `turn_on_AOMs()` no longer sets
  the box DDSes:
  `initialize_hardware()` has set each DDS to its binding frequency and
  power, and the stabilizer leaves its corrected amplitude in each DDS it
  adjusts (`fast_feedback_dds_list`: A1-A6 and the FORT on both nodes).
  Before, `turn_on_AOMs()` set every ticked DDS to its prepare-time power,
  so after feedback the correction was saved for the next run but undone for
  this one (inherited from the mixin). The one set() left is Node2's
  combined GRIN2 excitation tone. Compared with the version before over 589
  (Node1) and 824 (Node2) argument sets, both methods, with the test
  stabilizer and with one that changes the powers: identical apart from the
  removed set() calls, which rewrote the value the DDS already held except
  (a) the post-feedback overwrites this removes and (b) on Node2, with
  `D1_pumping_switch_ON` and `Node2_GRIN2_AOM_ON` both ticked, a ~9 ms rewrite
  of GRIN2's D1 tone over its excitation tone after feedback, also gone.
  Tests: `test_feedback_corrected_powers_are_what_comes_back_on`,
  `test_turn_on_only_switches_and_each_dds_keeps_its_binding`.
- Build: argument metadata and the suffixed superset (unconfigured Base).
  Prepare: `configure_execution('single_node', NODE)`, loads only that node's
  variables plus globals.
- `apply_manual_state()` calls `initialize_hardware(turn_off_dds_channels=True,
  turn_off_zotinos=True, initialize_spcms=False, reset_core=False)` so outputs
  established on the other node are preserved; then a second explicit
  OFF/0 V pass (`everything_off()`); then `turn_on_AOMs()` (requested DDS
  states, active-low gates, Node2 combined GRIN modes), `set_microwaves()`,
  `set_coils()`. The other
  node is never assigned a state, and AOMsCoils never calls
  `force_other_node_off`. Making a node safe means running that node's own
  experiment with controls unchecked. Caveat: `run()` always calls
  `run_feedback()`, and (when any waveplate box is ticked)
  `k10cr1_operations()`; both start with an unconditional `self.core.reset()`,
  which undercuts `reset_core=False` (O-11).
- Controls: FORT, cooling DP, A1-A6, GRIN1and2, the active-low optical gates
  (repump, pumping-repump, excitation0, GRIN1, GRIN2, microwave, and each
  node's dedicated D1 gate `ttl_D1_pumping_switch`, box `D1_pumping_switch_ON`,
  which on Node2 also switches the D1 DDS; Node1 D1 also opens the GRIN1 gate,
  `if self.GRIN1_AOM_switch_ON or self.D1_pumping_switch_ON`), microwave DDS,
  MW RF, coils, K10CR1,
  `run_laser_feedback`. DDS frequency/power are the authoritative `*_NodeX`
  variables, which `initialize_hardware()` writes through Base's resolver
  DDS bindings; `turn_on_AOMs()` only switches (except Node2's combined GRIN2
  mode, which sets `f/p_GRIN2_excitation_Node2`).
- Node1 D1 since 2026-10-19 (user): the D1 AOM runs from an external RF
  source switched by `ttl_D1_pumping_switch` (Node1 ttl15, active low), as on
  Node2. `D1_pumping_switch_ON` opens that gate (and GRIN1); Node1's Urukul D1
  channel (`dds_D1_pumping_DP`, urukul2_ch2) is no longer switched on by it
  and stays off. Test: `test_node1_d1_is_switched_by_its_external_source_ttl`.
  The Node1 sequences still switch D1 through that Urukul channel; with ttl15
  blocked at init they get no D1 light until they are changed.
- 2026-10-20 (user): the D1 box `D1_pumping_DP_ON` is renamed
  `D1_pumping_switch_ON` and moved from the DDS group to the "NodeX
  TTL-controlled AOM gates" group, in both files; what it switches is
  unchanged.
- Node2 combined `Node2_GRIN1_AOM_ON` / `Node2_GRIN2_AOM_ON` modes
  (`0a4ed7a`) are applied after the individual controls, as in standalone, so
  they are authoritative: their OFF branches override `GRIN1and2_ON`, the GRIN
  gates and the D1 DDS (`D1_pumping_switch_ON` on Node2 effectively drives only
  `ttl_D1_pumping_switch`). The GRIN2 mode reprograms `dds_D1_pumping_DP` to
  `f/p_GRIN2_excitation_Node2` before switching it on. `run_feedback()` re-applies
  them too after the stabilizer (2026-10-19); before that it re-applied only
  the individual controls, so a Node2 run with every stabilized AOM on could
  end with GRIN RF on or the GRIN1 gate open where standalone ends them off
  (`test_run_feedback_keeps_node2_combined_grin_state`).
- Coils (`71b4405`): one `disable_coils` checkbox (default False). Unchecked,
  the four coils go to the persistent `*_volts_MOT_NodeX` values read at
  prepare; checked, 0 V. Nothing is typed in the GUI or persisted.
  The channels are Base's `coil_channels` lists, as `Zotino.set_dac` needs
  (typed `list(int32)`; the old `COIL_CHANNELS` tuples had to be converted,
  `14d8bee`).
- Microwave/RF defaults OFF, is hard-gated and needs
  `yes_Im_sure_I_want_MW_or_RF_dds_ON`.
- K10CR1: thirteen booleans in the 780/852 groups drive this node's axes only
  (suffixed names); 852 target moves are gated on the 852 booleans (fixing the
  standalone copy-paste bug that gated them on the 780 ones).
- `run_laser_feedback` (`ff9e8ec`): only when A1-A6 and the cooling DP are
  all requested on, run the stabilizer (or `monitor()` when unticked), then
  re-apply the requested DDS and gate states without setting any stabilized
  DDS, so the corrected powers stay for the rest of the run; results persist
  node-suffixed.
- Excluded: old independent-Kasli networking TTLs, Rigol, other USB devices.
- Hardware: Node1 first ran 2026-09-09 with all controls off (RIDs
  38393-38401, after the three compile fixes of `14d8bee`); Node2 from
  2026-09-15; dozens of runs since. A documented per-channel polarity and
  coil-identity walk does not exist (O-19).

---

## 12. DRTIO and timing: constraints and measurements

### 12.1 Constraints

1. Access remote devices normally through unified database names; never
   invent destination arguments.
2. One canonical core; one unified database.
3. Generated gateware/database is the only source of encoded channels and
   destinations.
4. Wait for destination 1 before touching Node2 hardware
   (`initialize_hardware`, `force_other_node_off`).
5. Remote input, SPI, Sampler and TTL-timestamp operations are
   latency-sensitive (section 12.2).
6. RTIO underflow is a scheduling-slack problem, not a readiness or DRTIO
   latency problem (section 12.4).
7. Revalidate every DMA trace that involves satellite hardware (section
   12.7).
8. Inspect SED lane capacity before broad parallel RTIO (section 12.5).
9. AD9910 synchronization is calibrated for the final gateware and clocking
   (section 12.8).
10. Preserve edge counters on both nodes; never replace edge-counter
    channels with plain TTL inputs; include SPCM/counter checks in sanity
    tests.
11. Herald/mapping, photon-timestamp and cross-node optical timing tests are
    deferred until the two-node photon work starts.

### 12.2 Remote input and TTL-timestamp latency

(Cited by `tests/measure_drtio_ttl_skew.py`.) The retired PLAN.md (section
22, "Satellite Readiness") stated, verbatim:

> Remote:
>
> - SPI operations,
> - input operations,
> - Sampler operations,
> - TTL input/timestamp operations,
>
> must be treated as latency-sensitive where applicable.

The retired detail spec (section 20) condensed it to "Remote input/SPI
operations are latency-sensitive." This remains a binding rule.

Status 2026-10-19: the probe ran on 2026-10-18 (RIDs 39100 and 39102; the
numbers are in its RESULT block). The cables it used were removed on
2026-10-19 (section 2.7), so it cannot be rerun without re-cabling. Its leg 2
(Node1 ttl15 -> Node2 ttl26) was then removed from the probe, because ttl15
became Node1's D1 gate and the probe left it low (open); legs 1, 3, 4 and 5
keep their numbers. The 2026-10-14 status below is history.

Status 2026-10-14: the OUTPUT side is measured (section 12.3: no DRTIO
per-event cost). Remote INPUT latency and cross-node TTL skew have never been
measured. `tests/measure_drtio_ttl_skew.py` (`5383ee5`, 2026-10-07) exists but
its RESULT block is still a placeholder and no h5 from it exists, so
`t_Node2_rtio_offset_mu` is still the unmeasured default 0. The probe:

- T1 input skew (pulse master ttl12; difference of the satellite and master
  timestamps of one edge), T2 output skew (pulse satellite ttl31, the other
  way), T3 consistency, T4 counter aperture; each at four stability levels:
  within a kernel, across kernels, across `core.reset`, and (manually) across
  a link re-establish.
- Acceptance: |T1|, |T2| < 1 us and standard deviation < 100 ns.
- Cabling: the existing handshake coax plus two BNC T-adapters and two short
  coax: master ttl12 -> satellite ttl24 (existing) and -> master ttl3 (new);
  satellite ttl31 -> master ttl11 (existing) and -> satellite ttl27 (new). One
  output edge gives two timestamps, so the output path cancels and the
  difference is pure input-side skew. Cross-node cable length differences are
  < 3 ns (lab-measured), under half an RTIO period.
- Copy T2's median into `t_Node2_rtio_offset_mu` with its sign (section 7.5).
- It measures layer 1 (electrical) only; layer 2 needs photons.
- Decision logic (from the probe's docstring). Stability level 4 (across a
  link re-establish) is the risk: if levels 1-3 are sub-ns but level 4 moves
  by a clock period (8 ns), `t_Node2_rtio_offset_mu` is a per-session
  quantity, not a constant to bake into a `delay()`. That is irrelevant for
  loading gates (10 us-1 ms) but decisive for herald/photon timing. T4
  exists because the loaders gate EDGE COUNTERS, not timestamps; it must
  agree with T1. T3 must agree within jitter; a mismatch means a cable or
  channel assumption is wrong, not that DRTIO is bad.
- A failed acceptance gate means a structural fault (clocking, link, SED
  lanes). The probe said the two-node loader must not be written until this
  is understood; the loader and the alternating readout were written before
  the gate ran, so a failure reopens them.
- After the first run, fill in the RESULT block in
  `tests/measure_drtio_ttl_skew.py`, as `measure_drtio_event_cost.py` does.
- Free channels: master ttl2/ttl3/ttl10 are free counter-capable TTLInOut;
  satellite ttl27 (Node2's ttl11) is free.

### 12.3 Output event cost (measured)

RID 38429 (2026-09-10, `tests/measure_drtio_event_cost.py`; `ee0bdc2`,
`56f269a`): `spi_urukul1` (0x00002b, dest 0) against `spi_urukul4`
(0x01002b, dest 1), same peripheral, slot and local channel. Submission cost
713 mu/event (master) vs 711 mu/event (satellite). Minimum lead 1000 mu on
both (8/8 succeed at 1000 mu, 0/8 at 500 mu). `[DEST#1]` buffer space 128.
DRTIO costs nothing per event because events are buffered at the satellite.
Never explain an underflow by "it crosses DRTIO". The probe uses
`set_config_mu` (one output to the SPI config register, cs=0, no pin moves)
and avoids `cpld.sta_read()` (R21). It is a cabling-free control: keep it
separate from the skew probe, and re-run it after any gateware or link change.

### 12.4 Slack erosion

- Underflows on the feedback path are slack EROSION: host round trips
  (`set_dataset`, `append_to_dataset`, `print_async`, sampler reads) advance
  wall-clock time while the timeline cursor stands still, accumulating over
  loop iterations. The first such failure had -585808 mu of slack (about 500x
  the 1000 mu floor) on channel 65579 (`spi_urukul4`, dest 1); Node1 hit the
  same class on channel 39 (`ttl_urukul0_sw0`, dest 0). Not Node2-specific.
- Per-call cost through this master (2026-09-15, dashboard plus 9 applets
  subscribed): no RPC 0.30 us; blocking `set_dataset(broadcast)` 4.08 us;
  async `append_to_dataset` 2.73 us. Why standalone (with blocking appends)
  never hit this is unexplained; do not repeat the guess that the mixin made
  appends blocking.
- Fix pattern (R15): `core.break_realtime()` only where the preceding timing
  is irrelevant, outside any `which_node` guard, inside the shared helper.
  Per-site patching did not converge.
- Fix history: `edfb3c1` (2026-09-10, smooth FORT recycle loader feedback
  branch; its "crosses DRTIO" rationale was wrong and corrected by `56f269a`);
  `74ba1e5` (same entry in `load_MOT_and_FORT_until_atom` and `_recycle`; an
  audit of 54 MOT-set sites found 7 in loops); `6329a1a` (before
  `dds_FORT.sw.on()` in `run_feedback_and_record_FORT_MM_power`, -5.9 ms ch
  39); `dce7c80` (2026-09-15, before the post-loop `zotino0.set_dac` in three
  loaders, ch 33 `spi_zotino0`); `4549e4a` (2026-09-16, unconditional restore
  at the END of `run_feedback_and_record_FORT_MM_power`; the
  `record_power=False` path had none and died on ch 56 `dds_microwaves.sw`,
  -5.58 ms; covers about 44 of 51 host-heavy call sites); one restore after
  `tune_shims_for_atom_loading`; `5e8e82f` (2026-10-05, at the end of
  `end_measurement`, which makes 25 host RPCs and a blocking
  `measurements_progress` write; ch 36, -6.14 ms; covers all 47 call sites, so
  the hand-added `delay(5*ms)` lines are redundant but kept); `180f287`
  (2026-10-12, ch 65569 `spi_zotino1`); optimizer and two-node restores
  (sections 9.5, 11).
- These supersede the retired rule "add no speculative post-reset readiness
  wait or slack until hardware demonstrates a failure": hardware did
  demonstrate it. The rule still holds for readiness waits.

### 12.5 SED lanes

The SED has 8 FIFOs (lanes, `artiq/gateware/rtio/sed/__init__.py`). An event
with an equal timestamp takes a new lane; past 8 the event is DISCARDED with
`do_write` low. No exception: the only trace is "RTIO sequence error
involving channel 0x..." in the core log. Simultaneous events and timeline
rewinds (backward `at_mu`, each extra statement in `with parallel`) share
the budget. Evidence: RID 38676 (2026-10-07), the two-node loader switched 16
events (both repumps, both cooling DPs, all 12 fiber AOMs) on one timestamp;
events 9-16 were dropped (logged 0x01002f = `ttl_urukul4_sw1` =
`dds_AOM_A6_Node2`), so A3-A6 never switched on on either node and the MOT ran
on two of six beams (likely why the MOT check read about 12,700 c/s against
an 89,000 c/s threshold). The legacy code was safe only because two kernels
each issued 8.

Rule (lab owner, 2026-10-07): `delay(5*us)` after every FOUR events (2x FIFO
headroom). 5 us rather than 8 ns because each submission costs 0.3-0.7 us of
wall clock, so 5 us keeps slack positive. Space the edges, never inside gated
windows or physics-length intervals (`t_FORT_drop`, `t_PGC_after_loading`,
`t_alternating_RO_window`). Padding after every event is worse: inside
`at_mu`-anchored blocks the resync rewinds onto live timestamps. Two-node mode
is where this bites, because one timeline doubles every switching point.

### 12.6 `with parallel` policy

Use `with parallel` for stages whose durations differ per node
(`recooling_after_first_shot`, the loader's PGC block, the FORT drop;
`204e6c1`, `364e3c8`, 2026-10-07). It replaced hand-rolled
`now_mu()`/`t_longest`/`at_mu` resynchronization, which drifted once padding
was added and left the cursor on a live event's timestamp. Do NOT wrap
same-timestamp pairs (e.g. `ttl_repump_switch_Node1.off()` / `_Node2.off()`):
the gateware already fires them together, and `parallel` adds at_mu + now_mu
+ compare + select per statement (`artiq_ir_generator.py`), costing CPU and
slack. Do not make each node's six-AOM group one branch (8 events on the
entry timestamp, zero margin); interleaving N1/N2 is better. Wall clock
advances across branches even though the cursor resets. The same-timestamp
guard test models parallel blocks (every branch's first event sits on the
entry time) and forbids the `at_mu(now_mu())` no-op. Only
`two_node_alternating_shot`/`_windows` still place events with `at_mu`.

### 12.7 DMA

Existing DMA sequences were validated under standalone; a trace with
satellite devices is not assumed identical. Status: used on hardware without
a formal revalidation record: Node2 single-node DMA optical pumping
(MicrowaveScanOptimizer Node2 since 2026-09; `single_photon_experiment_9`
since 2026-09-17; `atom_photon_parity_11` 2026-10-12..14). The master DMA CRI
fans out to every destination, and parallel playback would force a
whole-trace SED rewind, hence one merged trace for both nodes (section 9.9).
Chopped readout makes DMA mandatory (section 18, O-31). Keep standalone DMA
implementations unless a master-satellite change is required.

### 12.8 Urukul SYNC window measurement

- `tests/measure_urukul_sync_windows.py` (`4f0b0d6`, 2026-09-10): 7 seeds
  across the 0-31 delay line, never writes EEPROM, per-channel init with
  `blind=True`, PLL lock reported separately.
- Read the WINDOW width (setup/hold margin), never the delay spread. Delays
  cluster about 11.8 taps apart, one SYSCLK period (1 ns at 1 GHz), giving a
  measured tap of about 85 ps (nominal 75 ps; full line about 2.3 ns). Every
  repeat catches the same SYNC_IN edge (62.5 MHz). `sync_delay_seed` is a TIME
  position 0-31 taps (15 is just the midpoint); `io_update_delay` is 0-3
  SYSCLK cycles; `tune_sync_delay` searches only ±6 taps (`search_span=13`,
  FIXME Urukul#16), so the tester's chosen value depends on its seed.
- Verdict rule (`a68d1f4`, 2026-09-15): 0 zero-window seeds = healthy; 1 =
  inconclusive, re-run; 2+ = MARGINAL. Rejected alternatives: "zeros >= half
  the seeds" missed two bad channels on RID 38431; "any zero" flagged two
  healthy ones on RID 38432 (urukul0_ch0 and the new urukul4_ch1). Never
  condemn or clear a card on one run; if a future run disagrees with the
  hardware picture, suspect the rule first.
- A `sinara_tester` pass is not proof of margin (on 2026-09-10 it passed 3 of
  the 5 channels the sweep flagged). It writes EEPROM as it goes and aborts on
  the first failure, which is why "re-run until it passes" used to work.
- History: RID 38430 (2026-09-10) found 20 of 24 channels with window 1-2 from
  every seed and `urukul4_ch2` plus all of `urukul5` with window 0 from 4-5 of
  7 seeds; RID 38431 (no power cycle) flagged the same five. `urukul3` in the
  same crate was perfect, exonerating crate, EEM and clock. Five other
  hypotheses (including DRTIO latency, clock distribution, supply droop) were
  disproved. The fault dated from a 2025 incident (a GVS over
  `t_microwave_pulse` entered in 335 MHz units, read as about 3877 days).
  Resolved by the 2026-09-15 card swap (section 2.10); RIDs 38432, 38436,
  38438, 38439.

### 12.9 Edge counters

`EdgeCounter` entries exist on both destinations (`ttl0_counter`,
`ttl16_counter` ...). `Base.EDGE_COUNTERS` binds them per node.
`artiq-master/drtio_passive_check.py` completed a `ttl0_counter` gate (master)
and a `ttl16_counter` gate (destination 1) on 2026-09-09. `gate_rising`
emits one input event per gate; `ififo_depth` is 64 per channel.

---

## 13. Host-controlled devices

- All Node1 and Node2 USB/host-controlled hardware connects to the MASTER-side
  host; this is final (no distributed USB design). Host devices are not in the
  per-node RTIO translation.
- K10CR1: one `k10cr1_ndsp` controller entry in the unified database (host
  `::1`, port 8080; appended 2026-08-31, `3b4dc43`, `a88f61e`) owns all eight
  rotators through one USB hub. `sn_list`: Node1 [55420984, 55422044,
  55000740, 55000759], Node2 [55422744, 55000915, 55105704, 55105694], with
  nicknames in the same order: `780_QWP_Node1`, `780_HWP_Node1`,
  `852_QWP_Node1`, `852_HWP_Node1`, `780_QWP_Node2`, `780_HWP_Node2`,
  `852_QWP_Node2`, `852_HWP_Node2`. `ndsp/k10cr1/launcher_multi_rotor.py`
  reads `sn_list`/`nickname_list` from that entry, so the launcher, driver and
  `subroutines/k10cr1_functions.py` are unchanged.
- Single-node proxy: Base publishes `k10cr1_ndsp` as `_NodeAxisK10CR1Proxy`.
  Bare legacy axis names (`780_HWP`, `780_QWP`, `852_HWP`, `852_QWP`) resolve
  to the selected node's nicknames; suffixed names pass through; the real
  sipyco client (which connects on creation) is fetched only on first call,
  so experiments that never rotate never connect. Connection failure raises
  `RuntimeError`. Master-satellite files use explicit suffixed names.
- No `best_effort` on the rotators (silent `None` returns would drop motion
  commands). Never propose per-node controller processes.
- After editing the deployed database, `artiq_client scan-devices`.
- In use: Node2 scans over `target_780_HWP` on 2026-10-13/14.
- Rigol and other host instruments remain deferred and must not block other
  validation; `Rigol_*` variables stay declared per node (historical).

---

## 14. Applets and dashboard

### 14.1 CCB auto-applets

`applets_master_satellite.py` (since `a15de11`, 2026-09-10) holds one table of
applet specs; the running experiment creates its set with
`ccb.issue('create_applet', ...)`. Dataset names are written legacy-unsuffixed
and resolved through `resolve_experiment_variable_target` then
`resolve_result_dataset_name`, so an applet cannot drift from what the run
writes. Applets are never deleted, only unchecked, so geometry survives. GVS
argument `create_applets` (default True) opts out; applet failure only logs a
warning.

REQUIRED: the applet dock's GLOBAL CCB policy must be "Create and
enable/disable applets" (`ccbp_global='enable'`). Under "Create applets"
entries appear unchecked and disabling is ignored. A per-group policy
overrides the global one and makes one node work while the other silently
does not.

### 14.2 Groups and specs

Top-level groups: `SINGLE_NODE_APPLET_GROUP = "Single node"` (one group for
whichever node runs, since `78a1b19`, 2026-10-13; only option variables such as
`single_atom_threshold_NodeX`, `t_SPCM_*_NodeX`, `t_atom_check_time_NodeX`
follow the node, and a node switch restarts only the applets whose command
changed); `TWO_NODE_APPLET_GROUP = "TwoNodes"`; `SHARED_APPLET_GROUP = "Both
nodes"` (node-pinned per-node applets, never retired). Functional subgroups:
"GVS and Cycler", "Optimization", "Monitor in each measurement", "Single
Photon Experiment", "K10CR1".

Spec sets: `APPLET_SPECS` (progress bar, retention and loading AllSPCMs,
RO1/RO2 histograms, Atom loading time); `GVS_APPLET_SPECS` (atom-check loading
histogram vs `single_atom_threshold_for_loading`); `SCAN_APPLET_SPECS` (Time
without atom); `TWO_NODE_APPLET_SPECS` (retention on the joint
`two_atom_threshold`, cutoff = `two_atom_threshold` x `t_SPCM_first_shot_Node1`,
420 counts at the live 42000 against 140 / 170 single-node for Node1 / Node2;
uses
`t_SPCM_first_shot_Node1` because the bare name resolves to a stale standalone
dataset); `OPTIONAL_APPLET_SPECS` (e.g. `n_excitation_cycles`);
`SPCM_MONITOR_APPLET_SPECS` (five specs: the four per-SPCM count rates plus
"All SPCMs count rate", MonitorSPCMinApplet only); `TWO_NODE_ALTERNATING_APPLET_SPECS` (gated by the empty
`TWO_NODE_ALTERNATING_EXPERIMENTS`); `TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS`
(4 histograms + 2 retention plots for the alternating retention experiment,
which retires the joint trio); `K10CR1_APPLET_SPECS` (waveplate trajectory);
`SHARED_APPLET_SPECS` (per node: Microwaves Health Check, feedback RF, FORT APD
normalized setpoint, FORT MM normalized to ref, FORT MM over APD normalized
via `plot_FORT_MM_norm_to_power.py`).

`MAGNETOMETER_APPLET_SPECS` (`bfca64d`, O-35): three applets under
"Monitor in each measurement", "Magnetometer (Zero)", "Magnetometer (OP)" and
"Magnetometer (MOT)", one per coil setting the readers measure, each plotting
`Magnetometer_<phase>_X/Y/Z` of the running node with `pts = n_measurements`
and the legend `magnetometer_axis_labels` (section 6.9), via
`plot_xy_multichannel.py`. MANUAL ONLY, by the user's choice
(2026-10-14): no experiment creates them, whatever
`monitor_magnetometer_in_end_measurement` is, and none closes them, so an
applet opened by hand stays open. Open them from the dashboard with the
command `build_applet_command` gives (the docstring above the tuple shows the
Node1 Zero one), or request the tuple explicitly from an experiment like any
opt-in set. With `monitors_for_atom_loading` also on, the loaders' samples
interleave with the end-of-measurement ones. No unit in the titles: Node1
scales 350 mG/V, Node2 `mGauss_per_Volt = 1` (volts).

Entry points: `create_applets_for` (single_node only, with retirement),
`create_two_node_applets_for` (two_nodes only), `create_shared_applets_for`
(either mode, no retirement). Callers: single GVS (APPLET + GVS, shared +
SCAN); two GVS (TWO_NODE + variants, shared + SCAN); MicrowaveScanOptimizer
(APPLET, shared + SCAN); FORT optimizer (shared + K10CR1); SamplerMOT
(shared); MonitorSPCMinApplet (APPLET + SPCM_MONITOR).

### 14.3 Retirement rules

A leftover applet from the other mode keeps plotting live data under the
wrong criterion and looks plausible (a single-node copy reads high during a
two-node run; a two-node copy reads about 0 during single-node). Retirement is
per applet, not per group (`204e6c1`): the conflicting set is the
intersection of the two spec sets, and `exclusive_applet_titles` derives the
single-experiment applets. Retirement never worked until `c1ddbd2`
(2026-10-12): `ccb_disable_applet` resolves its target with
`locate_applet(create_groups=False)`, which matches only a DIRECT child, so a
bare top-level group silently resolved to nothing. Always pass the FULL nested
path; a test compares every retired path against the spec-derived group. Old
Node1/Node2 groups are closed per applet by `_retire_legacy_node_applets`
(never `disable_applet_group`). `MAGNETOMETER_APPLET_SPECS` is deliberately
NOT in `EXCLUSIVE_APPLET_SPEC_SETS`: the user opens those applets by hand, and
retiring them would close a hand-opened applet on the next run.

### 14.4 Dashboard operations

- When `artiq_dashboard` crashes, its applet subprocesses survive; each
  restart adds a set and every orphan subscribes to port 3250. 80 orphans
  caused an endless "qasync Task exception was never retrieved" crash loop
  (2026-09-09/10). Before restarting, count python processes matching
  `applets|artiq\.applets`; if no dashboard runs, all are orphans and may be
  stopped.
- Dashboard state: `C:\Users\QC\AppData\Local\m-labs\artiq\7\`
  (`artiq_dashboard_..1_3250.pyon`; a `.bak-2026-09-10-crashloop` copy
  exists). Backups and a README (deleting a corrupt state file when the
  dashboard "flashes open and quits") are in OneDrive
  `Desktop/Eunji/Backup - ARTIQ Dashboard/`.
- "OSError WinError 6" from sipyco `pipe_ipc _autoclose` is expected, once per
  stopped applet; a burst after a node switch is normal.
- InfluxDB/Grafana telemetry has no code and stays future work (R24).

---

## 15. Analysis

### 15.1 Layout and naming

Split 2026-09-17 (`2a5796a`, `9af6737`, `2972733`): `Analysis/standalone/`
holds the pre-changeover tree; `Analysis/master_satellite/` holds
post-changeover notebooks; general material (`h5_data_utilities.py`,
LabMagnetometer, physics scratchpad, two plotting notebooks) stays at the
`Analysis/` top level. `Analysis/standalone/h5_data_utilities.py` is a
forwarding stub that loads the real module by file path under a private name
(a plain import would be a silent circular self-import).
`h5_data_utilities_node1.py` was merged and removed. Notebooks are named on
the node-count axis only (`2026_Single_Node_*`, future `Two_Nodes_*`), with
no folder-redundant prefix such as `MS_`; keep existing names when porting.
Current notebooks: `2026_Monitor_Magnetometer_Analysis` (no node prefix by
design), `2026_Single_Node_Atom_Loading`, `_Microwaves`,
`_Parity_Analysis_AllSPCM`, `_Single_Photon_Analysis`, and
`2026_Two_Nodes_Atom_Loading` (2026-10-19), ported from
`standalone/Two photon analysis/two_nodes_atom_loading.ipynb`. That
folder's two-photon notebook is not ported: there are no master-satellite
two-photon runs yet.

### 15.2 Master-satellite loaders (`Analysis/h5_data_utilities.py`)

- `ms_node_of(f, filename)`: node detection order: expid `selected_node` /
  `which_node` (EXACT match; 2403 standalone runs carry `alice`/`bob`), then
  the filename tag, then dataset key suffixes; raises if both suffixes are
  present.
- `ms_files(date_filters, node, name_filters)` returns `(file, node)` pairs;
  `MASTER_SATELLITE_NAME_FILTERS` cover GVS, MicrowaveScanOptimizer and
  AOMsCoils (not FORT, AtomLoading, MonitorMagnetometer, SamplerMOT,
  MonitorSPCMinApplet or HealthCheck files).
- `ms_archive_and_datasets_to_locals(f, parent_locals, quiet, node)`: strips
  the running node's suffix so notebooks use unsuffixed names; on a
  bare-plus-suffixed collision keeps the suffixed name and prints a note even
  when quiet (`d989540`); injects expid arguments; aliases
  `override_ExperimentVariables` to the running node's dictionary; re-applies
  overrides the file did not record (`_ms_unrecorded_overrides`, `ef55cae`:
  re-evaluates the submitted override text against the archive with ARTIQ
  unit names rebuilt, skipping names already recorded, scanned variables and
  other-node names; checked on 238 October files, only overridden names
  change); exposes `overrides_applied`; sets `node`. Two-node files are refused
  by design.
- Two-node files (2026-10-19), for `2026_Two_Nodes_*`: `ms_two_nodes_files`
  returns ONE file per RID (the plain file ARTIQ writes at the end of the
  run, else the `_TwoNodes_..._scan_over_...` copy the GVS rewrites after
  every scan point, read from its name and never opened), optionally
  filtered by `experiment_function`.
  `ms_two_nodes_archive_and_datasets_to_locals` keeps every suffix, applies
  unrecorded overrides under the name the run resolved them to, rebuilds
  the bare joint scalars (`MS_TWO_NODES_JOINT_SCALARS`, a copy of
  `TWO_NODE_AGREEING_SCALARS`) from their agreeing `_Node1`/`_Node2` pair,
  strips the seed from the four `SPCM*_OtherNode_RO*` series too, and
  injects expid arguments. The archived `*_two_node` timing globals are
  not what the sequence read (RID 39108: `t_atom_check_time_two_node`
  0.01, gate 0.005).
- `ms_get(source, key, default)`: datasets, then archive, then expid, with
  suffix resolution, returning DECODED values. Never wrap it in `str_from_h5`
  (`str_from_h5('dummy_variable')` gives `'mmy_variabl'`). GVS does not
  archive `scan_variable1_name`; `n_measurements` is in the archive.
  `ms_scan_label`, `ms_overrides`, `ms_print_overrides` (`8db3f6b`),
  `ms_show_node` (`6a89075`).
- A sweep of 226 runs since 2026-10-01 (`d989540`) found 171 loading
  cleanly, 7 newly loadable (`_Node2`, RID 38725 onward) and 46 two-node
  files raising by design. These sum to 224; the commit body does not
  account for the other two (possibly the unreadable RIDs 38816/38826).

### 15.3 History of "oops... something wrong with <file>"

(1) Master-satellite GVS files lack `scan_variable1_name` in datasets, and a
bare except reported the missing NAME as a corrupt FILE; fixed with `ms_get`.
(2) A node-selection knob `ms_filter_by_node` was added and reverted
(`0b59304`); `ms_files` still accepts `node=` but no notebook uses it. (3) A
live kernel served a pre-commit `h5_data_utilities` from `sys.modules`;
restarting the kernel fixed it. Disproved: Python/h5py version differences,
h5py file-mode conflicts. All bare excepts in the master_satellite notebooks
now print the exception (`4aed4c1`).

### 15.4 Known analysis issues

- `ms_node_of` can tag a two-node file as Node1: step 2 matches `_Node1_`
  anywhere in the filename, or a name ending `_Node1.h5`, before the
  both-suffix check, and the two-node GVS default scan variable is
  `f_FORT_Node1`, so its named copy ends `_scan_over_f_FORT_Node1.h5`.
  CONFIRMED on RID 38766's file with a renamed path; no current file is
  affected because two-node runs scan `dummy_variable` (O-27).
- `get_files_by_criteria` leaks one unclosed `h5py.File` per candidate and
  matches with `any()`, not `all()` (two-element filter: any = 71, all = 0);
  shared with standalone and left alone; the comment was corrected only in
  the master_satellite copies.
- `h5_archive_and_datasets_to_locals` strips the leading `[0]` seed from
  `SPCM0/1_RO*`, `AllSPCMs_RO*` and `BothSPCMs_RO*` but not from the four
  `SPCM*_OtherNode_RO*` series, so every notebook that histograms those
  per iteration (standalone and `2026_Single_Node_*`) pairs each shot with
  the previous one. Only the two-node loader strips them; the plain
  loader is left as it is.
- The `2026_Single_Node_*` name filter `GeneralVariableScan_master_satellite`
  also matches `..._two_nodes` files, which `ms_archive_and_datasets_to_locals`
  then refuses.
- In the single-photon notebook the decay-fit cells (48, 52) and efficiency
  cells (26, 53) still index `[0]` into per-iteration arrays, possibly on
  purpose.

### 15.5 Jupyter live-session caveats

Notebooks are usually open in a live Jupyter session. Editing one on disk
collides with the open tab, and saving from the stale tab overwrites the edit
(a moving `execution_count` is the tell): say explicitly that the file must be
closed and reopened. `from h5_data_utilities import *` is served from
`sys.modules`, so new module functions are invisible until a kernel restart.
Put logic in `Analysis/h5_data_utilities.py`. Check `git status` for
notebooks modified on disk. Never leave a bare except around unfamiliar code.
If EVERY file fails identically, the fault is in the cell's namespace (a stale
import), not the data ("found 33 files, 0 labelled, 33 oops", 2026-09-17).

---

## 16. Testing and verification

### 16.1 Unit tests (hardware-free)

- `unittest` by project convention (pytest is not used; it is installed in
  the `artiq7` env and Anaconda base, but not in `C:/Python39`). From the
  repository root with an artiq-free interpreter: `python -m unittest
  discover -s tests -p "test_*.py"`. Here `python` on the Git Bash PATH is
  `C:/Python39/python.exe` (3.9.5, no artiq); `C:/ProgramData/Anaconda3/python.exe`
  (Anaconda base, also artiq-free; an older memory recipe added
  `PYTHONPATH=<repo>;<repo>/tests`, which is unnecessary) gives the same
  result; add `-B` (or `PYTHONDONTWRITEBYTECODE=1`) to avoid
  `__pycache__`. Always use discover, not single modules (the shared artiq
  stub accumulates in `sys.modules`). Under the `artiq7` interpreter the real
  `EnvExperiment` makes many tests error; use an artiq-free interpreter.
- 2026-10-16: 154 tests OK (150 at `ef55cae`, +2 in `619d3ec`, +2 in
  `bfca64d`). Per file:
  `test_aoms_coils_master_satellite` 18, `test_atom_loading_optimizer_master_satellite`
  22, `test_base_experiment_master_satellite` 32,
  `test_device_aliases_master_satellite` 6, `test_experiment_functions_two_nodes`
  12, `test_experiment_variables_master_satellite` 10,
  `test_fort_polarization_optimizer_master_satellite` 6,
  `test_general_variable_scan_master_satellite` 31,
  `test_microwave_scan_optimizer_master_satellite` 17.
  `tests/test_ccb_applet.py` is a dashboard EnvExperiment, not a unit test.
- Coverage includes: every public build with `None` arguments and empty
  storage; strict missing-dataset failures; GUI `n_measurements` precedence;
  underflow retry; registry contents; per-node override merge; override
  record ordering; `n_measurements` broadcast+persist; SED same-timestamp
  guard; single-node wiring/amplitudes unchanged by two-node code;
  `test_global_initialization_has_no_node_side_effect` (a hardcoded global
  inventory, on purpose, so adding a global is a conscious act).
- Count history: 56 (2026-08-30), 58, 64-87 (09-02..09-17), 95-97 (10-01),
  119 (10-02), 147 (10-12/13), 148, 150 (`ef55cae`).

### 16.2 Offline kernel compile sweep

`C:/Networking Experiment/artiq codes/artiq-master/compile_check_master_satellite.py`
lives outside the repository (not version-controlled, so the repository scan
never imports it). It compiles every master-satellite kernel without hardware:
patches `Core._run_compiled` to a no-op and wraps `Core.compile` to record
every kernel; copies `dataset_db.pyon` to a temp file; builds with
cwd=`artiq-master` and chdirs to a temp results dir before prepare/run like
the worker; uses DummyScheduler/DummyCCB and `device_db.py`.

- Run after ANY master-satellite edit. Activate `artiq7` (or put
  `<artiq7 env>/Library/bin` on PATH): without `ld.lld` every case FAILs at
  link with "(0 kernels compiled)"; the script exits 2 if `ld.lld` is missing.
  `python compile_check_master_satellite.py [label-substring]` from
  `artiq-master`.
- Statuses: PASS (prepare+run OK); COMPILED (all kernels compiled, host error
  afterwards because kernels never execute, e.g. `retention_array[-1]`);
  FAIL.
- A green result means nothing unless the target kernel appears in the
  printed compiled-kernel list: flag-guarded kernels compile only in a case
  that switches the flag on (a masked failure once showed 2 kernels where 3
  were expected).
- `_seed_missing_variables` fills absent variables from code defaults in the
  COPY, so a green sweep cannot detect a missing persistent dataset. Read the
  end-of-run "ACTION NEEDED" block; if it is non-empty the change is not safe
  to deploy without seeding.
- 41 cases: AOMsCoils Node1 x4, Node2 x4; GVS single-node x9 (TTL pulse
  N1/N2, atom loading, single photon 9 N2, microwave Rabi 2 confirm ON,
  map01_map11 confirm ON, parity 11 N2 confirm ON, atom loading 2 N1/N2); GVS
  two-node x5 (namespace sanity, atom loading symmetric / confirm-by-first-shot
  / asymmetric per-node branches, alternating retention); MicrowaveScanOptimizer
  x3; FORT optimizer x2; AtomLoadingOptimizer x2; DRTIO event cost; Urukul
  sync sweep; DRTIO TTL skew probe; SamplerMOT x2; MonitorSPCMinApplet;
  MonitorMagnetometer x4; MOTBeamPowerCheck x2. No case covers parity 2-7 or most of the 31 wrapped
  functions (O-14).
- Last full runs: 2026-10-14 about 20:05-20:10, and again read-only about
  20:33, on the working tree including the optimizer fixes later committed
  as `619d3ec`:
  "All cases passed", 31 PASS, 6 COMPILED (MicrowaveScanOptimizer x3,
  AtomLoadingOptimizer x2, DRTIO event cost), 0 FAIL, nothing seeded. About
  5 minutes. `subroutines/experiment_functions.py` was edited after both
  (O-34, O-35). The five microwave cases (needle `icrowave`) were re-run
  after the magnetometer changes, about 21:25: both GVS microwave cases PASS,
  the three MicrowaveScanOptimizer cases COMPILED with unchanged kernel
  counts (3/3/25). Re-run the full sweep before committing.
- Stale text in it: the docstring ("87 hardware-free tests"), a CASES comment
  ("Only two functions are registered in two-node mode"; there are three) and
  the comment calling `tune_852_waveplates_to_target_in_experiment` "NOT a
  declared master-satellite variable" (section 7.7).

### 16.3 Positive controls and replay probes

- Pair every detector with a control that must fire: the reverse-transform
  proof was first shown to fail on three deliberate corruptions; `urukul3`
  served as the healthy in-crate control for the SYNC fault; the event-cost
  probe is the cabling-free control for the skew probe.
- Replay real data through real code before trusting a fix: RID 38903 was
  replayed through ARTIQ's `DatasetManager` to verify override recording; the
  `n_measurements` deletion was reproduced the same way (2026-10-14); the
  analysis override fallback was checked on 238 October files.
- For throughput, compute meas-per-load from the h5 (section 9.8).
- Suggested but not done: measure a kernel-generated TTL frequency on a
  counter referenced to the lab clock (locked reads exact; a free-running
  crystal reads tens of ppm off) (O-21).

### 16.4 Hardware diagnostics kept in `tests/`

`measure_drtio_event_cost.py` (cabling-free control), `measure_drtio_ttl_skew.py`
(needs T-adapters; unrun), `measure_urukul_sync_windows.py`,
`bringup_test_rtio_clock.py` (measures the RTIO clock against a
reference-locked TTL through an EdgeCounter; to be deleted once clock
verification is no longer needed), `MonitorSPCMinApplet.py`, and
`MOTBeamPowerCheck_master_satellite.py` (2026-10-16): one node's six MOT beam
photodiodes, one beam at a time exactly as feedback measures its series
channels, at the `p_AOM_A*_<node>` powers re-read at run start (feedback
leaves them there; nothing is written back). Sampler card and channel come from
the node's `feedback_channels.json` via `prepare_laser_stabilizer`. Raw volts
per round in `MOT_beam_PD_A1..A6_<node>`, applied dB in
`MOT_beam_p_AOM_dB_<node>`, one applet "MOT beam PDs <node>" (group "Both
nodes", only with `show_applet`), and a per-beam summary in `analyze()`. Each
point averages `aom_feedback_averages_<node>` readings; the rounds are not
averaged. Ends with cooling and the six MOT AOMs on, as feedback leaves them.
Compile-checker cases for both nodes.

### 16.5 Hardware sanity principle

Before any physics on new or repaired hardware, check independently and with
small tests: canonical core; destination readiness; DDS, CPLDs, TTLs, edge
counters/SPCM inputs, Zotino, Sampler and USB devices per node. Never debug
infrastructure with a large physics experiment. On failure report mode/node,
exception, device/operation, reset stage, destination readiness, resolved
alias, and whether slack is relevant, before proposing a targeted fix.

---

## 17. Status (2026-10-14)

### 17.1 Software timeline

| Date | Commits | Change |
|---|---|---|
| 2026-08-26 | `e86046f` | Legacy two_nodes node mode removed. |
| 2026-08-27/28 | `cc399dc`, `1f83f11` | Physical mapping; resolver; Base. |
| 2026-08-30 | `d44d528`, `a3eb7fe`, `f06f9be`, `078dc0e`, `522d878`, `8963734` | Master-satellite GVS; examination-safe build; combined AOMsCoils; plan documents. 56 tests. |
| 2026-08-31 | `025c027`, `83ae355`, `a7d8dff`, `71b4405`, `3b4dc43`, `ff9e8ec`, `0a4ed7a`, `614549e`, `9904f9b`, `1c428e4`, `c18a6f2`, `3d0ec82` | GVS per-mode split; EV rename; AOMsCoils per node, `disable_coils`, K10CR1, laser feedback, Node2 GRIN modes; FORT optimizer; underflow retry folded into GVS (both stacks); repository-scan safety. |
| 2026-09-01/02 | `1277477`, `bb72047`, `a88f61e` | Microwave optimizer; refine fix; K10CR1 proxy. |
| 2026-09-03 | `1a6b953`, `331b669`, `1a1bb6c`, `3257b1c`, `4e636da` | `standalone/` folder; `t_OP_coil_relax_time`; parity 11. |
| 2026-09-09 | `14d8bee`, `6f4001b` | First hardware run (AOMsCoils Node1); compile fixes; compile checker; GVS stabilizer/result-state requirement. |
| 2026-09-10 | `6329a1a`, `45e0848`, `a15de11`, `1382a49`, `ee0bdc2`, `56f269a`, `4f0b0d6`, `68d63ee`, `da1eca7`, `0ae281a`, `edfb3c1`, `74ba1e5` | Dataset seeding surface; CCB applets; deployed DB tracked; event-cost and SYNC probes; debug-note convention; slack fixes. |
| 2026-09-15/16 | `a68d1f4`, `dce7c80`, `bb5d687`, `a5e685b`, `4549e4a`, `0fea421` | SYNC verdict rule; diagnostics to `tests/`; SamplerMOT port; feedback-helper slack restore; MonitorSPCMinApplet port. |
| 2026-09-17 | `21f6392`, `5a0bc76`, `7ea1476`, `2a5796a`, `9af6737`, `2972733`, `0b59304`, `4aed4c1` | `force_other_node_off`; shared applets; global `two_atom_threshold`; `result_name_tag` and named saves; MW centring; Analysis split and loaders. |
| 2026-10-01/02 | `5c1f508`, `73a11e3`, `ecc83ca`, `3d99643` | Per-node override dictionaries; AtomLoadingOptimizer; SamplerMOT split. |
| 2026-10-05 | `82f5f6f`, `d41b9d6`, `0ff183e`, `369f805`, `5e8e82f`, `e79b0fd` | `n_measurements` persisted; MonitorMagnetometer; `end_measurement` slack restore. |
| 2026-10-06/07 | `6a89075`, `e768d50`, `af8b343`, `204e6c1`, `364e3c8`, `5383ee5`, `fd10950`, `acd63f5`, `1c46286`, `3a90e2d` | Bare override names; Base and feedback in two_nodes; cross-node timing rename; two-node atom loading; `with parallel`; skew probe; alternating diagnostic (reverted next day). |
| 2026-10-08/09 | `ae93327`, `30d7be1`, `f9a587a`, `19c289f`, `12c5314`, `50f9049`, `b6b4de1`, `c6546a1`, `d989540` | Revert; alternating rename and gate fix; parallel-feedback design note; alternating loader disabled. |
| 2026-10-11 | `94d0ef6` | PR #68 merged into `main`. |
| 2026-10-12 | `c1ddbd2`, `0e54c84`, `180f287`, `38cc8c5`, `2bbe4fb` | Alternating retention; one-trace OP for both nodes; advance flag; dead OP variables removed; upgrade notes. |
| 2026-10-13 | `43455f7`, `2a3592c`, `1cb0cc1`, `2fb2fbf`, `0edd7ce`, `78a1b19`, `f36a225`, `34bff1e`, `b770aa0`, `5e10655` | Confirm-by-first-shot rollout (31 sites); removals; `max_loading_rounds`; one "Single node" applet group. |
| 2026-10-14 | `5b2f618`, `ef55cae` | Analysis notebooks (`Atom_Loading`, `Parity_Analysis_AllSPCM`, `Single_Photon_Analysis`, `singe_photon_analysis_20260622`); override recording and loader fallback. The MW optimizer fixes, magnetometer work and plan consolidation were written this day and committed on 2026-10-16. |
| 2026-10-16 | `619d3ec`, `bfca64d`, this file's commit | MW optimizer: pre-scan calibration restored before the fit check, overrides recorded before any broadcast; magnetometer monitoring (optimizer seeding, `break_realtime` in the readers, the user's 2 ms coil-settle delays, manual-only applets, axis labels); the four plan files replaced by these two, citations repointed. |

An unmerged side branch `origin/20260831_testing_PGC_in_bias_B` holds
`c8607cd` (2026-09-03, "testing 1D pgc with bias field", +325 lines in
`experiment_functions.py`), on neither `main` nor this branch.

### 17.2 Hardware bring-up history

- 2026-09-07: Node1 flashed with `drtio_master` gateware. "RTIO PLL failed to
  lock" was real (`clock_alive_check.py`: RTIO counter frozen at 0 over 1e8
  reads). Root cause: the v1.0 board's GTX/MGTREFCLK path (Si5324 -> clk_gtp
  -> CPLL), used only by DRTIO gateware (standalone uses PLLE2_ADV).
  Exonerated: the gateware, JSONs and nix build (the same image locked on
  another v1.0 board at .131, and Node2 locked as master), the reference and
  `rtio_clock` (fails on `int_125` too). Node1 also failed as satellite ("RTIO
  clock failed"). 2026-09-08: the standalone image locked on Node1.
- 2026-09-09: Node1 board swapped for a v1.1 board with the unchanged master
  image; link up; `drtio_passive_check.py` (artiq_run, no outputs) completed a
  master `ttl0_counter` gate and a satellite `ttl16_counter` gate; first
  AOMsCoils Node1 run with all controls off (RIDs 38393-38401).
- 2026-09-10: first GVS physics, `atom_loading_2_experiment` on Node1 and
  Node2 (RID 38402 onward); event-cost probe (RID 38429); SYNC sweeps (RIDs
  38430, 38431); dashboard orphan-applet crash loop resolved.
- 2026-09-15: `urukul4/5` replaced; tester and sweep clean (RID 38432);
  AOMsCoils Node2 runs.
- 2026-09-16 onward: MicrowaveScanOptimizer on both nodes; 2026-09-17 onward:
  `single_photon_experiment_9` on Node2.
- 2026-10-01/02: AtomLoadingOptimizer; SamplerMOT.
- 2026-10-06/07: namespace sanity and two-node atom loading on hardware.
- 2026-10-08..12: alternating shot and retention; MonitorSPCMinApplet.
- 2026-10-12..14: parity 11 Node2, sp9 Node2, MW optimizer with the confirm
  loop.

### 17.3 Retired migration plan, phase by phase

| Phase (retired PLAN.md) | Status 2026-10-14 |
|---|---|
| 0 Preserve standalone | Done; standalone files under `standalone/`. |
| 1 Physical mapping | Done (`cc399dc`). |
| 2 Resolver | Done (`1f83f11`). |
| 3 Base | Done; two_nodes completed in `e768d50`. |
| 4 ExperimentVariables | Done (renamed files); the sequential/parallel selection exists only as the unread, reserved `parallel_AOM_feedback`. |
| 5 Master-satellite GVS | Done; split per mode. |
| 6 Feedback | Sequential single- and two-node feedback on hardware; SED inspected; PARALLEL not implemented. |
| 7 Compatibility tests | Master-satellite single-node on both nodes done; old GVS on standalone gateware never re-run (O-20). |
| 8 Timing validation | Event cost measured; AD9910 sync done; DMA in use without sign-off; remote input latency unmeasured. |
| 9 New physics | In progress: two-node loading and alternating retention; two-node photon/herald/atom-atom not started. |

The retired first-validation sequence was followed through manual control;
its magnetometer steps (single Node1 then Node2, `n_measurements = 1`) were
skipped in practice: no master-satellite magnetometer run exists (O-18). The
retired gateware gate was completed in substance (images archived, rollback
image located, link verified), except the full SD-card images and single
checksummed release bundle, which are not recorded.

### 17.4 What is in hardware use

Single-node on both nodes: GVS (`atom_loading_2`, microwave, single-photon
sp9, parity 11; Node2 `target_780_HWP` waveplate scans on 2026-10-13/14, the
latest of them, RIDs 38952-38960, recording overrides), MicrowaveScanOptimizer
(including the fixes of section 11.1, and with the magnetometer monitor on
by override on Node1: RIDs 38969, 38972-38973, 38975-38977, 38980, 38982),
AtomLoadingOptimizer,
SamplerMOT, AOMsCoils, FORT optimizer (both nodes; last run RID 38772 on
2026-10-12), MonitorSPCMinApplet. Two-node: namespace sanity, joint atom
loading (confirm flag: one 83-measurement run, RID 38785; no two-node run
since), alternating retention. Not yet run: MonitorMagnetometer
master-satellite, the TTL skew probe, any two-node photon experiment.

---

## 18. Open issues, known faults and next steps

Each item has an ID for citation. "Evidence" is where to verify it.

### 18.1 Physics and sequence faults

- **O-1 Alternating loader loads no atoms.**
  `load_until_atom_in_both_nodes_alternating_shot_recycle` reports success
  but readout sees none (RID 38718 vs RID 38710). Disabled. Suspicion untested
  (section 9.6).
- **O-2 Remote input latency / cross-node TTL skew unmeasured.** Needs the two
  T-adapters and the skew probe (section 12.2); `t_Node2_rtio_offset_mu` is the
  default 0 and the alternating readout consumes it.
- **O-3 `t_Node2_excitation_delay_mu` unmeasured** (189 default); needs a
  photon-coincidence (HOM) measurement; two-node single-photon/herald work is
  not ported.
- **O-4 Parallel AOM/FORT feedback not implemented** (section 10.5); measure
  convergence first.
- **O-5 Two-node FORT polarization optimizer not implemented** (waits on O-4).
- **O-6 Parity 10/11 attempt-loop underflow** at `n_excitation_attempts` 3
  (RID 38813 with the flag on; RID 38786 before parity 11 was wrapped: the
  flag is archived True there, but its tries series is seed-only). The
  `break_realtime` guard every sibling has inside `for excitation_attempt in
  range(self.n_excitation_attempts)` was commented out in `a2ff08d`
  (2026-08-17). Run flag-on checks at `n_excitation_attempts` 1 until settled.
- **O-7 Two-node confirm loop defeats recycling.**
  `Two_nodes_atom_loading_experiment` zeroes `AllSPCMs_RO2` on every attempt
  including the first; the single-node form clears only on retry. Hardware
  impact unmeasured: RID 38785 (83 measurements) is the only two-node run
  with the flag, and no two-node run followed through 2026-10-14.
- **O-8 Two-node loaders run PGC at the loading FORT depth**; single-node
  loaders drop to `stabilizer_FORT.amplitudes[1]` first.
- **O-9 Single-node `first_shot` A5/A6 asymmetry** (RO1 six beams, RO2 four
  when `PGC_and_RO_with_on_chip_beams` is True); fixed only in the two-node
  copy.
- **O-10 Unprotected slack.** Direct `laser_stabilizer.run()` calls without a
  restore in `MOT_loading_experiment`, `atom_loading_optimizer_experiment`,
  `single_photon_experiment_3_atom_heat_test`,
  `single_photon_experiment_4_atom_loading_advance`,
  `single_photon_experiment_5_atom_loading_advance_AllSPCM`,
  `FORT_monitoring_with_Luca_experiment`; the loops of
  `atom_loading_optimizer_experiment` and the legacy
  `load_until_atom_in_both_nodes_recycle` were deliberately left alone.
- **O-11 AOMsCoils core resets.** `run_feedback()` and `k10cr1_operations()`
  start with an unconditional `self.core.reset()`, undercutting
  `reset_core=False`; whether it disturbs the other node's TTL outputs is
  unverified on hardware.
- **O-12 GVS `microwave_Rabi_2_and_EXC` runs 38856-38863 show seed-only
  `n_first_shot_confirm_tries`** although the function is wrapped; cause not
  verified.
- **O-13 Pre-existing compile faults**: parity 2 and 3 play back
  `ba_dma_handle` without fetching it; parity 4-7 store float
  `target_780_HWP/_QWP` into int lists (`angle_780_HWP = [-1] * n`; parity
  9-11 fixed this with `[-1.0]`). They compiled wrapped only in scratch copies
  with the faults stepped around.
- **O-14 Compile-sweep coverage**: most of the 31 wrapped functions and all of
  parity 2-7 have no case; add one per family touched.
- **O-15 Stage 2b of confirm-by-first-shot (not done)**: delete the now
  redundant in-loop microwave arming (delay 1 ms; `ttl_microwave_switch.on()`;
  delay 20 us; `dds_microwaves.sw.on()`) at the 7 arming sites (parity 4, 5,
  6, 7, 9, 10, 11), as its own physics-visible commit after hardware runs;
  saves 4 RTIO events per measurement.
- **O-16 `single_atom_threshold_Node2` = 17000** sat below the whole
  first-shot distribution when the confirm loop was planned (median 22,700
  counts/s). Raising it trades against the recycle rate (same variable): scan
  it, do not guess.
- **O-17 MicrowaveScanOptimizer fixes not yet exercised** (section 11.1;
  committed 2026-10-16 as `619d3ec` at the user's request). The code has run
  on hardware since 2026-10-14 19:44 (RIDs 38942-38951, 38961-38965 and later;
  R26), but neither fixed path (an `n_measurements` override; a completed
  Frequency_11/m10/m11 scan whose fit check should now print the pre-scan
  calibration) has been exercised through 2026-10-16. Confirm from such a
  run. The magnetometer seeding in the same file is O-35.
- **O-33 `_wait_for_satellite` does not wait** (section 6.8): its 100 polls
  are separated only by timeline `delay()`, so they take microseconds of
  wall time and a down destination raises at once. A real wait needs a
  wall-clock pause (for example a host-side retry, or busy-waiting on
  `core.get_rtio_counter_mu()`). The skew probe has the same loop. Never
  exercised, because destination 1 has always been up.
- **O-34 The user's working-tree edits to `subroutines/experiment_functions.py`**
  (shared with standalone, R1), 2026-10-14 evening, committed in `bfca64d`
  on 2026-10-16: in
  `measure_Magnetometer` (Node1) the coil-settle delays before the Zero and MOT
  samples went from 200 ms to 2 ms and the one before OP from 1 ms to 2 ms;
  in `end_measurement`, three already-commented lines were re-commented. The
  version of 20:48 instead had a live `measure_Magnetometer(self)` call in the
  `which_node == "alice"` block of `end_measurement`, which RID 38965 ran
  (it then failed, O-35). The assistant's addition to the same file is the
  `break_realtime` of O-35. The unit tests and the microwave compile cases
  pass with both; the full compile sweep was last run before them. The two-node
  `end_measurement` is a separate copy and is unaffected.
- **O-35 Magnetometer monitoring** (written 2026-10-14, about 21:04-21:30;
  committed in `bfca64d` on 2026-10-16). Turning on
  `monitor_magnetometer_in_end_measurement` failed twice in
  MicrowaveScanOptimizer, and both causes are fixed:
  1. Missing datasets: the optimizer never created `Magnetometer_*_NodeX`
     (RID 38965, KeyError). Fixed by the seeding in section 11.1.
  2. Slack: RID 38968 (Node1 Frequency_11, override
     `monitor_magnetometer_in_end_measurement` True, 21:10) recorded six
     readings of all three phases, then raised `RTIOUnderflow` at
     `zotino0.set_dac` (channel 33, slack -549596 mu, i.e. -550 us), the
     readers' first output. `end_measurement` opens with
     `set_dataset(measurements_progress)`, a SYNCHRONOUS RPC (only
     `append_to_dataset` is async here), and its 4 x 1 ms delays did not
     always cover the host round trip; the next output used to be one that
     re-synchronised, so this never showed with the monitor off. Fixed by
     `self.core.break_realtime()` at the top of `measure_Magnetometer` and
     `measure_Magnetometer_Node2_X_Y_Z` (forward-only; nothing in them is
     timed against an earlier event), which covers `end_measurement` and the
     four loaders' `monitors_for_atom_loading` calls alike.
  Also added: three magnetometer applet specs, MANUAL ONLY -- an automatic
  version (created only while the flag was on, closed otherwise) was built and
  then removed at the user's request the same evening -- and
  `magnetometer_axis_labels` (sections 6.9, 14.2, 14.3), with tests
  `test_magnetometer_applets_are_manual_only` and
  `test_magnetometer_results_exist_for_every_experiment_call`; positive
  controls fail as expected. RID 38968's six readings, mG, Node1: Zero X 143,
  Y -213, Z -1125; OP X 76, Y -348, Z -1096; MOT X 1300, Y -208, Z -1629.
  Caveats: each reader call drives the coils through three settings (Zero,
  OP, MOT) every measurement; `measure_Magnetometer_Node2_X_Y_Z` is a
  coil-test variant whose first phase sets coil channel 0 to 3 V, and it
  scales by `mGauss_per_Volt = 1`, i.e. records volts. Live
  `monitor_magnetometer_in_end_measurement_NodeX` and
  `monitors_for_atom_loading_NodeX` remain False on both nodes; RID 38968 set
  the flag by override. On hardware after the slack fix (saved 21:13),
  Node1 health checks with the monitor on (by override) completed 100
  measurements with 100 readings in RIDs 38976 (2026-10-14), 38980 and 38982
  (2026-10-16); 38972, 38973 and 38975 stopped after 69-75 measurements and
  their h5 files do not record why; 38969 started a minute before the fix
  was saved. No Node2 run with the monitor on exists yet.
- **O-36 Feedback buffer index for a subset stabilizer** (found 2026-10-16).
  `AOMPowerStabilizer` sets `buffer_index = sampler_ch + 8*i`, with `i` the
  group's position in `feedback_channels.json`, but samples only the groups
  that hold one of its `dds_names`. A stabilizer built for channels that all
  sit in a LATER group therefore indexes past its buffer: Node2 with
  `dds_names=["dds_AOM_A1"]` (index 14, buffer 8) and, all along, Node1 with
  `["dds_FORT"]`. Not hit today: every experiment builds from
  `fast_feedback_dds_list_NodeX`, which includes `sampler0` channels on both
  nodes. Fix (shared code, user to decide): number by the cards actually
  sampled, `8 * (len(self.sampler_list) - 1)` after the append.
- **O-18 No `MonitorMagnetometer_master_satellite` run** (nor GVS
  `test_magnetometer_experiment`) exists through 2026-10-14. The first
  master-satellite magnetometer readings are from `end_measurement` in RID
  38968 (Node1, 2026-10-14 21:10, O-35).
- **O-19 Per-channel validation of AOMsCoils** (active-low polarity, coil
  identity/polarity) is undocumented.
- **O-20 Standalone rollback** with the old GVS on standalone gateware has not
  been exercised since the migration.
- **O-21 Clock check** (kernel TTL frequency against the lab reference) not
  recorded; the boot-time sinara_tester ritual should stay gone over several
  power cycles (only one boot with the new Urukuls was observed).
- **O-22 DMA traces with satellite hardware** are in use without a formal
  revalidation record. Separately, the alias-list guard test proposed in the
  D1 upgrade note (alias list == `ALIAS_MAP` keys == `DDS_DEFAULTS` keys for
  alice and bob; `utilities/DeviceAliases.py` silently skips missing aliases)
  does not exist.
- **O-23 No master-satellite HealthCheck files**; when written they must
  submit `MicrowaveScanOptimizer_master_satellite` with `which_node`.

### 18.3 Documentation and code hygiene

- **O-24 Repoint citations** to this file (section 1.6). DONE in the working
  tree on 2026-10-14 for the four code sites and both memory notes; commit
  with this file. Was: FORT optimizer
  docstring (TWO sites: "plan_codex_detail.md section 18" in the header of
  the two_nodes design block and the bare "once section 18 allows it" in its
  feedback bullet; both go to section 10, and the bare one would otherwise
  silently point at this file's section 18), `_install_wiring_metadata`
  docstring, skew-probe docstring,
  `standalone/README.md`, memory notes `gui-arguments-win-over-datasets.md`
  and `master-satellite-mixin-conventions.md`.
- **O-25 Stale comments and notes**: GVS mixin `prepare()` comment ("two-node
  mode must not", followed by "Both modes now"); `_create_node_applets`
  docstring (says `prepare_laser_stabilizer` is skipped in two-node); Base
  `_compute_derived_amplitudes` docstring ("prepare_laser_stabilizer refuses
  outside single_node"); `AtomLoadingOptimizer` module docstring (same claim);
  Base `_bind_other_node_dds_switches` comment ("AOMsCoils relies on that";
  AOMsCoils configures in prepare); `initialize_result_state` comment naming
  the removed `Two_nodes_atom_loading_2_experiment`; references to the deleted
  `Two_nodes_alternating_shot_experiment` in `experiment_functions_two_nodes.py`
  and `applets_master_satellite.py` (revival instructions point at a function
  that no longer exists); MicrowaveScanOptimizer docstring (`fit_model_dict`
  "shared"); compile-checker docstring counts; `master_satellite_test_files/README.md`
  misnamed DB; compile checker: docstring count, the CASES comment "Only two
  functions are registered in two-node mode" and the comment that
  `tune_852_waveplates_to_target_in_experiment` is not declared (section
  7.7); Base `_publish_broadcast_aliases` docstring (still describes bare
  switching aliases that act on BOTH crates, removed in `e768d50`); GVS mixin
  `_create_node_applets` docstring ("Two-node mode gets the SHARED applets
  only", though it also creates the TWO_NODE specs). Memory notes `keep-sequences-linear.md` (still describes the
  removed broadcast alias objects), `ccb-auto-applets.md` (says two-node
  raises NotImplementedError and skips the stabilizer),
  `master-satellite-mixin-conventions.md` (`initialize_single_node_result_state`,
  "aom_feedback stays identical to main", test counts) and
  `parallel-aom-feedback-grouping.md` (iteration cap 200 vs live 100) are
  stale. The root `README.md` describes only the standalone era.
- **O-26 Orphan dataset** `max_loading_rounds_two_node` in `dataset_db.pyon`
  (delete with `artiq_client del-dataset` when the user agrees). Of the seven
  declared but unread globals, the five `*_two_node` timings could be
  removed; keep `parallel_AOM_feedback` (future feedback selector) and
  `t_Node2_excitation_delay_mu` (placeholder) (section 7.7).
- **O-27 `ms_node_of` misclassification risk** (section 15.4).
- **O-28 `alternating_shot_threshold_for_loading_NodeX`**: live 50000/45000
  vs declared 62000/60000; irrelevant while O-1 stands.
- **O-29 Side branch `c8607cd`** (1D PGC with bias field) is unmerged.
- **O-30 Missing ARTIQ manual**: the retired PLAN.md named
  `documentation/ARTIQ manual 7.82.pdf` as the authoritative ARTIQ reference;
  it is not in the repository (`*.pdf` is gitignored) and was not found
  nearby. Use the ARTIQ release-7 source (e.g. `gateware/rtio/sed`, runtime
  `rtio_mgt.rs`, `coredevice/ad9910.py`, `master/worker_db.py`) and the
  official ARTIQ 7 manual instead of guessing.

### 18.4 Queued upgrades (documented, deliberately not started)

- **O-31 Chopped readout** (`documentation/upgrade_notes/chopped_RO_ideas.md`):
  low-loss readout chopping the FORT out of phase with the light (Saffman lab:
  50% duty at 1.25-1.42 MHz; target 0.5-1.5 MHz, f_chop >> radial trap
  frequency 90-105 kHz). As built, `record_chopped_RO` is played back by
  `first_shot`/`second_shot` behind `use_chopped_readout_NodeX` (False on both);
  Node1 values: FORT off = pulse = gate = 500 ns, pulse offset 100 ns, gate
  offset 900 ns, 4000 cycles. DMA is mandatory (32,000 events x 713 ns = 22.8
  ms to submit 7.6 ms, 300% duty). The 900 ns gate offset implies about 800
  ns of cooling-light path delay (slow DP AOM), capping duty near 26% at 526
  kHz; focusing at the DP AOM (100-150 um waist; check TeO2 slow-shear) should
  shrink `chop_RO_gate_offset` (falsifiable). Traps: it reads only 2 of 4
  SPCMs and AVERAGES them (about 4x smaller numbers, so thresholds do not
  transfer); it uses raw `self.ttl0/ttl1` (in two-node mode needs the `_NodeX`
  counters and all 4 SPCMs); re-record the trace per scan point and keep the
  optics block outside it.
- **O-32 D1/GRIN alias rename** (`D1_optical_pumping_merge_ideas.md`,
  2026-10-11): rename `GRIN1and2_dds -> GRIN1_dds` and `dds_D1_pumping_DP ->
  GRIN2_dds` in both stacks and the alias/feedback JSON keys (about 307
  occurrences in 21 files; zero electrical change because only alias KEYS
  change), verify Node2 on hardware, then give Node1 a real D1 TTL (now ttl15,
  section 2.7) and
  an external DDS. Prerequisite: the alias guard test (O-22). Latent bugs
  recorded there and left alone: `chopped_optical_pumping` has no
  Node2-correct path; several Node1-only unconditional sequences mis-drive
  Node2; the `aom_feedback` GRIN special case never fires in two-node mode;
  Bob's `feedback_channels.json` references the undeclared `set_point_D1_DP`;
  Node1 pre-programs `dds_D1_pumping_DP` at init and Node2 does not.
- **Shared blocker for O-31 and O-32**: is two-node STANDALONE synchronization
  retired? On Node1, ttl12 (`ttl_node1_output1`) is driven by the standalone
  handshake (about 10 sites) and ttl15 is the other handshake line. Retiring
  the handshake frees ttl12 (D1 TTL on both nodes) and ttl15 (chop gate);
  the alternative is moving the handshake to ttl8/9/10, which Node1 does not
  alias. User decision. Answered 2026-10-19: retired, cables removed (section
  2.7). Node1 ttl15 is now Node1's `ttl_D1_pumping_switch`. ttl12 is Node1's
  dummy output and Node2's ttl15 (ttl31) Node2's, so reusing either means
  moving that node's dummy first.

### 18.5 Next steps (suggested order, user decides)

1. Exercise the committed MicrowaveScanOptimizer fixes (O-17) with a run
   that overrides `n_measurements` and a completed Frequency_11 scan.
   (Magnetometer-monitored health checks already complete on Node1, O-35;
   GVS override recording after `ef55cae` is confirmed on hardware, section
   8.5.)
2. Fix the stale comments (O-25) in one reviewable change (the plan
   citations, O-24, are already repointed).
3. Run the TTL skew probe once the T-adapters are fitted (O-2).
4. Settle O-6 and O-7; decide O-15 and O-16.
5. Measure feedback convergence, then design parallel feedback (O-4).
6. Two-node physics beyond loading: herald/mapping timing, photon
   coincidence (O-3), a `Two_Nodes_` analysis notebook.

---

## 19. Decision log and history

Condensed and dated. "Superseded" entries record claims from the retired
files that are no longer true.

- **2026-08-26** The legacy standalone `two_nodes` node mode is removed for
  good (`e86046f`).
- **2026-08-27..30** Architecture chosen: separate stacks, one unified
  database, one explicit physical mapping, a resolver separate from Base,
  suffix naming, master-local canonical SPCMs, USB on the master host,
  examination-safe deferred configuration. Development happened on branch
  `20260826_master_satellite_codex` on a non-experiment computer; work now
  happens on the experiment PC (`artiq-master` has live `results/` and
  `dataset_db.pyon`).
- **2026-08-30** Explorer duplicate `GeneralVariableScan_master_satellite1`
  traced to a public import of the GVS class in the CatchError module; fixed
  by a private alias, later made moot by the `_mixin` convention.
- **2026-08-31** GVS split per mode (`025c027`); module renamed
  `experiment_functions_two_nodes.py`. CatchError GVS files restored
  (`02b5820`) and deleted the same day in favour of `enable_Catch_UnderFlow`
  (`9904f9b`). EV files renamed with `_master_satellite_` and `global`.
  AOMsCoils split per node; the original combined design (Base in two_nodes so
  an unselected node could be forced safe; selector `which_node = node1 |
  node2 | two_nodes`; run-local GUI coil voltages with examination-safe 0 V
  defaults; both Zotinos zeroed; feedback and K10CR1 excluded) is superseded
  by per-node files, `disable_coils`, `run_laser_feedback` and K10CR1.
- **2026-09-01** Standalone MicrowaveScanOptimizer is NOT compatible
  (standalone Base, unsuffixed targets, legacy feedback/result state, repeated
  `prepare()`, unvalidated DMA); a separate port was written and pulled ahead
  of the validation order by explicit decision.
- **2026-09-03** Standalone entry points moved to `standalone/`.
- **2026-09-07..09** Hardware gate executed; Node1 board replaced; link up.
  Superseded: "both physical nodes still run standalone gateware", "Node1
  carrier revision v1.0", "no deployable image, route table or proof of the
  master SFP cage", "readiness for controlled gateware integration, not
  proof".
- **2026-09-09** Superseded: "GVS-level feedback hardware validation remains
  outstanding" (feedback ran in GVS from 2026-09-10).
- **2026-09-10** Underflow rationale corrected: slack erosion, not DRTIO
  latency (`56f269a` corrects `edfb3c1`). Superseded: "add no speculative
  slack until hardware demonstrates a failure". Applets became experiment-driven
  through CCB; superseded: "applet restructuring is last priority, do not
  modify applets during the migration".
- **2026-09-15** Urukul SYNC: hardware fault, cards replaced; verdict rule
  fixed after two rejected rules. Superseded: "recalibrate AD9910 after final
  gateware" as an open step.
- **2026-09-17** Joint thresholds made global. `force_other_node_off`
  introduced; superseded: "single-node removes the unselected node from the
  lifecycle so hardware access is not broadened" and "there is no cross-node
  force-safe gate" (both now true only of AOMsCoils / init).
- **2026-10-01** Per-node override dictionaries; the shared field removed.
  Superseded: the GUI organization with one `override_experiment_variables`.
- **2026-10-05** `n_measurements` persisted by every run, against Base's own
  "broadcast, never persisted" rule, to protect applets and self-heal
  standalone de-persistence.
- **2026-10-06** `t_delay_in_bob_mu` split into `t_Node2_rtio_offset_mu`
  (electrical) and `t_Node2_excitation_delay_mu` (optical). Superseded: the
  global name `t_delay_in_bob_mu`.
- **2026-10-07** Two-node mode in Base and feedback (`e768d50`). Decided: no
  bare device aliases in two_nodes (the broadcast-wrapper design that let one
  bare object drive both crates was tried and removed after it drove Node1
  hardware from `laser_stabilizer_Node2`); per-crate sequential feedback via
  `aom_feedback` `node_suffix`. Superseded: "single-node feedback is
  implemented without modifying `aom_feedback.py`; shared feedback code stays
  byte-identical to main" and "two-node mode must not call
  `prepare_laser_stabilizer` or the result state; master-satellite feedback is
  single-node only". A decomposed two-node draft was rejected in favour of the
  linear legacy shape (R9). `with parallel` adopted for differing durations.
  SED discard found and the 4-events rule adopted. The alternating diagnostic
  (`3a90e2d`) broke AOMsCoils through new variables and was reverted next day;
  rule R12 dates from this.
- **2026-10-08** Parallel feedback grouping decided (optical constraint), not
  implemented. MonitorSPCMinApplet fixed to Node1.
- **2026-10-09** Per-node alternating atom check disabled (does not load).
- **2026-10-11** PR #68 merged.
- **2026-10-12** Node-resolved alternating retention replaces the alternating
  shot experiment; one atom-loading experiment with the advance flag replaces
  two; dead `t_coil_relaxation_time_OP_NodeX` removed (superseded: the retired
  spec's "confirmed values" listing them as 1 ms / 0.4 ms seeds). First
  all-at-once confirm-loop attempt crashed the compiler and was reverted.
- **2026-10-13** Confirm-by-first-shot rolled out one site first, then by
  family (R8). Parity 12/13 and `single_photon_experiment_10` removed after
  parity 12's RO2 reset showed a roughly 7x throughput loss. One "Single node"
  applet group. `max_loading_rounds` replaces `max_loading_rounds_two_node`.
- **2026-10-14** Overrides recorded in the h5; loader fallback for older
  files; first hardware use of the recording (RIDs 38952-38957). Recorded
  that the master serves the working tree (R26), after the uncommitted
  MicrowaveScanOptimizer fix was found to have run on hardware. Plan
  documents consolidated into `plan_claude_detail.md` /
  `plan_claude_concise.md`. Corrected during review of this file's draft:
  the `tune_852_waveplates_to_target_in_experiment_NodeX` history (section
  7.7), the RID 38710 loading-time figures (section 9.5) and the
  confirm-loop evidence labels (section 9.7).
- **2026-10-16** Committed `619d3ec` (MicrowaveScanOptimizer fixes),
  `bfca64d` (magnetometer monitoring) and the plan consolidation. The
  magnetometer applets were first built to open automatically while
  `monitor_magnetometer_in_end_measurement` was on and close otherwise; at
  the user's request they became manual-only the same evening, before
  anything was committed.

Other superseded claims from the retired files:

- File names: `GeneralVariableScan_master_satellite.py`,
  `GeneralVariableScan_CatchError_master_satellite.py`,
  `AOMsCoils_master_satellite.py`, `ExperimentVariables_Node1/_Node2.py`,
  `ExperimentVariables_master_satellite.py`,
  `subroutines/experiment_functions_master_satellite.py`, top-level
  `GeneralVariableScan.py`/`ExperimentVariables.py` (now in `standalone/`).
- The GVS selector `experiment_mode` as a GUI value and `which_node` as the
  node argument (now per-file `EXPERIMENT_MODE` and `selected_node`).
- Explorer lists of 6 and 7 entries; test counts 56 and 58.
- "Approved globals" limited to four names.
- The magnetometer as the first reused GVS target (validation went AOMsCoils,
  then atom loading).
- "K10CR1 is omitted from the unified database" (it is appended there).
- "Confirm standalone files were untouched" as a reporting step (shared and
  standalone files are now deliberately co-edited under R1).
- `PLAN.md` example `self.dds_FORT_Node1 -> urukul0_ch1` (wrong; section 5.5)
  and legacy SPCM alias names written as `SPCM0` (actual `ttl_SPCM0`).
- "ARTIQ manual 7.82.pdf in `documentation/` is the authoritative ARTIQ
  reference" (the file is absent; section 18, O-30).
- "Base.initialize_hardware resets the core once per run" (once per call; GVS
  calls it per scan point).

---

## 20. Reporting contract

### 20.1 Before editing

Before acting: read `plan_claude_concise.md` and the sections of this file
that cover the files you will touch; inspect the actual code (this file may
lag it); identify the smallest relevant set of files, their dependencies and
possible regressions; prefer the smallest safe change; do not blindly
normalize legacy naming or silently migrate historical logic; when hardware
semantics are uncertain, say so instead of guessing (R5). Remember that a
saved edit is live for the next submission (R26).

Report: files inspected, facts established (with evidence), ambiguities,
risks, and the exact bounded scope. For structural changes also give current
behaviour, relevant code paths, proposed behaviour, files to modify,
compatibility impact, alternatives (R6) and needed hardware validation. For
large cross-cutting changes propose an upgrade note instead (R7). For
repetitive changes name the pilot site (R8).

### 20.2 After editing

- List changed files and APIs.
- Run the unit tests (section 16.1) and, after any master-satellite code
  change, the compile sweep (section 16.2); report PASS/COMPILED/FAIL counts,
  the target kernels in the compiled list, and the ACTION NEEDED block.
- Run `git diff --check` and show `git diff --stat`.
- State which shared or standalone files changed and why the change is
  backward compatible (R1).
- State whether any ExperimentVariable was added (R12) and whether seeding is
  required before the next run.
- Identify what still needs hardware validation and how to run it (Explorer
  entry, arguments, expected result).
- Update this file (and the concise file) when the design moved.
- Never commit or push automatically (R2).

### 20.3 Hardware work records

For hardware work record: image hashes, routing/configuration, cabling, boot
logs, destination readiness, the exact operation tested, RIDs, and pass/fail
evidence. Decode RTIO channel numbers with `channel >> 16` = destination.
Local Python tests and the compile sweep never substitute for
experiment-system validation.
