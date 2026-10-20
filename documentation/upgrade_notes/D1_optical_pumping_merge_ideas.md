# D1 and optical pumping: merging Node1 and Node2

**Started:** 2026-10-11 **Status:** **Partly applied** — Node1's D1 hardware
change and its D1 switch are done (Status 2026-10-20, below); the renames
(stages 1-2), Node1's GRIN1 path and the stage 4 merge are not started.
**Scope:** both stacks (standalone and master-satellite) plus the per-node
alias config. Not the device_db.

---

## Status 2026-10-20

- Hardware (user, 2026-10-19): Node1's D1 AOM runs from an external RF source
  switched by `ttl_D1_pumping_switch` = Node1 `ttl15` (active low: `off()` = D1
  on; blocked at init). `dds_D1_pumping_DP` (`urukul2_ch2`) no longer drives it.
- Code (2026-10-20, the user's rule: a set of `dds_D1_pumping_DP` is commented
  out with a note, `sw.on()` becomes `ttl_D1_pumping_switch.off()`, `sw.off()`
  becomes `.on()`). Node1's D1 is switched by the TTL in
  `CW_optical_pumping_node1`, `record_CW_optical_pumping_node1`,
  `single_photon_experiment_4_atom_loading_advance` and
  `single_photon_experiment_5_atom_loading_advance_AllSPCM` (both Node1 only),
  `measure_GRIN1` and `record_chopped_optical_pumping` (treated as Node1 code),
  `record_CW_optical_pumping_both_nodes` (Node1 side), `standalone/AOMsCoils.py`
  (alice branch), the two `MOT_experiments` scripts and
  `tests/A2-SPCMGateTimeVsCounts.py`. Standalone alice init's explicit D1 tone
  set is commented out. AOMsCoils master-satellite Node1 did this on 2026-10-19.
- Left as they are:
  - the `which_node == "bob"` branches, and the Node2-written functions that
    use `dds_D1_pumping_DP` as GRIN2;
  - `dds_D1_pumping_DP.sw.off()` after the excitation loop in single photon
    6-9, and right after the chopped-pumping DMA in parity 2-5 (that DMA now
    ends every chop with D1 blocked). On Node1 these switch an unused channel;
    on Node2 they turn its GRIN2 RF off;
  - the init loops that program `urukul2_ch2` with its RF off: master-satellite
    `_initialize_dds_group` and standalone `named_devices.initialize()`;
  - the inactive D1 channel in alice's `feedback_channels*.json`;
  - `chopped_optical_pumping`, which nothing calls;
  - `tests/A2-ExcitationSPCMGateStartTime.py`, written for Node2: pending the
    user's decision.
- Not done (user: "D1 switch only for now"): Node1's D1 now passes the GRIN1
  AOM, but `CW_optical_pumping_node1`, `record_CW_optical_pumping_node1` and
  the two-node Node1 side neither open `ttl_GRIN1_switch` nor turn on
  `GRIN1and2_dds` (parity 2-5, single photon 4/5 and `measure_GRIN1` keep their
  own GRIN1 lines). Check on hardware whether D1 reaches the atoms.
- Later on 2026-10-20 the unused functions moved to
  `subroutines/experiment_functions_obsolete.py`. Of those named here, single
  photon 3_node2, 4, 5, 6, 7 and 8, parity 2-5 and 6_node2 and
  `test_excitation_rise_time_experiment` are there now; the pumping functions,
  `measure_GRIN1`, the chopped-pumping DMA, `chopped_optical_pumping` and single
  photon 9 stayed.

## The problem

*As written on 2026-10-11; Node1's D1 hardware has changed since (Status
2026-10-20, above).*

Both nodes run the same `subroutines/experiment_functions.py`. To keep that
possible, Node2 was given new hardware without new names: it reuses two Node1
alias names for physically different devices. The code says so itself, in the
docstring of `CW_optical_pumping_node2` at `subroutines/experiment_functions.py:2409-2414`:

```
Note: To avoid conflict with Node1 codes, I left the names of the dds channels
that are now used for GRIN1 and GRIN2 dds.

Name of the dds channels >>>>   what actually does here
GRIN1and2_dds            >>>>   GRIN1 dds
dds_D1_pumping_DP        >>>>   GRIN2 dds
```

| Alias | Node1 / alice | Node2 / bob |
|---|---|---|
| `GRIN1and2_dds` | `urukul0_ch3` — one DDS feeding **both** GRINs, carries excitation (`f_excitation`) | `urukul2_ch2`→`urukul5_ch2` — **GRIN1 only**, carries D1 pumping (`f_GRIN1_D1_pumping`) |
| `dds_D1_pumping_DP` | `urukul2_ch2` — a real D1 double-pass AOM at 373.075 MHz | `urukul2_ch0`→`urukul5_ch0` — **GRIN2**, carries excitation (`f_GRIN2_excitation`) |
| `ttl_D1_pumping_switch` | **dummy**, aliased to `ttl3` — "not used in node1" | real gate, `ttl12`, in series ahead of GRIN1 |
| `ttl_GRIN1_switch` | `ttl14` — routes the single shared tone to GRIN1 (D1) | `ttl14` — per-GRIN optical gate, own DDS behind it |
| `ttl_GRIN2_switch` | `ttl13` — routes the single shared tone to GRIN2 (EXC) | `ttl13` — per-GRIN optical gate, own DDS behind it |

**Node1 has one tone and two optical gates. Node2 has two tones plus a
dedicated D1 gate.** The names describe Node1.

The per-node amplitudes are the clearest evidence of the split on Node2:
`p_GRIN1_D1_pumping_Node2 = -8` (live) with `p_GRIN1_excitation_Node2 = -150`
(killed), and `p_GRIN2_excitation_Node2 = -6` (live) with
`p_GRIN2_D1_pumping_Node2 = -150` (killed). One function each: GRIN1 ≡ D1,
GRIN2 ≡ excitation. On Node1 all four sit at a neutral −8 because none is used
there.

## The goal

Both nodes use the same names **and** the same configuration, so that a
rollback to standalone gateware still works. Three stages:

1. **Rename** `GRIN1and2_dds` → `GRIN1_dds`, `dds_D1_pumping_DP` → `GRIN2_dds`.
   Makes the names true on Node2 and knowingly false on Node1.
2. **Verify Node2** still behaves identically — the rename should be invisible.
3. **Re-configure Node1** to match Node2, after which (4) Node1 can share
   Node2's OP/EXC code.

---

## Stage 1 — the rename

### What makes it safe

**No `device_db*.py` contains any logical alias name.** The device_db knows only
`urukul{0..5}_ch{0..3}`, `ttl{0..31}`, `sampler{0..5}`, `zotino{0,1}`.
Resolution is two JSON stages:

```
logical alias ──config/<alice|bob>/device_aliases.json ALIAS_MAP──> standalone name
standalone name ──config/master_satellite/device_aliases.json[node]──> unified device_db key
```

Only the **first** stage knows the names being renamed. Renaming an `ALIAS_MAP`
*key* while leaving its *value* untouched means **zero electrical change** —
every DDS keeps its exact urukul channel on both nodes.

`subroutines/experiment_functions_two_nodes.py` has **zero** occurrences, so
the two-node sequence and its AST guards are unaffected.

### Scale

307 occurrences across 21 files: `GRIN1and2_dds` ×199, `dds_D1_pumping_DP` ×108.
`GRIN1_dds` and `GRIN2_dds` appear **nowhere** today — the target names are
free. The only near-relatives are the dead commented-out `# 'dds_GRIN1'` /
`# 'dds_GRIN2'` at `utilities/BaseExperiment.py:694-695`, which this finally
makes real and which should be deleted as part of the change.

### What to change, by layer

**Alias and feedback config — keys only, defaults unchanged.**
`utilities/config/alice/device_aliases.json` and `.../bob/device_aliases.json`:
the two `ALIAS_MAP` keys and the two `DDS_DEFAULTS` keys. The defaults must keep
recording the per-node difference:

| | `GRIN1_dds` defaults | `GRIN2_dds` defaults |
|---|---|---|
| alice | `f_excitation` / `p_excitation` | `f_D1_pumping_DP` / `p_D1_pumping_DP` |
| bob | `f_GRIN1_D1_pumping` / `p_GRIN1_D1_pumping` | `f_GRIN2_D1_pumping` / `p_GRIN2_D1_pumping` |

Also the four feedback configs, which key on the same alias names:
`config/alice/feedback_channels.json`,
`config/alice/feedback_channels_FORT_APD.json`,
`config/alice/feedback_channels_FORT_uses_MM_fiber_monitor.json`,
`config/bob/feedback_channels.json`.

**Binding declarations.**
`utilities/BaseExperiment_master_satellite.py` — `DDS_ALIASES` at `:128`, `:131`.
`utilities/BaseExperiment.py` — both alias lists, alice `:592`/`:595`, bob
`:697`/`:698`, plus deleting the commented pair at `:694-695`.

**Shared code — imported by both stacks.**
`subroutines/aom_feedback.py`, 10 sites, including `ambient("GRIN1and2_dds")`
at `:357` and the two `if ch.name == 'GRIN1and2_dds'` compares at `:691`, `:730`.
`subroutines/experiment_functions.py`, 212 sites.

**Experiments.**
`AOMsCoils_master_satellite_Node1.py` and `AOMsCoils_master_satellite_Node2.py`
(self-contained since 2026-10-19, replacing the deleted shared mixin) — the
`GRIN1and2_ON` / `D1_pumping_switch_ON` boxes in `build()` and `turn_on_AOMs()`,
Node1's D1-opens-GRIN1 rule, Node2's combined GRIN blocks at the end of
`turn_on_AOMs()`. `standalone/AOMsCoils.py` (9 sites).
`MOT_experiments/MonitorMOTandExternalBeamPositions.py`,
`MOT_experiments/MOTDetuningAndGradBScan.py`.
`tests/A2-ExcitationSPCMGateStartTime.py` (34),
`tests/A2-SPCMGateTimeVsCounts.py` — bring-up scripts, not unittests.

**Tests that pin the old names.**
`tests/test_base_experiment_master_satellite.py:549`,
`tests/test_aoms_coils_master_satellite.py`:
`test_node2_d1_uses_authoritative_grin2_defaults`, and every test that sets the
`GRIN1and2_ON` / `D1_pumping_switch_ON` boxes by name (search for them; the
`combined=False` subtest of `test_run_feedback_keeps_node2_combined_grin_state`
expects everything off, so it keeps passing with stale names),
`tests/test_experiment_variables_master_satellite.py:376-387`.

**Explicitly NOT changed:** any `device_db*.py`,
`config/master_satellite/device_aliases.json` (physical names only), any
ExperimentVariable, `dataset_db.pyon`.

### Decisions already taken

| Question | Decision |
|---|---|
| Scope | Both stacks **and** the alias config. Standalone must keep working — the whole point is that a rollback to standalone gateware still runs. |
| ExperimentVariables | **Leave alone.** `GRIN2_dds` keeps defaulting to `f_D1_pumping_DP` on Node1 — a name mismatch one layer down, accepted because renaming variables makes every master-satellite experiment fail in `build()` until both node initializers are run. Node1's values change at stage 3 anyway, which is the natural time. |
| Latent bugs | **Record, don't fix.** Keeps the rename mechanical, so a token-count check proves it complete and any hardware regression is unambiguously a rename error. |
| GUI arg names | Rename — see the open question below. |

### Open question: the GUI argument names

The AOMsCoils checkboxes are named after the old devices:

| file | old | proposed |
|---|---|---|
| `AOMsCoils_master_satellite_Node1.py`, `_Node2.py` | `GRIN1and2_ON` | `GRIN1_DDS_ON` |
| " | `D1_pumping_switch_ON` | `GRIN2_DDS_ON` |
| `standalone/AOMsCoils.py` | `GRIN1and2_DDS_ON` | `GRIN1_DDS_ON` |
| " | `D1_pumping_DP_AOM_ON` | `GRIN2_DDS_ON` |

2026-10-20: the user renamed the master-satellite `D1_pumping_DP_ON` to
`D1_pumping_switch_ON` and moved it to the TTL-gates group (it opens
`ttl_D1_pumping_switch` on both nodes), so the table starts from that name.

The short forms `GRIN1_ON` / `GRIN2_ON` were the first choice. The `_DDS_` form
is suggested instead because `GRIN1_AOM_switch_ON` / `GRIN2_AOM_switch_ON`
already exist in the same GUI group for the *optical gates*, and the longer name
keeps RF and gate distinguishable. **Not settled.**

Either way this resets saved AOMsCoils argument state once — the checkboxes come
back at their defaults. `Node2_GRIN1_AOM_ON` / `Node2_GRIN2_AOM_ON` (the Node2
combined DDS+gate modes) stay as they are.

### Method

Word-boundary regex on exactly two tokens:

```
\bGRIN1and2_dds\b      ->  GRIN1_dds
\bdds_D1_pumping_DP\b  ->  GRIN2_dds
```

The `dds_` prefix is load-bearing: it stops the second pattern touching
`f_D1_pumping_DP`, `p_D1_pumping_DP` or standalone's `D1_pumping_DP_AOM_ON`,
which are renamed separately (GUI) or not at all (variables). The GUI renames
are done by hand.

### The standalone trap — add this guard first

`utilities/DeviceAliases.py:20` fails **soft**:

```python
for alias in device_aliases:
    if alias in self.experiment.alias_map.keys():
        ...
```

An alias present in the Python list but missing from `ALIAS_MAP` is **silently
skipped** — not even a printed warning — so the attribute is never set and the
failure surfaces much later as an apparently unrelated compile error. Standalone
has **no compile checker**, and `device_db.py` is symlinked to the
master-satellite DB, so a half-finished rename would be invisible on that path.

Before renaming, add a unittest asserting, for both alice and bob, that
`BaseExperiment.py`'s alias list == `ALIAS_MAP` keys == `DDS_DEFAULTS` keys.
Cheap, and it catches this entire class of error.

### Verification

1. **Completeness.** Zero remaining `GRIN1and2_dds` / `dds_D1_pumping_DP`; new
   counts equal the old ones (199 / 108); counts of `f_D1_pumping_DP`,
   `p_D1_pumping_DP` and `set_point_D1_SP` **unchanged**.
2. **Unit tests.** `python -m unittest discover -s tests -p "test_*.py"` with the
   base interpreter — 147 today, plus the new guard.
3. **Compile sweep.** `compile_check_master_satellite.py` from `artiq-master`
   with the artiq7 interpreter. This is the only check that resolves against the
   real `device_db.py`; it covers AOMsCoils, GVS, MicrowaveScanOptimizer and
   MonitorMagnetometer for both nodes and both modes, and a missing alias
   surfaces as `DeviceResolutionError`. Read the ACTION NEEDED block at the end.
4. **Hardware — this is stage 2.** Run `AOMsCoils_master_satellite_Node2`,
   toggle the two GRIN DDS controls, and confirm on the scope that the same
   channels respond as before the rename (`urukul5_ch2` and `urukul5_ch0` on
   Node2; `urukul0_ch3` and `urukul2_ch2` on Node1). Then run optical pumping on
   Node2 via single-node GVS (`CW_optical_pumping_node2`) and confirm OP
   performance is unchanged. The rename should be invisible.

---

## Stages 3 and 4 — Node1 to match Node2

Node1 gains a **real D1 TTL (`ttl15`, done 2026-10-19) plus an external
source**. D1 pumping light is delivered to GRIN1, so `GRIN1_dds` must also be
on to pass it — the same arrangement Node2 already uses.

Node1 then matches Node2 structurally, which resolves stage 4 to **yes**: once
Node1 has a real `ttl_D1_pumping_switch` and two independent GRIN tones,
`CW_optical_pumping_node1` and `record_CW_optical_pumping_node1`
(`experiment_functions.py:2333`, `:2475`) become redundant, and the roughly
twelve duplicated `if self.which_node == "alice":` excitation-arming blocks
collapse to the Node2 form. Those blocks are at `experiment_functions.py:5386`,
`:7197`, `:8586`, `:8902`, `:9183`, `:9450`, `:13006`, `:13457`, `:13756`,
`:15204`, `:15642`, `:15973`, with Node2-only teardowns at `:15473`, `:15852`,
`:16038`.

Node1's old dedicated D1 double-pass AOM (`urukul2_ch2`, 373.075 MHz) is free
since 2026-10-19. To make the configs genuinely identical, Node1 would take
Node2's channel assignment: `GRIN1_dds` → `urukul2_ch2`, `GRIN2_dds` →
`urukul2_ch0`. That is a physical repatch and has not been scoped.

### Blocking conflict — settle before any wiring

**Settled 2026-10-19 for master-satellite.** The handshake cables are
removed, and Node1's `ttl_D1_pumping_switch` is now `ttl15` (`ttl_D1_pumping_switch_Node1`
in two-node mode); ttl12 became Node1's dummy output instead
(`plan_claude_detail.md` section 2.7). The Node1 sequences drive it since
2026-10-20; init blocks it, like Node2's ttl12. Standalone alice matches (the
user's edit to `utilities/BaseExperiment.py`): `ttl_D1_pumping_switch` is ttl15
there too, and its init blocks D1. The text below is the analysis as it
stood before.

On Node1, **`ttl12` is currently `ttl_node1_output1`**
(`utilities/BaseExperiment.py:219`,
`utilities/BaseExperiment_master_satellite.py:171`), and standalone's inter-node
handshake drives it at roughly ten sites: `experiment_functions.py:14219`,
`:14285`, `:14361`, `:14442`, `:14516`, `:14523`, `:14549`, `:14663`, and the
`two_nodes_synchronization` / `_synchronization2` / `_trial` family from `:14199`.

`ttl11`, `ttl12` and `ttl15` are indeed unused in master-satellite — DRTIO gives
one timeline and made the handshake obsolete. But standalone still needs that
output line, and standalone has to keep working. `ttl15` is **not** a free
replacement: it is already `ttl_node1_output2`, pulsed at `:14270`.

Two options to weigh:

* move the standalone handshake to a TTL free on both nodes (Node1 does not
  alias `ttl8`, `ttl9` or `ttl10`), which is a one-line change to
  `NODE_TTL_ALIASES` and the alice block of `BaseExperiment.py` because the code
  refers to the alias, not the raw TTL; or
* accept that two-node **standalone** synchronization is retired, and only
  single-node standalone needs to survive a rollback. Defensible, and probably
  the real answer — but it should be a decision, not a discovery.

---

## Latent bugs found during reconnaissance — deliberately NOT fixed

Left alone so the rename stays mechanical and any hardware regression is
unambiguously a rename error. The rename actually makes these easier to see,
because the names will finally mean what they say.

* **`chopped_optical_pumping` has no Node2-correct path.**
  `experiment_functions.py:2248`, DMA body recorded at `:2175`. It chops
  `dds_D1_pumping_DP` — which on Node2 is GRIN2/excitation — and sets
  `GRIN1and2_dds` to `f_excitation` on *both* nodes (`:2286`, `:2289`), where
  Node2's D1 source needs `f_GRIN1_D1_pumping`.
* **Node1-only-correct and unconditional**, so these mis-drive Node2:
  `single_photon_experiment_4_atom_loading_advance:7662`,
  `atom_photon_parity_8_experiment:12831`,
  `test_excitation_rise_time_experiment:243`, `measure_GRIN1:2791`.
  `measure_GRIN1` is already flagged "not implemented in node2"
  (`BaseExperiment.py:353`) yet `GRIN1_sampler_ch = 4` is still populated for
  Node2 in both base layers.
* **`aom_feedback.py`'s GRIN special-case never fires in two-node mode.**
  `:691` and `:730` test `ch.name == 'GRIN1and2_dds'`, but `:390` sets
  `name = ch_name + suffix`, so in two-node mode `ch.name` is
  `'GRIN1and2_dds_Node1'` / `'_Node2'` and neither branch runs. Latent only
  because no GRIN channel sits in a live `fast_feedback_dds_list` — both nodes'
  lists hold only `dds_AOM_A1..A6` and `dds_FORT`.
* **`config/bob/feedback_channels.json:81` references `set_point_D1_DP`**, which
  is declared nowhere; only `set_point_D1_SP` exists, and
  `tests/test_experiment_variables_master_satellite.py:398-399` asserts
  `set_point_D1_DP` is deliberately absent. Latent for the same reason.
* **`f_GRIN2_excitation` is 150.60 MHz on Node1** and in the standalone defaults,
  but **150.50 MHz on Node2** — the only node where it is physically used.
* **Node1 pre-programs `dds_D1_pumping_DP` at init** (`BaseExperiment.py:1019-1020`);
  Node2 has no equivalent, so Node2's GRIN2 tone is only ever set inside the
  OP/EXC arming blocks. (2026-10-20: that explicit Node1 set is commented out;
  both nodes now program the channel only through `named_devices.initialize()`.)
