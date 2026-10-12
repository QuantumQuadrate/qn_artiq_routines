"""Declarative applet set for the master-satellite stack.

Applets are created at run time through the dashboard CCB instead of being
maintained by hand, so their dataset names always match the node that is
actually running. One table here replaces the node x mode multiplication of
applet definitions, and only the active node's applets are ever spawned --
which also keeps the dashboard light.

This module defines no experiment class, so ARTIQ Explorer never lists it.

REQUIRED DASHBOARD SETTING
--------------------------
The applet dock's CCB policy must be "Create and enable/disable applets"
(``ccbp_global = 'enable'``), not merely "Create applets". Under the weaker
``'create'`` policy the dashboard creates each applet entry but leaves it
UNCHECKED -- ``ccb_create_applet`` only ticks the box when the policy is
``'enable'`` -- and ``disable_applet_group`` is ignored outright, so the
other node's applets would never be retired. Set it once, from the applet
dock's context menu.

The table below was reconstructed from the working dashboard backup
(Backup - ARTIQ Dashboard/2026-09-10 Node2), keeping the applets that were
actually ENABLED there and the functional groups they lived in. Applets are
created under a nested group [<node>, <functional group>], so Node1's and
Node2's sets stay separate while the familiar organisation is preserved.

Adding an applet
----------------
Append an :class:`AppletSpec`. Give dataset names in their LEGACY (unsuffixed)
form exactly as the standalone applets use them; :func:`resolve_applet_dataset`
maps each to whatever the running node actually writes:

  * global variables (``n_measurements``)            -> unchanged
  * per-node variables (``single_atom_threshold``)   -> ``..._Node1``
  * redirected results (``FORT_MM_monitor``)         -> ``..._Node1``
  * display config (``photocount_bins``)             -> unchanged

That resolution goes through the same Base methods the experiment uses, so an
applet cannot drift away from the dataset the experiment writes.

``args`` are positional and MUST follow the order of the applet's own
``add_dataset`` declarations; ``options`` are its optional datasets. Values in
``options`` that are not dataset names (a literal like ``1000``) are passed
through untouched.
"""

from collections import namedtuple
from pathlib import Path

APPLET_DIRECTORY = Path(__file__).resolve().parent / "applets"

#: title    - shown in the dashboard applet list
#: script   - filename in applets/, or ("builtin", "<artiq.applets module>")
#: args     - legacy dataset names, positional, in the applet's declared order
#: options  - ((flag, dataset-or-literal), ...) for the applet's optional args
#: group    - functional group, nested under the node group
AppletSpec = namedtuple("AppletSpec", "title script args options group")


def _spec(title, script, args, options=(), group=None):
    return AppletSpec(title, script, tuple(args), tuple(options), group)


_GVS = "GVS and Cycler"
_OPT = "Optimization"
_MON = "Monitor in each measurement"
_PHOTON = "Single Photon Experiment"
_K10 = "K10CR1"

# What every experiment gets, and nothing more. Originally the nine applets
# Node1 actually had ENABLED, transcribed from its live dashboard state
# (AppData/Local/m-labs/artiq/7). Four of those moved out on 2026-09-17 so
# they stop coming up on every scan:
#
#   SPCM count rates (x5)     -> SPCM_MONITOR_APPLET_SPECS, asked for only by
#                                MonitorSPCMinApplet
#   n_excitation_cycles       -> OPTIONAL_APPLET_SPECS, single-photon work only
#   FORT APD / FORT MM        -> SHARED_APPLET_SPECS, now one applet per NODE
#                                and always up, since the FORT is used
#                                everywhere
#
# Keeping this set small is deliberate: 65 applets were defined in that
# dashboard state but only nine running, and applet count is what makes the
# dashboard struggle. Put anything that only one experiment cares about in its
# own tuple, as SPCM_MONITOR_APPLET_SPECS does, or in OPTIONAL_APPLET_SPECS --
# do not grow this list.
APPLET_SPECS = (
    _spec("measurements progress", ("builtin", "progress_bar"),
          ("measurements_progress",)),

    _spec("retention and loading AllSPCMs", "plot_retention_and_loading.py",
          ("AllSPCMs_RO1", "AllSPCMs_RO2", "n_measurements",
           "single_atom_threshold", "t_SPCM_first_shot"),
          (("--scan_vars", "scan_variables"),
           ("--scan_sequence1", "scan_sequence1")),
          group=_GVS),
    _spec("All SPCMs RO1 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_RO1_current_iteration",),
          (("--x", "photocount_bins"), ("--iteration", "iteration"),
           ("--t_exposure", "t_SPCM_first_shot")),
          group=_GVS),
    _spec("All SPCMs RO2 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_RO2_current_iteration",),
          (("--x", "photocount_bins"),
           ("--color", "second_shot_hist_color"),
           ("--iteration", "iteration"),
           ("--t_exposure", "t_SPCM_second_shot")),
          group=_GVS),

    _spec("Atom loading time (s)", "plot_xy_multichannel.py",
          ("Atom_loading_time", "n_measurements"), group=_MON),
)

#: Asked for by GeneralVariableScan only:
#:
#:     create_applets_for(self, self.base,
#:                        specs=APPLET_SPECS + GVS_APPLET_SPECS)
#:
#: Deliberately NOT in APPLET_SPECS. AllSPCMs_atom_check_in_loading is seeded
#: by initialize_result_state on every experiment, so putting it
#: there would open an empty histogram on every run that never loads an atom.
GVS_APPLET_SPECS = (
    # Written by load_MOT_and_FORT_until_atom and its relatives, two values
    # per loaded atom: the count that crossed the threshold, and the first
    # try's count, which did not. So the histogram shows BOTH populations --
    # which is the point, since the experiment function's own comment says it
    # exists "to find a good single_atom_threshold_for_loading".
    #
    # The counts are over t_atom_check_time (10 ms by default), summed across
    # all four SPCMs, and the threshold they are compared against is
    # single_atom_threshold_for_loading (44000 c/s on Node1) -- NOT
    # single_atom_threshold, which discriminates the readout shots.
    #
    # No --iteration: unlike the RO1/RO2 histograms this dataset accumulates
    # over the whole run rather than per scan point, so a per-iteration title
    # would be misleading. plot_hist_autosize.py treats it as optional.
    _spec("AllSPCMs atom check in loading histogram",
          "plot_hist_autosize.py",
          ("AllSPCMs_atom_check_in_loading",),
          (("--x", "photocount_bins"),
           ("--t_exposure", "t_atom_check_time")),
          group=_GVS),
)

#: Asked for by GeneralVariableScan AND MicrowaveScanOptimizer, as
#: shared_specs rather than per-node specs:
#:
#:     create_applets_for(self, self.base,
#:                        shared_specs=SHARED_APPLET_SPECS + SCAN_APPLET_SPECS)
#:
#: Passed through the SHARED path on purpose. time_without_atom is NOT in the
#: dataset-redirect set, so both nodes write that one unsuffixed name: there is
#: a single dataset, so a single applet is correct, and SHARED_APPLET_GROUP is
#: never retired. Routing it through `specs` instead would put it under the
#: running node's group for GVS and under the shared group for the optimizer
#: -- and applet identity is (name, group), so the dashboard would grow two
#: applets with the same title showing the same number.
#:
#: Not in SHARED_APPLET_SPECS itself, which the FORT optimizer and
#: SamplerMOTCoilAndBeamBalanceTune also create: neither loads atoms, so the
#: number would sit there stale.
SCAN_APPLET_SPECS = (
    _spec("Time without atom (s)", ("builtin", "big_number"),
          ("time_without_atom",), group=_GVS),
)

#: Two-node results. Loading is a JOINT measurement -- all four SPCMs see both
#: traps through the beamsplitter fan-out -- so these read the same unsuffixed
#: AllSPCMs_RO1/RO2 the sequence writes, and the discriminator is the JOINT
#: two_atom_threshold, not single_atom_threshold. That is the whole difference
#: from the single-node retention applet: "loaded" here means an atom in BOTH
#: traps, and "retained" means both survived.
#:
#: t_SPCM_first_shot is named WITH Node1's suffix on purpose. The bare name
#: resolves to a leftover STANDALONE dataset in two-node mode (it is 0.01 in
#: dataset_db today, which merely happens to match), and reading a stale
#: standalone value to scale a live threshold is the n_measurements trap over
#: again. Node1's copy is provably equal to Node2's: t_SPCM_first_shot is in
#: TWO_NODE_AGREEING_SCALARS, so prepare and every scan point RAISE if the two
#: ever diverge.
TWO_NODE_APPLET_SPECS = (
    _spec("retention and loading AllSPCMs", "plot_retention_and_loading.py",
          ("AllSPCMs_RO1", "AllSPCMs_RO2", "n_measurements",
           "two_atom_threshold", "t_SPCM_first_shot_Node1"),
          (("--scan_vars", "scan_variables"),
           ("--scan_sequence1", "scan_sequence1")),
          group=_GVS),
)

#: Applets that were defined but not enabled on Node1 (several were enabled
#: on Node2). Pass them explicitly when you want them:
#:     create_applets_for(self, self.base,
#:                        specs=APPLET_SPECS + OPTIONAL_APPLET_SPECS)
OPTIONAL_APPLET_SPECS = (
    # The unsuffixed "Microwaves Health Check" that used to live here is
    # superseded by the two per-node entries in SHARED_APPLET_SPECS below,
    # which are always created and show both nodes at once.
    #
    # n_excitation_cycles moved down here 2026-09-17: only the single-photon
    # work looks at it, so it should not come up on every run.
    _spec("n_excitation_cycles", "plot_xyline.py", ("n_excitation_cycles",),
          (("--pts", "applet_plot_points_short"),), group=_PHOTON),
    # "Time without atom (s)" moved to SCAN_APPLET_SPECS, which GVS and
    # MicrowaveScanOptimizer both ask for, so it is no longer opt-in.
    _spec("optimization cost", ("builtin", "plot_xy"), ("cost",), group=_OPT),
    # The unsuffixed "feedback RF" that used to live here is superseded by the
    # two per-node entries in SHARED_APPLET_SPECS below, which are always
    # created and show both nodes at once.
    # "waveplate trajectory" moved to K10CR1_APPLET_SPECS, which the FORT
    # polarization optimizer asks for: one per node, with the suffixed dataset
    # names those angles are actually written under.
    # "FORT MM (normalized and accounted for power)" moved to
    # SHARED_APPLET_SPECS on 2026-10-08 as "FORT MM over APD normalized
    # NodeX": same script and same arithmetic, but one applet per NODE with
    # suffixed dataset names, and always up alongside the other two FORT
    # monitors instead of opt-in under _K10. The unsuffixed version plotted a
    # stale standalone dataset in two-node mode.
)

#: The SPCM count-rate applets. Deliberately NOT in APPLET_SPECS: they would
#: otherwise be created on every scan, and MonitorSPCMinApplet is the one
#: experiment that exists to watch them, so it asks for them explicitly:
#:
#:     create_applets_for(self, self.base,
#:                        specs=APPLET_SPECS + SPCM_MONITOR_APPLET_SPECS)
#:
#: Titles carry the canonical detector names (SPCM_H1 = SPCM0, SPCM_V1 =
#: SPCM1, SPCM_H2 = SPCM0_OtherNode, SPCM_V2 = SPCM1_OtherNode) while the
#: dataset names stay LEGACY, because that is what the experiments write.
#: Not node-suffixed: the detectors are master-local, so both nodes write the
#: same five datasets. Commands copied verbatim from the dashboard backup.
SPCM_MONITOR_APPLET_SPECS = (
    _spec("SPCM_H1 count rate", "plot_xyline.py",
          ("SPCM0_counts_per_s",),
          (("--pts", "applet_plot_points_short"),), group=_OPT),
    _spec("SPCM_V1 count rate", "plot_xyline.py",
          ("SPCM1_counts_per_s",),
          (("--pts", "applet_plot_points_short"),), group=_OPT),
    _spec("SPCM_H2 count rate", "plot_xyline.py",
          ("SPCM0_OtherNode_counts_per_s",),
          (("--pts", "applet_plot_points_short"),), group=_OPT),
    _spec("SPCM_V2 count rate", "plot_xyline.py",
          ("SPCM1_OtherNode_counts_per_s",),
          (("--pts", "applet_plot_points_short"),), group=_OPT),
    _spec("All SPCMs count rate", "plot_xyline.py",
          ("AllSPCMs_counts_per_s",),
          (("--pts", "applet_plot_points_medium"),), group=_OPT),
)


#: Top-level group for applets that belong to NEITHER node. It is never one of
#: base.VALID_NODES, so create_applets_for's disable_applet_group pass leaves
#: it alone -- which is what keeps these up no matter which node is running.
SHARED_APPLET_GROUP = "Both nodes"

#: Top-level group for applets that only make sense during a TWO-NODE run.
TWO_NODE_APPLET_GROUP = "TwoNodes"


#: Asked for ONLY when the selected experiment_function uses the alternating
#: atom check. Deliberately not in TWO_NODE_APPLET_SPECS, which is created on
#: every two-node run: the datasets are seeded on every run, so a default
#: entry would open two empty histograms whenever the joint loader is used.
#: Same reasoning as GVS_APPLET_SPECS, which exists for exactly this.
#:
#: These ARE the calibration. load_until_atom_in_both_nodes_alternating_shot
#: _recycle appends two values per loaded round per node -- the count that
#: crossed and the first try's count that did not -- so each histogram is
#: bimodal and alternating_shot_threshold_for_loading_NodeX is read off the
#: valley between the two populations.
#:
#: --t_exposure is t_alternating_RO_exposure, which Base derives and
#: broadcasts. NOT t_atom_check_time: that is the joint gate's duration, and
#: in two-node mode it resolves to a bare name this stack never declares as a
#: dataset, so the applet would be answered by a leftover STANDALONE value of
#: 20 ms -- twice the real exposure.
#:
#: No --iteration: like the joint atom-check histogram, these accumulate over
#: the whole run and are never reset per scan point, so a per-iteration title
#: would be misleading.
TWO_NODE_ALTERNATING_APPLET_SPECS = (
    _spec("Node1 atom check in loading histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO_atom_check_in_loading_Node1",),
          (("--x", "photocount_bins"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
    _spec("Node2 atom check in loading histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO_atom_check_in_loading_Node2",),
          (("--x", "photocount_bins"),
           ("--color", "second_shot_hist_color_Node2"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
)

#: Two-node experiment_functions whose loader writes the per-node atom-check
#: datasets. Keyed by name because the applet set is chosen in prepare(),
#: where the selected function is already known.
#:
#: EMPTY as of 2026-10-08. The only member was
#: Two_nodes_alternating_shot_experiment, which has been swapped back to the
#: joint loader because load_until_atom_in_both_nodes_alternating_shot_recycle
#: does not load atoms (see its docstring). The joint loader never writes
#: AllSPCMs_alternating_RO_atom_check_in_loading_Node1/_Node2, so leaving the name here would
#: open two permanently empty histograms on every alternating run -- the exact
#: thing this gate exists to prevent. Put the name back when the per-node
#: loader works; nothing else needs changing.
TWO_NODE_ALTERNATING_EXPERIMENTS = ()

#: Asked for ONLY by the experiments whose SHOTS are node-resolved, for the
#: same reason as the tuple above: the joint experiments never write these, so
#: a default entry would open six permanently empty plots.
#:
#: FOUR HISTOGRAMS -- two shots x two nodes -- plus the two retention plots
#: they feed. The histograms are where alternating_shot_threshold_NodeX is read
#: off: each should be bimodal, with the threshold in the valley.
#:
#: --t_exposure is t_alternating_RO_exposure, which Base derives and broadcasts
#: as n_alternating_RO_windows_per_node * t_alternating_RO_window. NOT
#: t_SPCM_first_shot: that is the JOINT gate's duration, and these windows are
#: a different length, so it would scale every rate wrongly.
TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS = (
    _spec("Node1 alternating RO1 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO1_Node1_current_iteration",),
          (("--x", "photocount_bins"), ("--iteration", "iteration"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
    _spec("Node2 alternating RO1 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO1_Node2_current_iteration",),
          (("--x", "photocount_bins"), ("--iteration", "iteration"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
    _spec("Node1 alternating RO2 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO2_Node1_current_iteration",),
          (("--x", "photocount_bins"),
           ("--color", "second_shot_hist_color_Node1"),
           ("--iteration", "iteration"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
    _spec("Node2 alternating RO2 histogram", "plot_hist_autosize.py",
          ("AllSPCMs_alternating_RO2_Node2_current_iteration",),
          (("--x", "photocount_bins"),
           ("--color", "second_shot_hist_color_Node2"),
           ("--iteration", "iteration"),
           ("--t_exposure", "t_alternating_RO_exposure")),
          group=_GVS),
    # The point of the experiment: retention PER NODE, each against its own
    # threshold. Deliberately NOT titled "retention and loading AllSPCMs" --
    # that title is the one conflicting_applet_titles() retires across modes,
    # and these are additional plots, not another copy of it.
    _spec("Node1 alternating retention", "plot_retention_and_loading.py",
          ("AllSPCMs_alternating_RO1_Node1",
           "AllSPCMs_alternating_RO2_Node1",
           "n_measurements",
           "alternating_shot_threshold_Node1",
           "t_alternating_RO_exposure"),
          (("--scan_vars", "scan_variables"),
           ("--scan_sequence1", "scan_sequence1")),
          group=_GVS),
    _spec("Node2 alternating retention", "plot_retention_and_loading.py",
          ("AllSPCMs_alternating_RO1_Node2",
           "AllSPCMs_alternating_RO2_Node2",
           "n_measurements",
           "alternating_shot_threshold_Node2",
           "t_alternating_RO_exposure"),
          (("--scan_vars", "scan_variables"),
           ("--scan_sequence1", "scan_sequence1")),
          group=_GVS),
)

#: Two-node experiment_functions whose SHOTS are node-resolved.
TWO_NODE_PER_NODE_RETENTION_EXPERIMENTS = (
    "Two_nodes_alternating_retention_experiment",
)

#: Titles to CLOSE when those experiments run, because they write none of the
#: datasets behind them: the node-resolved shots never call first_shot or
#: second_shot, so AllSPCMs_RO1/RO2 and their _current_iteration histograms
#: stay at their seeded zero for the whole run. Left open, each keeps showing
#: the previous run's data -- a retention plot reading 0% and two histograms
#: frozen at the last joint run, sitting beside the six correct per-node plots.
TWO_NODE_PER_NODE_RETENTION_RETIRES = (
    "retention and loading AllSPCMs",
    "All SPCMs RO1 histogram",
    "All SPCMs RO2 histogram",
)


def execution_applet_groups(base):
    """Every top-level group that belongs to ONE execution mode."""
    return tuple(base.VALID_NODES) + (TWO_NODE_APPLET_GROUP,)


def conflicting_applet_titles():
    """Titles that exist in BOTH spec sets meaning DIFFERENT things.

    Derived, not listed: the intersection of the single-node and two-node spec
    titles. Today that is "retention and loading AllSPCMs", which appears in
    each set with a different threshold -- single_atom_threshold_NodeX for one
    atom, two_atom_threshold for both. Same title, same unsuffixed
    AllSPCMs_RO1/RO2, different criterion.

    Only THESE need retiring when the mode changes. The rest of the idle
    group's applets stay up on purpose: readout histograms and the atom
    loading time are just as meaningful during a two-node run, and the owner
    wants to keep watching them. Closing a whole group to fix one applet
    throws those away.
    """
    single = {spec.title for spec in APPLET_SPECS}
    two_node = {spec.title for spec in TWO_NODE_APPLET_SPECS}
    return single & two_node


def applet_paths_for_titles(titles):
    """(title, functional group) for each named title, taken from the SPECS.

    The functional group is what a bare title cannot carry and what every
    retirement below cannot do without: ccb_disable_applet resolves its group
    argument with locate_applet(create_groups=False), which looks for the
    applet as a DIRECT CHILD of the last group named and silently does nothing
    when it is one level deeper. Read from the spec rather than assumed, so a
    title that later moves functional group is still reached.

    EVERY spec set is searched, and a title present in more than one yields
    one entry per functional group it occupies -- so a disagreement about
    where it lives becomes two retirements rather than one silent miss.
    """
    wanted = set(titles)
    return {
        (spec.title, spec.group)
        for spec in (APPLET_SPECS + TWO_NODE_APPLET_SPECS + GVS_APPLET_SPECS
                     + TWO_NODE_ALTERNATING_APPLET_SPECS
                     + TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS)
        if spec.title in wanted
    }


def conflicting_applet_paths():
    """(title, functional group) for each applet that means two things."""
    return applet_paths_for_titles(conflicting_applet_titles())


#: Spec sets whose datasets are written by exactly ONE experiment_function.
#: Any run that is not that experiment must close them -- a plot of a dataset
#: nobody is writing does not go blank, it keeps showing the previous run's
#: data, which is indistinguishable from a result.
#:
#: GVS_APPLET_SPECS is deliberately NOT here. Its atom-check histogram plots
#: AllSPCMs_atom_check_in_loading, which the two-node JOINT loader writes too,
#: so a two-node run must leave it alone rather than close it.
EXCLUSIVE_APPLET_SPEC_SETS = (
    TWO_NODE_ALTERNATING_APPLET_SPECS,
    TWO_NODE_PER_NODE_RETENTION_APPLET_SPECS,
)


def exclusive_applet_titles(specs_in_use):
    """Titles of single-experiment applets this run is NOT creating.

    Derived from EXCLUSIVE_APPLET_SPEC_SETS rather than listed, so a new
    opt-in applet is covered the moment it joins one of those tuples -- which
    is the failure this closes: the six per-node retention plots have titles
    that appear in NEITHER spec-set intersection, so conflicting_applet_titles
    could never see them and a single-node run left all six open, still
    plotting the last two-node run's AllSPCMs_alternating_RO1/RO2_NodeX.
    """
    in_use = {spec.title for spec in specs_in_use}
    exclusive = set()
    for spec_set in EXCLUSIVE_APPLET_SPEC_SETS:
        exclusive |= {spec.title for spec in spec_set}
    return exclusive - in_use


def _retire_conflicting_applets(ccb, live_group, base):
    """Close the other mode's copy of each conflicting applet, and only that.

    A leftover copy does not go blank, which is what makes it dangerous: it
    keeps plotting live data under the wrong criterion. A single-node copy
    left open during a two-node run judges joint counts against
    single_atom_threshold and reads high; a two-node copy left open during a
    single-node run judges one-atom counts against two_atom_threshold and
    reads ~0. Either sits beside the correct plot looking plausible.

    Per-applet rather than per-group because disable_applet_group would also
    retire the histograms, the loading-time plot and everything else in the
    idle node's group, which stay useful in either mode.

    THE GROUP MUST BE THE FULL NESTED PATH, and getting this wrong is silent.
    create_applet places an applet at [top_level, functional_group], so
    disable_applet has to name both: ccb_disable_applet resolves its group
    argument with locate_applet(..., create_groups=False), which looks for the
    applet as a DIRECT CHILD of the last group named and returns None -- doing
    nothing at all -- when it is one level deeper. Passing only the top-level
    group made this function a no-op for its entire existence. It was caught on
    2026-10-11 from the live dashboard state, where the Node1 AND Node2
    retention applets were both enabled at once, each plotting the same bare
    AllSPCMs_RO1/RO2 against a different threshold.
    """
    for group in execution_applet_groups(base):
        if group == live_group:
            continue
        for title, functional_group in conflicting_applet_paths():
            path = [group] if functional_group is None else [
                group, functional_group
            ]
            ccb.issue("disable_applet", title, path)


def _retire_named_applets(ccb, base, titles):
    """Close specific applets in EVERY execution group, the live one included.

    For an experiment that writes none of their datasets. A plot whose dataset
    nobody is writing does not go blank -- it keeps showing the previous run's
    data, which is indistinguishable from a result, and is how a 0% retention
    or a stale histogram ends up sitting beside the correct plots.

    Unlike _retire_conflicting_applets this does NOT skip the live group,
    because the reason is the experiment rather than the mode: the joint
    AllSPCMs_RO1/RO2 plots are equally meaningless under TwoNodes as under
    Node1 when the shots are node-resolved.
    """
    for group in execution_applet_groups(base):
        for title, functional_group in applet_paths_for_titles(titles):
            path = [group] if functional_group is None else [
                group, functional_group
            ]
            ccb.issue("disable_applet", title, path)

_HEALTH_CHECK_DATASETS = (
    "health_check_uw_freq00",
    "health_check_uw_freq01",
    "health_check_uw_freq11",
    "health_check_uw_freqm10",
    "health_check_uw_freqm11",
)

_FEEDBACK_RF_ARGS = ("p_AOM_A1_history", "MOT_beam_monitor_points")
_FEEDBACK_RF_OPTIONS = (
    ("--y2", "p_AOM_A2_history"),
    ("--y3", "p_AOM_A3_history"),
    ("--y4", "p_AOM_A4_history"),
    ("--y5", "p_AOM_A5_history"),
    ("--y6", "p_AOM_A6_history"),
    ("--y7", "p_FORT_loading_history"),
    ("--labels", "feedbackchannels"),
)


def _per_node_specs(title, script, args, options=(), group=None):
    """One spec per node, with every dataset name pinned to that node.

    Every entry in ``args`` and every VALUE in ``options`` must be a dataset
    name rather than a literal, since each simply gets the node suffix
    appended.
    """
    return tuple(
        _spec(
            f"{title} {node}",
            script,
            tuple(f"{name}_{node}" for name in args),
            tuple((flag, f"{value}_{node}") for flag, value in options),
            group=group,
        )
        for node in ("Node1", "Node2")
    )


#: Applets that show BOTH nodes and stay up regardless of which one is
#: running. Unlike every spec above, their dataset names are written out
#: node-suffixed rather than left legacy-unsuffixed: resolve_applet_dataset
#: would otherwise rewrite them to whichever node happens to be running, and
#: the whole point here is to pin one applet per node. Suffixed names survive
#: that resolution in both directions -- the running node's are returned
#: unchanged, and the other node's raise ValueError and fall through to
#: resolve_result_dataset_name, which passes them through untouched.
#:
#: Both families keep showing the idle node's data rather than going blank:
#: the microwave fidelities are persistent per-node ExperimentVariables, and
#: the feedback histories are per-node result datasets. Note the suffix goes
#: on the END of the whole legacy name -- p_AOM_A1_history_Node1, not
#: p_AOM_A1_Node1_history -- which is what feedback_dataset_map produces.
#: Asked for by FORT_Polarization_Optimizer_master_satellite, as shared_specs:
#:
#:     create_shared_applets_for(
#:         self, self.base,
#:         shared_specs=SHARED_APPLET_SPECS + K10CR1_APPLET_SPECS)
#:
#: One per node, because HWP_angle and QWP_angle are in
#: POLARIZATION_RESULT_DATASETS and so are written node-suffixed: there really
#: are two trajectories, and a single unsuffixed applet could only ever show
#: whichever node ran last.
#:
#: Not in SHARED_APPLET_SPECS itself, which GVS, both optimizers and
#: SamplerMOT all create: this is the only experiment that moves waveplates
#: and writes those datasets, so anywhere else the plot would sit stale.
K10CR1_APPLET_SPECS = _per_node_specs(
    "waveplate trajectory",
    ("builtin", "plot_xy"),
    ("HWP_angle",),
    (("--x", "QWP_angle"),),
    group=_K10,
)

SHARED_APPLET_SPECS = (
    _per_node_specs(
        "Microwaves Health Check",
        "bar_plot_microwaves_health_check.py",
        _HEALTH_CHECK_DATASETS,
    )
    + _per_node_specs(
        "feedback RF",
        "plot_xy_multichannel.py",
        _FEEDBACK_RF_ARGS,
        _FEEDBACK_RF_OPTIONS,
    )
    # The FORT monitors are shared for the same reason: the FORT is used
    # everywhere, so both nodes' traces should stay on screen whichever node
    # is running. Unlike the two families above, the plotted y-series
    # (FORT_APD_monitor_NodeX, FORT_MM_monitor_NodeX) are broadcast-only
    # transients rather than persistent datasets -- they live in the running
    # master's memory, so the APPLET persists across runs but its data starts
    # empty after a master restart until feedback republishes it. The --y0
    # reference values are persistent per-node variables.
    + _per_node_specs(
        "FORT APD normalized setpoint",
        "plot_xyline_relative_y.py",
        ("FORT_APD_monitor",),
        (("--y0", "set_point_FORT_APD_loading"),),
    )
    + _per_node_specs(
        "FORT MM normalized to ref",
        "plot_xyline_relative_y.py",
        ("FORT_MM_monitor",),
        (("--y0", "best_852_power_ref"),),
    )
    # The ratio of the two applets above: (MM / best_852_power_ref) divided by
    # (APD / set_point_FORT_APD_loading). It divides out the FORT power the APD
    # actually saw, so what is left is mode-matching alone -- a drop here is
    # the fiber mode degrading rather than the 852 laser drifting, which the
    # two separate traces cannot distinguish on their own.
    #
    # plot_FORT_MM_norm_to_power.py already does exactly this arithmetic
    # (y /= y0; y1 /= y2; y = y/y1), so there is no new applet script. What was
    # missing is that it was registered ONCE with legacy-unsuffixed names in
    # the opt-in _K10 group, so it never came up per node -- and in two-node
    # mode a bare per-node name resolves ambiguously and falls through to a
    # stale STANDALONE dataset. That entry is gone; this is its replacement.
    #
    # Do NOT add --pts here. The script truncates x and y to the last pts
    # points but not y1, and the zip then pairs the LAST pts MM values against
    # the FIRST pts APD values. None of the FORT specs pass it, so the bug is
    # unreachable; passing it would silently plot nonsense.
    + _per_node_specs(
        "FORT MM over APD normalized",
        "plot_FORT_MM_norm_to_power.py",
        ("FORT_MM_monitor",),
        (("--y0", "best_852_power_ref"),
         ("--y1", "FORT_APD_monitor"),
         ("--y2", "set_point_FORT_APD_loading")),
    )
)


def resolve_applet_dataset(base, name):
    """Map a legacy dataset name to what the running node actually writes."""
    try:
        # Globals return unchanged; per-node variables gain the node suffix.
        return base.resolve_experiment_variable_target(name)
    except ValueError:
        # Not an experiment variable: a result dataset. The redirect suffixes
        # the node-specific ones and passes everything else through.
        return base.resolve_result_dataset_name(name)


def _resolve_option_value(base, value):
    """Options may carry a literal (e.g. a point count) instead of a name."""
    if not isinstance(value, str) or not value or value[0].isdigit():
        return str(value)
    return resolve_applet_dataset(base, value)


def build_applet_command(spec, base):
    """Return the full command line for one applet, resolved for this node."""
    if isinstance(spec.script, tuple):
        # ${artiq_applet} is the dashboard's own substitution for the
        # "<interpreter> -m artiq.applets." prefix; this is the form the
        # working dashboard configuration used.
        parts = ["${artiq_applet}" + spec.script[1]]
    else:
        applet_path = APPLET_DIRECTORY / spec.script
        if not applet_path.exists():
            raise FileNotFoundError(
                f"Applet {spec.script!r} for {spec.title!r} is missing from "
                f"{APPLET_DIRECTORY}."
            )
        parts = ['python "{}"'.format(applet_path)]

    parts += [resolve_applet_dataset(base, name) for name in spec.args]
    for flag, value in spec.options:
        parts += [flag, _resolve_option_value(base, value)]
    return " ".join(parts)


def applet_group_for(base):
    """Top-level dashboard group for the currently configured execution."""
    if base.experiment_mode == "single_node":
        return base.which_node
    return TWO_NODE_APPLET_GROUP


def _dashboard_ccb(experiment):
    """Fetch the dashboard-supplied CCB virtual device.

    Experiments need not bind it, so fetch rather than assume an attribute.
    """
    ccb = getattr(experiment, "ccb", None)
    if ccb is None:
        ccb = experiment.get_device("ccb")
    return ccb


def create_shared_applets_for(experiment, base,
                              shared_specs=SHARED_APPLET_SPECS):
    """Create the node-independent applets and nothing else.

    These live under SHARED_APPLET_GROUP, which is never one of
    base.VALID_NODES, so no disable_applet_group pass can retire them.

    Use this instead of create_applets_for from an experiment that should
    keep the shared applets up WITHOUT taking over the per-node applet set
    and without retiring the idle node's group. Passing specs=() to
    create_applets_for would not do: it still runs the disable pass.

    Works in either execution mode. The shared specs name their datasets
    node-suffixed, and suffixed names resolve unchanged in two_nodes too.

    Returns the (title, command) pairs issued.
    """
    ccb = _dashboard_ccb(experiment)
    issued = []
    for spec in shared_specs:
        command = build_applet_command(spec, base)
        group = (
            [SHARED_APPLET_GROUP] if spec.group is None
            else [SHARED_APPLET_GROUP, spec.group]
        )
        ccb.issue("create_applet", spec.title, command, group=group)
        issued.append((spec.title, command))
    return issued


def create_two_node_applets_for(experiment, base,
                                specs=TWO_NODE_APPLET_SPECS,
                                shared_specs=SHARED_APPLET_SPECS,
                                retire_titles=()):
    """Create the two-node applets, the shared ones, and retire both nodes.

    The two-node counterpart of create_applets_for. It exists rather than
    being a mode branch inside that function because the two SPEC SETS are
    different -- APPLET_SPECS names per-node result datasets that two-node
    mode does not suffix -- but the RETIREMENT rule is identical and shared.

    Retiring Node1 and Node2 is the point. "retention and loading AllSPCMs"
    exists in both spec sets under the same title, reading the same
    unsuffixed AllSPCMs_RO1/RO2, and differing only in the threshold. Without
    this pass a single-node run's copy stays open through a two-node run and
    plots a confidently wrong retention, judging two-atom counts against
    single_atom_threshold, right beside the correct one.
    """
    if base.experiment_mode == "single_node":
        raise NotImplementedError(
            "create_two_node_applets_for is for two-node runs; single-node "
            "runs go through create_applets_for."
        )

    ccb = _dashboard_ccb(experiment)
    group_name = applet_group_for(base)
    issued = []
    for spec in specs:
        command = build_applet_command(spec, base)
        group = ([group_name] if spec.group is None
                 else [group_name, spec.group])
        ccb.issue("create_applet", spec.title, command, group=group)
        issued.append((spec.title, command))

    issued.extend(create_shared_applets_for(experiment, base, shared_specs))

    _retire_conflicting_applets(ccb, group_name, base)
    # Titles this particular experiment writes nothing behind. Retired in every
    # group, live one included -- see _retire_named_applets. Two sources: the
    # caller's explicit list (the per-node retention experiment closing the
    # JOINT plots) and the derived one (any run closing the single-experiment
    # plots it is not creating, which is the reverse direction).
    _retire_named_applets(
        ccb, base, set(retire_titles) | exclusive_applet_titles(specs)
    )
    return issued


def create_applets_for(experiment, base, specs=APPLET_SPECS,
                       shared_specs=SHARED_APPLET_SPECS):
    """Create this node's applets and retire the other node's group.

    ``specs`` are per-node: they are created under the running node's group
    and the other node's group is disabled. ``shared_specs`` are created
    under SHARED_APPLET_GROUP instead, which is never disabled, so they stay
    up whichever node runs. Pass ``shared_specs=()`` to skip them.

    Returns the (title, command) pairs issued. Safe to call on every run:
    create_applet replaces the spec of an existing applet with the same name
    in the same group and restarts it, rather than duplicating it.
    """
    if base.experiment_mode != "single_node":
        raise NotImplementedError(
            "create_applets_for is single-node only: APPLET_SPECS names "
            "PER-NODE result datasets, and in two-node mode "
            "resolve_result_dataset_name does not suffix at all, so both "
            "nodes write the same names and a per-node applet has nothing to "
            "point at. Two-node runs go through create_shared_applets_for "
            "with SHARED_APPLET_SPECS + SCAN_APPLET_SPECS + "
            "TWO_NODE_APPLET_SPECS, whose datasets resolve unchanged in "
            "either mode."
        )

    ccb = _dashboard_ccb(experiment)

    node_group = applet_group_for(base)
    issued = []
    for spec in specs:
        command = build_applet_command(spec, base)
        group = [node_group] if spec.group is None else [node_group, spec.group]
        ccb.issue("create_applet", spec.title, command, group=group)
        issued.append((spec.title, command))

    # Node-independent applets go in their own top-level group, so the
    # disable pass below cannot retire them along with the idle node.
    issued.extend(create_shared_applets_for(experiment, base, shared_specs))

    # Close ONLY the other mode's conflicting applets -- the retention plot
    # that shares this one's title and reads the same unsuffixed
    # AllSPCMs_RO1/RO2 under a different threshold. Everything else in the
    # idle group stays up; histograms and loading time are useful either way.
    _retire_conflicting_applets(ccb, node_group, base)
    # ... and the single-experiment applets this run writes nothing behind. A
    # single-node run writes no AllSPCMs_alternating_RO1/RO2_NodeX, so the six
    # per-node retention plots would otherwise stay up showing the last
    # two-node run's data.
    _retire_named_applets(ccb, base, exclusive_applet_titles(specs))
    return issued
