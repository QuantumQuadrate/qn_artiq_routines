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

# The nine applets Node1 actually had ENABLED, transcribed from its live
# dashboard state (AppData/Local/m-labs/artiq/7), plus the four per-SPCM
# count-rate applets restored 2026-09-16 for MonitorSPCMinApplet -- they were
# defined in that dashboard state but not enabled, and their commands here are
# copied from it verbatim. Keeping the default set small is deliberate: 65
# applets were defined there but only nine running, and applet count is what
# makes the dashboard struggle. Move anything rarely watched down to
# OPTIONAL_APPLET_SPECS rather than growing this list further.
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

    # The four master-local SPCMs, written by MonitorSPCMinApplet. The dataset
    # names stay in their LEGACY form because that is what the experiments
    # actually write; only the titles carry the canonical detector names
    # (SPCM_H1 = SPCM0, SPCM_V1 = SPCM1, SPCM_H2 = SPCM0_OtherNode,
    # SPCM_V2 = SPCM1_OtherNode). These are not node-suffixed: the detectors
    # are master-local, so both nodes write the same five datasets.
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

    _spec("Atom loading time (s)", "plot_xy_multichannel.py",
          ("Atom_loading_time", "n_measurements"), group=_MON),

    _spec("n_excitation_cycles", "plot_xyline.py", ("n_excitation_cycles",),
          (("--pts", "applet_plot_points_short"),), group=_PHOTON),

    _spec("FORT APD normalized setpoint", "plot_xyline_relative_y.py",
          ("FORT_APD_monitor",),
          (("--y0", "set_point_FORT_APD_loading"),), group=_K10),
    _spec("FORT MM normalized to ref", "plot_xyline_relative_y.py",
          ("FORT_MM_monitor",),
          (("--y0", "best_852_power_ref"),), group=_K10),
)

#: Applets that were defined but not enabled on Node1 (several were enabled
#: on Node2). Pass them explicitly when you want them:
#:     create_applets_for(self, self.base,
#:                        specs=APPLET_SPECS + OPTIONAL_APPLET_SPECS)
OPTIONAL_APPLET_SPECS = (
    _spec("Microwaves Health Check", "bar_plot_microwaves_health_check.py",
          ("health_check_uw_freq00", "health_check_uw_freq01",
           "health_check_uw_freq11", "health_check_uw_freqm10",
           "health_check_uw_freqm11")),
    _spec("Time without atom (s)", ("builtin", "big_number"),
          ("time_without_atom",), group=_GVS),
    _spec("optimization cost", ("builtin", "plot_xy"), ("cost",), group=_OPT),
    _spec("feedback RF", "plot_xy_multichannel.py",
          ("p_AOM_A1_history", "MOT_beam_monitor_points"),
          (("--y2", "p_AOM_A2_history"), ("--y3", "p_AOM_A3_history"),
           ("--y4", "p_AOM_A4_history"), ("--y5", "p_AOM_A5_history"),
           ("--y6", "p_AOM_A6_history"), ("--y7", "p_FORT_loading_history"),
           ("--labels", "feedbackchannels")),
          group=_MON),
    _spec("waveplate trajectory", ("builtin", "plot_xy"), ("HWP_angle",),
          (("--x", "QWP_angle"),), group=_K10),
    _spec("FORT MM (normalized and accounted for power)",
          "plot_FORT_MM_norm_to_power.py", ("FORT_MM_monitor",),
          (("--y0", "best_852_power_ref"), ("--y1", "FORT_APD_monitor"),
           ("--y2", "set_point_FORT_APD_loading")), group=_K10),
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
    return "TwoNodes"


def create_applets_for(experiment, base, specs=APPLET_SPECS):
    """Create this node's applets and retire the other node's group.

    Returns the (title, command) pairs issued. Safe to call on every run:
    create_applet replaces the spec of an existing applet with the same name
    in the same group and restarts it, rather than duplicating it.
    """
    if base.experiment_mode != "single_node":
        raise NotImplementedError(
            "Applet creation is single-node only. In two-node mode "
            "resolve_result_dataset_name does not suffix at all, so both "
            "nodes write the same result datasets and there is nothing for "
            "per-node applets to point at. Decide how two-node results are "
            "stored first. (The old dashboard's two-node retention applet "
            "differed only by using two_atom_threshold.)"
        )

    # "ccb" is a virtual device supplied by the dashboard/worker; experiments
    # need not bind it, so fetch it rather than assuming an attribute.
    ccb = getattr(experiment, "ccb", None)
    if ccb is None:
        ccb = experiment.get_device("ccb")

    node_group = applet_group_for(base)
    issued = []
    for spec in specs:
        command = build_applet_command(spec, base)
        group = [node_group] if spec.group is None else [node_group, spec.group]
        ccb.issue("create_applet", spec.title, command, group=group)
        issued.append((spec.title, command))

    # Only one node runs at a time, so retire the other node's applets rather
    # than leaving them subscribed to datasets nothing is updating.
    for other in base.VALID_NODES:
        if other != node_group:
            ccb.issue("disable_applet_group", other)
    return issued
