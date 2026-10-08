from collections import namedtuple

from artiq.experiment import *


Variable = namedtuple("Variable", ["name", "value", "value_type", "kwargs", "group"])


MASTER_SATELLITE_VARIABLES = (
    Variable(
        "n_measurements",
        100,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "general",
    ),
    # Cross-node timing comes in TWO independent layers. Conflating them is
    # what made the pre-DRTIO code hard to reason about, so they are two
    # variables with two names and two units of meaning.
    #
    # LAYER 1, electrical. Do the two crates' output PINS move together for
    # one scheduled timestamp? Everything with a gate wider than a microsecond
    # -- atom loading, the alternating readout -- cares only about this.
    # Measured by tests/measure_drtio_ttl_skew.py; expected to be ~0 because
    # DRTIO buffers events at the satellite ahead of time (measured
    # 2026-09-10, RID 38429: no per-event penalty, identical 1000 mu minimum
    # lead on both destinations), but it must be a MEASURED zero.
    #
    # SIGN CONVENTION, and it matters: this is POSITIVE when Node2's output pin
    # moves LATER than Node1's for the same scheduled timestamp -- i.e. set it
    # to the T2 output skew that measure_drtio_ttl_skew.py reports, with its
    # sign unchanged. The sequence therefore schedules Node2's events EARLIER
    # by this amount, as at_mu(t - t_Node2_rtio_offset_mu), to cancel it.
    # Getting the sign backwards doubles the skew instead of removing it, and
    # the only symptom is one node's readout light leaking into the other's
    # window -- which reads as a physics problem, not a timing one. At the
    # default 0 the sign is unobservable, so it is written down here rather
    # than left to be rediscovered.
    #
    # This replaces alice_extra_offset_mu (subroutines/experiment_functions.py
    # :14308, ":14382", and the inline +200 at ":14280", all "Calibrated this
    # on the scope"). That constant characterised the old two-crate TTL
    # handshake path and carries no meaning now that both nodes share one
    # timeline, so the concept survives the migration but the value does not.
    Variable(
        "t_Node2_rtio_offset_mu",
        0,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "master-satellite timing",
    ),
    # LAYER 2, optical. Firing both excitation pulses at the same now_mu does
    # NOT make the light reach the two atoms together, nor the two photons
    # reach the detector together: fibre and free-space path lengths differ,
    # AOM turn-on differs, and each atom's emission is referenced to its own
    # local pulse. This is the fine delay that brings the two photons into
    # coincidence AT THE SPCM.
    #
    # Renamed from t_delay_in_bob_mu, which collided: the same dataset name
    # was declared 189 here and 0 in standalone/ExperimentVariables.py:77, so
    # whichever initializer ran last won and the live value was 0. That is the
    # same hazard documented for n_measurements in
    # utilities/BaseExperiment_master_satellite.py. The rename fixes the
    # collision and the name together.
    #
    # Not measurable by TTL loopback -- it needs a photon coincidence
    # measurement (HOM dip / cross-correlation), which is the deferred herald
    # work. Read today only by the single-photon experiments in
    # subroutines/experiment_functions.py (":15291", ":15559", ":15681",
    # ":16011"), which are reachable from neither master-satellite registry;
    # they also hardcode 189 as the old zero-point, so when they are ported
    # this variable should mean the delay directly.
    Variable(
        "t_Node2_excitation_delay_mu",
        189,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "two-node optical timing",
    ),
    Variable(
        "parallel_AOM_feedback",
        True,
        BooleanValue,
        {},
        "master-satellite feedback",
    ),
    # A scratch scan target with no physical meaning, let alone a per-node one:
    # the point of it is to scan SOMETHING without having to declare a new
    # variable first. It was per-node until 2026-10-07, which made
    # scan_variable1_name = "dummy_variable" fail in two_nodes mode with
    # "Ambiguous node-specific experiment-variable target" -- asking which
    # node's dummy you meant, a question with no answer. Both copies held 0.0
    # and every reader wants the bare name (subroutines/experiment_functions.py
    # reads self.dummy_variable).
    #
    # Shares its name with standalone/ExperimentVariables.py, deliberately and
    # harmlessly: both default to 0.0, nothing ever writes the dataset (it is
    # only ever read, or overridden for one run), so the same scan string works
    # in either stack. That is unlike n_measurements, where the shared name
    # carries a real value AND the standalone base rewrites it without
    # persist=True.
    Variable(
        "dummy_variable",
        0.0,
        NumberValue,
        {"type": "float"},
        "debugging",
    ),
    # Two-atom thresholds are GLOBAL on purpose: they describe the joint
    # two-node readout, not one node's detector, so both nodes must agree on
    # them. They were per-node until 2026-09-17 and carried identical values
    # on both, which is the duplication this removes. Being global also makes
    # the bare name exist in two_nodes mode -- the per-node declarations only
    # ever produced two_atom_threshold_NodeX, while the two-node experiment
    # functions read self.two_atom_threshold.
    Variable(
        "two_atom_threshold",
        39000.0,
        NumberValue,
        {"type": "float"},
        "Thresholds and cut-offs",
    ),
    Variable(
        "two_atom_threshold_for_loading",
        89000.0,
        NumberValue,
        {"type": "float"},
        "Thresholds and cut-offs",
    ),
    # Joint two-node timing. These are GLOBAL for the same reason the two-atom
    # thresholds are: a two-node run opens ONE gate across all four
    # master-local SPCM counters, so there is one duration, not one per node.
    # The thresholds above are RATES, compared against counts/t, so the gate
    # duration and the threshold have to be the same joint quantity or the
    # load criterion is silently wrong by the ratio of the two nodes'
    # settings. Per-node durations stay per-node (coil volts, PGC time, FORT
    # drop); only what the joint gate spans lives here.
    Variable(
        "t_atom_check_time_two_node",
        0.01,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 3},
        "two-node joint timing",
    ),
    Variable(
        "t_SPCM_first_shot_two_node",
        0.01,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 3},
        "two-node joint timing",
    ),
    Variable(
        "t_SPCM_second_shot_two_node",
        0.01,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 3},
        "two-node joint timing",
    ),
    Variable(
        "t_delay_between_shots_two_node",
        0.0,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 3},
        "two-node joint timing",
    ),
    Variable(
        "t_MOT_dissipation_two_node",
        0.05,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 3},
        "two-node joint timing",
    ),
    # Replaces the literal 100 at subroutines/experiment_functions.py:14793.
    Variable(
        "max_atom_check_tries_two_node",
        100,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "two-node joint timing",
    ),
    # NEW, and not optional. The legacy loader's outer `while True` has no
    # bound, and its try_n counter is reset ONLY inside the laser-feedback
    # branch (":14882"). Two-node mode runs open-loop, so without a bound the
    # loader would spin forever once the inner loop stopped checking. On
    # reaching this many rounds the loader gives up, records time_without_atom
    # and returns without an atom, so the caller skips a measurement instead
    # of hanging the dashboard.
    Variable(
        "max_loading_rounds_two_node",
        20,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "two-node joint timing",
    ),
    # Alternating-readout geometry. Replaces the literals at
    # subroutines/experiment_functions.py:14983 and ":14985", and the three
    # repeated 0.1 ms pads.
    Variable(
        "n_alternating_RO_windows_per_node",
        10,
        NumberValue,
        {"type": "int", "ndecimals": 0, "step": 1, "scale": 1},
        "two-node alternating readout",
    ),
    # SECONDS, like every other time variable in this stack -- the unit kwarg
    # is display scaling only. 0.001 is the legacy 1.0 * ms readout window and
    # 0.0001 is its 0.1 ms settle pad. Getting this wrong by 1000x would make
    # one "window" a full second.
    Variable(
        "t_alternating_RO_window",
        0.001,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 4},
        "two-node alternating readout",
    ),
    Variable(
        "t_alternating_RO_pad",
        0.0001,
        NumberValue,
        {"type": "float", "unit": "ms", "ndecimals": 4},
        "two-node alternating readout",
    ),
)


class ExperimentVariablesMasterSatelliteGlobal(EnvExperiment):
    """ExperimentVariables_master_satellite_global

    Initialize or intentionally update shared master-satellite variables.
    """

    def build(self):
        self.vars_list = list(MASTER_SATELLITE_VARIABLES)
        for variable in self.vars_list:
            try:
                value = self.get_dataset(variable.name)
            except KeyError:
                value = variable.value
            self.setattr_argument(
                variable.name,
                variable.value_type(value, **variable.kwargs),
                variable.group,
            )

    def run(self):
        for variable in self.vars_list:
            self.set_dataset(
                variable.name,
                getattr(self, variable.name),
                broadcast=True,
                persist=True,
            )
