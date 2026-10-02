"""Hardware-free tests for AtomLoadingOptimizer_load_until_atom_master_satellite.

The physics sequence is byte-for-byte the standalone one, so these tests cover
only what the port changed: node selection, the two GUI arguments that shadow
persistent variables, and the authoritative names the optimizer persists at
finish. That last one is the reason this file exists -- see
test_persisted_names_resolve_to_the_selected_node.
"""
import json
import sys
import types
import unittest
from pathlib import Path

import numpy as np


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


class _StubParallel:
    def __enter__(self):
        return self

    def __exit__(self, *exception):
        return False


_STUB_MARKER = "_qn_hardware_free_stub"


def _artiq_experiment_stub_or_none():
    """Return the artiq.experiment stub to (re)configure, or None.

    The ARTIQ repository scan imports this file inside a live worker, where a
    real artiq installation must never be modified or shadowed. A stub is
    installed (or extended) only when artiq is genuinely absent - the
    hardware-free test environment - or when an earlier hardware-free test
    module already installed the marked stub.
    """
    existing_artiq = sys.modules.get("artiq")
    if existing_artiq is not None:
        if getattr(existing_artiq, _STUB_MARKER, False):
            return sys.modules["artiq.experiment"]
        return None
    try:
        import artiq  # noqa: F401
    except ImportError:
        artiq_module = types.ModuleType("artiq")
        setattr(artiq_module, _STUB_MARKER, True)
        experiment_module = types.ModuleType("artiq.experiment")
        artiq_module.experiment = experiment_module
        sys.modules["artiq"] = artiq_module
        sys.modules["artiq.experiment"] = experiment_module
        return experiment_module
    return None


_experiment_stub = _artiq_experiment_stub_or_none()
if _experiment_stub is not None:
    _stub_exports = {
        "EnvExperiment": object,
        "kernel": _identity_decorator,
        "rpc": _identity_decorator,
        "delay": lambda duration: None,
        "delay_mu": lambda duration: None,
        "now_mu": lambda: 0,
        "parallel": _StubParallel(),
        "sequential": _StubParallel(),
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value=False, **kwargs: value,
        "StringValue": lambda value, **kwargs: value,
        "EnumerationValue": lambda values, **kwargs: tuple(values)[0],
        "TBool": bool,
        "TFloat": float,
        "TInt32": int,
        "TInt64": int,
        "TStr": str,
        "TArray": lambda *args, **kwargs: list,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
        "MHz": 1e6,
        "kHz": 1e3,
        "dB": 1.0,
        "V": 1.0,
        "s": 1.0,
    }
    for _name, _value in _stub_exports.items():
        if not hasattr(_experiment_stub, _name):
            setattr(_experiment_stub, _name, _value)
    if hasattr(_experiment_stub, "__all__"):
        _experiment_stub.__all__ = sorted(
            set(_experiment_stub.__all__) | set(_stub_exports)
        )

    # The optimizer now imports run_feedback_and_record_FORT_MM_power from
    # subroutines.experiment_functions, which pulls artiq.coredevice.ad9910
    # and .urukul at module level. Same block the GVS and microwave optimizer
    # test modules already carry; setdefault so whichever test module loads
    # first wins and the others extend rather than clobber.
    _coredevice = types.ModuleType("artiq.coredevice")
    _exceptions = types.ModuleType("artiq.coredevice.exceptions")

    class _RTIOUnderflow(Exception):
        pass

    _exceptions.RTIOUnderflow = _RTIOUnderflow
    _ad9910 = types.ModuleType("artiq.coredevice.ad9910")
    for _constant in (
        "PHASE_MODE_ABSOLUTE",
        "PHASE_MODE_CONTINUOUS",
        "PHASE_MODE_TRACKING",
        "RAM_DEST_ASF",
        "RAM_MODE_RAMPUP",
    ):
        setattr(_ad9910, _constant, 0)
    _urukul = types.ModuleType("artiq.coredevice.urukul")
    _urukul.CFG_MASK_NU = 0
    _language = types.ModuleType("artiq.language")
    _language.us = 1e-6
    _language.ns = 1e-9
    _language.MHz = 1e6
    sys.modules.setdefault("artiq.coredevice", _coredevice)
    sys.modules.setdefault("artiq.coredevice.exceptions", _exceptions)
    sys.modules.setdefault("artiq.coredevice.ad9910", _ad9910)
    sys.modules.setdefault("artiq.coredevice.urukul", _urukul)
    sys.modules.setdefault("artiq.language", _language)
    sys.modules.setdefault("pyvisa", types.ModuleType("pyvisa"))


def _install_mloop_stub_if_absent():
    """Stub M-LOOP only when it is genuinely missing.

    Same rule as the artiq stub above and for the same reason: the ARTIQ
    repository scan imports this file in a live worker, where the real mloop
    must not be shadowed. The real package pulls in TensorFlow, so the
    hardware-free environment deliberately does without it.
    """
    try:
        import mloop.interfaces  # noqa: F401
        import mloop.controllers  # noqa: F401
        import mloop.visualizations  # noqa: F401
        return
    except ImportError:
        pass

    mloop_module = types.ModuleType("mloop")
    setattr(mloop_module, _STUB_MARKER, True)

    interfaces = types.ModuleType("mloop.interfaces")

    class Interface:
        def __init__(self, *args, **kwargs):
            pass

    interfaces.Interface = Interface

    controllers = types.ModuleType("mloop.controllers")

    class _StubController:
        def __init__(self, interface, **kwargs):
            self.interface = interface
            self.settings = kwargs
            self.best_params = [0.0] * int(kwargs.get("num_params", 0))
            self.optimize_calls = 0

        def optimize(self):
            self.optimize_calls += 1

    controllers.create_controller = (
        lambda interface, **kwargs: _StubController(interface, **kwargs)
    )

    visualizations = types.ModuleType("mloop.visualizations")
    visualizations.show_all_default_visualizations = lambda controller: None

    mloop_module.interfaces = interfaces
    mloop_module.controllers = controllers
    mloop_module.visualizations = visualizations
    sys.modules["mloop"] = mloop_module
    sys.modules["mloop.interfaces"] = interfaces
    sys.modules["mloop.controllers"] = controllers
    sys.modules["mloop.visualizations"] = visualizations


_install_mloop_stub_if_absent()

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from AtomLoadingOptimizer_load_until_atom_master_satellite import (  # noqa: E402
    AtomLoadingOptimizer_load_until_atom_master_satellite,
)
from ExperimentVariables_master_satellite_Node1 import NODE1_VARIABLES  # noqa: E402
from ExperimentVariables_master_satellite_Node2 import NODE2_VARIABLES  # noqa: E402
from ExperimentVariables_master_satellite_global import (  # noqa: E402
    MASTER_SATELLITE_VARIABLES,
)


class FakeDevice:
    def __init__(self, name, log):
        self.name = name
        self.log = log
        self.sw = self
        self.state = None
        self.dac_values = {}

    def _record(self, operation, value=None):
        self.log.append((self.name, operation, value))

    def init(self): self._record("init")
    def input(self): self._record("input")
    def output(self): self._record("output")
    def load(self): self._record("load")
    def set_att(self, value): self._record("set_att", value)
    def set_gain_mu(self, channel, gain): self._record("gain", (channel, gain))
    def on(self): self.state = True; self._record("on")
    def off(self): self.state = False; self._record("off")
    def set(self, **kwargs): self._record("set", kwargs)

    def write_dac(self, channel, value):
        self.dac_values[channel] = value
        self._record("write_dac", (channel, value))

    def set_dac(self, values, channels):
        for channel, value in zip(channels, values):
            self.dac_values[channel] = value
        self._record("set_dac", (tuple(channels), tuple(values)))

    def sample(self, buffer):
        for channel in range(len(buffer)):
            buffer[channel] = 0.0


class FakeCore(FakeDevice):
    def reset(self): self._record("reset")
    def break_realtime(self): self._record("break_realtime")
    def get_rtio_destination_status(self, destination): return True
    def mu_to_seconds(self, duration): return float(duration) * 1e-9


class FakeStabilizerChannel:
    def __init__(self, dataset, dB_dataset):
        self.dataset = dataset
        self.dB_dataset = dB_dataset
        self.dB_history_dataset = dB_dataset + "_history"
        self.set_points = [0.0]
        self.amplitude = 0.05

    def run(self, setpoint_index=0):
        pass


class FakeStabilizer:
    """Stands in for AOMPowerStabilizer.

    The real one reads a config file resolved from the process working
    directory, which is meaningless under the test runner.
    """

    def __init__(self, experiment, dds_names, iterations, averages, **kwargs):
        self.exp = experiment
        self.dds_names = dds_names
        self.all_channels = [
            FakeStabilizerChannel("MOT1_monitor", "p_AOM_A1"),
        ]
        experiment.stabilizer_FORT = FakeStabilizerChannel(
            "FORT_monitor", "p_FORT_loading"
        )
        for index in range(1, 7):
            setattr(
                experiment,
                f"stabilizer_AOM_A{index}",
                FakeStabilizerChannel(
                    f"MOT{index}_monitor", f"p_AOM_A{index}"
                ),
            )

    def run(self):
        pass


class OptimizerEnvironment:
    def initialize_environment(self, devices, datasets, log):
        self.devices = devices
        self.datasets = datasets
        self.log = log
        self.dataset_reads = []
        self.dataset_writes = []
        self.dataset_appends = []
        self.broadcast_sets = {}

    def setattr_argument(self, name, value, *args, **kwargs):
        setattr(self, name, value)

    def setattr_device(self, name):
        setattr(self, name, self.devices[name])

    def get_device(self, name):
        return self.devices[name]

    def get_dataset(self, name, *args, **kwargs):
        self.dataset_reads.append(name)
        return self.datasets[name]

    def set_dataset(self, name, value, **kwargs):
        self.dataset_writes.append((name, value, kwargs))
        self.broadcast_sets[name] = value

    def append_to_dataset(self, name, value):
        self.dataset_appends.append((name, value))
        if name not in self.broadcast_sets:
            raise KeyError(name)
        self.broadcast_sets[name] = list(self.broadcast_sets[name]) + [value]


# The harness sits after the experiment class in the MRO so the dataset
# redirection layer resolves names before the harness records them.
class OptimizerExperiment(
    AtomLoadingOptimizer_load_until_atom_master_satellite, OptimizerEnvironment
):
    _stabilizer_factory = FakeStabilizer


class _OptimizerHarness(unittest.TestCase):
    """Shared fixture. Holds no tests of its own."""

    def setUp(self):
        self.log = []
        self.devices = {"core": FakeCore("core", self.log)}
        config = (
            Path(__file__).resolve().parent.parent
            / "utilities" / "config" / "master_satellite"
            / "device_aliases.json"
        )
        with config.open() as file:
            mapping = json.load(file)
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                self.devices.setdefault(
                    unified_name, FakeDevice(unified_name, self.log)
                )
        for extra in ("core_dma", "scheduler", "k10cr1_ndsp"):
            self.devices.setdefault(extra, FakeDevice(extra, self.log))
        self.datasets = {
            variable.name: variable.value
            for variable in (
                NODE1_VARIABLES + NODE2_VARIABLES + MASTER_SATELLITE_VARIABLES
            )
        }

    def make_experiment(self, updates=None):
        experiment = OptimizerExperiment()
        experiment.initialize_environment(self.devices, self.datasets, self.log)
        experiment.build()
        for name, value in (updates or {}).items():
            setattr(experiment, name, value)
        experiment.prepare()
        return experiment


class AtomLoadingOptimizerMasterSatelliteTests(_OptimizerHarness):
    def test_examination_build_with_none_arguments_reads_no_datasets(self):
        """ARTIQ's repository scan supplies None for every argument."""
        experiment = OptimizerExperiment()
        experiment.initialize_environment(self.devices, {}, self.log)
        experiment.setattr_argument = lambda name, *args, **kwargs: setattr(
            experiment, name, None
        )
        experiment.build()
        self.assertEqual(experiment.dataset_reads, [])

    def test_node_selection_publishes_the_legacy_presentation(self):
        """which_node is consumed and replaced by alice/bob.

        The reused single-node code reads self.which_node expecting the legacy
        spelling; aom_feedback.py builds
        utilities/config/<which_node>/feedback_channels.json from it, and there
        is no config/Node1 directory, so leaving 'Node1' on the attribute makes
        the stabilizer raise FileNotFoundError.
        """
        for node, legacy in (("Node1", "alice"), ("Node2", "bob")):
            experiment = self.make_experiment({"which_node": node})
            self.assertEqual(experiment.which_node, legacy)
            self.assertEqual(experiment._selected_node, node)
            self.assertEqual(experiment.base.which_node, node)
            self.assertEqual(experiment.base.experiment_mode, "single_node")

    def test_persisted_names_resolve_to_the_selected_node(self):
        """The names written with persist=True must be node-suffixed.

        Both lists are persisted by set_experiment_variables_to_best_params.
        Only the coil volts are covered by the dataset redirect (they are
        exactly COIL_CALIBRATION_DATASETS); the six set points are NOT, since
        the feedback map holds each channel's monitor and power datasets rather
        than its set point. Unsuffixed, the beam-power half of every
        optimization would be written where nothing reads it.
        """
        for node in ("Node1", "Node2"):
            experiment = self.make_experiment({"which_node": node})
            self.assertEqual(
                experiment.volt_datasets,
                [f"AZ_bottom_volts_MOT_{node}", f"AZ_top_volts_MOT_{node}",
                 f"AX_volts_MOT_{node}", f"AY_volts_MOT_{node}"],
            )
            self.assertEqual(
                experiment.setpoint_datasets,
                [f"set_point_PD{i}_AOM_A{i}_{node}" for i in range(1, 7)],
            )
            for name in experiment.volt_datasets + experiment.setpoint_datasets:
                self.assertIn(
                    name, self.datasets,
                    f"{name} is not an authoritative master-satellite variable",
                )

    def test_fort_power_monitor_datasets_are_seeded(self):
        """record_FORT_MM/APD_power append without creating; prepare must seed.

        Without this, the first append raised
        KeyError: Cannot mutate nonexistent dataset 'FORT_MM_monitor_Node1'
        on hardware, after the run had already done its feedback. The harness
        mirrors ARTIQ here: append_to_dataset raises for an unseeded name.
        """
        for node in ("Node1", "Node2"):
            experiment = self.make_experiment({"which_node": node})
            for name in (f"FORT_MM_monitor_{node}", f"FORT_APD_monitor_{node}"):
                self.assertIn(
                    name, experiment.broadcast_sets,
                    f"{name} was not seeded, so the first append will raise",
                )
                self.assertEqual(experiment.broadcast_sets[name], [])
            # And appending now works, which is what the kernel path does.
            experiment.append_to_dataset("FORT_MM_monitor", 0.25)
            self.assertEqual(
                experiment.broadcast_sets[f"FORT_MM_monitor_{node}"], [0.25]
            )

    def test_default_setpoints_come_from_the_selected_node(self):
        node1 = self.make_experiment({"which_node": "Node1"})
        node2 = self.make_experiment({"which_node": "Node2"})
        for index in range(6):
            name = f"set_point_PD{index + 1}_AOM_A{index + 1}"
            self.assertEqual(
                node1.default_setpoints[index],
                self.datasets[f"{name}_Node1"],
            )
            self.assertEqual(
                node2.default_setpoints[index],
                self.datasets[f"{name}_Node2"],
            )

    def test_submitted_gui_values_survive_the_dataset_load(self):
        """t_SPCM_exposure and n_measurements shadow persistent variables.

        The standalone stack relies on base.set_datasets_from_gui_args(), which
        has no master-satellite equivalent, so prepare() captures these before
        configure_execution and re-asserts them after.
        """
        experiment = self.make_experiment({
            "which_node": "Node1",
            "t_SPCM_exposure": 0.123,
            "n_measurements": 3,
            "atom_counts_per_s_threshold": 20000.0,
        })
        self.assertEqual(experiment.t_SPCM_exposure, 0.123)
        self.assertEqual(experiment.n_measurements, 3)
        self.assertNotEqual(
            experiment.t_SPCM_exposure, self.datasets["t_SPCM_exposure_Node1"]
        )
        # atom_counts_threshold is derived from the submitted exposure.
        self.assertAlmostEqual(
            experiment.atom_counts_threshold, 20000.0 * 0.123
        )
        self.assertEqual(len(experiment.atom_loading_time_list), 3)

    def test_unsupported_node_fails_clearly(self):
        experiment = OptimizerExperiment()
        experiment.initialize_environment(self.devices, self.datasets, self.log)
        experiment.build()
        experiment.which_node = "two_nodes"
        with self.assertRaisesRegex(ValueError, "Unsupported which_node"):
            experiment.prepare()

    def test_differential_bounds_are_centred_on_the_selected_node(self):
        """Differential boundaries must bracket that node's own coil volts."""
        for node in ("Node1", "Node2"):
            experiment = self.make_experiment({
                "which_node": node,
                "use_differential_boundaries": True,
                "dV_AZ_bottom": 0.05,
                "dV_AZ_top": 0.05,
                "dV_AX": 0.05,
                "dV_AY": 0.05,
            })
            centre = self.datasets[f"AZ_bottom_volts_MOT_{node}"]
            self.assertAlmostEqual(experiment.V_AZ_bottom_min, centre - 0.05)
            self.assertAlmostEqual(experiment.V_AZ_bottom_max, centre + 0.05)

    def test_tuning_mode_selects_the_optimizer_parameters(self):
        both = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "coils and beam powers",
            "disable_z_beam_tuning": False,
        })
        self.assertTrue(both.tune_coils)
        self.assertTrue(both.tune_beams)
        self.assertEqual(len(both.best_params), 10)

        coils = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "coils only",
        })
        self.assertTrue(coils.tune_coils)
        self.assertFalse(coils.tune_beams)
        self.assertEqual(len(coils.best_params), 4)

        beams = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "beam powers only",
            "disable_z_beam_tuning": True,
        })
        self.assertFalse(beams.tune_coils)
        self.assertTrue(beams.tune_beams)
        self.assertEqual(len(beams.best_params), 4)


class ComparePreviousAndOptimizedTests(_OptimizerHarness):
    """The post-optimization check that guards against a worse 'optimum'."""

    @staticmethod
    def _script_costs(experiment, costs):
        """Return scripted costs from optimization_routine, recording params."""
        calls = []
        remaining = list(costs)

        def fake_routine(params):
            calls.append(np.array(params, dtype=float))
            return remaining.pop(0)

        experiment.optimization_routine = fake_routine
        return calls

    def test_previous_params_describe_the_pre_optimization_settings(self):
        """Coil volts enter absolute; set points enter as a 1.0 multiplier."""
        experiment = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "coils and beam powers",
            "disable_z_beam_tuning": False,
        })
        expected = [
            self.datasets["AZ_bottom_volts_MOT_Node1"],
            self.datasets["AZ_top_volts_MOT_Node1"],
            self.datasets["AX_volts_MOT_Node1"],
            self.datasets["AY_volts_MOT_Node1"],
        ] + [1.0] * 6
        self.assertEqual(list(experiment.previous_params), expected)
        # Same length as the vector M-LOOP will return.
        self.assertEqual(
            len(experiment.previous_params), len(experiment.best_params)
        )

        coils_only = self.make_experiment({
            "which_node": "Node2", "what_to_tune": "coils only",
        })
        self.assertEqual(len(coils_only.previous_params), 4)
        self.assertEqual(
            coils_only.previous_params[0],
            self.datasets["AZ_bottom_volts_MOT_Node2"],
        )

        beams_only = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "beam powers only",
            "disable_z_beam_tuning": True,
        })
        self.assertEqual(list(beams_only.previous_params), [1.0] * 4)

    def test_previous_params_survive_the_search_overwriting_coil_values(self):
        """optimization_routine reassigns self.coil_values; the snapshot must not."""
        experiment = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils only",
        })
        snapshot = list(experiment.previous_params)
        # What the search does to self.coil_values on every iteration.
        experiment.coil_values = np.array([9.9, 9.9, 9.9, 9.9])
        self.assertEqual(list(experiment.previous_params), snapshot)

    def test_improved_comparison_returns_the_optimized_vector(self):
        """Accepted only when it beats the previous values AND target_cost."""
        experiment = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils and beam powers",
        })
        optimized = np.array([4.1, 4.9, 0.11, -0.22,
                              1.05, 0.95, 1.02, 0.98, 1.01, 0.99])
        # Cost is -1000/t averaged, so more negative is better. -4000 beats
        # both the previous -2000 and the default target_cost of -3000.
        self.assertEqual(experiment.target_cost, -3000.0)
        calls = self._script_costs(experiment, [-2000, -4000])

        kept = experiment.compare_previous_and_optimized(optimized)

        self.assertEqual(list(kept), list(optimized))
        self.assertEqual(len(calls), 2, "both settings must be re-measured")
        self.assertEqual(list(calls[0]), list(experiment.previous_params),
                         "the previous values are measured first")
        self.assertEqual(list(calls[1]), list(optimized))

    def test_better_than_previous_but_short_of_target_keeps_previous(self):
        """The case RID 38552 hit: improved, but nowhere near target_cost.

        M-LOOP stops at max_runs whether or not the target was reached, so
        beating the starting point is not sufficient to overwrite tuned
        settings.
        """
        experiment = self.make_experiment({
            "which_node": "Node2", "what_to_tune": "coils only",
        })
        optimized = np.array([4.1, 4.9, 0.11, -0.22])
        # -2000 is better than the previous -1000 but does not reach -3000.
        self.assertEqual(experiment.target_cost, -3000.0)
        self._script_costs(experiment, [-1000, -2000])

        kept = experiment.compare_previous_and_optimized(optimized)
        self.assertEqual(list(kept), list(experiment.previous_params))

    def test_target_met_but_not_better_than_previous_keeps_previous(self):
        """Both bars are required, not either one."""
        experiment = self.make_experiment({
            "which_node": "Node2", "what_to_tune": "coils only",
        })
        optimized = np.array([4.1, 4.9, 0.11, -0.22])
        # -5000 clears target -3000, but the previous values measured better.
        self._script_costs(experiment, [-6000, -5000])

        kept = experiment.compare_previous_and_optimized(optimized)
        self.assertEqual(list(kept), list(experiment.previous_params))

    def test_unimproved_comparison_keeps_the_previous_vector(self):
        experiment = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils and beam powers",
        })
        optimized = np.array([4.1, 4.9, 0.11, -0.22,
                              1.05, 0.95, 1.02, 0.98, 1.01, 0.99])
        self._script_costs(experiment, [-3000, -2000])

        kept = experiment.compare_previous_and_optimized(optimized)
        self.assertEqual(list(kept), list(experiment.previous_params))

    def test_tied_comparison_keeps_the_previous_vector(self):
        """No gain is not an improvement: do not overwrite working settings."""
        experiment = self.make_experiment({
            "which_node": "Node2", "what_to_tune": "coils only",
        })
        optimized = np.array([4.1, 4.9, 0.11, -0.22])
        self._script_costs(experiment, [-2500, -2500])

        kept = experiment.compare_previous_and_optimized(optimized)
        self.assertEqual(list(kept), list(experiment.previous_params))

    def test_non_finite_optimized_parameters_skip_the_measurement(self):
        """M-LOOP reports nan when its learner fails; never drive coils with it."""
        experiment = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils only",
        })
        calls = self._script_costs(experiment, [])
        kept = experiment.compare_previous_and_optimized(
            np.array([float("nan")] * 4)
        )
        self.assertEqual(list(kept), list(experiment.previous_params))
        self.assertEqual(calls, [], "nothing should have been measured")

    def test_mismatched_vector_lengths_fail_clearly(self):
        experiment = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils only",
        })
        self._script_costs(experiment, [])
        with self.assertRaisesRegex(ValueError, "must match"):
            experiment.compare_previous_and_optimized(np.zeros(7))

    def test_keeping_previous_values_writes_back_the_original_numbers(self):
        """The 'keep' path still goes through the normal persist call.

        Passing the previous vector to set_experiment_variables_to_best_params
        writes the original numbers back, so the datasets end where they began.
        """
        for node in ("Node1", "Node2"):
            experiment = self.make_experiment({
                "which_node": node,
                "what_to_tune": "coils and beam powers",
                "disable_z_beam_tuning": False,
                "set_best_parameters_at_finish": True,
            })
            self._script_costs(experiment, [-3000, -1000])
            kept = experiment.compare_previous_and_optimized(
                np.array([9.9, 9.9, 9.9, 9.9] + [1.5] * 6)
            )

            before = len(experiment.dataset_writes)
            experiment.set_experiment_variables_to_best_params(kept)
            written = dict(
                (name, value)
                for name, value, _ in experiment.dataset_writes[before:]
            )
            for name in experiment.volt_datasets + experiment.setpoint_datasets:
                self.assertIn(name, written)
                self.assertAlmostEqual(
                    written[name], self.datasets[name],
                    msg=f"{name} was not written back to its original value",
                )

    def test_describe_parameters_reports_absolute_set_points(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "coils and beam powers",
            "disable_z_beam_tuning": False,
        })
        described = experiment.describe_parameters(
            np.array([4.1, 4.9, 0.11, -0.22] + [1.5] * 6)
        )
        names = [name for name, _ in described]
        self.assertEqual(
            names, experiment.volt_datasets + experiment.setpoint_datasets
        )
        self.assertAlmostEqual(described[0][1], 4.1)
        # Set points are reported as default * multiplier, matching what is
        # persisted, not as the raw multiplier.
        self.assertAlmostEqual(
            described[4][1], self.datasets["set_point_PD1_AOM_A1_Node1"] * 1.5
        )

    def test_comparison_is_skipped_when_the_boolean_is_off(self):
        """The default is on; unticking it must restore the old behaviour."""
        experiment = self.make_experiment({
            "which_node": "Node1",
            "what_to_tune": "coils only",
            "compare_prev_and_optimized_values": False,
        })
        self.assertFalse(experiment.compare_prev_and_optimized_values)
        default = self.make_experiment({
            "which_node": "Node1", "what_to_tune": "coils only",
        })
        self.assertTrue(default.compare_prev_and_optimized_values)


class FeedbackSlackGuardTests(unittest.TestCase):
    """Every feedback call must re-arm the timeline first.

    subroutines/aom_feedback.py manages no slack of its own: it advances the
    timeline by 0.2 ms per iteration while spending real time on sampler reads
    and host RPCs. The master-satellite feedback variables are much heavier
    than the standalone ones -- aom_feedback_averages 20 vs 10, and
    aom_feedback_iterations capped at 200 (Node1) / 100 (Node2) vs 15 -- so
    consecutive runs erode slack until the next RTIO event underflows. That is
    what happened on hardware: -4.9 ms of slack on the repump switch, inside
    the warm_up loop.

    Checked at source level on purpose. The failure appears only on real
    hardware, so nothing a hardware-free suite can execute would notice an
    unguarded call being added back.
    """

    @staticmethod
    def _next_executable_line(lines, index):
        for candidate in lines[index + 1:]:
            stripped = candidate.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
        return ""

    def test_every_feedback_call_is_followed_by_a_slack_rearm(self):
        """AFTER, not before: the run returns with the cursor already behind.

        Re-arming before the run was tried on hardware and only moved the
        underflow from inside the loop to the first RTIO event after it
        (dds_FORT.sw.off(), -5.8 ms). What has to be guarded is everything
        downstream of a run, which re-arming after covers -- including the next
        iteration of a loop.
        """
        source = (
            Path(__file__).resolve().parent.parent
            / "AtomLoadingOptimizer_load_until_atom_master_satellite.py"
        )
        lines = source.read_text(encoding="utf-8").splitlines()

        found = 0
        unguarded = []
        for index, line in enumerate(lines):
            if "self.laser_stabilizer.run()" not in line:
                continue
            if line.strip().startswith("#"):
                continue
            found += 1
            following = self._next_executable_line(lines, index)
            if "break_realtime()" not in following:
                unguarded.append(f"line {index + 1}: followed by {following!r}")

        self.assertGreaterEqual(
            found, 5, "the feedback call sites were not found; has the file "
                      "been restructured?"
        )
        self.assertEqual(
            unguarded, [],
            "these laser_stabilizer.run() calls are not immediately followed "
            f"by core.break_realtime(): {unguarded}",
        )


if __name__ == "__main__":
    unittest.main()
