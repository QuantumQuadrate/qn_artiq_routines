import ast
import json
import sys
import types
import unittest
from pathlib import Path
from unittest import mock


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


class _StubParallel:
    def __enter__(self):
        return None

    def __exit__(self, *exc):
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
        "kernel": _identity_decorator,
        "rpc": _identity_decorator,
        "delay": lambda duration: None,
        "parallel": _StubParallel(),
        "EnvExperiment": object,
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value=False, **kwargs: value,
        "StringValue": lambda value, **kwargs: value,
        "EnumerationValue": lambda values, **kwargs: tuple(values)[0],
        "MHz": 1e6,
        "kHz": 1e3,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
        "s": 1.0,
        "V": 1.0,
        "TFloat": float,
        "TStr": str,
        "TInt32": int,
        "TInt64": int,
        "TBool": bool,
    }
    for name, value in _stub_exports.items():
        if not hasattr(_experiment_stub, name):
            setattr(_experiment_stub, name, value)
    if hasattr(_experiment_stub, "__all__"):
        _experiment_stub.__all__ = sorted(
            set(_experiment_stub.__all__) | set(_stub_exports)
        )

    # Heavier third-party modules the imported production code pulls in;
    # only ever stubbed alongside our marked artiq stub, never for real.
    language = types.ModuleType("artiq.language")
    language.us = 1e-6
    language.ns = 1e-9
    language.MHz = 1e6
    sys.modules.setdefault("artiq.language", language)
    sys.modules.setdefault("pyvisa", types.ModuleType("pyvisa"))

    coredevice = types.ModuleType("artiq.coredevice")
    exceptions = types.ModuleType("artiq.coredevice.exceptions")

    class _StubRTIOUnderflow(Exception):
        pass

    exceptions.RTIOUnderflow = _StubRTIOUnderflow
    ad9910 = types.ModuleType("artiq.coredevice.ad9910")
    for constant in (
        "PHASE_MODE_ABSOLUTE",
        "PHASE_MODE_CONTINUOUS",
        "PHASE_MODE_TRACKING",
        "RAM_DEST_ASF",
        "RAM_MODE_RAMPUP",
    ):
        setattr(ad9910, constant, 0)
    urukul = types.ModuleType("artiq.coredevice.urukul")
    urukul.CFG_MASK_NU = 0
    sys.modules.setdefault("artiq.coredevice", coredevice)
    sys.modules.setdefault("artiq.coredevice.exceptions", exceptions)
    sys.modules.setdefault("artiq.coredevice.ad9910", ad9910)
    sys.modules.setdefault("artiq.coredevice.urukul", urukul)

    def _stub_module(name, **attrs):
        module = types.ModuleType(name)
        for key, value in attrs.items():
            setattr(module, key, value)
        sys.modules.setdefault(name, module)
        return sys.modules[name]

    _noop = lambda *args, **kwargs: None
    scipy = _stub_module("scipy")
    scipy.optimize = _stub_module(
        "scipy.optimize", curve_fit=_noop, minimize_scalar=_noop
    )
    scipy.signal = _stub_module("scipy.signal", lombscargle=_noop)
    scipy.special = _stub_module(
        "scipy.special", eval_genlaguerre=_noop, factorial=_noop
    )


from MicrowaveScanOptimizer_master_satellite import (  # noqa: E402
    MicrowaveScanOptimizer_master_satellite,
    scan_dict,
    scan_options,
)
import subroutines.experiment_functions as exp_functions  # noqa: E402
from utilities.conversions import dB_to_V  # noqa: E402
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

    def set_dac(self, values, channels):
        for channel, value in zip(channels, values):
            self.dac_values[channel] = value


class FakeCore(FakeDevice):
    def reset(self): self._record("reset")
    def break_realtime(self): self._record("break_realtime")
    def get_rtio_destination_status(self, destination): return True


class FakeScheduler:
    def __init__(self, rid=100):
        self.rid = rid
        self.submitted = []

    def get_status(self):
        return {}

    def submit(self, pipeline_name, expid, priority, due_date, flush):
        self.submitted.append((pipeline_name, expid, priority))

    def request_termination(self, rid):
        self.terminated = rid


class FakeStabilizerChannel:
    def __init__(self, dataset, dB_dataset):
        self.dataset = dataset
        self.dB_dataset = dB_dataset
        self.dB_history_dataset = dB_dataset + "_history"


class FakeStabilizer:
    def __init__(self, experiment, dds_names, iterations, averages, **kwargs):
        self.exp = experiment
        self.dds_names = dds_names
        self.all_channels = [
            FakeStabilizerChannel("MOT1_monitor", "p_AOM_A1"),
        ]


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

    def get_dataset(self, name, *args, **kwargs):
        self.dataset_reads.append(name)
        if name in self.broadcast_sets:
            return self.broadcast_sets[name]
        return self.datasets[name]

    def set_dataset(self, name, value, **kwargs):
        self.dataset_writes.append((name, value, kwargs))
        self.broadcast_sets[name] = value

    def append_to_dataset(self, name, value):
        self.dataset_appends.append((name, value))
        if name not in self.broadcast_sets:
            raise KeyError(name)
        self.broadcast_sets[name] = list(self.broadcast_sets[name]) + [value]

    def write_results(self, kwargs={}):
        self.log.append(("write_results", kwargs.get("name"), None))


# The harness sits after the experiment class in the MRO so the dataset
# redirection layer resolves names before the harness records them.
class OptimizerExperiment(
    MicrowaveScanOptimizer_master_satellite, OptimizerEnvironment
):
    _stabilizer_factory = FakeStabilizer


class MicrowaveScanOptimizerMasterSatelliteTests(unittest.TestCase):
    def setUp(self):
        self.log = []
        self.devices = {"core": FakeCore("core", self.log)}
        with Path(
            "utilities/config/master_satellite/device_aliases.json"
        ).open() as file:
            mapping = json.load(file)
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                self.devices.setdefault(
                    unified_name, FakeDevice(unified_name, self.log)
                )
        self.devices["core_dma"] = FakeDevice("core_dma", self.log)
        self.scheduler = FakeScheduler()
        self.devices["scheduler"] = self.scheduler
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
        # Base.build installs the real write_results wrapper as an instance
        # attribute; replace it with the recording stub for tests.
        experiment.write_results = lambda kwargs={}: self.log.append(
            ("write_results", kwargs.get("name"), None)
        )
        for name, value in (updates or {}).items():
            setattr(experiment, name, value)
        experiment.prepare()
        return experiment

    def test_scan_definitions_are_loaded_from_the_standalone_source(self):
        self.assertEqual(len(scan_options), 13)
        self.assertEqual(set(scan_dict), set(scan_options))
        for definition in scan_dict.values():
            self.assertIn("experiment_name", definition)
        # No class import from the standalone module, so ARTIQ cannot
        # rediscover the standalone experiment through this file.
        source = Path("MicrowaveScanOptimizer_master_satellite.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("from MicrowaveScanOptimizer import", source)
        tree = ast.parse(source)
        public = {
            node.name for node in tree.body
            if isinstance(node, ast.ClassDef) and not node.name.startswith("_")
        }
        self.assertEqual(public, {"MicrowaveScanOptimizer_master_satellite"})

    def test_examination_build_with_none_arguments_reads_no_datasets(self):
        experiment = OptimizerExperiment()
        experiment.initialize_environment(self.devices, {}, self.log)
        experiment.setattr_argument = lambda name, *args, **kwargs: setattr(
            experiment, name, None
        )
        experiment.build()
        self.assertFalse(experiment.base._execution_configured)
        self.assertEqual(experiment.dataset_reads, [])
        with self.assertRaisesRegex(ValueError, "Unsupported which_node"):
            experiment.prepare()

    def test_prepare_maps_node_and_resolves_scan_targets(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 7,
        })
        self.assertEqual(experiment.base.which_node, "Node1")
        self.assertEqual(experiment.which_node, "alice")
        self.assertEqual(experiment._selected_node, "Node1")
        self.assertEqual(experiment.n_measurements, 7)

        self.assertEqual(experiment.scan_type, "Frequency_00_Scan")
        self.assertEqual(experiment.scan_variable1_name, "f_microwaves_dds")
        self.assertEqual(experiment.scan_variable1, "f_microwaves_dds_Node1")
        # override_items: t_microwave_pulse takes the 00 pi-pulse value.
        self.assertEqual(
            experiment.t_microwave_pulse_Node1,
            self.datasets["t_microwave_00_pulse_Node1"],
        )
        self.assertEqual(
            experiment.t_microwave_pulse, experiment.t_microwave_pulse_Node1
        )
        self.assertGreater(len(experiment.scan_sequence1), 0)
        self.assertIs(
            experiment._selected_experiment_function,
            getattr(exp_functions, "microwave_Rabi_2_experiment"),
        )

    def test_derived_amplitudes_follow_authoritative_values(self):
        experiment = self.make_experiment({
            "which_node": "Node2",
            "Time_00_Scan": True,
            "override_ExperimentVariables_Node2": "{'p_cooling_DP_MOT': -10.0}",
        })
        self.assertAlmostEqual(
            experiment.ampl_cooling_DP_MOT, dB_to_V(-10.0)
        )
        self.assertAlmostEqual(
            experiment.ampl_cooling_DP_RO,
            dB_to_V(-10.0) * experiment.p_cooling_DP_RO,
        )
        self.assertAlmostEqual(
            experiment.ampl_AOM_A1,
            dB_to_V(self.datasets["p_AOM_A1_Node2"]),
        )
        self.assertEqual(
            experiment.base.resolve_result_dataset_name(
                "health_check_uw_freq00"
            ),
            "health_check_uw_freq00_Node2",
        )

    def test_result_state_covers_the_microwave_chain(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
        })
        experiment.initialize_datasets()
        for dataset_name in (
            "AllSPCMs_RO1", "AllSPCMs_RO2",
            "SPCM0_RO1", "SPCM1_OtherNode_RO2",
            "SPCM0_RO1_in_health_check", "AllSPCMs_RO2_in_health_check",
            "n_feedback_per_iteration", "n_atom_loaded_per_iteration",
            "atom_loading_wall_clock", "Atom_loading_time",
            "time_without_atom", "GRIN1_D1_monitor", "GRIN1_EXC_monitor",
            "Magnetometer_TEST_X", "Sampler0_test", "photocount_bins",
            "AllSPCMs_atom_check_in_loading",
        ):
            self.assertIn(dataset_name, experiment.broadcast_sets)
        self.assertEqual(len(experiment.AllSPCMs_RO1_list), 4)
        self.assertEqual(len(experiment.atom_loading_time_list), 4)

    def test_health_check_passes_and_stores_suffixed_fidelity(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
            "target_fidelity": 0.8,
        })
        experiment.initialize_datasets()

        cutoff = (
            experiment.single_atom_threshold * experiment.t_SPCM_first_shot
        )
        loaded = cutoff * 10

        def fake_scan_point(inner_self):
            # One iteration: every atom loads, none is retained -> for a
            # resonance dip the fidelity is 1.0.
            shots1 = [loaded] * 4
            shots2 = [0.0] * 4
            inner_self.set_dataset(
                "AllSPCMs_RO1", [0] + shots1, broadcast=True
            )
            inner_self.set_dataset(
                "AllSPCMs_RO2", [0] + shots2, broadcast=True
            )

        with mock.patch.object(
            exp_functions, "microwave_Rabi_2_experiment", fake_scan_point
        ):
            passed = experiment.health_check_general()

        self.assertTrue(passed)
        self.assertEqual(
            experiment.broadcast_sets["health_check_uw_freq00_Node1"], 1.0
        )
        self.assertNotIn("health_check_uw_freq00", experiment.broadcast_sets)

    def test_failed_health_check_reports_low_fidelity(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
            "target_fidelity": 0.8,
        })
        experiment.initialize_datasets()

        cutoff = (
            experiment.single_atom_threshold * experiment.t_SPCM_first_shot
        )
        loaded = cutoff * 10

        def fake_scan_point(inner_self):
            # Full retention -> dip fidelity 0.0 -> health check fails.
            inner_self.set_dataset(
                "AllSPCMs_RO1", [0] + [loaded] * 4, broadcast=True
            )
            inner_self.set_dataset(
                "AllSPCMs_RO2", [0] + [loaded] * 4, broadcast=True
            )

        with mock.patch.object(
            exp_functions, "microwave_Rabi_2_experiment", fake_scan_point
        ):
            passed = experiment.health_check_general()

        self.assertFalse(passed)
        self.assertEqual(
            experiment.broadcast_sets["health_check_uw_freq00_Node1"], 0.0
        )

    def test_resubmission_targets_master_satellite_class_with_which_node(self):
        experiment = self.make_experiment({
            "which_node": "Node2",
            "Frequency_01_Scan": True,
            "n_measurements": 6,
            # Both dictionaries populated. The resubmitted expid is the only
            # carrier of the idle node's, so dropping it here would silently
            # empty that field for every later run of the resubmitted chain.
            "override_ExperimentVariables_Node1": "{'f_FORT_Node1': 241.0}",
            "override_ExperimentVariables_Node2": "{'p_cooling_DP_MOT': -9.0}",
        })
        experiment.submit_opt_exp_general()

        self.assertEqual(len(self.scheduler.submitted), 1)
        _, expid, priority = self.scheduler.submitted[0]
        self.assertEqual(priority, 99)
        self.assertEqual(
            expid["file"],
            "qn_artiq_routines\\MicrowaveScanOptimizer_master_satellite.py",
        )
        self.assertEqual(
            expid["class_name"], "MicrowaveScanOptimizer_master_satellite"
        )
        self.assertEqual(expid["arguments"]["which_node"], "Node2")
        self.assertTrue(expid["arguments"]["Frequency_01_Scan"])
        self.assertFalse(expid["arguments"]["run_health_check_and_optimize"])
        self.assertEqual(expid["arguments"]["n_measurements"], 6)
        self.assertEqual(
            expid["arguments"]["override_ExperimentVariables_Node1"],
            "{'f_FORT_Node1': 241.0}",
        )
        self.assertEqual(
            expid["arguments"]["override_ExperimentVariables_Node2"],
            "{'p_cooling_DP_MOT': -9.0}",
        )
        self.assertNotIn(
            "override_ExperimentVariables", expid["arguments"],
            "the removed node-agnostic field is still being submitted",
        )

    def _adaptive_offsets(self, node, scan, geometric, pulse_name, t_pi):
        self.datasets[f"{pulse_name}_{node}"] = t_pi
        experiment = self.make_experiment({
            "which_node": node,
            scan: True,
            "enable_node_specific_adaptive_scan": True,
            "enable_geometric_frequency_scan": geometric,
        })
        offsets = [value - experiment.center for value in experiment.scan_sequence1]
        return experiment, offsets

    def test_adaptive_frequency_scan_follows_each_nodes_pi_time(self):
        """Node1's 01 line (17 us pulse, ~47 kHz wide) and Node2's (5.4 us,
        ~148 kHz) get scans sized to their own linewidth."""
        for node, t_pi in (("Node1", 17e-6), ("Node2", 5.4e-6)):
            linewidth = 0.80 / t_pi
            experiment, offsets = self._adaptive_offsets(
                node, "Frequency_01_Scan", False, "t_microwave_01_pulse", t_pi)
            self.assertEqual(len(offsets), 13)
            self.assertAlmostEqual(offsets[0], -2.0 * linewidth, delta=1.0)
            self.assertAlmostEqual(offsets[-1], 2.0 * linewidth, delta=1.0)
            self.assertAlmostEqual(offsets[1] - offsets[0], linewidth / 3.0, delta=1.0)
            self.assertAlmostEqual(
                experiment.initialise["sigma"], linewidth / 2.3548, delta=1.0)
            # The broader search and re-centering read the same run-local values.
            self.assertAlmostEqual(
                experiment.freq_scan_step_size_kHz * 1e3, linewidth / 3.0, delta=1.0)
            self.assertAlmostEqual(
                experiment.freq_scan_min_step_size_kHz * 1e3, 0.2 * linewidth, delta=1.0)
            # The shared scan definition is not modified.
            self.assertEqual(scan_dict["Frequency_01_Scan"]["initialise"], {"sigma": 50000.0})

    def test_adaptive_geometric_list_is_in_linewidth_units(self):
        linewidth = 0.80 / 4.4e-6
        _, offsets = self._adaptive_offsets(
            "Node1", "Frequency_11_Scan", True, "t_microwave_11_pulse", 4.4e-6)
        self.assertAlmostEqual(max(abs(o) for o in offsets), 2.0 * linewidth, delta=1.0)
        self.assertEqual(len(offsets), 9)
        self.assertIn(0.0, offsets)

    def test_adaptive_time_scan_runs_from_zero_to_two_pi(self):
        self.datasets["t_microwave_11_pulse_Node2"] = 9e-6
        experiment = self.make_experiment({
            "which_node": "Node2",
            "Time_11_Scan": True,
            "enable_node_specific_adaptive_scan": True,
        })
        sequence = list(experiment.scan_sequence1)
        self.assertEqual(len(sequence), 11)
        self.assertEqual(sequence[0], 0.0)
        self.assertAlmostEqual(sequence[-1], 18e-6, places=12)

    def test_unchecked_adaptive_scan_keeps_the_gui_ranges(self):
        self.datasets["t_microwave_01_pulse_Node1"] = 17e-6
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_01_Scan": True,
            "enable_geometric_frequency_scan": False,
        })
        offsets = [value - experiment.center for value in experiment.scan_sequence1]
        self.assertAlmostEqual(offsets[0], -125e3, delta=1.0)
        self.assertAlmostEqual(offsets[1] - offsets[0], 20e3, delta=1.0)
        self.assertEqual(experiment.initialise, {"sigma": 50000.0})

    def test_resubmission_carries_the_adaptive_scan_setting(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_11_Scan": True,
            "enable_node_specific_adaptive_scan": True,
            "adaptive_half_range_linewidths": 2.5,
        })
        experiment.submit_opt_exp_general()
        _, expid, _ = self.scheduler.submitted[-1]
        self.assertTrue(expid["arguments"]["enable_node_specific_adaptive_scan"])
        self.assertEqual(expid["arguments"]["adaptive_half_range_linewidths"], 2.5)
        for name in ("adaptive_points_per_linewidth", "adaptive_time_scan_max_pi",
                     "adaptive_time_scan_n_steps"):
            self.assertIn(name, expid["arguments"])

    def test_idle_node_overrides_do_not_block_the_running_node(self):
        """The optimizer had the same failure GeneralVariableScan did.

        Node1 entries left in place must not break a Node2 run. The Node1
        dictionary here reads self.p_cooling_DP_MOT_Node1, which does not
        exist in single_node Node2 mode, so the run only succeeds if that
        dictionary is never evaluated.
        """
        experiment = self.make_experiment({
            "which_node": "Node2",
            "Time_00_Scan": True,
            "override_ExperimentVariables_Node1": (
                "{'p_cooling_DP_MOT_Node1': self.p_cooling_DP_MOT_Node1}"
            ),
            "override_ExperimentVariables_Node2": "{'p_cooling_DP_MOT': -10.0}",
        })
        self.assertEqual(
            experiment._user_authoritative_overrides,
            {"p_cooling_DP_MOT_Node2": -10.0},
        )
        self.assertAlmostEqual(
            experiment.ampl_cooling_DP_MOT, dB_to_V(-10.0)
        )

    def test_no_node_agnostic_override_field_remains(self):
        """The shared field was removed: with per-node ones it is never read."""
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Time_00_Scan": True,
        })
        self.assertFalse(
            hasattr(experiment, "override_ExperimentVariables"),
            "the removed node-agnostic override field is still declared",
        )
        for node in ("Node1", "Node2"):
            self.assertTrue(
                hasattr(experiment, f"override_ExperimentVariables_{node}"),
                f"override_ExperimentVariables_{node} is not declared",
            )

    def test_min_step_displacement_uses_khz_units_in_both_sources(self):
        import re

        for filename in (
            "standalone/MicrowaveScanOptimizer.py",
            "MicrowaveScanOptimizer_master_satellite.py",
        ):
            source = Path(filename).read_text(encoding="utf-8")
            converted = re.findall(
                r"new_center = f\d [+-] "
                r"self\.freq_scan_min_step_size_kHz \* kHz",
                source,
            )
            self.assertEqual(len(converted), 4, filename)
            bare = re.findall(
                r"new_center = f\d [+-] "
                r"self\.freq_scan_min_step_size_kHz(?! \* kHz)",
                source,
            )
            self.assertEqual(bare, [], filename)
            # Difference-based offsets are already in Hz and stay unconverted.
            self.assertIn("step = abs(f2 - f0)", source)

    def test_adaptive_refine_displaces_center_by_ten_kilohertz(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
            "run_health_check_and_optimize": False,
            "enable_geometric_frequency_scan": True,
            "enable_fitting": False,
            "freq_scan_min_step_size_kHz": 10.0,
        })

        cutoff = (
            experiment.single_atom_threshold * experiment.t_SPCM_first_shot
        )
        loaded = cutoff * 10
        # Retentions per scanned point: the first point (f0 = center - x) is
        # the lowest and far below 0.2, the left slope is positive, so the
        # adaptive logic must take Case 2.1 and shift the center left by
        # exactly freq_scan_min_step_size_kHz expressed in Hz.
        retained_counts = iter([0, 3, 3, 3] + [3] * 50)

        def fake_scan_point(inner_self):
            shots1 = [loaded] * 4
            retained = next(retained_counts)
            shots2 = [loaded] * retained + [0.0] * (4 - retained)
            ro1 = list(inner_self.get_dataset("AllSPCMs_RO1")) + shots1
            ro2 = list(inner_self.get_dataset("AllSPCMs_RO2")) + shots2
            inner_self.set_dataset("AllSPCMs_RO1", ro1, broadcast=True)
            inner_self.set_dataset("AllSPCMs_RO2", ro2, broadcast=True)

        # prepare() bound the real function by reference; substitute the
        # fake for the scan loop itself.
        experiment._selected_experiment_function = fake_scan_point
        experiment.experiment_function = lambda: fake_scan_point(experiment)
        experiment.run()

        scanned = list(experiment.broadcast_sets["scan_sequence1"])
        self.assertGreater(len(scanned), 4)
        f0 = scanned[0]
        # Pair-mode scan lists end on their center point, so the last refined
        # point IS the new center: f0 - 10 kHz, not f0 - 10 Hz.
        refined_center = scanned[-1]
        self.assertAlmostEqual(f0 - refined_center, 10_000.0)

    def test_fit_check_sees_the_calibration_from_before_the_scan(self):
        """The fit check's "original" is the calibration before the scan.

        _apply_scan_point writes every point INTO the scanned attribute, and
        for the 11, m10 and m11 scans that attribute is the calibration
        itself: scan_dict names their scan variable and their center /
        pi_pulse alike. RID 38934, a Frequency_11 scan of a 339.213 MHz
        calibration, printed "Original f0 value: 339.073 MHz" -- the last
        point it scanned, its refined center. Both scan branches (fixed and
        adaptive) and a time scan are covered.
        """
        cases = (
            ("Frequency_11_Scan", False, "f_microwaves_11_dds_Node1",
             {"f0": 339.2e6, "sigma": 50e3, "C": 0.9}),
            ("Frequency_11_Scan", True, "f_microwaves_11_dds_Node1",
             {"f0": 339.2e6, "sigma": 50e3, "C": 0.9}),
            ("Time_11_Scan", False, "t_microwave_11_pulse_Node1",
             {"t_pi": 5e-6}),
        )

        class FakeFit:
            def __init__(self, parameters):
                self.parameters = parameters

            def fit(self, t, y, **kwargs):
                return (dict(self.parameters),
                        {name: 0.0 for name in self.parameters})

        for scan_flag, adaptive, calibration, parameters in cases:
            with self.subTest(scan=scan_flag, adaptive=adaptive):
                experiment = self.make_experiment({
                    "which_node": "Node1",
                    scan_flag: True,
                    "n_measurements": 4,
                    "run_health_check_and_optimize": False,
                    "enable_geometric_frequency_scan": adaptive,
                    "enable_fitting": True,
                    "freq_scan_min_step_size_kHz": 10.0,
                    "time_scan_sequence": "np.arange(0,10,1)*us",
                })
                before = getattr(experiment, calibration)

                cutoff = (
                    experiment.single_atom_threshold
                    * experiment.t_SPCM_first_shot
                )
                loaded = cutoff * 10
                # Uneven first four points, so the adaptive branch refines
                # (Frequency_11 fits a resonance PEAK, so this moves the
                # center right) and the scan ends away from `before`.
                retained_counts = iter([0, 3, 3, 3] + [3] * 50)

                def fake_scan_point(inner_self):
                    retained = next(retained_counts)
                    ro1 = (list(inner_self.get_dataset("AllSPCMs_RO1"))
                           + [loaded] * 4)
                    ro2 = (list(inner_self.get_dataset("AllSPCMs_RO2"))
                           + [loaded] * retained + [0.0] * (4 - retained))
                    inner_self.set_dataset("AllSPCMs_RO1", ro1, broadcast=True)
                    inner_self.set_dataset("AllSPCMs_RO2", ro2, broadcast=True)

                experiment._selected_experiment_function = fake_scan_point
                experiment.experiment_function = (
                    lambda: fake_scan_point(experiment)
                )
                experiment.which_fit_model = FakeFit(parameters)
                seen = []

                def recording_health_check(fit_check=False):
                    seen.append(getattr(experiment, calibration))
                    return False

                experiment.health_check_general = recording_health_check

                experiment.run()

                scanned = list(experiment.broadcast_sets["scan_sequence1"])
                self.assertNotEqual(
                    scanned[-1], before,
                    "the scan must end away from the calibration, or this "
                    "test cannot tell a restore from no change",
                )
                self.assertEqual(seen, [before])

    def test_no_local_write_follows_a_broadcast_of_the_same_name(self):
        """No dataset is written locally after this run broadcast it.

        ARTIQ's DatasetManager.set deletes a key from the master DB when a
        local write follows a broadcast of the same key within one run
        (artiq/master/worker_db.py, the `elif key in
        self._broadcaster.raw_view` branch). initialize_datasets used to
        record the overrides AFTER initialize_result_state, which broadcasts
        and persists n_measurements, so an override of n_measurements
        deleted it from the DB and every master-satellite experiment then
        failed in prepare(). Checked over the whole run, build to end, so a
        new local write anywhere is caught too.
        """
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
            "run_health_check_and_optimize": False,
            "enable_geometric_frequency_scan": False,
            "enable_fitting": False,
            "override_ExperimentVariables_Node1": "{'n_measurements': 4}",
        })
        experiment.experiment_function = lambda: None
        experiment.run()

        broadcast = set()
        offenders = []
        local = []
        for name, _, kwargs in experiment.dataset_writes:
            if kwargs.get("broadcast") or kwargs.get("persist"):
                broadcast.add(name)
            else:
                local.append(name)
                if name in broadcast:
                    offenders.append(name)
        self.assertIn(
            "n_measurements", local,
            "the harness must see the override record, or this test is moot",
        )
        self.assertEqual(offenders, [])

    def test_magnetometer_results_exist_for_every_experiment_call(self):
        """measure_Magnetometer can append on every path the optimizer runs.

        Base creates Magnetometer_*_NodeX in initialize_result_datasets, once
        per run; GVS calls it and this optimizer did not, so with
        monitor_magnetometer_in_end_measurement on the first append raised
        "Cannot mutate nonexistent dataset 'Magnetometer_Zero_X_Node1'" --
        RID 38965, in the health check. The harness's append_to_dataset raises
        the same KeyError for a dataset that was never set, as ARTIQ does.
        Covers the health-check path and the scan path.
        """
        def appending_scan_point(inner_self):
            # What measure_Magnetometer does, through the same redirect.
            for axis in ("X", "Y", "Z"):
                inner_self.append_to_dataset(f"Magnetometer_Zero_{axis}", 1.0)

        for node in ("Node1", "Node2"):
            with self.subTest(node=node, path="health check"):
                experiment = self.make_experiment({
                    "which_node": node,
                    "Frequency_00_Scan": True,
                    "n_measurements": 4,
                    "target_fidelity": 0.8,
                })
                experiment.initialize_datasets()
                with mock.patch.object(
                    exp_functions, "microwave_Rabi_2_experiment",
                    appending_scan_point,
                ):
                    try:
                        experiment.health_check_general()
                    except IndexError:
                        # No SPCM data from this fake, so the fidelity step
                        # has nothing to index; the appends already ran.
                        pass
                self.assertEqual(
                    experiment.broadcast_sets[f"Magnetometer_Zero_X_{node}"],
                    [0.0, 1.0],
                )

            with self.subTest(node=node, path="scan"):
                experiment = self.make_experiment({
                    "which_node": node,
                    "Frequency_00_Scan": True,
                    "n_measurements": 4,
                    "run_health_check_and_optimize": False,
                    "enable_geometric_frequency_scan": False,
                    "enable_fitting": False,
                })
                calls = []

                def counting_scan_point():
                    calls.append(None)
                    appending_scan_point(experiment)

                experiment.experiment_function = counting_scan_point
                experiment.run()
                # Seeded once per run and never cleared between scan points
                # (user, 2026-10-23), so one sample per experiment call.
                self.assertGreater(len(calls), 1)
                self.assertEqual(
                    experiment.broadcast_sets[f"Magnetometer_Zero_Z_{node}"],
                    [0.0] + [1.0] * len(calls),
                )

    def test_retention_slicing_skips_the_dataset_seed(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
            "n_measurements": 4,
        })
        loaded = 200.0
        cutoff = 100.0
        # One iteration of four shots behind the [0] seed: three atoms load
        # (indices 0, 1, 3), two of those are retained (indices 0 and 3).
        shot1 = [0] + [loaded, loaded, 0.0, loaded]
        shot2 = [0] + [loaded, 0.0, 0.0, loaded]
        retention = experiment.get_retention(shot1, shot2, 4, 1, cutoff)
        self.assertAlmostEqual(retention[0], 2 / 3)

        retention_array, loading_rate, n_loaded = (
            experiment.get_loading_and_retention(shot1, shot2, 4, 1, cutoff)
        )
        self.assertAlmostEqual(retention_array[0], 2 / 3)
        self.assertAlmostEqual(loading_rate[0], 3 / 4)
        self.assertEqual(n_loaded[0], 3)

    def test_scan_point_updates_authoritative_and_compatibility_values(self):
        experiment = self.make_experiment({
            "which_node": "Node1",
            "Frequency_00_Scan": True,
        })
        experiment._apply_scan_point(123456.0)
        self.assertEqual(experiment.f_microwaves_dds_Node1, 123456.0)
        self.assertEqual(experiment.f_microwaves_dds, 123456.0)


if __name__ == "__main__":
    unittest.main()
