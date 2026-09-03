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
            "override_ExperimentVariables": "{'p_cooling_DP_MOT': -10.0}",
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
