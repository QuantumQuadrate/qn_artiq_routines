import json
import math
import sys
import types
import unittest
from pathlib import Path


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


from FORT_Polarization_Optimizer_master_satellite import (  # noqa: E402
    FORT_Polarization_Optimizer_master_satellite,
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
        self.last_set = None
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
    def set(self, **kwargs): self.last_set = kwargs; self._record("set", kwargs)

    def write_dac(self, channel, value):
        self.dac_values[channel] = value
        self._record("write_dac", (channel, value))

    def set_dac(self, values, channels):
        for channel, value in zip(channels, values):
            self.dac_values[channel] = value
        self._record("set_dac", (tuple(channels), tuple(values)))


class FakeCore(FakeDevice):
    def reset(self): self._record("reset")
    def break_realtime(self): self._record("break_realtime")
    def get_rtio_destination_status(self, destination): return True


class FakeRotator(FakeDevice):
    """Single k10cr1_ndsp controller; node identity lives in the axis name."""

    def __init__(self, name, log):
        super().__init__(name, log)
        self.positions = {}

    def move_to(self, position, axis):
        self.positions[axis] = position
        self._record("move_to", (axis, position))

    def home(self, axis): self._record("home", axis)
    def is_homed(self, axis): return True
    def is_moving(self, axis): return False
    def get_position(self, axis): return self.positions.get(axis, 0)


class FakeSampler(FakeDevice):
    """Channel values may be constants or callables evaluated per sample."""

    def __init__(self, name, log, channel_values=None):
        super().__init__(name, log)
        self.channel_values = channel_values or {}
        self.sample_calls = 0

    def sample(self, buffer):
        self.sample_calls += 1
        for channel in range(8):
            value = self.channel_values.get(channel, 0.0)
            buffer[channel] = value() if callable(value) else value


class FakeStabilizerFORTChannel:
    def __init__(self):
        self.amplitude = 0.05
        self.run_calls = 0

    def run(self, setpoint_index=0):
        self.run_calls += 1


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
        experiment.stabilizer_FORT = FakeStabilizerFORTChannel()


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
    FORT_Polarization_Optimizer_master_satellite, OptimizerEnvironment
):
    _stabilizer_factory = FakeStabilizer


DEG_TO_POS = 136533


class FortPolarizationOptimizerMasterSatelliteTests(unittest.TestCase):
    def setUp(self):
        self.log = []
        self.rotator = FakeRotator("k10cr1_ndsp", self.log)
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
        self.devices["scheduler"] = FakeDevice("scheduler", self.log)
        self.devices["k10cr1_ndsp"] = self.rotator
        self.datasets = {
            variable.name: variable.value
            for variable in (
                NODE1_VARIABLES + NODE2_VARIABLES + MASTER_SATELLITE_VARIABLES
            )
        }

    def install_peaked_node1_samplers(self, peak_hwp, peak_qwp):
        """Node1 MM and APD both live on unified sampler1 (ch 7 and ch 6)."""

        def mm_power():
            hwp = self.rotator.positions.get("852_HWP_Node1", 0) / DEG_TO_POS
            qwp = self.rotator.positions.get("852_QWP_Node1", 0) / DEG_TO_POS
            return math.exp(
                -((hwp - peak_hwp) ** 2 + (qwp - peak_qwp) ** 2) / 200.0
            )

        self.devices["sampler1"] = FakeSampler(
            "sampler1", self.log, {7: mm_power, 6: 0.58}
        )

    def make_experiment(self, updates=None):
        experiment = OptimizerExperiment()
        experiment.initialize_environment(self.devices, self.datasets, self.log)
        experiment.build()
        for name, value in (updates or {}).items():
            setattr(experiment, name, value)
        experiment.prepare()
        return experiment

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

    def test_node_selection_maps_to_single_node_modes(self):
        node1 = self.make_experiment({"which_node": "node1"})
        self.assertEqual(node1.base.experiment_mode, "single_node")
        self.assertEqual(node1.base.which_node, "Node1")
        self.assertEqual(node1.which_node, "alice")
        self.assertEqual(node1._axis_852_HWP, "852_HWP_Node1")
        # The Base-published proxy forwards to the single controller with
        # node-resolved axis names (bare names get this node's suffix).
        node1.k10cr1_ndsp.home("852_HWP")
        self.assertIn(("k10cr1_ndsp", "home", "852_HWP_Node1"), self.log)

        node2 = self.make_experiment({"which_node": "node2"})
        self.assertEqual(node2.base.which_node, "Node2")
        self.assertEqual(node2.which_node, "bob")
        self.assertEqual(node2._axis_852_QWP, "852_QWP_Node2")
        self.assertEqual(
            node2.base.resolve_result_dataset_name("best_852_power_ref"),
            "best_852_power_ref_Node2",
        )

    def test_gui_scan_parameters_win_over_persistent_datasets(self):
        experiment = self.make_experiment({
            "which_node": "node1",
            "tolerance_deg": 3.0,
            "full_range": 16.0,
            "sample_pts": 5,
        })
        self.assertEqual(experiment.tolerance_deg, 3.0)
        self.assertEqual(experiment.full_range, 16.0)
        self.assertEqual(experiment.sample_pts, 5)
        # The persistent values differ, proving the re-assert mattered.
        self.assertNotEqual(self.datasets["tolerance_deg_Node1"], 3.0)
        self.assertNotEqual(self.datasets["full_range_Node1"], 16.0)

    def test_optimization_persists_suffixed_best_and_moves_to_it(self):
        self.install_peaked_node1_samplers(peak_hwp=-4.0, peak_qwp=-14.0)
        experiment = self.make_experiment({
            "which_node": "node1",
            "tolerance_deg": 3.0,
            "full_range": 16.0,
            "sample_pts": 5,
            "search_start_from_scratch": False,
        })
        experiment.initialize_datasets()
        experiment.optimization_routine_zigzag_power_normalized()

        hwp_angles = experiment.broadcast_sets["HWP_angle_Node1"]
        qwp_angles = experiment.broadcast_sets["QWP_angle_Node1"]
        mm_powers = experiment.broadcast_sets["FORT_MM_monitor_Node1"]
        apd_powers = experiment.broadcast_sets["FORT_APD_monitor_Node1"]
        self.assertEqual(len(hwp_angles), len(mm_powers))

        set_point = self.datasets["set_point_FORT_APD_loading_Node1"]
        normalized = [
            mm / (apd / set_point) for mm, apd in zip(mm_powers, apd_powers)
        ]
        best_index = max(
            range(len(normalized)), key=lambda i: normalized[i]
        )

        best_hwp = experiment.broadcast_sets["best_852HWP_to_max_Node1"]
        best_qwp = experiment.broadcast_sets["best_852QWP_to_max_Node1"]
        best_power = experiment.broadcast_sets["best_852_power_Node1"]
        self.assertEqual(best_hwp, hwp_angles[best_index])
        self.assertEqual(best_qwp, qwp_angles[best_index])
        self.assertAlmostEqual(best_power, normalized[best_index])
        persisted = {
            name: kwargs
            for name, _, kwargs in experiment.dataset_writes
            if name.startswith("best_852")
        }
        self.assertTrue(persisted["best_852HWP_to_max_Node1"]["persist"])

        # The routine leaves the plates at the best point it found.
        self.assertEqual(
            self.rotator.positions["852_HWP_Node1"],
            int(best_hwp * DEG_TO_POS),
        )
        self.assertEqual(
            self.rotator.positions["852_QWP_Node1"],
            int(best_qwp * DEG_TO_POS),
        )
        # Feedback ran through the stabilizer FORT channel.
        self.assertGreaterEqual(experiment.stabilizer_FORT.run_calls, 1)

    def test_ref_power_step_persists_suffixed_reference(self):
        self.install_peaked_node1_samplers(peak_hwp=0.0, peak_qwp=0.0)
        experiment = self.make_experiment({"which_node": "node1"})
        experiment.initialize_datasets()
        experiment.run_feedback_and_record_ref_power()
        self.assertIn("best_852_power_ref_Node1", experiment.broadcast_sets)
        self.assertGreaterEqual(experiment.stabilizer_FORT.run_calls, 1)

    def test_node2_apd_reads_its_own_sampler0(self):
        from subroutines.k10cr1_functions import record_FORT_APD_power

        apd_node2 = FakeSampler("sampler3", self.log, {6: 0.44})
        self.devices["sampler3"] = apd_node2
        experiment = self.make_experiment({"which_node": "node2"})
        experiment.initialize_datasets()

        power = record_FORT_APD_power(experiment)
        self.assertAlmostEqual(power, 0.44)
        self.assertGreater(apd_node2.sample_calls, 0)
        self.assertEqual(
            experiment.broadcast_sets["FORT_APD_monitor_Node2"][-1], power
        )


if __name__ == "__main__":
    unittest.main()
