import ast
import json
import sys
import types
import unittest
from pathlib import Path


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


if "artiq.experiment" not in sys.modules:
    artiq_module = types.ModuleType("artiq")
    experiment_module = types.ModuleType("artiq.experiment")
    experiment_module.kernel = _identity_decorator
    experiment_module.rpc = _identity_decorator
    experiment_module.delay = lambda duration: None
    experiment_module.EnvExperiment = object
    experiment_module.NumberValue = lambda value, **kwargs: value
    experiment_module.BooleanValue = lambda value=False, **kwargs: value
    experiment_module.StringValue = lambda value, **kwargs: value
    experiment_module.EnumerationValue = lambda values, **kwargs: tuple(values)[0]
    experiment_module.MHz = 1e6
    experiment_module.kHz = 1e3
    experiment_module.ms = 1e-3
    experiment_module.us = 1e-6
    experiment_module.ns = 1e-9
    experiment_module.s = 1.0
    experiment_module.TFloat = float
    experiment_module.TStr = str
    experiment_module.TInt32 = int
    experiment_module.TInt64 = int
    experiment_module.TBool = bool
    artiq_module.experiment = experiment_module
    sys.modules["artiq"] = artiq_module
    sys.modules["artiq.experiment"] = experiment_module
else:
    experiment_module = sys.modules["artiq.experiment"]
    for name, value in {
        "rpc": _identity_decorator,
        "kernel": _identity_decorator,
        "TStr": str,
        "TInt32": int,
        "TInt64": int,
        "TBool": bool,
        "s": 1.0,
    }.items():
        if name in {"rpc", "kernel"} or not hasattr(experiment_module, name):
            setattr(experiment_module, name, value)
    if hasattr(experiment_module, "__all__"):
        experiment_module.__all__ = sorted(
            set(experiment_module.__all__)
            | {"rpc", "kernel", "TStr", "TInt32", "TInt64", "TBool", "s"}
        )


from AOMsCoils_master_satellite_mixin import _AOMsCoilsMasterSatelliteMixin  # noqa: E402
from AOMsCoils_master_satellite_Node1 import (  # noqa: E402
    AOMsCoils_master_satellite_Node1,
)
from AOMsCoils_master_satellite_Node2 import (  # noqa: E402
    AOMsCoils_master_satellite_Node2,
)
from ExperimentVariables_master_satellite_Node1 import NODE1_VARIABLES  # noqa: E402
from ExperimentVariables_master_satellite_Node2 import NODE2_VARIABLES  # noqa: E402
from ExperimentVariables_master_satellite_global import MASTER_SATELLITE_VARIABLES  # noqa: E402


class FakeDevice:
    def __init__(self, name, log):
        self.name = name
        self.log = log
        self.sw = self
        self.state = None
        self.last_set = None
        self.dac_values = {}

    def _record(self, operation, value=None): self.log.append((self.name, operation, value))
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


def make_shared_hardware(log):
    devices = {"core": FakeCore("core", log)}
    with Path("utilities/config/master_satellite/device_aliases.json").open() as file:
        mapping = json.load(file)
    for node_mapping in mapping.values():
        for unified_name in node_mapping.values():
            devices.setdefault(unified_name, FakeDevice(unified_name, log))
    devices["core_dma"] = FakeDevice("core_dma", log)
    devices["scheduler"] = FakeDevice("scheduler", log)
    devices["k10cr1_ndsp"] = FakeDevice("k10cr1_ndsp", log)
    return devices


class ManualExperimentEnvironment:
    def initialize_environment(self, devices, datasets, log):
        self.devices = devices
        self.datasets = datasets
        self.log = log
        self.dataset_reads = []
        self.dataset_writes = []
    def setattr_argument(self, name, value, *args, **kwargs): setattr(self, name, value)
    def setattr_device(self, name): setattr(self, name, self.devices[name])
    def get_dataset(self, name): self.dataset_reads.append(name); return self.datasets[name]
    def set_dataset(self, name, value, **kwargs): self.dataset_writes.append((name, value, kwargs))


class Node1ManualExperiment(ManualExperimentEnvironment, AOMsCoils_master_satellite_Node1):
    pass


class Node2ManualExperiment(ManualExperimentEnvironment, AOMsCoils_master_satellite_Node2):
    pass


class AOMsCoilsMasterSatelliteTests(unittest.TestCase):
    def setUp(self):
        self.log = []
        self.devices = make_shared_hardware(self.log)
        self.datasets = {
            variable.name: variable.value
            for variable in NODE1_VARIABLES + NODE2_VARIABLES + MASTER_SATELLITE_VARIABLES
        }

    def make_experiment(self, experiment_class, updates=None):
        experiment = experiment_class()
        experiment.initialize_environment(self.devices, self.datasets, self.log)
        experiment.build()
        for name, value in (updates or {}).items(): setattr(experiment, name, value)
        experiment.prepare()
        return experiment

    def test_fixed_nodes_unsuffixed_gui_and_coil_channels(self):
        node1 = self.make_experiment(Node1ManualExperiment)
        node2 = self.make_experiment(Node2ManualExperiment)
        self.assertEqual(node1.base.which_node, "Node1")
        self.assertEqual(node2.base.which_node, "Node2")
        self.assertTrue(hasattr(node1, "FORT_AOM_ON"))
        self.assertFalse(hasattr(node1, "FORT_AOM_ON_Node1"))
        self.assertEqual(node1.COIL_CHANNELS, (0, 1, 13, 14))
        self.assertEqual(node2.COIL_CHANNELS, (0, 1, 2, 3))

    def test_node1_then_node2_preserves_node1_outputs(self):
        node1 = self.make_experiment(Node1ManualExperiment, {
            "FORT_AOM_ON": True, "Repump_AOM_switch_ON": True,
        })
        node1.apply_manual_state()
        self.assertFalse(any(operation == "reset" for _, operation, _ in self.log))
        fort1, repump1, zotino1 = (
            self.devices["urukul0_ch0"], self.devices["ttl5"], self.devices["zotino0"]
        )
        state1 = (fort1.state, repump1.state, dict(zotino1.dac_values))
        self.assertEqual(state1, (True, False, {**state1[2]}))
        self.assertEqual(
            zotino1.dac_values[14], self.datasets["AY_volts_MOT_Node1"]
        )

        node2 = self.make_experiment(Node2ManualExperiment, {
            "AOM_A1_ON": True,
        })
        node2_log_start = len(self.log)
        node2.apply_manual_state()
        self.assertFalse(any(operation == "reset" for _, operation, _ in self.log))
        self.assertEqual((fort1.state, repump1.state, dict(zotino1.dac_values)), state1)
        self.assertTrue(self.devices["urukul3_ch0"].state)
        self.assertEqual(
            self.devices["zotino1"].dac_values[2],
            self.datasets["AX_volts_MOT_Node2"],
        )
        with Path("utilities/config/master_satellite/device_aliases.json").open() as file:
            node1_devices = set(json.load(file)["Node1"].values())
        self.assertTrue(
            all(device not in node1_devices for device, _, _ in self.log[node2_log_start:])
        )

    def test_node2_then_node1_preserves_node2_outputs(self):
        node2 = self.make_experiment(Node2ManualExperiment, {
            "FORT_AOM_ON": True,
        })
        node2.apply_manual_state()
        fort2, zotino2 = self.devices["urukul4_ch3"], self.devices["zotino1"]
        self.assertEqual(
            zotino2.dac_values[1], self.datasets["AZ_top_volts_MOT_Node2"]
        )
        state2 = (fort2.state, dict(zotino2.dac_values))
        node1 = self.make_experiment(Node1ManualExperiment, {"AOM_A6_ON": True})
        node1_log_start = len(self.log)
        node1.apply_manual_state()
        self.assertEqual((fort2.state, dict(zotino2.dac_values)), state2)
        with Path("utilities/config/master_satellite/device_aliases.json").open() as file:
            node2_devices = set(json.load(file)["Node2"].values())
        self.assertTrue(
            all(device not in node2_devices for device, _, _ in self.log[node1_log_start:])
        )

    def test_disable_coils_forces_all_four_coils_to_zero(self):
        node1 = self.make_experiment(
            Node1ManualExperiment, {"disable_coils": True}
        )
        node1.apply_manual_state()
        zotino = self.devices["zotino0"]
        self.assertEqual(
            [zotino.dac_values[channel] for channel in (0, 1, 13, 14)],
            [0.0, 0.0, 0.0, 0.0],
        )

    def test_mot_coil_voltages_come_from_persistent_calibration(self):
        node1 = self.make_experiment(Node1ManualExperiment)
        self.assertFalse(node1.disable_coils)
        node1.apply_manual_state()
        zotino = self.devices["zotino0"]
        self.assertEqual(
            [zotino.dac_values[channel] for channel in (0, 1, 13, 14)],
            [
                self.datasets["AZ_bottom_volts_MOT_Node1"],
                self.datasets["AZ_top_volts_MOT_Node1"],
                self.datasets["AX_volts_MOT_Node1"],
                self.datasets["AY_volts_MOT_Node1"],
            ],
        )

    def test_only_selected_node_datasets_are_loaded(self):
        node1 = self.make_experiment(Node1ManualExperiment)
        self.assertTrue({v.name for v in NODE2_VARIABLES}.isdisjoint(node1.dataset_reads))
        node2 = self.make_experiment(Node2ManualExperiment)
        self.assertTrue({v.name for v in NODE1_VARIABLES}.isdisjoint(node2.dataset_reads))

    def test_node2_d1_uses_authoritative_grin2_defaults(self):
        node2 = self.make_experiment(Node2ManualExperiment)
        binding = next(item for item in node2.base.node_resolvers["Node2"].dds_bindings
                       if item["logical_alias"] == "dds_D1_pumping_DP")
        self.assertEqual(binding["frequency_attribute"], "f_GRIN2_D1_pumping_Node2")
        self.assertEqual(binding["power_attribute"], "p_GRIN2_D1_pumping_Node2")

    def test_microwave_defaults_off_and_no_calibration_is_persisted(self):
        node1 = self.make_experiment(Node1ManualExperiment, {
            "microwave_dds_ON": True,
            "yes_Im_sure_I_want_MW_or_RF_dds_ON": False,
        })
        node1.apply_manual_state()
        self.assertFalse(node1.dds_microwaves.state)
        self.assertTrue(node1.ttl_microwave_switch.state)
        self.assertEqual(node1.dataset_writes, [])

    def test_k10cr1_binds_only_when_a_waveplate_action_is_selected(self):
        node1 = self.make_experiment(Node1ManualExperiment)
        self.assertFalse(node1._k10cr1_requested)
        self.assertFalse(hasattr(node1, "k10cr1_ndsp"))

        node2 = self.make_experiment(
            Node2ManualExperiment, {"go_to_home_852QWP": True}
        )
        self.assertTrue(node2._k10cr1_requested)
        self.assertIs(node2.k10cr1_ndsp, self.devices["k10cr1_ndsp"])
        self.assertEqual(node2._axis_780_HWP, "780_HWP_Node2")
        self.assertEqual(node2._axis_852_QWP, "852_QWP_Node2")

    def test_852_target_moves_are_gated_on_852_booleans(self):
        import AOMsCoils_master_satellite_mixin as mixin_module
        from unittest import mock

        calls = []
        with mock.patch.object(
            mixin_module, "go_to_home",
            side_effect=lambda experiment, name: calls.append(("home", name)),
        ), mock.patch.object(
            mixin_module, "move_to_target_deg",
            side_effect=lambda experiment, name, target_deg: calls.append(
                ("target", name, target_deg)
            ),
        ), mock.patch.object(
            mixin_module, "move_by_deg",
            side_effect=lambda experiment, name, target_deg: calls.append(
                ("by", name, target_deg)
            ),
        ):
            node1 = self.make_experiment(Node1ManualExperiment, {
                "go_to_target_780HWP": True,
                "go_to_target_852QWP": True,
            })
            node1.run()

        self.assertEqual(calls, [
            ("target", "780_HWP_Node1", self.datasets["target_780_HWP_Node1"]),
            ("target", "852_QWP_Node1", self.datasets["target_852_QWP_Node1"]),
        ])

    def test_common_module_is_not_an_explorer_experiment(self):
        self.assertFalse(hasattr(_AOMsCoilsMasterSatelliteMixin, "build"))
        expected = {
            "AOMsCoils_master_satellite_mixin.py": set(),
            "AOMsCoils_master_satellite_Node1.py": {"AOMsCoils_master_satellite_Node1"},
            "AOMsCoils_master_satellite_Node2.py": {"AOMsCoils_master_satellite_Node2"},
        }
        for filename, classes in expected.items():
            tree = ast.parse(Path(filename).read_text())
            public = {node.name for node in tree.body if isinstance(node, ast.ClassDef)
                      and not node.name.startswith("_")
                      and any(isinstance(base, ast.Name) and base.id == "EnvExperiment"
                              for base in node.bases)}
            self.assertEqual(public, classes)


if __name__ == "__main__":
    unittest.main()
