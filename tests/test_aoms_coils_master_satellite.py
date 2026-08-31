import json
import sys
import types
import unittest
from pathlib import Path


if "artiq.experiment" not in sys.modules:
    artiq_module = types.ModuleType("artiq")
    experiment_module = types.ModuleType("artiq.experiment")
    experiment_module.kernel = lambda function: function
    experiment_module.delay = lambda duration: None
    experiment_module.EnvExperiment = object
    experiment_module.NumberValue = lambda value=None, default=None, **kwargs: (
        value if default is None else default
    )
    experiment_module.BooleanValue = lambda default=False, **kwargs: default
    experiment_module.StringValue = lambda value=None, default=None, **kwargs: (
        value if default is None else default
    )
    experiment_module.EnumerationValue = (
        lambda choices, default=None, **kwargs: default or choices[0]
    )
    experiment_module.MHz = 1e6
    experiment_module.kHz = 1e3
    experiment_module.ms = 1e-3
    experiment_module.us = 1e-6
    experiment_module.ns = 1e-9
    experiment_module.TFloat = float
    artiq_module.experiment = experiment_module
    sys.modules["artiq"] = artiq_module
    sys.modules["artiq.experiment"] = experiment_module
else:
    experiment_module = sys.modules["artiq.experiment"]
    if not hasattr(experiment_module, "EnumerationValue"):
        experiment_module.EnumerationValue = (
            lambda choices, default=None, **kwargs: default or choices[0]
        )


from AOMsCoils_master_satellite import AOMsCoils_master_satellite  # noqa: E402
from ExperimentVariables_Node1 import NODE1_VARIABLES  # noqa: E402
from ExperimentVariables_Node2 import NODE2_VARIABLES  # noqa: E402
from ExperimentVariables_master_satellite import (  # noqa: E402
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

    def on(self):
        self.state = True
        self._record("on")

    def off(self):
        self.state = False
        self._record("off")

    def set(self, **kwargs):
        self.last_set = kwargs
        self._record("set", kwargs)

    def write_dac(self, channel, value):
        self.dac_values[channel] = value
        self._record("write_dac", (channel, value))

    def set_dac(self, values, channels):
        for channel, value in zip(channels, values):
            self.dac_values[channel] = value
        self._record("set_dac", (tuple(channels), tuple(values)))


class FakeCore(FakeDevice):
    def __init__(self, log):
        super().__init__("core", log)
        self.reset_calls = 0

    def reset(self):
        self.reset_calls += 1
        self._record("reset")

    def break_realtime(self): self._record("break_realtime")
    def get_rtio_destination_status(self, destination): return True


class ManualExperiment(AOMsCoils_master_satellite):
    def __init__(self):
        self.log = []
        self.dataset_writes = []
        self.datasets = {
            variable.name: variable.value
            for variable in (
                NODE1_VARIABLES
                + NODE2_VARIABLES
                + MASTER_SATELLITE_VARIABLES
            )
        }
        self.core = FakeCore(self.log)
        self.devices = {"core": self.core}
        mapping_path = Path(
            "utilities/config/master_satellite/device_aliases.json"
        )
        with mapping_path.open() as config_file:
            mapping = json.load(config_file)
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                self.devices.setdefault(
                    unified_name, FakeDevice(unified_name, self.log)
                )
        self.devices["core_dma"] = FakeDevice("core_dma", self.log)
        self.devices["scheduler"] = FakeDevice("scheduler", self.log)

    def setattr_argument(self, name, value, *args, **kwargs):
        setattr(self, name, value)

    def setattr_device(self, name):
        setattr(self, name, self.devices[name])

    def get_dataset(self, name):
        return self.datasets[name]

    def set_dataset(self, name, value, **kwargs):
        self.dataset_writes.append((name, value, kwargs))


class AOMsCoilsMasterSatelliteTests(unittest.TestCase):
    def make_experiment(self, which_node, updates=None):
        experiment = ManualExperiment()
        experiment.build()
        experiment.which_node = which_node
        for name, value in (updates or {}).items():
            setattr(experiment, name, value)
        experiment.prepare()
        return experiment

    def test_always_uses_two_node_base_and_explicit_coil_mappings(self):
        experiment = self.make_experiment("node1")
        self.assertEqual(experiment.base.experiment_mode, "two_nodes")
        self.assertEqual(experiment.COIL_CHANNELS_NODE1, (0, 1, 13, 14))
        self.assertEqual(experiment.COIL_CHANNELS_NODE2, (0, 1, 2, 3))
        self.assertIs(experiment.zotino0_Node1, experiment.devices["zotino0"])
        self.assertIs(experiment.zotino0_Node2, experiment.devices["zotino1"])

    def test_base_switches_dds_off_before_programming(self):
        experiment = self.make_experiment("node1")
        experiment.apply_manual_state()
        operations = [
            operation
            for device, operation, _ in experiment.log
            if device == "urukul0_ch0"
        ]
        self.assertLess(operations.index("off"), operations.index("set_att"))
        self.assertLess(operations.index("off"), operations.index("set"))

    def test_node1_selection_hard_gates_all_node2_requests(self):
        updates = {
            "FORT_AOM_ON_Node1": True,
            "AOM_A1_ON_Node1": False,
            "AZ_bottom_coil_ON_Node1": True,
            "AZ_bottom_voltage_Node1": 1.25,
        }
        for control, _ in AOMsCoils_master_satellite.DDS_CONTROLS:
            updates[f"{control}_Node2"] = True
        for coil in AOMsCoils_master_satellite.COIL_NAMES:
            updates[f"{coil}_coil_ON_Node2"] = True
            updates[f"{coil}_voltage_Node2"] = 3.0
        experiment = self.make_experiment("node1", updates)
        experiment.apply_manual_state()

        self.assertTrue(experiment.dds_FORT_Node1.state)
        self.assertFalse(experiment.dds_AOM_A1_Node1.state)
        for device in experiment._manual_dds_node2:
            self.assertFalse(device.state)
        self.assertEqual(
            [experiment.zotino0_Node1.dac_values[channel]
             for channel in experiment.COIL_CHANNELS_NODE1],
            [1.25, 0.0, 0.0, 0.0],
        )
        self.assertEqual(
            [experiment.zotino0_Node2.dac_values[channel]
             for channel in experiment.COIL_CHANNELS_NODE2],
            [0.0, 0.0, 0.0, 0.0],
        )

    def test_node2_selection_hard_gates_node1_and_respects_one_coil(self):
        experiment = self.make_experiment(
            "node2",
            {
                "FORT_AOM_ON_Node1": True,
                "FORT_AOM_ON_Node2": True,
                "AY_coil_ON_Node2": True,
                "AY_voltage_Node2": -1.5,
                "AX_coil_ON_Node1": True,
                "AX_voltage_Node1": 2.0,
            },
        )
        experiment.apply_manual_state()
        self.assertFalse(experiment.dds_FORT_Node1.state)
        self.assertTrue(experiment.dds_FORT_Node2.state)
        self.assertEqual(
            [experiment.zotino0_Node1.dac_values[channel]
             for channel in experiment.COIL_CHANNELS_NODE1],
            [0.0, 0.0, 0.0, 0.0],
        )
        self.assertEqual(
            [experiment.zotino0_Node2.dac_values[channel]
             for channel in experiment.COIL_CHANNELS_NODE2],
            [0.0, 0.0, 0.0, -1.5],
        )

    def test_two_nodes_respects_independent_states(self):
        experiment = self.make_experiment(
            "two_nodes",
            {
                "FORT_AOM_ON_Node1": True,
                "FORT_AOM_ON_Node2": False,
                "AY_coil_ON_Node1": False,
                "AY_coil_ON_Node2": True,
                "AY_voltage_Node2": 0.75,
            },
        )
        experiment.apply_manual_state()
        self.assertTrue(experiment.dds_FORT_Node1.state)
        self.assertFalse(experiment.dds_FORT_Node2.state)
        self.assertEqual(experiment.zotino0_Node1.dac_values[14], 0.0)
        self.assertEqual(experiment.zotino0_Node2.dac_values[3], 0.75)

    def test_dds_defaults_and_physical_resolution_are_authoritative(self):
        experiment = self.make_experiment("two_nodes")
        expected = {
            "Node1": {
                "dds_FORT": "urukul0_ch0",
                "dds_cooling_DP": "urukul0_ch1",
                "dds_AOM_A1": "urukul1_ch2",
                "dds_AOM_A2": "urukul1_ch0",
                "dds_AOM_A3": "urukul1_ch1",
                "dds_AOM_A4": "urukul2_ch0",
                "dds_AOM_A5": "urukul2_ch1",
                "dds_AOM_A6": "urukul1_ch3",
                "dds_D1_pumping_DP": "urukul2_ch2",
                "GRIN1and2_dds": "urukul0_ch3",
            },
            "Node2": {
                "dds_FORT": "urukul4_ch3",
                "dds_cooling_DP": "urukul4_ch2",
                "dds_AOM_A1": "urukul3_ch0",
                "dds_AOM_A2": "urukul3_ch1",
                "dds_AOM_A3": "urukul3_ch2",
                "dds_AOM_A4": "urukul3_ch3",
                "dds_AOM_A5": "urukul4_ch0",
                "dds_AOM_A6": "urukul4_ch1",
                "dds_D1_pumping_DP": "urukul5_ch0",
                "GRIN1and2_dds": "urukul5_ch2",
            },
        }
        for node, mappings in expected.items():
            resolver = experiment.base.node_resolvers[node]
            for alias, physical in mappings.items():
                self.assertEqual(
                    resolver.resolve_dds_physical_name(alias)[1], physical
                )

        node2_d1 = next(
            binding
            for binding in experiment.base.node_resolvers["Node2"].dds_bindings
            if binding["logical_alias"] == "dds_D1_pumping_DP"
        )
        self.assertEqual(
            node2_d1["frequency_attribute"], "f_GRIN2_D1_pumping_Node2"
        )
        self.assertEqual(
            node2_d1["power_attribute"], "p_GRIN2_D1_pumping_Node2"
        )

    def test_no_calibration_persistence_and_run_local_coil_values(self):
        experiment = self.make_experiment(
            "node1", {"AZ_bottom_voltage_Node1": 2.345}
        )
        experiment.apply_manual_state()
        self.assertEqual(experiment.dataset_writes, [])
        self.assertNotEqual(
            experiment.AZ_bottom_voltage_Node1,
            experiment.AZ_bottom_volts_MOT_Node1,
        )

    def test_microwave_and_rf_default_off_and_require_confirmation(self):
        experiment = self.make_experiment(
            "two_nodes",
            {
                "microwave_dds_ON_Node1": True,
                "MW_RF_dds_ON_Node2": True,
                "yes_Im_sure_I_want_MW_or_RF_dds_ON": False,
            },
        )
        experiment.apply_manual_state()
        self.assertFalse(experiment.dds_microwaves_Node1.state)
        self.assertFalse(experiment.dds_MW_RF_Node2.state)
        self.assertTrue(experiment.ttl_microwave_switch_Node1.state)
        self.assertTrue(experiment.ttl_microwave_switch_Node2.state)

    def test_old_architecture_features_are_absent(self):
        experiment = self.make_experiment("node1")
        for name in (
            "Atom_signal_to_other_node_ON",
            "Other_node_atom_signal_checker",
            "run_laser_feedback",
            "laser_stabilizer",
            "go_to_home_780HWP",
            "rigol",
        ):
            self.assertFalse(hasattr(experiment, name), name)


if __name__ == "__main__":
    unittest.main()
