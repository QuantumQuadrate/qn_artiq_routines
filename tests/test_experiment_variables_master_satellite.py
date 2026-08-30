import importlib
import json
import sys
import types
import unittest
from pathlib import Path


class FakeArgumentValue:
    def __init__(self, default, **kwargs):
        self.default = default
        self.kwargs = kwargs


class FakeEnvExperiment:
    def __init__(self):
        self.datasets = {}
        self.dataset_writes = []

    def get_dataset(self, name):
        if name not in self.datasets:
            raise KeyError(name)
        return self.datasets[name]

    def setattr_argument(self, name, argument, group=None):
        setattr(self, name, getattr(argument, "default", argument))

    def set_dataset(self, name, value, broadcast=False, persist=False):
        self.datasets[name] = value
        self.dataset_writes.append((name, value, broadcast, persist))


if "artiq.experiment" not in sys.modules:
    artiq_module = types.ModuleType("artiq")
    experiment_module = types.ModuleType("artiq.experiment")
    experiment_module.EnvExperiment = FakeEnvExperiment
    experiment_module.NumberValue = FakeArgumentValue
    experiment_module.BooleanValue = FakeArgumentValue
    experiment_module.StringValue = FakeArgumentValue
    experiment_module.kernel = lambda function: function
    experiment_module.delay = lambda duration: None
    experiment_module.MHz = 1e6
    experiment_module.kHz = 1e3
    experiment_module.ms = 1e-3
    experiment_module.us = 1e-6
    experiment_module.ns = 1e-9
    experiment_module.__all__ = [
        "EnvExperiment",
        "NumberValue",
        "BooleanValue",
        "StringValue",
        "kernel",
        "delay",
        "MHz",
        "kHz",
        "ms",
        "us",
        "ns",
    ]
    artiq_module.experiment = experiment_module
    sys.modules["artiq"] = artiq_module
    sys.modules["artiq.experiment"] = experiment_module

# Other hardware-free tests may have installed a smaller ARTIQ stub first.
experiment_module = sys.modules["artiq.experiment"]
for name, value in {
    "EnvExperiment": FakeEnvExperiment,
    "NumberValue": FakeArgumentValue,
    "BooleanValue": FakeArgumentValue,
    "StringValue": FakeArgumentValue,
    "MHz": 1e6,
    "kHz": 1e3,
    "ms": 1e-3,
    "us": 1e-6,
    "ns": 1e-9,
}.items():
    if not hasattr(experiment_module, name):
        setattr(experiment_module, name, value)

existing_exports = set(getattr(experiment_module, "__all__", []))
existing_exports.update(
    {
        "EnvExperiment",
        "NumberValue",
        "BooleanValue",
        "StringValue",
        "MHz",
        "kHz",
        "ms",
        "us",
        "ns",
    }
)
experiment_module.__all__ = sorted(existing_exports)


node1_module = importlib.import_module("ExperimentVariables_Node1")
node2_module = importlib.import_module("ExperimentVariables_Node2")
global_module = importlib.import_module(
    "ExperimentVariables_master_satellite"
)


class ExperimentVariablesMasterSatelliteTests(unittest.TestCase):
    def initialize(self, experiment_class, existing=None):
        experiment = experiment_class()
        FakeEnvExperiment.__init__(experiment)
        experiment.get_dataset = types.MethodType(
            FakeEnvExperiment.get_dataset, experiment
        )
        experiment.setattr_argument = types.MethodType(
            FakeEnvExperiment.setattr_argument, experiment
        )
        experiment.set_dataset = types.MethodType(
            FakeEnvExperiment.set_dataset, experiment
        )
        experiment.datasets.update(existing or {})
        experiment.build()
        experiment.run()
        return experiment

    def test_fresh_node1_initialization_is_suffixed(self):
        experiment = self.initialize(
            node1_module.ExperimentVariablesNode1
        )
        self.assertEqual(experiment.datasets["f_FORT_Node1"], 245e6)
        self.assertEqual(
            experiment.datasets["t_coil_relaxation_time_OP_Node1"],
            1e-3,
        )
        self.assertTrue(
            all(
                variable.name.endswith("_Node1")
                for variable in node1_module.NODE1_VARIABLES
            )
        )

    def test_fresh_node2_initialization_is_suffixed(self):
        experiment = self.initialize(
            node2_module.ExperimentVariablesNode2
        )
        self.assertEqual(experiment.datasets["f_FORT_Node2"], 240e6)
        self.assertEqual(
            experiment.datasets["t_coil_relaxation_time_OP_Node2"],
            0.4e-3,
        )
        self.assertTrue(
            all(
                variable.name.endswith("_Node2")
                for variable in node2_module.NODE2_VARIABLES
            )
        )

    def test_existing_persistent_value_is_used_instead_of_default(self):
        experiment = self.initialize(
            node1_module.ExperimentVariablesNode1,
            {"f_FORT_Node1": 251e6},
        )
        self.assertEqual(experiment.f_FORT_Node1, 251e6)
        self.assertEqual(experiment.datasets["f_FORT_Node1"], 251e6)

        globals_experiment = self.initialize(
            global_module.ExperimentVariablesMasterSatellite,
            {"t_delay_in_bob_mu": 211},
        )
        self.assertEqual(globals_experiment.t_delay_in_bob_mu, 211)
        self.assertEqual(
            globals_experiment.datasets["t_delay_in_bob_mu"], 211
        )

    def test_node_values_are_independent(self):
        node1 = self.initialize(node1_module.ExperimentVariablesNode1)
        node2 = self.initialize(node2_module.ExperimentVariablesNode2)
        self.assertNotEqual(
            node1.datasets["p_FORT_loading_Node1"],
            node2.datasets["p_FORT_loading_Node2"],
        )
        self.assertNotEqual(
            node1.datasets["single_atom_threshold_Node1"],
            node2.datasets["single_atom_threshold_Node2"],
        )

    def test_global_initialization_has_no_node_side_effect(self):
        experiment = self.initialize(
            global_module.ExperimentVariablesMasterSatellite
        )
        self.assertEqual(
            experiment.datasets,
            {
                "n_measurements": 100,
                "t_delay_in_bob_mu": 189,
                "parallel_AOM_feedback": True,
            },
        )
        self.assertFalse(
            any(
                name.endswith(("_Node1", "_Node2"))
                for name in experiment.datasets
            )
        )

    def test_node_initializers_create_no_unsuffixed_calibrations(self):
        node1 = self.initialize(node1_module.ExperimentVariablesNode1)
        node2 = self.initialize(node2_module.ExperimentVariablesNode2)
        forbidden = {
            "f_FORT",
            "p_FORT_loading",
            "f_cooling_DP_MOT",
            "p_cooling_DP_MOT",
            "set_point_FORT_APD_loading",
        }
        self.assertTrue(forbidden.isdisjoint(node1.datasets))
        self.assertTrue(forbidden.isdisjoint(node2.datasets))

    def test_all_dds_defaults_have_authoritative_node_variables(self):
        names_by_node = {
            "Node1": {
                variable.name
                for variable in node1_module.NODE1_VARIABLES
            },
            "Node2": {
                variable.name
                for variable in node2_module.NODE2_VARIABLES
            },
        }
        for node, standalone_config in (
            ("Node1", "alice"),
            ("Node2", "bob"),
        ):
            path = (
                Path("utilities/config")
                / standalone_config
                / "device_aliases.json"
            )
            with path.open() as config_file:
                defaults = json.load(config_file)["DDS_DEFAULTS"]
            for logical_alias, default_names in defaults.items():
                with self.subTest(node=node, alias=logical_alias):
                    self.assertIn(
                        default_names["frequency"] + "_" + node,
                        names_by_node[node],
                    )
                    self.assertIn(
                        default_names["power"] + "_" + node,
                        names_by_node[node],
                    )

        bob_defaults_path = Path(
            "utilities/config/bob/device_aliases.json"
        )
        with bob_defaults_path.open() as config_file:
            bob_defaults = json.load(config_file)["DDS_DEFAULTS"]
        self.assertEqual(
            bob_defaults["dds_D1_pumping_DP"],
            {
                "frequency": "f_GRIN2_D1_pumping",
                "power": "p_GRIN2_D1_pumping",
            },
        )

    def test_historical_and_result_datasets_are_not_migrated(self):
        all_names = {
            variable.name
            for variable in (
                node1_module.NODE1_VARIABLES
                + node2_module.NODE2_VARIABLES
                + global_module.MASTER_SATELLITE_VARIABLES
            )
        }
        self.assertNotIn("set_point_D1_DP_Node1", all_names)
        self.assertNotIn("set_point_D1_DP_Node2", all_names)
        self.assertNotIn("feedbackchannels", all_names)
        self.assertNotIn("parent_rid", all_names)
        self.assertNotIn("child_rid", all_names)
        self.assertFalse(any(name.startswith("fit_parameter_") for name in all_names))
        self.assertNotIn("which_node", all_names)
        self.assertNotIn("experiment_mode", all_names)


if __name__ == "__main__":
    unittest.main()
