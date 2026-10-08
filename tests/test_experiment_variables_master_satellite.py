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
        "EnvExperiment": FakeEnvExperiment,
        "NumberValue": FakeArgumentValue,
        "BooleanValue": FakeArgumentValue,
        "StringValue": FakeArgumentValue,
        "kernel": lambda function: function,
        "delay": lambda duration: None,
        "MHz": 1e6,
        "kHz": 1e3,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
    }
    for name, value in _stub_exports.items():
        if not hasattr(_experiment_stub, name):
            setattr(_experiment_stub, name, value)
    _existing_exports = set(getattr(_experiment_stub, "__all__", []))
    _existing_exports.update(_stub_exports)
    _experiment_stub.__all__ = sorted(_existing_exports)


node1_module = importlib.import_module("ExperimentVariables_master_satellite_Node1")
node2_module = importlib.import_module("ExperimentVariables_master_satellite_Node2")
global_module = importlib.import_module(
    "ExperimentVariables_master_satellite_global"
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
            node1_module.ExperimentVariablesMasterSatelliteNode1
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
            node2_module.ExperimentVariablesMasterSatelliteNode2
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
            node1_module.ExperimentVariablesMasterSatelliteNode1,
            {"f_FORT_Node1": 251e6},
        )
        self.assertEqual(experiment.f_FORT_Node1, 251e6)
        self.assertEqual(experiment.datasets["f_FORT_Node1"], 251e6)

        globals_experiment = self.initialize(
            global_module.ExperimentVariablesMasterSatelliteGlobal,
            {"t_Node2_excitation_delay_mu": 211},
        )
        self.assertEqual(
            globals_experiment.t_Node2_excitation_delay_mu, 211
        )
        self.assertEqual(
            globals_experiment.datasets["t_Node2_excitation_delay_mu"], 211
        )

    def test_node_values_are_independent(self):
        node1 = self.initialize(node1_module.ExperimentVariablesMasterSatelliteNode1)
        node2 = self.initialize(node2_module.ExperimentVariablesMasterSatelliteNode2)
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
            global_module.ExperimentVariablesMasterSatelliteGlobal
        )
        # Deliberately an exact inventory rather than something derived from
        # MASTER_SATELLITE_VARIABLES: adding a global should be a conscious
        # decision, and this is what forces it to be one. The two-atom
        # thresholds joined it on 2026-09-17, moving out of the per-node
        # files because they describe the joint two-node readout. The
        # two-node joint timing joined on 2026-10-06 for the same reason: a
        # two-node run opens ONE gate across all four master-local SPCM
        # counters, so there is one duration, and the thresholds above are
        # rates compared against counts/t, so gate and threshold must be the
        # same joint quantity.
        self.assertEqual(
            experiment.datasets,
            {
                "n_measurements": 100,
                # Layer 1 (electrical) and Layer 2 (optical) cross-node
                # timing. t_Node2_excitation_delay_mu was t_delay_in_bob_mu,
                # renamed because that name collided with the standalone
                # declaration of the same dataset.
                "t_Node2_rtio_offset_mu": 0,
                "t_Node2_excitation_delay_mu": 189,
                "parallel_AOM_feedback": True,
                # A scratch scan target: no physical meaning, so no per-node
                # meaning either. Was dummy_variable_Node1/_Node2 until
                # 2026-10-07, which made scanning it in two_nodes mode raise
                # "Ambiguous node-specific experiment-variable target".
                "dummy_variable": 0.0,
                "two_atom_threshold": 39000.0,
                "two_atom_threshold_for_loading": 89000.0,
                "t_atom_check_time_two_node": 0.01,
                "t_SPCM_first_shot_two_node": 0.01,
                "t_SPCM_second_shot_two_node": 0.01,
                "t_delay_between_shots_two_node": 0.0,
                "t_MOT_dissipation_two_node": 0.05,
                "max_atom_check_tries_two_node": 100,
                "max_loading_rounds_two_node": 20,
                "n_alternating_RO_windows_per_node": 10,
                # Seconds, not milliseconds: the unit kwarg is display only.
                "t_alternating_RO_window": 0.001,
                "t_alternating_RO_pad": 0.0001,
            },
        )
        # The two Node2 timing globals name the node they apply TO, in the
        # middle of the name; they are not node-OWNED variables and must not
        # pick up the per-node suffix treatment. Both end in _mu, so this
        # still holds -- and it is the assertion that would catch a future
        # global accidentally named with a trailing _Node2.
        self.assertFalse(
            any(
                name.endswith(("_Node1", "_Node2"))
                for name in experiment.datasets
            )
        )

    def test_node_initializers_create_no_unsuffixed_calibrations(self):
        node1 = self.initialize(node1_module.ExperimentVariablesMasterSatelliteNode1)
        node2 = self.initialize(node2_module.ExperimentVariablesMasterSatelliteNode2)
        forbidden = {
            "f_FORT",
            "p_FORT_loading",
            "f_cooling_DP_MOT",
            "p_cooling_DP_MOT",
            "set_point_FORT_APD_loading",
        }
        self.assertTrue(forbidden.isdisjoint(node1.datasets))
        self.assertTrue(forbidden.isdisjoint(node2.datasets))

    def test_two_atom_thresholds_are_global_and_never_suffixed(self):
        """Both nodes' fluorescence lands on the SAME four SPCMs.

        All four detectors see both traps through the beamsplitter fan-out, so
        the two-atom criterion is one joint number, not one per node. These
        were per-node until 2026-09-17 and carried identical values on both --
        the duplication this removed. Being global is also what makes the bare
        name exist in two-node mode, which is how the sequence reads it.

        Pinned because a per-node copy would be invisible until the two nodes'
        values drifted apart, at which point the atom-check criterion would
        depend on which node's copy happened to win.
        """
        thresholds = {"two_atom_threshold", "two_atom_threshold_for_loading"}

        globals_written = self.initialize(
            global_module.ExperimentVariablesMasterSatelliteGlobal
        ).datasets
        for name in thresholds:
            self.assertIn(
                name, globals_written,
                f"{name} must be declared in the GLOBAL variable file",
            )

        node1 = self.initialize(
            node1_module.ExperimentVariablesMasterSatelliteNode1
        ).datasets
        node2 = self.initialize(
            node2_module.ExperimentVariablesMasterSatelliteNode2
        ).datasets
        for name in thresholds:
            for node, written in (("Node1", node1), ("Node2", node2)):
                self.assertNotIn(
                    name, written,
                    f"{name} must not be declared per node",
                )
                self.assertNotIn(
                    f"{name}_{node}", written,
                    f"{name}_{node} must not exist -- the criterion is joint",
                )

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
