import sys
import types
import unittest
import json
from pathlib import Path


if "artiq.experiment" not in sys.modules:
    artiq_module = types.ModuleType("artiq")
    experiment_module = types.ModuleType("artiq.experiment")
    experiment_module.kernel = lambda function: function
    experiment_module.delay = lambda duration: None
    experiment_module.EnvExperiment = object
    experiment_module.NumberValue = lambda value, **kwargs: value
    experiment_module.BooleanValue = lambda value, **kwargs: value
    experiment_module.StringValue = lambda value, **kwargs: value
    experiment_module.MHz = 1e6
    experiment_module.kHz = 1e3
    experiment_module.ms = 1e-3
    experiment_module.us = 1e-6
    experiment_module.ns = 1e-9
    artiq_module.experiment = experiment_module
    sys.modules["artiq"] = artiq_module
    sys.modules["artiq.experiment"] = experiment_module


from utilities.BaseExperiment_master_satellite import (  # noqa: E402
    BaseExperimentMasterSatellite,
)
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

    def _record(self, operation):
        self.log.append((self.name, operation))

    def init(self): self._record("init")
    def input(self): self._record("input")
    def output(self): self._record("output")
    def on(self): self._record("on")
    def off(self): self._record("off")
    def load(self): self._record("load")
    def set_att(self, value): self._record("set_att")
    def set(self, **kwargs): self._record("set")
    def set_gain_mu(self, channel, gain): self._record("set_gain_mu")
    def write_dac(self, channel, value): self._record("write_dac")


class FakeCore:
    def __init__(self, log):
        self.log = log
        self.reset_calls = 0
        self.destination_status_calls = 0

    def reset(self):
        self.reset_calls += 1
        self.log.append(("core", "reset"))

    def break_realtime(self):
        self.log.append(("core", "break_realtime"))

    def get_rtio_destination_status(self, destination):
        self.destination_status_calls += 1
        self.log.append(("core", f"destination_{destination}"))
        return True


class FakeExperiment:
    def __init__(self):
        self.log = []
        self.dataset_reads = []
        self.dataset_writes = []
        self.setattr_device_calls = []
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

        with Path(
            "utilities/config/master_satellite/device_aliases.json"
        ).open() as config_file:
            mapping = json.load(config_file)
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                self.devices.setdefault(
                    unified_name, FakeDevice(unified_name, self.log)
                )
        self.devices["core_dma"] = FakeDevice("core_dma", self.log)
        self.devices["scheduler"] = FakeDevice("scheduler", self.log)

    def get_dataset(self, name):
        self.dataset_reads.append(name)
        if name not in self.datasets:
            raise KeyError(name)
        return self.datasets[name]

    def setattr_device(self, name):
        self.setattr_device_calls.append(name)
        if name not in self.devices:
            raise KeyError(name)
        setattr(self, name, self.devices[name])

    def set_dataset(self, name, value, **kwargs):
        self.dataset_writes.append((name, value, kwargs))
        self.datasets[name] = value


class BaseExperimentMasterSatelliteTests(unittest.TestCase):
    def build_and_prepare(self, mode, node=None):
        experiment = FakeExperiment()
        base = BaseExperimentMasterSatellite(experiment, mode, node)
        base.build()
        base.prepare()
        return experiment, base

    def test_deferred_build_binds_superset_then_loads_only_selected_node(self):
        experiment = FakeExperiment()
        base = BaseExperimentMasterSatellite(experiment)
        base.build()

        self.assertFalse(base._execution_configured)
        self.assertEqual(experiment.dataset_reads, [])
        self.assertEqual(experiment.sampler0_Node1.name, "sampler0")
        self.assertEqual(experiment.sampler0_Node2.name, "sampler3")

        base.configure_execution("single_node", "Node2")
        self.assertTrue(base._execution_configured)
        self.assertEqual(base.active_nodes, ("Node2",))
        self.assertEqual(base._cplds_node1, [])
        self.assertEqual(base._samplers_node1, [])
        self.assertEqual(base._zotinos_node1, [])
        self.assertEqual(base._ttl_outputs_node1, [])
        self.assertEqual(experiment.sampler0.name, "sampler3")
        self.assertEqual(experiment.f_FORT, experiment.f_FORT_Node2)
        self.assertFalse(hasattr(experiment, "f_FORT_Node1"))
        self.assertEqual(
            set(experiment.dataset_reads),
            {
                variable.name
                for variable in NODE2_VARIABLES + MASTER_SATELLITE_VARIABLES
            },
        )
        base.prepare()
        self.assertEqual(experiment.dds_FORT.name, "urukul4_ch3")

    def test_deferred_configuration_retains_strict_runtime_validation(self):
        for mode, node, message in (
            (None, None, "Unsupported master-satellite experiment_mode"),
            ("invalid", None, "Unsupported master-satellite experiment_mode"),
            ("single_node", None, "requires which_node"),
            ("single_node", "invalid", "requires which_node"),
            ("two_nodes", "Node1", "must be omitted"),
        ):
            with self.subTest(mode=mode, node=node):
                base = BaseExperimentMasterSatellite(FakeExperiment())
                with self.assertRaisesRegex(ValueError, message):
                    base.configure_execution(mode, node)

    def test_single_node1_bindings(self):
        experiment, base = self.build_and_prepare("single_node", "Node1")
        self.assertEqual(experiment.dds_FORT.name, "urukul0_ch0")
        self.assertEqual(experiment.sampler0.name, "sampler0")
        self.assertEqual(experiment.zotino0.name, "zotino0")
        self.assertEqual(experiment.coil_channels, [0, 1, 13, 14])
        self.assertEqual(experiment.Magnetometer_X_ch, 1)
        self.assertEqual(experiment.Magnetometer_Y_ch, 2)
        self.assertEqual(experiment.Magnetometer_Z_ch, 3)
        self.assertEqual(experiment.ttl_SPCM0.name, "ttl0")
        self.assertEqual(experiment.f_FORT, experiment.f_FORT_Node1)
        self.assertEqual(
            experiment.t_MOT_loading, experiment.t_MOT_loading_Node1
        )
        self.assertFalse(hasattr(experiment, "f_FORT_Node2"))
        self.assertEqual(experiment.n_measurements, 100)
        self.assertEqual(experiment.t_delay_in_bob_mu, 189)
        self.assertTrue(experiment.parallel_AOM_feedback)
        self.assertEqual(
            base.compatibility_variable_map["f_FORT"], "f_FORT_Node1"
        )

    def test_single_node2_uses_node2_devices_and_master_spcms(self):
        experiment, _ = self.build_and_prepare("single_node", "Node2")
        self.assertEqual(experiment.dds_FORT.name, "urukul4_ch3")
        self.assertEqual(experiment.sampler0.name, "sampler3")
        self.assertEqual(experiment.zotino0.name, "zotino1")
        self.assertEqual(experiment.coil_channels, [0, 1, 2, 3])
        self.assertEqual(experiment.Magnetometer_X_ch, 1)
        self.assertEqual(experiment.Magnetometer_Y_ch, 2)
        self.assertEqual(experiment.Magnetometer_Z_ch, 3)
        self.assertEqual(experiment.ttl0.name, "ttl16")
        self.assertEqual(experiment.ttl_SPCM0.name, "ttl0")
        self.assertEqual(experiment.ttl_SPCM0_counter.name, "ttl0_counter")
        self.assertEqual(experiment.f_FORT, experiment.f_FORT_Node2)
        self.assertEqual(
            experiment.t_MOT_loading, experiment.t_MOT_loading_Node2
        )
        self.assertFalse(hasattr(experiment, "f_FORT_Node1"))

    def test_single_node_wiring_metadata_is_complete_and_fixed(self):
        expected_common = {
            "coil_names": ["AZ bottom", "AZ top", "AX", "AY"],
            "AZ_bottom_Zotino_channel": 0,
            "AZ_top_Zotino_channel": 1,
            "UV_trig_channel": [8],
            "Osc_trig_channel": [10],
            "FORT_MM_sampler_ch": 7,
            "GRIN1_sampler_ch": 4,
            "Magnetometer_X_ch": 1,
            "Magnetometer_Y_ch": 2,
            "Magnetometer_Z_ch": 3,
        }
        for node, ax_channel, ay_channel in (
            ("Node1", 13, 14),
            ("Node2", 2, 3),
        ):
            experiment, _ = self.build_and_prepare("single_node", node)
            for name, value in expected_common.items():
                self.assertEqual(getattr(experiment, name), value)
            self.assertEqual(experiment.AX_Zotino_channel, ax_channel)
            self.assertEqual(experiment.AY_Zotino_channel, ay_channel)
            self.assertEqual(
                experiment.coil_channels,
                [0, 1, ax_channel, ay_channel],
            )
            self.assertEqual(
                experiment.measurements_progress, "measurements_progress"
            )

    def test_magnetometer_result_reset_is_transient_and_side_effect_free(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        reads_before = list(experiment.dataset_reads)
        device_calls_before = list(experiment.setattr_device_calls)
        reset_calls_before = experiment.core.reset_calls

        base.initialize_result_datasets()
        experiment.datasets["Magnetometer_MOT_X_Node2"] = [123.0]
        base.reset_result_state_for_scan_point()

        expected_names = {
            "measurements_progress",
            *(
                f"{name}_Node2"
                for name in base.MAGNETOMETER_RESULT_DATASETS
            ),
        }
        written_names = {
            name for name, _, _ in experiment.dataset_writes
        }
        self.assertEqual(written_names, expected_names)
        self.assertEqual(experiment.datasets["measurements_progress"], 0.0)
        for name in base.MAGNETOMETER_RESULT_DATASETS:
            self.assertEqual(experiment.datasets[f"{name}_Node2"], [0.0])
            self.assertNotIn(name, experiment.datasets)
        self.assertTrue(
            all(
                kwargs == {"broadcast": True, "persist": False}
                for _, _, kwargs in experiment.dataset_writes
            )
        )
        self.assertFalse(
            any(
                name.startswith(("f_", "p_"))
                for name, _, _ in experiment.dataset_writes
            )
        )
        self.assertEqual(experiment.dataset_reads, reads_before)
        self.assertEqual(experiment.setattr_device_calls, device_calls_before)
        self.assertEqual(experiment.core.reset_calls, reset_calls_before)

    def test_magnetometer_result_names_follow_selected_node(self):
        for node in ("Node1", "Node2"):
            experiment, base = self.build_and_prepare("single_node", node)
            for legacy_name in base.MAGNETOMETER_RESULT_DATASETS:
                self.assertEqual(
                    base.resolve_result_dataset_name(legacy_name),
                    f"{legacy_name}_{node}",
                )
            self.assertEqual(
                base.resolve_result_dataset_name("measurements_progress"),
                "measurements_progress",
            )

            base.initialize_result_datasets()
            written_names = {
                name for name, _, _ in experiment.dataset_writes
            }
            self.assertIn("measurements_progress", written_names)
            self.assertFalse(
                set(base.MAGNETOMETER_RESULT_DATASETS) & written_names
            )
            self.assertTrue(
                {
                    f"{name}_{node}"
                    for name in base.MAGNETOMETER_RESULT_DATASETS
                }.issubset(written_names)
            )

    def test_two_node_bindings(self):
        experiment, _ = self.build_and_prepare("two_nodes")
        self.assertEqual(experiment.dds_FORT_Node1.name, "urukul0_ch0")
        self.assertEqual(experiment.dds_FORT_Node2.name, "urukul4_ch3")
        self.assertEqual(experiment.sampler0_Node1.name, "sampler0")
        self.assertEqual(experiment.sampler0_Node2.name, "sampler3")
        self.assertEqual(experiment.SPCM_H1.name, "ttl0")
        self.assertFalse(hasattr(experiment, "ttl_SPCM0"))
        self.assertTrue(hasattr(experiment, "f_FORT_Node1"))
        self.assertTrue(hasattr(experiment, "f_FORT_Node2"))
        self.assertFalse(hasattr(experiment, "f_FORT"))

    def test_single_node_refresh_updates_only_in_memory_projection(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        reads_before = list(experiment.dataset_reads)
        experiment.f_FORT_Node2 = 249e6

        base.refresh_compatibility_variables()

        self.assertEqual(experiment.f_FORT, 249e6)
        self.assertEqual(experiment.dataset_reads, reads_before)
        self.assertEqual(experiment.datasets["f_FORT_Node2"], 240e6)

    def test_two_node_mode_does_not_publish_legacy_dataset_value(self):
        experiment = FakeExperiment()
        experiment.datasets["f_FORT"] = 999e6
        base = BaseExperimentMasterSatellite(experiment, "two_nodes")
        base.build()
        base.prepare()

        self.assertFalse(hasattr(experiment, "f_FORT"))
        self.assertEqual(experiment.f_FORT_Node1, 245e6)
        self.assertEqual(experiment.f_FORT_Node2, 240e6)

    def test_missing_dataset_identifies_owner_and_name(self):
        experiment = FakeExperiment()
        del experiment.datasets["f_FORT_Node2"]
        base = BaseExperimentMasterSatellite(
            experiment, "single_node", "Node2"
        )

        with self.assertRaisesRegex(
            RuntimeError,
            r"ExperimentVariables_Node2\.py: f_FORT_Node2",
        ):
            base.build()

    def test_resolves_single_node_targets_and_globals(self):
        _, node1_base = self.build_and_prepare("single_node", "Node1")
        self.assertEqual(
            node1_base.resolve_experiment_variable_target("f_FORT"),
            "f_FORT_Node1",
        )
        self.assertEqual(
            node1_base.resolve_experiment_variable_target("f_FORT_Node1"),
            "f_FORT_Node1",
        )
        self.assertEqual(
            node1_base.resolve_experiment_variable_target("n_measurements"),
            "n_measurements",
        )

        _, node2_base = self.build_and_prepare("single_node", "Node2")
        self.assertEqual(
            node2_base.resolve_experiment_variable_target("f_FORT"),
            "f_FORT_Node2",
        )
        self.assertEqual(
            node2_base.resolve_experiment_variable_target("f_FORT_Node2"),
            "f_FORT_Node2",
        )
        self.assertEqual(
            node2_base.resolve_experiment_variable_target("n_measurements"),
            "n_measurements",
        )

    def test_single_node_rejects_other_node_and_unknown_targets(self):
        _, base = self.build_and_prepare("single_node", "Node2")
        with self.assertRaisesRegex(ValueError, r"belongs to Node1"):
            base.resolve_experiment_variable_target("f_FORT_Node1")
        with self.assertRaisesRegex(ValueError, r"Unknown.*not_a_variable"):
            base.resolve_experiment_variable_target("not_a_variable")

    def test_resolves_two_node_targets_and_rejects_ambiguity(self):
        _, base = self.build_and_prepare("two_nodes")
        self.assertEqual(
            base.resolve_experiment_variable_target("f_FORT_Node1"),
            "f_FORT_Node1",
        )
        self.assertEqual(
            base.resolve_experiment_variable_target("f_FORT_Node2"),
            "f_FORT_Node2",
        )
        self.assertEqual(
            base.resolve_experiment_variable_target("n_measurements"),
            "n_measurements",
        )
        with self.assertRaisesRegex(
            ValueError, r"Ambiguous.*f_FORT.*f_FORT_Node1.*f_FORT_Node2"
        ):
            base.resolve_experiment_variable_target("f_FORT")

    def test_refresh_variable_dependent_state_updates_node2_dds_cache(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        resolver = base.node_resolvers["Node2"]
        fort_index = next(
            index
            for index, binding in enumerate(resolver.dds_bindings)
            if binding["logical_alias"] == "dds_FORT"
        )
        experiment.f_FORT_Node2 = 249e6
        experiment.p_FORT_loading_Node2 = -9.5
        base.refresh_compatibility_variables()
        reads_before = list(experiment.dataset_reads)
        writes_before = list(experiment.dataset_writes)

        base.refresh_variable_dependent_state()

        self.assertEqual(experiment.f_FORT, 249e6)
        self.assertEqual(base._dds_frequencies_node2[fort_index], 249e6)
        self.assertEqual(resolver.dds_frequencies[fort_index], 249e6)
        self.assertEqual(base._dds_powers_node2[fort_index], -9.5)
        self.assertEqual(resolver.dds_powers[fort_index], -9.5)
        self.assertEqual(experiment.dataset_reads, reads_before)
        self.assertEqual(experiment.dataset_writes, writes_before)
        self.assertEqual(experiment.core.reset_calls, 0)

    def test_reload_updates_attributes_and_cache_without_side_effects(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        resolver = base.node_resolvers["Node2"]
        fort_index = next(
            index
            for index, binding in enumerate(resolver.dds_bindings)
            if binding["logical_alias"] == "dds_FORT"
        )
        device_calls_before = list(experiment.setattr_device_calls)
        dataset_writes_before = list(experiment.dataset_writes)
        experiment.datasets["f_FORT_Node2"] = 247e6

        base.reload_experiment_variables()

        self.assertEqual(experiment.f_FORT_Node2, 247e6)
        self.assertEqual(experiment.f_FORT, 247e6)
        self.assertEqual(base._dds_frequencies_node2[fort_index], 247e6)
        self.assertEqual(experiment.core.reset_calls, 0)
        self.assertEqual(experiment.setattr_device_calls, device_calls_before)
        self.assertEqual(experiment.dataset_writes, dataset_writes_before)

    def test_reload_fails_if_required_dataset_disappears(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        del experiment.datasets["f_FORT_Node2"]

        with self.assertRaisesRegex(
            RuntimeError,
            r"ExperimentVariables_Node2\.py: f_FORT_Node2",
        ):
            base.reload_experiment_variables()

    def test_node2_initialization_waits_after_one_reset(self):
        experiment, base = self.build_and_prepare("single_node", "Node2")
        base.initialize_hardware()
        self.assertEqual(experiment.core.reset_calls, 1)
        self.assertEqual(experiment.core.destination_status_calls, 1)
        reset_index = experiment.log.index(("core", "reset"))
        ready_index = experiment.log.index(("core", "destination_1"))
        first_hardware_index = next(
            i
            for i, event in enumerate(experiment.log)
            if event[0] != "core"
        )
        self.assertLess(reset_index, ready_index)
        self.assertLess(ready_index, first_hardware_index)


if __name__ == "__main__":
    unittest.main()
