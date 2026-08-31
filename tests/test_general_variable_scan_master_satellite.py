import sys
import types
import unittest
import json
from pathlib import Path


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


if "artiq.experiment" not in sys.modules:
    artiq = types.ModuleType("artiq")
    experiment = types.ModuleType("artiq.experiment")
    experiment.EnvExperiment = object
    experiment.kernel = _identity_decorator
    experiment.rpc = _identity_decorator
    experiment.delay = lambda duration: None
    experiment.NumberValue = lambda value, **kwargs: value
    experiment.BooleanValue = lambda value, **kwargs: value
    experiment.StringValue = lambda value, **kwargs: value
    experiment.EnumerationValue = lambda values, **kwargs: tuple(values)[0]
    experiment.TBool = bool
    experiment.TFloat = float
    experiment.TInt32 = int
    experiment.TInt64 = int
    experiment.TStr = str
    experiment.MHz = 1e6
    experiment.kHz = 1e3
    experiment.ms = 1e-3
    experiment.us = 1e-6
    experiment.ns = 1e-9
    experiment.s = 1.0
    experiment.V = 1.0
    artiq.experiment = experiment
    sys.modules["artiq"] = artiq
    sys.modules["artiq.experiment"] = experiment
else:
    experiment = sys.modules["artiq.experiment"]
    supplemental_exports = {
        "kernel": _identity_decorator,
        "rpc": _identity_decorator,
        "delay": lambda duration: None,
        "EnvExperiment": object,
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value, **kwargs: value,
        "StringValue": lambda value, **kwargs: value,
        "EnumerationValue": lambda values, **kwargs: tuple(values)[0],
        "TBool": bool,
        "TFloat": float,
        "TInt32": int,
        "TInt64": int,
        "TStr": str,
        "MHz": 1e6,
        "kHz": 1e3,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
        "s": 1.0,
        "V": 1.0,
    }
    for name, value in supplemental_exports.items():
        if name in {"kernel", "rpc"} or not hasattr(experiment, name):
            setattr(experiment, name, value)
    if hasattr(experiment, "__all__"):
        experiment.__all__ = sorted(
            set(experiment.__all__) | set(supplemental_exports)
        )


coredevice = types.ModuleType("artiq.coredevice")
exceptions = types.ModuleType("artiq.coredevice.exceptions")


class RTIOUnderflow(Exception):
    pass


exceptions.RTIOUnderflow = RTIOUnderflow
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
language = types.ModuleType("artiq.language")
language.us = experiment.us
language.ns = experiment.ns
language.MHz = experiment.MHz
sys.modules.setdefault("artiq.coredevice", coredevice)
sys.modules.setdefault("artiq.coredevice.exceptions", exceptions)
sys.modules.setdefault("artiq.coredevice.ad9910", ad9910)
sys.modules.setdefault("artiq.coredevice.urukul", urukul)
sys.modules.setdefault("artiq.language", language)
sys.modules.setdefault("pyvisa", types.ModuleType("pyvisa"))


from GeneralVariableScan_master_satellite import (  # noqa: E402
    GeneralVariableScanMasterSatellite,
    HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS,
    build_experiment_function_registry,
    build_single_node_function_registry,
    build_two_node_function_registry,
)
from GeneralVariableScan_CatchError_master_satellite import (  # noqa: E402
    GeneralVariableScan_CatchError_master_satellite,
)
from subroutines.experiment_functions_master_satellite import (  # noqa: E402
    MASTER_SATELLITE_SANITY_ATTRIBUTES,
    master_satellite_namespace_sanity_experiment,
)
from AOMsCoils_master_satellite import AOMsCoils_master_satellite  # noqa: E402
from ExperimentVariables_Node1 import (  # noqa: E402
    ExperimentVariablesNode1,
    NODE1_VARIABLES,
)
from ExperimentVariables_Node2 import (  # noqa: E402
    ExperimentVariablesNode2,
    NODE2_VARIABLES,
)
from ExperimentVariables_master_satellite import (  # noqa: E402
    ExperimentVariablesMasterSatellite,
    MASTER_SATELLITE_VARIABLES,
)


class FakeBase:
    def __init__(self, mode="single_node", node="Node2"):
        self.mode = mode
        self.node = node
        self.resolve_calls = []
        self.compatibility_refreshes = 0
        self.dependent_refreshes = 0
        self.reload_calls = 0
        self.build_calls = 0
        self.prepare_calls = 0
        self.result_initializations = 0
        self.result_resets = 0

    def resolve_experiment_variable_target(self, name):
        self.resolve_calls.append(name)
        globals_ = {"n_measurements", "t_delay_in_bob_mu", "parallel_AOM_feedback"}
        if name in globals_:
            return name
        if self.mode == "single_node":
            other = "Node1" if self.node == "Node2" else "Node2"
            if name.endswith(f"_{other}"):
                raise ValueError("other node")
            if name.endswith(f"_{self.node}"):
                return name
            if name == "f_FORT":
                return f"f_FORT_{self.node}"
            raise ValueError("unknown")
        if name in {"f_FORT_Node1", "f_FORT_Node2"}:
            return name
        if name == "f_FORT":
            raise ValueError("ambiguous")
        raise ValueError("unknown")

    def resolve_result_dataset_name(self, name):
        magnetometer_names = {
            f"Magnetometer_{phase}_{axis}"
            for phase in ("Zero", "OP", "MOT")
            for axis in ("X", "Y", "Z")
        }
        if self.mode == "single_node" and name in magnetometer_names:
            return f"{name}_{self.node}"
        return name

    def refresh_compatibility_variables(self):
        self.compatibility_refreshes += 1
        experiment = self.experiment
        if self.mode == "single_node":
            experiment.f_FORT = getattr(experiment, f"f_FORT_{self.node}")

    def refresh_variable_dependent_state(self):
        self.dependent_refreshes += 1
        self.dds_frequency_cache = getattr(
            self.experiment, f"f_FORT_{self.node}"
        )

    def reload_experiment_variables(self):
        self.reload_calls += 1
        setattr(self.experiment, f"f_FORT_{self.node}", 250.0)
        self.refresh_compatibility_variables()
        self.refresh_variable_dependent_state()

    def initialize_result_datasets(self):
        self.result_initializations += 1

    def reset_result_state_for_scan_point(self):
        self.result_resets += 1


class GeneralVariableScanMasterSatelliteTests(unittest.TestCase):
    @staticmethod
    def _make_repository_examination_experiment(experiment_class):
        experiment_instance = experiment_class()
        experiment_instance.dataset_reads = []
        experiment_instance.datasets = {
            variable.name: variable.value
            for variable in (
                NODE1_VARIABLES
                + NODE2_VARIABLES
                + MASTER_SATELLITE_VARIABLES
            )
        }
        with Path(
            "utilities/config/master_satellite/device_aliases.json"
        ).open() as config_file:
            mapping = json.load(config_file)
        devices = {"core": object(), "core_dma": object(), "scheduler": object()}
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                devices.setdefault(unified_name, object())

        def get_dataset(name):
            experiment_instance.dataset_reads.append(name)
            if name not in experiment_instance.datasets:
                raise KeyError(name)
            return experiment_instance.datasets[name]

        experiment_instance.get_dataset = get_dataset
        experiment_instance.setattr_device = lambda name: setattr(
            experiment_instance, name, devices[name]
        )
        # This is the repository examiner behavior that triggered the bug:
        # argument metadata is registered, but submitted values are all None.
        experiment_instance.setattr_argument = lambda name, *args, **kwargs: setattr(
            experiment_instance, name, None
        )
        return experiment_instance

    def test_all_public_master_satellite_experiments_build_with_none_arguments(self):
        public_experiments = (
            ExperimentVariablesNode1,
            ExperimentVariablesNode2,
            ExperimentVariablesMasterSatellite,
            GeneralVariableScanMasterSatellite,
            GeneralVariableScan_CatchError_master_satellite,
            AOMsCoils_master_satellite,
        )
        for experiment_class in public_experiments:
            with self.subTest(experiment=experiment_class.__name__):
                instance = self._make_repository_examination_experiment(
                    experiment_class
                )
                instance.build()

        scan = self._make_repository_examination_experiment(
            GeneralVariableScanMasterSatellite
        )
        scan.build()
        self.assertFalse(scan.base._execution_configured)
        self.assertEqual(scan.dataset_reads, ["n_measurements"])

    def test_single_node_registry_uses_current_functions_and_exclusions(self):
        # Other hardware-free test modules may narrow the wildcard-exported
        # names on the shared ARTIQ stub during unittest discovery.
        required_exports = {
            "kernel": _identity_decorator,
            "rpc": _identity_decorator,
            "TBool": bool,
            "TFloat": float,
            "TInt32": int,
            "TInt64": int,
            "TStr": str,
            "delay": lambda duration: None,
            "MHz": 1e6,
            "kHz": 1e3,
            "ms": 1e-3,
            "us": 1e-6,
            "ns": 1e-9,
            "s": 1.0,
        }
        artiq_experiment = sys.modules["artiq.experiment"]
        for name, value in required_exports.items():
            if not hasattr(artiq_experiment, name):
                setattr(artiq_experiment, name, value)
        if hasattr(artiq_experiment, "__all__"):
            artiq_experiment.__all__ = sorted(
                set(artiq_experiment.__all__) | set(required_exports)
            )

        registry = build_single_node_function_registry()
        self.assertIn("atom_loading_experiment", registry)
        self.assertIn("test_ttl_pulse_experiment", registry)
        self.assertEqual(len(registry), 50)
        for name in HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS:
            self.assertNotIn(name, registry)
        self.assertNotIn("master_satellite_namespace_sanity_experiment", registry)

    def test_registry_excludes_imported_function(self):
        module = types.ModuleType("test_registry_module")

        def local_experiment():
            pass

        def imported_experiment():
            pass

        local_experiment.__module__ = module.__name__
        imported_experiment.__module__ = "another_module"
        module.local_experiment = local_experiment
        module.imported_experiment = imported_experiment

        self.assertEqual(
            build_experiment_function_registry(module),
            {"local_experiment": local_experiment},
        )

    def test_native_registry_contains_only_native_sanity_function(self):
        registry = build_two_node_function_registry()
        self.assertEqual(
            registry,
            {
                "master_satellite_namespace_sanity_experiment":
                    master_satellite_namespace_sanity_experiment
            },
        )
        for name in HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS:
            self.assertNotIn(name, registry)

    def test_wrong_mode_or_unknown_function_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "not available"):
            GeneralVariableScanMasterSatellite._select_experiment_function(
                build_two_node_function_registry(),
                "atom_loading_experiment",
            )

    def test_selected_node_has_separate_legacy_presentation(self):
        scan = GeneralVariableScanMasterSatellite()
        scan.experiment_mode = "single_node"
        scan.selected_node = "Node1"
        scan._publish_legacy_node_compatibility()
        self.assertEqual(scan.selected_node, "Node1")
        self.assertEqual(scan.which_node, "alice")

        scan.selected_node = "Node2"
        scan._publish_legacy_node_compatibility()
        self.assertEqual(scan.selected_node, "Node2")
        self.assertEqual(scan.which_node, "bob")

    def test_scan_target_resolution_is_delegated_to_base(self):
        for node, expected in (
            ("Node1", "f_FORT_Node1"),
            ("Node2", "f_FORT_Node2"),
        ):
            base = FakeBase("single_node", node)
            self.assertEqual(
                base.resolve_experiment_variable_target("f_FORT"), expected
            )
            self.assertEqual(
                base.resolve_experiment_variable_target("n_measurements"),
                "n_measurements",
            )

        base = FakeBase("two_nodes", None)
        self.assertEqual(
            base.resolve_experiment_variable_target("f_FORT_Node1"),
            "f_FORT_Node1",
        )
        self.assertEqual(
            base.resolve_experiment_variable_target("f_FORT_Node2"),
            "f_FORT_Node2",
        )
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            base.resolve_experiment_variable_target("f_FORT")

    def test_scan_and_override_refresh_authoritative_state_only(self):
        scan = GeneralVariableScanMasterSatellite()
        scan.f_FORT_Node2 = 240.0
        scan.f_FORT = 240.0
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.authoritative_overrides = {"f_FORT_Node2": 245.0}
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None

        scan._apply_run_wide_overrides()
        self.assertEqual(scan.f_FORT_Node2, 245.0)
        self.assertEqual(scan.f_FORT, 245.0)
        self.assertEqual(scan.base.dds_frequency_cache, 245.0)

        scan._apply_scan_point(247.0)
        self.assertEqual(scan.f_FORT_Node2, 247.0)
        self.assertEqual(scan.f_FORT, 247.0)
        self.assertEqual(scan.base.dds_frequency_cache, 247.0)

    def test_magnetometer_append_names_resolve_to_selected_node(self):
        for node in ("Node1", "Node2"):
            base = FakeBase("single_node", node)
            self.assertEqual(
                base.resolve_result_dataset_name("Magnetometer_MOT_X"),
                f"Magnetometer_MOT_X_{node}",
            )
            self.assertEqual(
                base.resolve_result_dataset_name("measurements_progress"),
                "measurements_progress",
            )

    def test_run_uses_reload_and_never_rebuilds_prepares_or_persists_scan(self):
        scan = GeneralVariableScanMasterSatellite()
        scan.f_FORT_Node2 = 240.0
        scan.f_FORT = 240.0
        scan.n_measurements = 12
        scan.execution_n_measurements = 12
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.needs_experiment_variable_reload = True
        scan.authoritative_overrides = {"f_FORT_Node2": 245.0}
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None
        scan.scan_sequence1 = [246.0, 247.0]
        scan.scan_sequence2 = [0.0]
        scan.hardware_initializations = 0
        scan.function_calls = 0
        scan.dataset_writes = []
        scan.initialize_hardware = lambda: setattr(
            scan,
            "hardware_initializations",
            scan.hardware_initializations + 1,
        )
        scan._selected_experiment_function = lambda experiment: setattr(
            experiment,
            "function_calls",
            experiment.function_calls + 1,
        )
        scan.set_dataset = lambda name, value, **kwargs: scan.dataset_writes.append(
            (name, value, kwargs)
        )

        scan.run()

        self.assertEqual(scan.base.reload_calls, 1)
        self.assertEqual(scan.base.build_calls, 0)
        self.assertEqual(scan.base.prepare_calls, 0)
        self.assertEqual(scan.hardware_initializations, 2)
        self.assertEqual(scan.function_calls, 2)
        self.assertEqual(scan.base.result_initializations, 1)
        self.assertEqual(scan.base.result_resets, 2)
        self.assertEqual(scan.f_FORT_Node2, 247.0)
        self.assertEqual(scan.f_FORT, 247.0)
        self.assertEqual(scan.base.dds_frequency_cache, 247.0)
        self.assertEqual(scan.n_measurements, 12)
        self.assertTrue(all(name == "iteration" for name, _, _ in scan.dataset_writes))
        self.assertTrue(
            all(not kwargs.get("persist", False) for _, _, kwargs in scan.dataset_writes)
        )

    def test_catch_error_retries_only_failed_scan_point(self):
        scan = GeneralVariableScan_CatchError_master_satellite()
        scan.scan_sequence1 = [1.0, 2.0]
        scan.scan_sequence2 = [0.0]
        scan.underflow_max_retries = 3
        scan.underflow_backoff_ms = 0.0
        scan.skip_only_that_iteration_if_exhausted = False
        scan.dataset_writes = []
        scan.attempts = []
        scan._initialize_run_state = lambda: None
        scan.set_dataset = lambda name, value, **kwargs: scan.dataset_writes.append(
            (name, value, kwargs)
        )

        def run_point(value1, value2, iteration):
            scan.attempts.append((value1, iteration))
            if value1 == 1.0 and len(scan.attempts) == 1:
                raise RTIOUnderflow("transient")

        scan._run_scan_point = run_point
        scan.run()

        self.assertEqual(scan.attempts, [(1.0, 0), (1.0, 0), (2.0, 1)])
        self.assertTrue(
            all(not kwargs.get("persist", False) for _, _, kwargs in scan.dataset_writes)
        )

    def test_catch_error_does_not_hide_non_underflow_errors(self):
        scan = GeneralVariableScan_CatchError_master_satellite()
        scan.scan_sequence1 = [1.0]
        scan.scan_sequence2 = [0.0]
        scan.underflow_max_retries = 3
        scan.underflow_backoff_ms = 0.0
        scan.skip_only_that_iteration_if_exhausted = True
        scan._initialize_run_state = lambda: None
        scan.set_dataset = lambda *args, **kwargs: None
        scan._run_scan_point = lambda *args: (_ for _ in ()).throw(
            RuntimeError("DRTIO unavailable")
        )

        with self.assertRaisesRegex(RuntimeError, "DRTIO unavailable"):
            scan.run()

    def test_namespace_sanity_is_non_destructive(self):
        experiment_object = types.SimpleNamespace()
        devices = []
        for name in MASTER_SATELLITE_SANITY_ATTRIBUTES:
            device = object()
            devices.append(device)
            setattr(experiment_object, name, device)

        self.assertTrue(
            master_satellite_namespace_sanity_experiment(experiment_object)
        )
        self.assertEqual(
            [getattr(experiment_object, name) for name in MASTER_SATELLITE_SANITY_ATTRIBUTES],
            devices,
        )

        del experiment_object.dds_FORT_Node2
        del experiment_object.SPCM_H2_counter
        with self.assertRaisesRegex(
            RuntimeError, r"dds_FORT_Node2.*SPCM_H2_counter"
        ):
            master_satellite_namespace_sanity_experiment(experiment_object)


if __name__ == "__main__":
    unittest.main()
