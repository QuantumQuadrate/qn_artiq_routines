import sys
import types
import unittest
import json
import ast
from pathlib import Path


def _identity_decorator(function=None, **kwargs):
    if function is not None:
        return function
    return lambda decorated: decorated


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
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value=False, **kwargs: value,
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
    for name, value in _stub_exports.items():
        if not hasattr(_experiment_stub, name):
            setattr(_experiment_stub, name, value)
    if hasattr(_experiment_stub, "__all__"):
        _experiment_stub.__all__ = sorted(
            set(_experiment_stub.__all__) | set(_stub_exports)
        )

    class RTIOUnderflow(Exception):
        pass

    coredevice = types.ModuleType("artiq.coredevice")
    exceptions = types.ModuleType("artiq.coredevice.exceptions")
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
    language.us = 1e-6
    language.ns = 1e-9
    language.MHz = 1e6
    sys.modules.setdefault("artiq.coredevice", coredevice)
    sys.modules.setdefault("artiq.coredevice.exceptions", exceptions)
    sys.modules.setdefault("artiq.coredevice.ad9910", ad9910)
    sys.modules.setdefault("artiq.coredevice.urukul", urukul)
    sys.modules.setdefault("artiq.language", language)
    sys.modules.setdefault("pyvisa", types.ModuleType("pyvisa"))
    RTIOUnderflow = sys.modules["artiq.coredevice.exceptions"].RTIOUnderflow
else:
    from artiq.coredevice.exceptions import RTIOUnderflow  # noqa: F401


from GeneralVariableScan_master_satellite_mixin import (  # noqa: E402
    _GeneralVariableScanMasterSatelliteMixin,
    HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS,
    build_experiment_function_registry,
    build_single_node_function_registry,
    build_two_node_function_registry,
)
from GeneralVariableScan_master_satellite_single_node import (  # noqa: E402
    GeneralVariableScan_master_satellite_single_node,
)
from GeneralVariableScan_master_satellite_two_nodes import (  # noqa: E402
    GeneralVariableScan_master_satellite_two_nodes,
)
from subroutines.experiment_functions_two_nodes import (  # noqa: E402
    MASTER_SATELLITE_SANITY_ATTRIBUTES,
    master_satellite_namespace_sanity_experiment,
)
from AOMsCoils_master_satellite_Node1 import (  # noqa: E402
    AOMsCoils_master_satellite_Node1,
)
from AOMsCoils_master_satellite_Node2 import (  # noqa: E402
    AOMsCoils_master_satellite_Node2,
)
from ExperimentVariables_master_satellite_Node1 import (  # noqa: E402
    ExperimentVariablesMasterSatelliteNode1,
    NODE1_VARIABLES,
)
from ExperimentVariables_master_satellite_Node2 import (  # noqa: E402
    ExperimentVariablesMasterSatelliteNode2,
    NODE2_VARIABLES,
)
from ExperimentVariables_master_satellite_global import (  # noqa: E402
    ExperimentVariablesMasterSatelliteGlobal,
    MASTER_SATELLITE_VARIABLES,
)


class _FakeStabilizerChannel:
    def __init__(self, dataset, dB_dataset):
        self.dataset = dataset
        self.dB_dataset = dB_dataset
        self.dB_history_dataset = dB_dataset + "_history"


class _FakeStabilizer:
    """Stand-in for AOMPowerStabilizer.

    The real one reads utilities/config/<node>/feedback_channels.json via a
    path built from the process cwd, which is meaningless under the test
    runner; only the attribute surface Base touches is reproduced here.
    """

    def __init__(self, experiment, dds_names, iterations, averages, **kwargs):
        self.exp = experiment
        self.dds_names = dds_names
        self.iterations = iterations
        self.averages = averages
        self.kwargs = kwargs
        self.all_channels = [
            _FakeStabilizerChannel("MOT1_monitor", "p_AOM_A1"),
            _FakeStabilizerChannel("FORT_monitor", "p_FORT_loading"),
        ]

    def run(self):
        pass

    def monitor(self):
        pass


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
        self.single_node_result_initializations = 0
        self.stabilizer_preparations = 0
        self.result_resets = 0
        self.force_off_calls = 0

    def resolve_experiment_variable_target(self, name, node=None):
        """Mirror the real resolver, including its node-aware override path.

        node is the node that OWNS the name -- the per-node override
        dictionary it was written in -- and makes a bare node-specific name
        unambiguous. A scan variable passes no node, so a bare name stays
        ambiguous in two_nodes mode there.
        """
        self.resolve_calls.append(name)
        globals_ = {"n_measurements", "t_delay_in_bob_mu", "parallel_AOM_feedback"}
        if name in globals_:
            return name
        if node is not None:
            if name.endswith(f"_{node}"):
                return name
            if name == "f_FORT":
                return f"f_FORT_{node}"
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

    def initialize_single_node_result_state(self):
        # Single-node GVS also builds the full atom-physics result surface,
        # not just the magnetometer datasets.
        self.single_node_result_initializations += 1

    def prepare_laser_stabilizer(self, stabilizer_factory=None):
        # Single-node GVS builds the stabilizer so the reused experiment
        # functions find self.stabilizer_<channel>.
        self.stabilizer_preparations += 1

    def reset_result_state_for_scan_point(self):
        self.result_resets += 1

    def force_other_node_off(self):
        # Counted so a test can assert the unselected node is driven safe
        # once per RUN rather than once per scan point.
        self.force_off_calls += 1

    def result_name_tag(self):
        # The mixin builds the write_results name dict before calling, so
        # this runs even when write_results itself is stubbed out.
        if self.mode == "single_node":
            return self.node
        return "TwoNodes"


class GeneralVariableScanMasterSatelliteTests(unittest.TestCase):
    @staticmethod
    def _make_repository_examination_experiment(experiment_class):
        experiment_instance = experiment_class()
        experiment_instance.dataset_reads = []
        experiment_instance.datasets = {}
        with Path(
            "utilities/config/master_satellite/device_aliases.json"
        ).open() as config_file:
            mapping = json.load(config_file)
        devices = {"core": object(), "core_dma": object(), "scheduler": object()}
        for node_mapping in mapping.values():
            for unified_name in node_mapping.values():
                devices.setdefault(unified_name, object())

        def get_dataset(name, *args, **kwargs):
            experiment_instance.dataset_reads.append(name)
            if name not in experiment_instance.datasets:
                raise KeyError(name)
            return experiment_instance.datasets[name]

        def set_dataset(name, value, *args, **kwargs):
            experiment_instance.datasets[name] = value

        def append_to_dataset(name, value):
            experiment_instance.datasets.setdefault(name, []).append(value)

        experiment_instance.get_dataset = get_dataset
        # prepare() now builds the laser stabilizer in single-node mode, which
        # writes feedbackchannels and seeds the per-channel dB histories.
        experiment_instance.set_dataset = set_dataset
        experiment_instance.append_to_dataset = append_to_dataset
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
            ExperimentVariablesMasterSatelliteNode1,
            ExperimentVariablesMasterSatelliteNode2,
            ExperimentVariablesMasterSatelliteGlobal,
            GeneralVariableScan_master_satellite_single_node,
            GeneralVariableScan_master_satellite_two_nodes,
            AOMsCoils_master_satellite_Node1,
            AOMsCoils_master_satellite_Node2,
        )
        for experiment_class in public_experiments:
            with self.subTest(experiment=experiment_class.__name__):
                instance = self._make_repository_examination_experiment(
                    experiment_class
                )
                instance.build()

        scan = self._make_repository_examination_experiment(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.build()
        self.assertFalse(scan.base._execution_configured)
        self.assertEqual(scan.dataset_reads, [])

        manual = self._make_repository_examination_experiment(
            AOMsCoils_master_satellite_Node1
        )
        manual.build()
        self.assertFalse(manual.base._execution_configured)
        self.assertEqual(manual.dataset_reads, [])

    def test_empty_examination_does_not_weaken_runtime_dataset_validation(self):
        for experiment_class, configure in (
            (
                GeneralVariableScan_master_satellite_single_node,
                lambda instance: (
                    setattr(instance, "selected_node", "Node1"),
                ),
            ),
            (
                GeneralVariableScan_master_satellite_two_nodes,
                lambda instance: None,
            ),
            (
                AOMsCoils_master_satellite_Node1,
                lambda instance: None,
            ),
            (
                AOMsCoils_master_satellite_Node2,
                lambda instance: None,
            ),
        ):
            with self.subTest(experiment=experiment_class.__name__):
                instance = self._make_repository_examination_experiment(
                    experiment_class
                )
                instance.build()
                configure(instance)
                with self.assertRaisesRegex(
                    RuntimeError,
                    "Missing required master-satellite persistent datasets",
                ):
                    instance.prepare()

    def test_public_gvs_classes_have_fixed_distinct_modes(self):
        self.assertEqual(
            GeneralVariableScan_master_satellite_single_node.EXPERIMENT_MODE,
            "single_node",
        )
        self.assertEqual(
            GeneralVariableScan_master_satellite_two_nodes.EXPERIMENT_MODE,
            "two_nodes",
        )
        self.assertFalse(hasattr(_GeneralVariableScanMasterSatelliteMixin, "build"))

        single_registry = (
            GeneralVariableScan_master_satellite_single_node()
            ._active_function_registry()
        )
        two_registry = (
            GeneralVariableScan_master_satellite_two_nodes()
            ._active_function_registry()
        )
        self.assertIn("atom_loading_experiment", single_registry)
        self.assertNotIn(
            "master_satellite_namespace_sanity_experiment", single_registry
        )
        self.assertEqual(
            set(two_registry),
            {"master_satellite_namespace_sanity_experiment"},
        )

    def test_public_gvs_envexperiment_classes_match_the_split_files(self):
        expected = {
            "GeneralVariableScan_master_satellite_mixin.py": set(),
            "GeneralVariableScan_master_satellite_single_node.py": {
                "GeneralVariableScan_master_satellite_single_node"
            },
            "GeneralVariableScan_master_satellite_two_nodes.py": {
                "GeneralVariableScan_master_satellite_two_nodes"
            },
        }
        for filename, expected_classes in expected.items():
            tree = ast.parse(Path(filename).read_text())
            public_experiment_classes = {
                node.name
                for node in tree.body
                if isinstance(node, ast.ClassDef)
                and not node.name.startswith("_")
                and any(
                    isinstance(base, ast.Name) and base.id == "EnvExperiment"
                    for base in node.bases
                )
            }
            self.assertEqual(public_experiment_classes, expected_classes)

        for retired_filename in (
            "GeneralVariableScan_master_satellite.py",
            "GeneralVariableScan_CatchError_master_satellite.py",
            "GeneralVariableScan_CatchError_master_satellite_single_node.py",
            "GeneralVariableScan_CatchError_master_satellite_two_nodes.py",
            "AOMsCoils_master_satellite.py",
            "subroutines/experiment_functions_master_satellite.py",
        ):
            self.assertFalse(Path(retired_filename).exists())

    def test_single_node_registry_uses_current_functions_and_exclusions(self):
        # Other hardware-free test modules may narrow the wildcard-exported
        # names on the shared ARTIQ stub during unittest discovery. Only ever
        # touch our own marked stub, never a real artiq installation.
        if getattr(sys.modules.get("artiq"), _STUB_MARKER, False):
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
        self.assertIn("atom_photon_parity_11_AllSPCM_experiment", registry)
        self.assertEqual(len(registry), 51)
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
            _GeneralVariableScanMasterSatelliteMixin._select_experiment_function(
                build_two_node_function_registry(),
                "atom_loading_experiment",
            )

    def test_selected_node_has_separate_legacy_presentation(self):
        scan = GeneralVariableScan_master_satellite_single_node()
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
        scan = GeneralVariableScan_master_satellite_single_node()
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

    def _make_override_scan(self, node):
        """A single-node scan wired up just enough to collect overrides."""
        scan = GeneralVariableScan_master_satellite_single_node()
        scan.selected_node = node
        scan.base = FakeBase("single_node", node)
        scan.base.experiment = scan
        scan.override_ExperimentVariables_Node1 = "{}"
        scan.override_ExperimentVariables_Node2 = "{}"
        return scan

    def test_idle_node_overrides_do_not_block_the_running_node(self):
        """Both per-node dictionaries stay populated; only one takes effect.

        This is the reported failure: Node2 entries left in override_
        ExperimentVariables made a Node1 run fail, so the dictionary had to be
        emptied by hand to switch nodes.
        """
        populated = {
            "override_ExperimentVariables_Node1": "{'f_FORT_Node1': 241.0}",
            "override_ExperimentVariables_Node2": "{'f_FORT_Node2': 242.0}",
        }
        for node, expected in (
            ("Node1", {"f_FORT_Node1": 241.0}),
            ("Node2", {"f_FORT_Node2": 242.0}),
        ):
            scan = self._make_override_scan(node)
            for name, value in populated.items():
                setattr(scan, name, value)
            self.assertEqual(
                scan._collect_authoritative_overrides(), expected,
                f"{node} run did not apply exactly its own overrides",
            )

    def test_idle_node_overrides_are_never_evaluated(self):
        """Skipped before eval, not merely excused after a failed resolve.

        In single_node mode the idle node's attributes are never loaded, so an
        override that reads one cannot even be evaluated. The positive control
        below runs the same text on the node that owns it and does fail, which
        is what makes the negative case evidence of a skipped branch rather
        than of a harmless expression.
        """
        reads_node2 = "{'f_FORT_Node2': self.f_FORT_Node2}"

        scan = self._make_override_scan("Node1")
        scan.override_ExperimentVariables_Node2 = reads_node2
        self.assertEqual(scan._collect_authoritative_overrides(), {})
        # Nothing from the idle dictionary reached the resolver either.
        self.assertEqual(scan.base.resolve_calls, [])

        scan = self._make_override_scan("Node2")
        scan.override_ExperimentVariables_Node2 = reads_node2
        with self.assertRaisesRegex(ValueError, "Could not evaluate"):
            scan._collect_authoritative_overrides()

    def test_globals_and_bare_names_resolve_inside_the_node_dictionary(self):
        """There is no shared field, so globals go in the node's dictionary.

        A global resolves to itself and a bare name to the running node's
        suffixed attribute, both from the same dictionary.
        """
        scan = self._make_override_scan("Node2")
        scan.override_ExperimentVariables_Node2 = (
            "{'n_measurements': 5, 'f_FORT': 243.0}"
        )
        self.assertEqual(
            scan._collect_authoritative_overrides(),
            {"n_measurements": 5, "f_FORT_Node2": 243.0},
        )

    def test_wrong_node_name_inside_a_node_dictionary_still_raises(self):
        """The exemption is per-field, not blanket.

        The Node1 field IS read on a Node1 run, so a Node2 name written there
        is a genuine mistake and must not be silently dropped.
        """
        scan = self._make_override_scan("Node1")
        scan.override_ExperimentVariables_Node1 = "{'f_FORT_Node2': 242.0}"
        with self.assertRaisesRegex(ValueError, "other node"):
            scan._collect_authoritative_overrides()

    def test_one_variable_overridden_twice_raises(self):
        """Key order must not quietly decide which value wins.

        Two spellings of one variable resolve to the same target: on a Node2
        run 'f_FORT' and 'f_FORT_Node2' are the same attribute.
        """
        scan = self._make_override_scan("Node2")
        scan.override_ExperimentVariables_Node2 = (
            "{'f_FORT': 243.0, 'f_FORT_Node2': 244.0}"
        )
        with self.assertRaisesRegex(ValueError, "overridden twice"):
            scan._collect_authoritative_overrides()

    def test_blank_override_field_means_no_overrides(self):
        scan = self._make_override_scan("Node1")
        scan.override_ExperimentVariables_Node1 = "   "
        self.assertEqual(scan._collect_authoritative_overrides(), {})

    def test_two_nodes_mode_applies_both_per_node_dictionaries(self):
        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = FakeBase("two_nodes", None)
        scan.base.experiment = scan
        scan.override_ExperimentVariables_Node1 = "{'f_FORT_Node1': 241.0}"
        scan.override_ExperimentVariables_Node2 = "{'f_FORT_Node2': 242.0}"
        self.assertEqual(
            scan._collect_authoritative_overrides(),
            {"f_FORT_Node1": 241.0, "f_FORT_Node2": 242.0},
        )

    def test_two_nodes_mode_suffixes_bare_names_per_dictionary(self):
        """A bare name in two_nodes mode takes the suffix of its own field.

        The field says which node the entry belongs to, so nothing is
        ambiguous -- unlike a scan variable, which carries no node and is
        still rejected (see the resolver's two_nodes branch). This is what
        lets the same dictionary text be moved between the two fields, and
        between single_node and two_nodes, without editing every key.
        """
        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = FakeBase("two_nodes", None)
        scan.base.experiment = scan
        # Bare on one side, already spelled out on the other, plus a global.
        scan.override_ExperimentVariables_Node1 = (
            "{'f_FORT': 241.0, 'n_measurements': 7}"
        )
        scan.override_ExperimentVariables_Node2 = "{'f_FORT_Node2': 242.0}"
        self.assertEqual(
            scan._collect_authoritative_overrides(),
            {"f_FORT_Node1": 241.0, "f_FORT_Node2": 242.0, "n_measurements": 7},
        )

    def test_two_nodes_mode_still_catches_one_target_named_twice(self):
        """Bare and suffixed spellings must collide, not silently merge."""
        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = FakeBase("two_nodes", None)
        scan.base.experiment = scan
        scan.override_ExperimentVariables_Node1 = "{'f_FORT': 241.0}"
        scan.override_ExperimentVariables_Node2 = "{'f_FORT_Node1': 242.0}"
        with self.assertRaisesRegex(ValueError, "overridden twice"):
            scan._collect_authoritative_overrides()

    def test_both_public_gvs_classes_declare_the_per_node_fields(self):
        for experiment_class in (
            GeneralVariableScan_master_satellite_single_node,
            GeneralVariableScan_master_satellite_two_nodes,
        ):
            scan = self._make_repository_examination_experiment(
                experiment_class
            )
            scan.build()
            for node in ("Node1", "Node2"):
                self.assertTrue(
                    hasattr(scan, f"override_ExperimentVariables_{node}"),
                    f"{experiment_class.__name__} does not declare "
                    f"override_ExperimentVariables_{node}",
                )
            # The node-agnostic field was removed deliberately: with per-node
            # dictionaries it is never read, and leaving it on the dashboard
            # invites overrides that silently do nothing.
            self.assertFalse(
                hasattr(scan, "override_ExperimentVariables"),
                f"{experiment_class.__name__} still declares the removed "
                "node-agnostic override_ExperimentVariables field",
            )

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
        scan = GeneralVariableScan_master_satellite_single_node()
        scan.f_FORT_Node2 = 240.0
        scan.f_FORT = 240.0
        scan.n_measurements = 12
        scan.execution_n_measurements = 12
        scan.enable_Catch_UnderFlow = False
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.needs_experiment_variable_reload = True
        scan.authoritative_overrides = {"f_FORT_Node2": 245.0}
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None
        scan.scan_sequence1 = [246.0, 247.0]
        scan.scan_sequence2 = [0.0]
        # Submitted argument names and the dataset names Base publishes; the
        # run now also broadcasts the scan labels the retention applet reads.
        scan.scan_variable1_name = "f_FORT"
        scan.scan_variable2_name = ""
        scan.scan_var_dataset = "scan_variables"
        scan.scan_sequence1_dataset = "scan_sequence1"
        scan.scan_sequence2_dataset = "scan_sequence2"
        # Applet creation is a dashboard convenience; off for this unit test.
        scan.create_applets = False
        scan.hardware_initializations = 0
        scan.function_calls = 0
        scan.dataset_writes = []
        # prepare() normally sets these two; this harness builds the scan by
        # hand, so supply them for the per-scan-point write_results name.
        scan.experiment_name = "atom_loading_experiment"
        scan.scan_var_filesuffix = "f_FORT"
        scan.result_names = []
        scan.write_results = lambda kwargs={}: scan.result_names.append(
            kwargs.get("name")
        )
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
        # The other node is driven safe ONCE for the whole run, not once per
        # scan point: zeroing its Zotino is an init plus sixteen DAC writes
        # over remote SPI, and two scan points ran above.
        self.assertEqual(scan.base.force_off_calls, 1)
        # One named save per scan point, each carrying the node FIRST in the
        # name. ARTIQ's own <rid>-<class>.h5 has no node in it, so this is
        # the only thing that makes results sortable by node in Analysis.
        self.assertEqual(
            scan.result_names,
            ["Node2_atom_loading_scan_over_f_FORT"] * 2,
        )
        self.assertEqual(scan.function_calls, 2)
        self.assertEqual(scan.base.result_initializations, 1)
        self.assertEqual(scan.base.result_resets, 2)
        self.assertEqual(scan.f_FORT_Node2, 247.0)
        self.assertEqual(scan.f_FORT, 247.0)
        self.assertEqual(scan.base.dds_frequency_cache, 247.0)
        self.assertEqual(scan.n_measurements, 12)
        # run() may broadcast progress and the scan labels the applets read,
        # and nothing else -- in particular never an authoritative variable.
        written = [name for name, _, _ in scan.dataset_writes]
        self.assertTrue(
            set(written) <= {"iteration", "scan_variables", "scan_sequence1",
                             "scan_sequence2"},
            f"unexpected dataset writes during run: {sorted(set(written))}",
        )
        self.assertNotIn("f_FORT_Node2", written)
        self.assertNotIn("f_FORT", written)
        self.assertTrue(
            all(not kwargs.get("persist", False) for _, _, kwargs in scan.dataset_writes)
        )

    def test_submitted_n_measurements_survives_prepare_dataset_load(self):
        scan = self._make_repository_examination_experiment(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.build()
        for variable in (
            *NODE1_VARIABLES,
            *NODE2_VARIABLES,
            *MASTER_SATELLITE_VARIABLES,
        ):
            scan.datasets[variable.name] = variable.value
        self.assertEqual(scan.datasets["n_measurements"], 100)

        # Submitted GUI values; n_measurements differs from the dataset.
        scan.selected_node = "Node1"
        scan.n_measurements = 7
        scan.scan_variable1_name = "f_FORT"
        scan.scan_sequence1 = "np.array([self.f_FORT])"
        scan.scan_variable2_name = ""
        scan.scan_sequence2 = "np.zeros(1)"
        scan.override_ExperimentVariables_Node1 = "{}"
        scan.override_ExperimentVariables_Node2 = "{}"
        scan.experiment_function = "atom_loading_experiment"
        scan.scheduler = types.SimpleNamespace(get_status=lambda: {}, rid=0)
        # Single-node prepare() now builds the laser stabilizer so the reused
        # experiment functions find self.stabilizer_<channel>. Inject a
        # stand-in: the real AOMPowerStabilizer reads a config file resolved
        # from the process cwd, which is meaningless under the test runner.
        scan._stabilizer_factory = _FakeStabilizer
        # This harness replaces get_dataset with an instance attribute, which
        # shadows the redirect mixin, so the stabilizer's channel datasets are
        # read under their unsuffixed legacy names here.
        scan.datasets["p_AOM_A1"] = 0.0
        scan.datasets["p_FORT_loading"] = 0.0

        scan.prepare()

        self.assertEqual(scan.n_measurements, 7)
        self.assertEqual(scan.execution_n_measurements, 7)
        self.assertEqual(scan.scan_variable1, "f_FORT_Node1")
        self.assertEqual(scan.f_FORT, scan.datasets["f_FORT_Node1"])

    def _make_catch_error_scan(self, scan_class):
        scan = scan_class()
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.f_FORT_Node2 = 240.0
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None
        scan.scan_sequence1 = [1.0, 2.0]
        scan.scan_sequence2 = [0.0]
        scan.enable_Catch_UnderFlow = True
        scan.underflow_max_retries = 3
        scan.underflow_backoff_ms = 0.0
        scan.skip_only_that_iteration_if_exhausted = False
        scan.dataset_writes = []
        scan.attempts = []
        scan._initialize_run_state = lambda: None
        scan.initialize_hardware = lambda: None
        # These tests are about underflow retry, not filenames, but the
        # per-scan-point save still runs; give it what it needs and swallow it.
        scan.experiment_name = "atom_loading_experiment"
        scan.scan_var_filesuffix = "f_FORT"
        scan.write_results = lambda kwargs={}: None
        scan.set_dataset = lambda name, value, **kwargs: scan.dataset_writes.append(
            (name, value, kwargs)
        )
        return scan

    def test_disabled_catch_underflow_propagates_immediately(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.enable_Catch_UnderFlow = False

        def selected_function(experiment):
            experiment.attempts.append(experiment.f_FORT_Node2)
            raise RTIOUnderflow("unmanaged")

        scan._selected_experiment_function = selected_function
        with self.assertRaises(RTIOUnderflow):
            scan.run()
        self.assertEqual(scan.attempts, [1.0])

    def test_catch_error_retries_only_failed_scan_point(self):
        for scan_class in (
            GeneralVariableScan_master_satellite_single_node,
            GeneralVariableScan_master_satellite_two_nodes,
        ):
            with self.subTest(experiment=scan_class.__name__):
                scan = self._make_catch_error_scan(scan_class)

                def selected_function(experiment):
                    experiment.attempts.append(experiment.f_FORT_Node2)
                    if len(experiment.attempts) == 1:
                        raise RTIOUnderflow("transient")

                scan._selected_experiment_function = selected_function
                scan.run()

                self.assertEqual(scan.attempts, [1.0, 1.0, 2.0])
                iteration_writes = [
                    value
                    for name, value, _ in scan.dataset_writes
                    if name == "iteration"
                ]
                self.assertEqual(iteration_writes, [0, 0, 0, 1])
                self.assertEqual(scan.base.result_resets, 3)
                self.assertTrue(
                    all(
                        not kwargs.get("persist", False)
                        for _, _, kwargs in scan.dataset_writes
                    )
                )

    def test_catch_error_does_not_hide_non_underflow_errors(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.skip_only_that_iteration_if_exhausted = True

        def selected_function(experiment):
            raise RuntimeError("DRTIO unavailable")

        scan._selected_experiment_function = selected_function
        with self.assertRaisesRegex(RuntimeError, "DRTIO unavailable"):
            scan.run()

    def test_catch_error_exhaustion_respects_skip_configuration(self):
        def always_underflow_first_point(experiment):
            experiment.attempts.append(experiment.f_FORT_Node2)
            if experiment.f_FORT_Node2 == 1.0:
                raise RTIOUnderflow("persistent")

        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.underflow_max_retries = 2
        scan.skip_only_that_iteration_if_exhausted = True
        scan._selected_experiment_function = always_underflow_first_point
        scan.run()
        self.assertEqual(scan.attempts, [1.0, 1.0, 2.0])

        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_single_node
        )
        scan.underflow_max_retries = 2
        scan.skip_only_that_iteration_if_exhausted = False
        scan._selected_experiment_function = always_underflow_first_point
        with self.assertRaises(RTIOUnderflow):
            scan.run()
        self.assertEqual(scan.attempts, [1.0, 1.0])

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
