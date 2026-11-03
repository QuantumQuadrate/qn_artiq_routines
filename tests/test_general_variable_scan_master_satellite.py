import sys
import types
import unittest
import json
import ast
from pathlib import Path

import numpy as np


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
        # the GVS health check catches it around scheduler.pause()
        "TerminationRequested": type(
            "TerminationRequested", (Exception,), {}),
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
    parse_underflow_message,
    repository_frames,
    roll_back_list_datasets,
    rtio_channel_names,
    snapshot_list_datasets,
)
from GeneralVariableScan_master_satellite_single_node import (  # noqa: E402
    GeneralVariableScan_master_satellite_single_node,
)
from GeneralVariableScan_master_satellite_two_nodes import (  # noqa: E402
    GeneralVariableScan_master_satellite_two_nodes,
)
from subroutines.experiment_functions_two_nodes import (  # noqa: E402
    MASTER_SATELLITE_SANITY_ATTRIBUTES,
    Two_nodes_atom_loading_experiment,
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


class _FakeBroadcaster:
    """DatasetManager's broadcaster as ARTIQ 7 uses it: raw_view holds
    (persist, value) per broadcast key, and every top-level write is
    published to the master -- recorded here, with a copy of the value as
    it was sent."""

    def __init__(self):
        self.raw_view = {}
        self.published = []

    def __setitem__(self, key, value):
        self.raw_view[key] = value
        persist, data = value
        self.published.append(
            ("setitem", key, persist,
             list(data) if isinstance(data, list) else data)
        )

    def __delitem__(self, key):
        del self.raw_view[key]
        self.published.append(("delitem", key))


class FakeDatasetManager:
    """artiq.master.worker_db.DatasetManager (ARTIQ 7), the parts a scan
    touches, with the same semantics: set() and append_to() are copied."""

    def __init__(self):
        self._broadcaster = _FakeBroadcaster()
        self.local = {}

    def set(self, key, value, broadcast=False, persist=False, archive=True):
        if persist:
            broadcast = True
        if broadcast:
            self._broadcaster[key] = persist, value
        elif key in self._broadcaster.raw_view:
            del self._broadcaster[key]
        if archive:
            self.local[key] = value
        elif key in self.local:
            del self.local[key]

    def append_to(self, key, value):
        target = self.local.get(key, None)
        if key in self._broadcaster.raw_view:
            if target is not None:
                assert target is self._broadcaster.raw_view[key][1]
            target = self._broadcaster.raw_view[key][1]
        if target is None:
            raise KeyError(f"Cannot mutate nonexistent dataset '{key}'")
        target.append(value)


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
        self.global_variable_names = frozenset({
            "n_measurements", "t_Node2_excitation_delay_mu",
            "t_Node2_rtio_offset_mu", "parallel_AOM_feedback",
            "t_pumping_two_nodes"})

    def resolve_experiment_variable_target(self, name, node=None):
        """Mirror the real resolver, including its node-aware override path.

        node is the node that OWNS the name -- the per-node override
        dictionary it was written in -- and makes a bare node-specific name
        unambiguous. A scan variable passes no node, so a bare name stays
        ambiguous in two_nodes mode there.
        """
        self.resolve_calls.append(name)
        if name in self.global_variable_names:
            return name
        # a two-node sequence setting: the bare name is the shared global
        if self.mode == "two_nodes" and name == "t_pumping":
            return "t_pumping_two_nodes"
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

    def initialize_result_state(self):
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
            {
                "master_satellite_namespace_sanity_experiment",
                # Keeps the legacy name: the two-node migration preserves the
                # original sequence shape rather than renaming it.
                "Two_nodes_atom_loading_experiment",
                # The alternating readout is its own entry point, not a flag:
                # two_node_alternating_shot is reachable from nothing else, and
                # ARTIQ only type-checks what an entry point reaches, so
                # without this it would never be compiled by anything.
                # Both shots node-resolved, for per-node retention.
                "Two_nodes_alternating_retention_experiment",
                # Port of the standalone Two_node_single_photon_2_experiment.
                "Two_nodes_HOM_experiment",
                # Two-photon herald, then each node's parity readout.
                "Two_nodes_atom_atom_parity_experiment",
            },
        )
        # The two registries must stay disjoint: the single-node one is built
        # from experiment_functions.py and the two-node one from
        # experiment_functions_two_nodes.py, so an overlap means something
        # landed in the wrong file.
        self.assertEqual(set(single_registry) & set(two_registry), set())

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
        # 23 since 2026-10-20: the functions no longer in use moved to
        # subroutines/experiment_functions_obsolete.py, the four historical
        # two-node experiments among them (51 before; 54 until 2026-10-13).
        self.assertEqual(len(registry), 23)
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

    def test_native_registry_contains_exactly_the_native_experiments(self):
        """An exact inventory, so adding a two-node entry point is deliberate.

        The registry is built by name matching -- every locally defined function
        whose name CONTAINS "experiment" -- so a helper that picked up that
        substring would silently appear in the dashboard dropdown as something
        runnable. This is what catches that.
        """
        registry = build_two_node_function_registry()
        self.assertEqual(
            set(registry),
            {
                "master_satellite_namespace_sanity_experiment",
                # Keeps the legacy name: the two-node migration preserves the
                # original sequence shape rather than renaming it.
                "Two_nodes_atom_loading_experiment",
                # The alternating readout is its own entry point, not a flag:
                # two_node_alternating_shot is reachable from nothing else, and
                # ARTIQ only type-checks what an entry point reaches, so
                # without this it would never be compiled by anything.
                # Both shots node-resolved, for per-node retention.
                "Two_nodes_alternating_retention_experiment",
                # Port of the standalone Two_node_single_photon_2_experiment.
                "Two_nodes_HOM_experiment",
                # Two-photon herald, then each node's parity readout.
                "Two_nodes_atom_atom_parity_experiment",
            },
        )
        self.assertIs(
            registry["master_satellite_namespace_sanity_experiment"],
            master_satellite_namespace_sanity_experiment,
        )
        # HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS keeps the old
        # two-independent-Kasli functions out of the SINGLE-node registry (see
        # test_single_node_registry_uses_current_functions_and_exclusions);
        # it is not a blanket ban on the names. The master-satellite rewrite of
        # the atom loading sequence deliberately keeps the legacy name, so the
        # same string is both excluded from single-node mode -- where the old
        # body in experiment_functions.py would otherwise be reachable and
        # would drive one crate only -- and present here, where it resolves to
        # the new body in experiment_functions_two_nodes.py.
        self.assertIs(
            registry["Two_nodes_atom_loading_experiment"],
            Two_nodes_atom_loading_experiment,
        )
        for name in HISTORICAL_INDEPENDENT_TWO_NODE_EXPERIMENTS:
            if name == "Two_nodes_atom_loading_experiment":
                continue
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

    def test_overrides_are_recorded_locally_like_standalone(self):
        """Every applied override is written as a LOCAL dataset for the h5.

        Standalone GeneralVariableScan has always done this
        (initialize_datasets: self.set_dataset(var, val)), which is how its
        analysis sees what a run actually used: the h5 keeps the persistent
        value under `archive` and the override under `datasets`, and the
        loader reads datasets last. Master-satellite GVS only used setattr, so
        RID 38903 -- t_MW_RF_pulse overridden to 0, Z parity -- was analysed
        as X parity from the archived 69.6 us.

        Pinned: the name written is the one the run read (suffixed per-node,
        bare global), never the bare per-node name; the write is local only,
        never broadcast or persisted, so the master's dataset DB is untouched;
        and an override of the scanned variable is not recorded, since the
        scan sets that value per point.
        """
        scan = GeneralVariableScan_master_satellite_single_node()
        scan.f_FORT_Node2 = 240.0
        scan.f_FORT = 240.0
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.authoritative_overrides = {
            "f_FORT_Node2": 245.0,
            "t_MW_RF_pulse_Node2": 0.0,
            "advance_if_atom_in_first_shot": True,
            "target_780_HWP_Node2": 10.0,
        }
        scan.scan_variable1 = "target_780_HWP_Node2"
        scan.scan_variable2 = None
        scan.dataset_writes = []
        scan.set_dataset = lambda name, value, **kwargs: scan.dataset_writes.append(
            (name, value, kwargs)
        )

        scan._apply_run_wide_overrides()

        written = {name: (value, kwargs)
                   for name, value, kwargs in scan.dataset_writes}
        self.assertEqual(
            set(written),
            {"f_FORT_Node2", "t_MW_RF_pulse_Node2",
             "advance_if_atom_in_first_shot"},
            "the scanned variable must not be recorded, and nothing bare "
            "may be written for a per-node variable",
        )
        self.assertEqual(written["f_FORT_Node2"][0], 245.0)
        self.assertEqual(written["t_MW_RF_pulse_Node2"][0], 0.0)
        self.assertIs(written["advance_if_atom_in_first_shot"][0], True)
        for name, (_, kwargs) in written.items():
            self.assertIs(kwargs.get("broadcast"), False, name)
            self.assertIs(kwargs.get("persist"), False, name)
        # The attribute the sequence reads still follows the override.
        self.assertEqual(scan.f_FORT_Node2, 245.0)
        self.assertEqual(scan.f_FORT, 245.0)

    def test_override_records_precede_every_broadcast(self):
        """The local override records are written before ANY broadcast.

        ARTIQ's DatasetManager.set deletes a key from the master DB when a
        local write follows a broadcast of the same key within one run
        (artiq/master/worker_db.py, the `elif key in
        self._broadcaster.raw_view` branch). initialize_result_state
        broadcasts and persists n_measurements, a valid override target, so
        if the override pass ever moved after it, overriding n_measurements
        would delete that dataset and every master-satellite experiment would
        then fail in prepare() with "Missing required master-satellite
        persistent datasets". The pass must stay first.
        """
        scan = GeneralVariableScan_master_satellite_single_node()
        scan.f_FORT_Node2 = 240.0
        scan.f_FORT = 240.0
        scan.n_measurements = 12
        scan.execution_n_measurements = 12
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        scan.needs_experiment_variable_reload = True
        scan.authoritative_overrides = {
            "n_measurements": 7,
            "t_MW_RF_pulse_Node2": 0.0,
        }
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None
        scan.scan_sequence1 = [246.0]
        scan.scan_sequence2 = [0.0]
        scan.scan_variable1_name = "f_FORT"
        scan.scan_variable2_name = ""
        scan.scan_var_dataset = "scan_variables"
        scan.scan_sequence1_dataset = "scan_sequence1"
        scan.scan_sequence2_dataset = "scan_sequence2"
        scan.create_applets = False

        events = []
        scan.set_dataset = lambda name, value, **kwargs: events.append(
            ("broadcast" if kwargs.get("broadcast") or kwargs.get("persist")
             else "local", name)
        )
        for method in ("initialize_result_datasets", "initialize_result_state"):
            original = getattr(scan.base, method)

            def marked(original=original, method=method):
                events.append(("broadcasts", method))
                return original()

            setattr(scan.base, method, marked)

        scan._initialize_run_state()

        local = [i for i, (kind, _) in enumerate(events) if kind == "local"]
        others = [i for i, (kind, _) in enumerate(events) if kind != "local"]
        self.assertEqual(
            {events[i][1] for i in local},
            {"n_measurements", "t_MW_RF_pulse_Node2"},
        )
        self.assertTrue(others, "the harness must see the later broadcasts")
        self.assertLess(
            max(local), min(others),
            f"an override record came after a broadcast: {events}",
        )

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
        """A global in a node's dictionary still works, with a warning.

        Submissions saved before override_ExperimentVariables_global existed
        put globals there. A global resolves to itself and a bare name to the
        running node's suffixed attribute, both from the same dictionary.
        """
        scan = self._make_override_scan("Node2")
        scan.override_ExperimentVariables_Node2 = (
            "{'n_measurements': 5, 'f_FORT': 243.0}"
        )
        with self.assertLogs(level="WARNING") as logs:
            self.assertEqual(
                scan._collect_authoritative_overrides(),
                {"n_measurements": 5, "f_FORT_Node2": 243.0},
            )
        self.assertIn("override_ExperimentVariables_global", logs.output[0])
        self.assertIn("'n_measurements' in override_ExperimentVariables_Node2",
                      logs.output[0])
        self.assertNotIn("f_FORT", logs.output[0])

    def test_global_field_takes_the_shared_variables_in_both_modes(self):
        """A global in either mode; in two_nodes also a two-node sequence
        setting written bare. No warning: this is where they belong."""
        import logging
        scan = self._make_override_scan("Node1")
        scan.override_ExperimentVariables_global = "{'n_measurements': 5}"
        scan.override_ExperimentVariables_Node1 = "{'f_FORT': 241.0}"
        ### Python 3.9 has no assertNoLogs: the sentinel must be the only line
        with self.assertLogs(level="WARNING") as logs:
            self.assertEqual(
                scan._collect_authoritative_overrides(),
                {"n_measurements": 5, "f_FORT_Node1": 241.0},
            )
            logging.getLogger().warning("sentinel")
        self.assertEqual(logs.output, ["WARNING:root:sentinel"])

        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = FakeBase("two_nodes", None)
        scan.base.experiment = scan
        scan.override_ExperimentVariables_global = (
            "{'t_pumping': 2e-4, 'n_measurements': 7}"
        )
        scan.override_ExperimentVariables_Node1 = "{'f_FORT': 241.0}"
        scan.override_ExperimentVariables_Node2 = "{'f_FORT_Node2': 242.0}"
        with self.assertLogs(level="WARNING") as logs:
            self.assertEqual(
                scan._collect_authoritative_overrides(),
                {"t_pumping_two_nodes": 2e-4, "n_measurements": 7,
                 "f_FORT_Node1": 241.0, "f_FORT_Node2": 242.0},
            )
            logging.getLogger().warning("sentinel")
        self.assertEqual(logs.output, ["WARNING:root:sentinel"])

    def test_two_node_runs_check_the_MOT_coil_set_their_loader_reads(self):
        """Base refuses an override or scan of the unread MOT coil set.

        Which set is read is the flag AS THE RUN SETS IT, so an override of
        override_coils_with_two_nodes_setting counts; a scan over the flag
        uses both sets, so there is nothing to check.
        """
        calls = []
        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = types.SimpleNamespace(
            check_two_node_MOT_coil_targets=(
                lambda targets, use: calls.append((set(targets), use))))
        scan.override_coils_with_two_nodes_setting = False  # stored value
        scan.scan_variable1 = "AX_volts_MOT_for_two_nodes_Node1"
        scan.scan_variable2 = None

        scan.authoritative_overrides = {"f_FORT_Node1": 241.0}
        scan._check_two_node_MOT_coil_targets()
        scan.authoritative_overrides = {
            "override_coils_with_two_nodes_setting": True}
        scan._check_two_node_MOT_coil_targets()
        scan.scan_variable1 = "override_coils_with_two_nodes_setting"
        scan._check_two_node_MOT_coil_targets()

        self.assertEqual(calls, [
            ({"f_FORT_Node1", "AX_volts_MOT_for_two_nodes_Node1"}, False),
            ({"override_coils_with_two_nodes_setting",
              "AX_volts_MOT_for_two_nodes_Node1"}, True),
        ])

    def test_global_field_refuses_per_node_variables(self):
        """A per-node name there would silently pick a node (single_node)
        or be ambiguous (two_nodes): it raises and names the node fields."""
        scan = self._make_override_scan("Node1")
        scan.override_ExperimentVariables_global = "{'f_FORT': 241.0}"
        with self.assertRaisesRegex(
                ValueError, "per-node variable 'f_FORT_Node1'.*"
                            "override_ExperimentVariables_Node1"):
            scan._collect_authoritative_overrides()

        for text in ("{'f_FORT': 241.0}", "{'f_FORT_Node2': 242.0}"):
            scan = GeneralVariableScan_master_satellite_two_nodes()
            scan.base = FakeBase("two_nodes", None)
            scan.base.experiment = scan
            scan.override_ExperimentVariables_global = text
            with self.assertRaisesRegex(
                    ValueError, "override_ExperimentVariables_global.*"
                                "Per-node variables go in "
                                "override_ExperimentVariables_Node1 or _Node2"):
                scan._collect_authoritative_overrides()

    def test_global_named_in_global_and_node_field_raises(self):
        """One variable in both fields is overridden twice, as between the
        two node fields."""
        scan = GeneralVariableScan_master_satellite_two_nodes()
        scan.base = FakeBase("two_nodes", None)
        scan.base.experiment = scan
        scan.override_ExperimentVariables_global = "{'t_pumping': 2e-4}"
        scan.override_ExperimentVariables_Node1 = "{'t_pumping': 3e-4}"
        scan.override_ExperimentVariables_Node2 = "{}"
        with self.assertRaisesRegex(
                ValueError, "'t_pumping_two_nodes' is overridden twice: once "
                            "via 't_pumping' in "
                            "override_ExperimentVariables_global and again via "
                            "'t_pumping' in override_ExperimentVariables_Node1"):
            scan._collect_authoritative_overrides()

    def test_missing_or_blank_global_field_means_no_overrides(self):
        """A submission saved before the field existed has no attribute."""
        scan = self._make_override_scan("Node1")
        self.assertFalse(hasattr(scan, "override_ExperimentVariables_global"))
        self.assertEqual(scan._collect_authoritative_overrides(), {})
        scan.override_ExperimentVariables_global = "  "
        self.assertEqual(scan._collect_authoritative_overrides(), {})

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
            for field in ("Node1", "Node2", "global"):
                self.assertTrue(
                    hasattr(scan, f"override_ExperimentVariables_{field}"),
                    f"{experiment_class.__name__} does not declare "
                    f"override_ExperimentVariables_{field}",
                )
            # The node-agnostic field was removed deliberately: with per-node
            # dictionaries it is never read, and leaving it on the dashboard
            # invites overrides that silently do nothing.
            self.assertFalse(
                hasattr(scan, "override_ExperimentVariables"),
                f"{experiment_class.__name__} still declares the removed "
                "node-agnostic override_ExperimentVariables field",
            )

    def test_first_shot_alive_applet_only_with_the_confirm_on(self):
        """Two-node runs ask for it only with advance_if_atom_in_first_shot on
        AND an experiment whose confirm loop writes the time."""
        import applets_master_satellite as applets

        title = applets.FIRST_SHOT_ALIVE_APPLET_SPECS[0].title
        cases = (
            ("Two_nodes_HOM_experiment", True, True),
            ("Two_nodes_atom_loading_experiment", True, True),
            ("Two_nodes_HOM_experiment", False, False),
            ("Two_nodes_alternating_retention_experiment", True, False),
            ("Two_nodes_atom_atom_parity_experiment", True, True),
            ("Two_nodes_atom_atom_parity_experiment", False, False),
        )
        original = applets.create_two_node_applets_for
        try:
            for experiment_name, flag, expected in cases:
                calls = []
                applets.create_two_node_applets_for = (
                    lambda experiment, base, **kwargs:
                    calls.append(kwargs) or []
                )
                scan = GeneralVariableScan_master_satellite_two_nodes()
                scan.create_applets = True
                scan.base = types.SimpleNamespace(which_node=None)
                scan.experiment_name = experiment_name
                scan.advance_if_atom_in_first_shot = flag
                scan._create_node_applets()
                self.assertEqual(len(calls), 1, experiment_name)
                titles = {spec.title for spec in calls[0]["specs"]}
                self.assertEqual(title in titles, expected,
                                 (experiment_name, flag))
        finally:
            applets.create_two_node_applets_for = original

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
        scan._HasEnvironment__dataset_mgr = FakeDatasetManager()

        def set_dataset(name, value, **kwargs):
            scan.dataset_writes.append((name, value, kwargs))
            scan._HasEnvironment__dataset_mgr.set(name, value, **kwargs)

        scan.set_dataset = set_dataset
        return scan

    def _seed_result_series(self, scan):
        """Run-wide series as initialize_result_state leaves them: broadcast
        and archived, archived only, and broadcast-only and persistent."""
        datasets = scan._HasEnvironment__dataset_mgr
        datasets.set("AllSPCMs_RO1", [0], broadcast=True)
        datasets.set("Atom_loading_time", [0.0])
        datasets.set("p_AOM_A1_history_Node2", [-3.0], broadcast=True,
                     persist=True, archive=False)
        return datasets

    @staticmethod
    def _measure(experiment, n, fail_after=None):
        """A stand-in experiment function: n measurements, each appending to
        the series; a per-point display list replaced first, as Base's reset
        does; raises RTIOUnderflow after fail_after measurements."""
        datasets = experiment._HasEnvironment__dataset_mgr
        experiment.attempts.append(experiment.f_FORT_Node2)
        datasets.set("AllSPCMs_RO1_current_iteration", [0], broadcast=True)
        for measurement in range(n):
            if measurement == fail_after:
                raise RTIOUnderflow("transient")
            value = experiment.f_FORT_Node2
            datasets.append_to("AllSPCMs_RO1", value)
            datasets.append_to("AllSPCMs_RO1_current_iteration", value)
            datasets.append_to("Atom_loading_time", value)
            datasets.append_to("p_AOM_A1_history_Node2", value)

    def test_catch_error_discards_the_failed_attempts_data(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes
        )
        datasets = self._seed_result_series(scan)
        seeded = {key: datasets.local.get(key, datasets._broadcaster.raw_view.get(key, (None, None))[1])
                  for key in ("AllSPCMs_RO1", "Atom_loading_time",
                              "p_AOM_A1_history_Node2")}

        def selected_function(experiment):
            first = len(experiment.attempts) == 0
            self._measure(experiment, 3, fail_after=2 if first else None)

        scan._selected_experiment_function = selected_function
        scan.run()

        self.assertEqual(scan.attempts, [1.0, 1.0, 2.0])
        whole_points = [1.0] * 3 + [2.0] * 3
        self.assertEqual(datasets.local["AllSPCMs_RO1"], [0] + whole_points)
        self.assertEqual(datasets.local["Atom_loading_time"], [0.0] + whole_points)
        self.assertEqual(datasets._broadcaster.raw_view["p_AOM_A1_history_Node2"],
                         (True, [-3.0] + whole_points))
        self.assertEqual(datasets.local["AllSPCMs_RO1_current_iteration"],
                         [0, 2.0, 2.0, 2.0])
        # same objects, flags unchanged: broadcast stays broadcast (a local
        # write would delete it from the master DB), archive-only stays
        # local, and the persistent broadcast-only one is still both
        self.assertIs(datasets.local["AllSPCMs_RO1"], seeded["AllSPCMs_RO1"])
        self.assertEqual(datasets._broadcaster.raw_view["AllSPCMs_RO1"],
                         (False, datasets.local["AllSPCMs_RO1"]))
        self.assertNotIn("Atom_loading_time", datasets._broadcaster.raw_view)
        self.assertNotIn("p_AOM_A1_history_Node2", datasets.local)
        # the master (and so the dashboard) got the cut series: the seed,
        # then the rollback (appends are not top-level writes)
        def published(key):
            return [m[2:] for m in datasets._broadcaster.published
                    if m[0] == "setitem" and m[1] == key]
        self.assertEqual(published("AllSPCMs_RO1"), [(False, [0])] * 2)
        self.assertEqual(published("p_AOM_A1_history_Node2"),
                         [(True, [-3.0])] * 2)
        self.assertNotIn("delitem",
                         [m[0] for m in datasets._broadcaster.published])
        self.assertEqual(datasets.local["underflow_retries_per_scan_point"], [1, 0])
        self.assertEqual(datasets.local["underflow_skipped_scan_points"], [])

    def test_catch_error_exhausted_point_leaves_no_data(self):
        def always_underflow_first_point(experiment):
            self._measure(experiment, 3,
                          fail_after=1 if experiment.f_FORT_Node2 == 1.0 else None)

        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes
        )
        datasets = self._seed_result_series(scan)
        scan.underflow_max_retries = 2
        scan.skip_only_that_iteration_if_exhausted = True
        scan._selected_experiment_function = always_underflow_first_point
        saves = []
        scan.write_results = lambda kwargs={}: saves.append(
            list(datasets.local["AllSPCMs_RO1"]))
        scan.run()
        self.assertEqual(datasets.local["AllSPCMs_RO1"], [0, 2.0, 2.0, 2.0])
        self.assertEqual(datasets.local["underflow_retries_per_scan_point"], [2, 0])
        self.assertEqual(datasets.local["underflow_skipped_scan_points"], [0])
        # the skipped point is saved too, without its attempts' data
        self.assertEqual(saves, [[0], [0, 2.0, 2.0, 2.0]])

        # skip off: the run stops, and the h5 written then holds no part of
        # the point
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes
        )
        datasets = self._seed_result_series(scan)
        scan.underflow_max_retries = 2
        scan.skip_only_that_iteration_if_exhausted = False
        scan._selected_experiment_function = always_underflow_first_point
        with self.assertRaises(RTIOUnderflow):
            scan.run()
        self.assertEqual(datasets.local["AllSPCMs_RO1"], [0])
        self.assertEqual(datasets.local["Atom_loading_time"], [0.0])
        self.assertEqual(datasets.local["underflow_retries_per_scan_point"], [2])
        self.assertEqual(datasets.local["underflow_skipped_scan_points"], [])

    def test_paired_scan_steps_both_variables_together(self):
        seen = []

        def selected_function(experiment):
            seen.append((experiment.f_FORT_Node2, experiment.dummy_variable))

        for paired, expected in (
            (True, [(1.0, 10.0), (2.0, 20.0)]),
            (False, [(1.0, 10.0), (1.0, 20.0), (2.0, 10.0), (2.0, 20.0)]),
        ):
            with self.subTest(paired=paired):
                del seen[:]
                scan = self._make_catch_error_scan(
                    GeneralVariableScan_master_satellite_two_nodes
                )
                scan.scan_variable2 = "dummy_variable"
                scan.scan_sequence2 = [10.0, 20.0]
                scan.pair_scan_sequences = paired
                scan._selected_experiment_function = selected_function
                scan.run()
                self.assertEqual(seen, expected)
                iteration_writes = [value for name, value, _ in
                                    scan.dataset_writes if name == "iteration"]
                self.assertEqual(iteration_writes,
                                 [0] + list(range(len(expected))))

    def test_paired_scan_validation(self):
        scan = GeneralVariableScan_master_satellite_two_nodes()
        self.assertIs(scan.pair_scan_sequences, False)
        scan._validate_paired_scan()  # off: nothing to check

        scan.pair_scan_sequences = True
        scan.scan_sequence1 = [0.0, 10.0]
        scan.scan_variable2 = None
        scan.scan_sequence2 = [0.0]
        with self.assertRaisesRegex(ValueError, "scan_variable2_name"):
            scan._validate_paired_scan()
        scan.scan_variable2 = "target_780_HWP_Node2"
        with self.assertRaisesRegex(ValueError, "equally long"):
            scan._validate_paired_scan()
        scan.scan_sequence2 = [5.0, 15.0]
        scan._validate_paired_scan()

    def test_parse_underflow_message(self):
        self.assertEqual(
            parse_underflow_message(
                "RTIO underflow at 123456789 mu, channel 65540, slack -1250 mu"),
            (123456789, 65540, -1250))
        # a DMA playback reports no slack
        self.assertEqual(
            parse_underflow_message("RTIO underflow at 42 mu, channel 7"),
            (42, 7, None))
        self.assertIsNone(parse_underflow_message("transient"))

    CORE_TRACEBACK = (
        'Core Device Traceback:\n'
        'Traceback (most recent call first):\n'
        '  File "ksupport/lib.rs", line 417, column 21, in (Rust function)\n'
        '    <unknown>\n'
        '  File "<artiq>/coredevice/dma.py", line 112, column 9, in '
        'artiq.coredevice.dma.CoreDMA.playback_handle (RA=+0x5c)\n'
        '    dma_playback(now_mu(), ptr)\n'
        '  File "C:\\Networking Experiment\\artiq codes\\artiq-master\\repository'
        '\\qn_artiq_routines\\subroutines\\experiment_functions_two_nodes.py", '
        'line 2862, column 17, in '
        'subroutines.experiment_functions_two_nodes.Two_nodes_HOM_experiment '
        '(RA=+0x3a0)\n'
        '    self.core_dma.playback_handle(recooling_handle)\n'
        'RTIOUnderflow(0): RTIO underflow at 3250000 mu, channel 65540\n'
        '\nEnd of Core Device Traceback\n')

    def test_repository_frames_and_channel_names(self):
        self.assertEqual(
            repository_frames(self.CORE_TRACEBACK),
            ["experiment_functions_two_nodes.py:2862 Two_nodes_HOM_experiment"])
        self.assertEqual(repository_frames(""), [])
        device_db = {
            "ttl14_Node2": {"type": "local", "module": "artiq.coredevice.ttl",
                            "class": "TTLOut", "arguments": {"channel": 0x010004}},
            "ttl_repump_switch_Node2": "ttl14_Node2",
            "repump": "ttl_repump_switch_Node2",
            "core": {"type": "local", "module": "artiq.coredevice.core",
                     "class": "Core", "arguments": {"host": "x"}},
            "loop_a": "loop_b",
            "loop_b": "loop_a",
        }
        self.assertEqual(rtio_channel_names(device_db),
                         {0x010004: "ttl14_Node2/ttl_repump_switch_Node2/repump"})

    def test_catch_error_records_where_each_underflow_happened(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes
        )
        datasets = self._seed_result_series(scan)
        datasets.set("n_first_shot_confirm_tries", [0])
        datasets.set("atom_loading_wall_clock", [0])
        scan.get_device_db = lambda: {
            "ttl14_Node2": {"type": "local", "module": "artiq.coredevice.ttl",
                            "class": "TTLOut", "arguments": {"channel": 0x010004}},
            "ttl_repump_switch_Node2": "ttl14_Node2",
        }
        core_traceback = self.CORE_TRACEBACK

        class CoreException:
            def __str__(self):
                return core_traceback

        def selected_function(experiment):
            # three measurements; the first attempt fails in the third, after
            # its load and first-shot confirm but before its readouts
            first = len(experiment.attempts) == 0
            experiment.attempts.append(experiment.f_FORT_Node2)
            for measurement in range(3):
                datasets.append_to("atom_loading_wall_clock",
                                   1000000 * (measurement + 1))
                datasets.append_to("n_first_shot_confirm_tries", 1)
                if first and measurement == 2:
                    error = RTIOUnderflow("RTIO underflow at 3250000 mu, "
                                          "channel 65540, slack -300 mu")
                    error.artiq_core_exception = CoreException()
                    raise error
                datasets.append_to("AllSPCMs_RO1", experiment.f_FORT_Node2)

        scan._selected_experiment_function = selected_function
        scan.run()

        self.assertEqual(datasets.local["AllSPCMs_RO1"],
                         [0] + [1.0] * 3 + [2.0] * 3)
        self.assertEqual(len(datasets.local["n_first_shot_confirm_tries"]), 7)
        self.assertEqual(datasets.local["underflow_retries_per_scan_point"], [1, 0])
        expected = {
            "underflow_iteration": [0],
            "underflow_retry": [1],
            "underflow_channel": [0x010004],
            "underflow_timestamp_mu": [3250000],
            "underflow_slack_mu": [-300.0],
            "underflow_measurements_done": [2],
            "underflow_confirmed_measurements": [3],
            "underflow_mu_after_last_load": [250000.0],
        }
        for key, value in expected.items():
            self.assertEqual(datasets.local[key], value, key)
        self.assertEqual(
            datasets.local["underflow_sites"],
            "iteration 0 retry 1: CPU write, slack -300 mu, channel 0x10004 = "
            "ttl14_Node2/ttl_repump_switch_Node2; at "
            "experiment_functions_two_nodes.py:2862 Two_nodes_HOM_experiment; "
            "the attempt had finished 2 measurement(s), 3 past the first-shot "
            "confirm, 0.250 ms after its last load stamp")
        self.assertIn("slack -300 mu", datasets.local["underflow_messages"])
        self.assertIn("End of Core Device Traceback",
                      datasets.local["underflow_messages"])
        # archived, never broadcast or persisted
        for name, _, kwargs in scan.dataset_writes:
            if name.startswith("underflow_"):
                self.assertFalse(kwargs.get("broadcast", False), name)
                self.assertFalse(kwargs.get("persist", False), name)

    def test_catch_error_records_an_unparsed_underflow(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_single_node
        )
        datasets = self._seed_result_series(scan)

        def selected_function(experiment):
            first = len(experiment.attempts) == 0
            self._measure(experiment, 3, fail_after=1 if first else None)

        scan._selected_experiment_function = selected_function
        scan.run()
        self.assertEqual(datasets.local["underflow_channel"], [-1])
        self.assertEqual(datasets.local["underflow_timestamp_mu"], [-1])
        self.assertTrue(np.isnan(datasets.local["underflow_slack_mu"][0]))
        self.assertTrue(np.isnan(datasets.local["underflow_mu_after_last_load"][0]))
        self.assertEqual(datasets.local["underflow_measurements_done"], [1])
        self.assertEqual(datasets.local["underflow_confirmed_measurements"], [0])
        self.assertEqual(
            datasets.local["underflow_sites"],
            "iteration 0 retry 1: unparsed message; at unknown line; the "
            "attempt had finished 1 measurement(s), 0 past the first-shot "
            "confirm")
        self.assertEqual(datasets.local["underflow_messages"],
                         "iteration 0 retry 1: transient\n")

    def test_disabled_catch_underflow_keeps_data_and_records_nothing(self):
        scan = self._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes
        )
        datasets = self._seed_result_series(scan)
        scan.enable_Catch_UnderFlow = False
        scan._selected_experiment_function = (
            lambda experiment: self._measure(experiment, 3, fail_after=2))
        with self.assertRaises(RTIOUnderflow):
            scan.run()
        self.assertEqual(datasets.local["AllSPCMs_RO1"], [0, 1.0, 1.0])
        self.assertNotIn("underflow_retries_per_scan_point", datasets.local)
        self.assertNotIn("underflow_skipped_scan_points", datasets.local)

    def test_roll_back_cuts_only_lists_grown_in_place(self):
        datasets = FakeDatasetManager()
        datasets.set("grown", [0, 1], broadcast=True)
        datasets.set("replaced", [0, 1])
        datasets.set("shrunk", [0, 1, 2])
        datasets.set("scalar", 5, broadcast=True)
        datasets.set("array", np.zeros(3))
        snapshot = snapshot_list_datasets(datasets)
        self.assertEqual(set(snapshot), {"grown", "replaced", "shrunk"})

        datasets.append_to("grown", 2)
        datasets.set("replaced", [0, 9, 9, 9])
        del datasets.local["shrunk"][1:]
        datasets.set("created", [0, 7])
        datasets.set("scalar", 6, broadcast=True)
        removed = roll_back_list_datasets(datasets, snapshot)

        self.assertEqual(removed, {"grown": 1})
        self.assertEqual(datasets.local["grown"], [0, 1])
        self.assertEqual(datasets.local["replaced"], [0, 9, 9, 9])
        self.assertEqual(datasets.local["shrunk"], [0])
        self.assertEqual(datasets.local["created"], [0, 7])
        self.assertEqual(datasets.local["scalar"], 6)
        self.assertEqual(roll_back_list_datasets(datasets, snapshot), {})

    def test_catch_underflow_defaults(self):
        """On by default since 2026-10-28, and a point that keeps
        underflowing stops the run rather than leaving a gap that would
        shift every later point's label."""
        defaults = {}
        scan = self._make_repository_examination_experiment(
            GeneralVariableScan_master_satellite_two_nodes
        )
        scan.setattr_argument = lambda name, processor=None, *args, **kwargs: (
            defaults.__setitem__(
                name, getattr(processor, "default_value", processor)),
            setattr(scan, name, None),
        )
        scan.build()
        self.assertIs(defaults["enable_Catch_UnderFlow"], True)
        self.assertIs(defaults["skip_only_that_iteration_if_exhausted"], False)

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


import copy  # noqa: E402

import GeneralVariableScan_master_satellite_mixin as gvs_mixin  # noqa: E402
from GeneralVariableScan_master_satellite_mixin import (  # noqa: E402
    HEALTH_CHECK_SCANS,
    format_override_dictionary,
    republish_run_datasets,
)
from utilities.BaseExperiment_master_satellite import (  # noqa: E402
    BaseExperimentMasterSatellite,
)

# The class the mixin catches, whichever test module installed the stub.
TerminationRequested = gvs_mixin.TerminationRequested

_OPTIMIZER_CLASS = "MicrowaveScanOptimizer_master_satellite"
_CENTERS = {flag: center for _, flag, center in HEALTH_CHECK_SCANS}


class _FakeComm:
    def __init__(self, log):
        self.log = log

    def close(self):
        self.log.append(("close",))


class HealthCheckScheduler:
    """ARTIQ 7's scheduler as a running GVS sees it (artiq/master/scheduler.py
    with the worker's Scheduler on top), driving a scripted optimizer.

    Every get_status() moves each run one step: pending -> preparing ->
    prepare_done, and a finished run from analyzing to gone. check_pause() is
    True when a prepared run has a higher priority than the GVS, or when
    termination was requested; pause() then runs the prepared runs, highest
    priority first, as the RunStage does while the GVS sits paused, and only
    THEN raises TerminationRequested if Terminate was pressed (worker_impl
    raises it when the master resumes the run). terminate_when_prepared=rid
    presses Terminate as that run reaches prepare_done.

    *plan* scripts the optimizer, one entry per run, check runs and the scans
    they submit alike: "pass", "fail" or "crash" (no result at all). A failed
    CHECK submits its scan at 99 with parent_rid carried on and
    n_measurements = n_measurements_optimization_scan when above 0, as
    submit_opt_exp_general does; a passing SCAN writes the new center.
    delete() only records, so a test can assert that nothing was deleted.
    """

    def __init__(self, master, log, plan=(), rid=500, priority=0,
                 pipeline_name="main", new_center=None,
                 prepare_scans_at_once=False, terminate_when_prepared=None):
        self.master = master
        self.log = log
        self.plan = list(plan)
        self.rid = rid
        self.priority = priority
        self.pipeline_name = pipeline_name
        self.new_center = new_center or {}
        self.prepare_scans_at_once = prepare_scans_at_once
        self.next_rid = rid + 1
        self.runs = {}
        self.submitted = []
        self.deleted = []
        self.pauses = 0
        self.terminate_at_pause = None
        self.terminate_when_prepared = terminate_when_prepared
        self.termination_requested = False

    def submit(self, pipeline_name=None, expid=None, priority=None,
               due_date=None, flush=False):
        rid = self.next_rid
        self.next_rid += 1
        active = [r for r, run in self.runs.items()
                  if run["status"] not in ("analyzing",)]
        self.log.append(("submit", rid, expid["arguments"]["which_node"],
                         [f for f in _CENTERS if expid["arguments"].get(f)],
                         expid["arguments"]["run_health_check_and_optimize"],
                         tuple(active)))
        self.runs[rid] = {"expid": expid, "status": "pending",
                          "priority": priority, "pipeline": pipeline_name}
        self.submitted.append((rid, pipeline_name, expid, priority))
        return rid

    def get_status(self):
        for rid in list(self.runs):
            run = self.runs[rid]
            if run["status"] == "analyzing":
                run["status"] = "gone"
            elif run["status"] == "pending":
                run["status"] = "preparing"
            elif run["status"] == "preparing":
                run["status"] = "prepare_done"
                if rid == self.terminate_when_prepared:
                    self.termination_requested = True
        for rid in [r for r, run in self.runs.items()
                    if run["status"] == "gone"]:
            del self.runs[rid]
        return {rid: {"expid": run["expid"], "status": run["status"],
                      "priority": run["priority"],
                      "pipeline": run["pipeline"]}
                for rid, run in self.runs.items()}

    def _prepared(self):
        return [rid for rid, run in self.runs.items()
                if run["status"] == "prepare_done"
                and run["priority"] > self.priority]

    def check_pause(self, rid=None):
        return self.termination_requested or bool(self._prepared())

    def check_termination(self, rid=None):
        return self.termination_requested

    def pause(self):
        self.pauses += 1
        self.log.append(("pause",))
        if self.terminate_at_pause == self.pauses:
            self.termination_requested = True
        while self._prepared():
            rid = max(self._prepared(),
                      key=lambda r: (self.runs[r]["priority"], -r))
            self.runs[rid]["status"] = "running"
            self._run_optimizer(rid, self.runs[rid]["expid"])
            self.runs[rid]["status"] = "analyzing"
            if self.prepare_scans_at_once:
                for run in self.runs.values():
                    if run["status"] == "pending":
                        run["status"] = "prepare_done"
        if self.prepare_scans_at_once:
            # finished and analysed before the GVS looks again
            for rid in [r for r, run in self.runs.items()
                        if run["status"] == "analyzing"]:
                del self.runs[rid]
        if self.termination_requested:
            raise TerminationRequested

    def delete(self, rid):
        self.deleted.append(rid)
        self.runs.pop(rid, None)

    def _run_optimizer(self, rid, expid):
        arguments = expid["arguments"]
        outcome = self.plan.pop(0)
        self.log.append(("optimizer", rid, outcome))
        if outcome == "crash":
            return
        flag = [f for f in _CENTERS if arguments[f]][0]
        node = arguments["which_node"]
        passed = outcome == "pass"
        if passed and not arguments["run_health_check_and_optimize"]:
            self.master[f"{_CENTERS[flag]}_{node}"] = self.new_center.get(
                (node, flag), 1.0e6 + rid)
        # what the optimizer run leaves in the master DB besides its result
        self.master["n_measurements"] = arguments["n_measurements"]
        fidelity = 0.9 if passed else 0.3
        for name, value in (("rid", rid),
                            ("parent_rid", arguments["parent_rid"]),
                            ("node", node), ("scan", flag),
                            ("fidelity", fidelity),
                            ("passed", passed)):
            self.master[f"health_check_result_{name}"] = value
        parent = int(arguments["parent_rid"])
        if parent != rid:
            self.master[f"health_check_result_for_parent_{parent}"] = [
                rid, node, flag, fidelity, int(passed)]
        if arguments["run_health_check_and_optimize"] and not passed:
            child = copy.deepcopy(expid)
            child["arguments"]["run_health_check_and_optimize"] = False
            scan_n = int(arguments.get("n_measurements_optimization_scan", 0))
            if scan_n > 0:
                child["arguments"]["n_measurements"] = scan_n
            self.submit("main", child, 99, None, False)


class GeneralVariableScanHealthCheckTests(unittest.TestCase):
    """The MW health check (run_health_check, 2026-10-29): pause, delegate to
    MicrowaveScanOptimizer_master_satellite, wait, resubmit until passed."""

    @staticmethod
    def _names_from_real_base(base):
        real = BaseExperimentMasterSatellite(experiment=None)
        base.node_variable_names = real.node_variable_names
        base.two_node_shared_variable_map = real.two_node_shared_variable_map
        base.global_variable_names = real.global_variable_names

    def _make_scan(self, scan_class=GeneralVariableScan_master_satellite_two_nodes,
                   selection="Node1", flags=("Frequency_11_Scan",), every_n=1,
                   sequence=(1.0, 2.0, 3.0), plan=(), overrides=None, **scheduler):
        scan = scan_class()
        scan.base = FakeBase("single_node", "Node2")
        scan.base.experiment = scan
        self._names_from_real_base(scan.base)
        if scan_class.EXPERIMENT_MODE == "single_node":
            scan.selected_node = "Node2"
        for variable in (*MASTER_SATELLITE_VARIABLES, *NODE1_VARIABLES,
                         *NODE2_VARIABLES):
            setattr(scan, variable.name, variable.value)
        scan.f_FORT_Node2 = 240.0
        scan.scan_variable1 = "f_FORT_Node2"
        scan.scan_variable2 = None
        scan.scan_sequence1 = list(sequence)
        scan.scan_sequence2 = [0.0]
        scan.authoritative_overrides = dict(overrides or {})
        scan.enable_Catch_UnderFlow = True
        scan.underflow_max_retries = 3
        scan.underflow_backoff_ms = 0.0
        scan.skip_only_that_iteration_if_exhausted = False
        scan.create_applets = False
        scan.experiment_name = "atom_loading_experiment"
        scan.scan_var_filesuffix = "f_FORT"
        scan._initialize_run_state = lambda: None
        scan.initialize_hardware = lambda: None
        scan.health_check_poll_interval = 0.0
        scan.health_check_crash_backoff_s = 0

        scan.log = []
        scan.saves = 0

        def write_results(kwargs={}):
            scan.saves += 1
            scan.log.append(("save",))

        scan.write_results = write_results
        scan.status_applets = 0
        scan._create_health_check_status_applet = lambda: setattr(
            scan, "status_applets", scan.status_applets + 1)
        scan.core = types.SimpleNamespace(comm=_FakeComm(scan.log))
        scan.master = {}
        scan.scheduler = HealthCheckScheduler(scan.master, scan.log, plan,
                                              **scheduler)
        dataset_mgr = FakeDatasetManager()
        scan._HasEnvironment__dataset_mgr = dataset_mgr
        scan.dataset_writes = []
        scan.dataset_reads = []

        def set_dataset(name, value, **kwargs):
            scan.dataset_writes.append((name, value, kwargs))
            dataset_mgr.set(name, value, **kwargs)

        def get_dataset(name, *args, **kwargs):
            scan.dataset_reads.append((name, kwargs.get("archive", True)))
            if name in dataset_mgr.local:
                return dataset_mgr.local[name]
            if name in scan.master:
                return scan.master[name]
            if args:
                return args[0]
            raise KeyError(name)

        scan.set_dataset = set_dataset
        scan.get_dataset = get_dataset

        def selected_function(experiment):
            experiment.log.append(("point", experiment.f_FORT_Node2))

        scan._selected_experiment_function = selected_function

        scan.run_health_check = True
        scan.health_check_every_n_ite = every_n
        scan.health_check_node_selection = selection
        for _, flag, _ in HEALTH_CHECK_SCANS:
            setattr(scan, flag, flag in flags)
        scan.target_fidelity = 0.8
        # neither the GUI defaults (200, 100) nor the optimizer's (100, 0)
        scan.health_check_n_measurements = 50
        scan.health_check_optimization_n_measurements = 20
        # neither the GUI defaults (1.0, 4.0) nor the optimizer's (2.0, 3.0)
        scan.health_check_adaptive_half_range_linewidths = 1.5
        scan.health_check_adaptive_points_per_linewidth = 5.0
        return scan

    def _prepared_scan(self, *args, **kwargs):
        scan = self._make_scan(*args, **kwargs)
        scan._prepare_health_check()
        return scan

    @staticmethod
    def _submissions(scan):
        return [entry for entry in scan.log if entry[0] == "submit"]

    @staticmethod
    def _points_and_checks(scan):
        """('point', value) and ('check', node, flag) in the order they ran."""
        order = []
        for entry in scan.log:
            if entry[0] == "point":
                order.append(entry)
            elif entry[0] == "submit" and entry[4]:
                order.append(("check", entry[2], entry[3][0]))
        return order

    def test_build_declares_the_health_check_inert_by_default(self):
        from unittest import mock
        scan = GeneralVariableScan_master_satellite_two_nodes()
        defaults = {}
        scan.setattr_argument = lambda name, value, *a, **k: defaults.__setitem__(
            name, value)
        enumerations = []

        def enumeration(values, **kwargs):
            enumerations.append((list(values), kwargs))
            return kwargs.get("default", tuple(values)[0])

        with mock.patch.object(gvs_mixin, "EnumerationValue", enumeration):
            scan._build_health_check_arguments()
        # an explicit default, or ARTIQ raises DefaultMissing for every expid
        # submitted before this argument existed (saved, artiq_client)
        self.assertEqual(enumerations, [(["Node1", "Node2", "both_nodes"],
                                         {"default": "Node1"})])
        self.assertIs(defaults["run_health_check"], False)
        self.assertIs(GeneralVariableScan_master_satellite_two_nodes.run_health_check,
                      False)
        self.assertEqual(defaults["health_check_every_n_ite"], 1)
        self.assertEqual(defaults["health_check_node_selection"], "Node1")
        for _, flag, _ in HEALTH_CHECK_SCANS:
            self.assertIs(defaults[flag], False, flag)
        self.assertEqual(defaults["target_fidelity"], 0.80)
        self.assertEqual(defaults["health_check_n_measurements"], 200)
        self.assertEqual(defaults["health_check_optimization_n_measurements"], 100)
        self.assertEqual(
            GeneralVariableScan_master_satellite_two_nodes.health_check_crash_backoff_s,
            60)
        self.assertEqual(defaults["health_check_adaptive_half_range_linewidths"], 1.0)
        self.assertEqual(defaults["health_check_adaptive_points_per_linewidth"], 4.0)
        # no ExperimentVariable shares a name with a new argument
        variables = {v.name for v in (*MASTER_SATELLITE_VARIABLES,
                                      *NODE1_VARIABLES, *NODE2_VARIABLES)}
        self.assertEqual(set(defaults) & variables, set())

    def test_off_by_default_changes_nothing(self):
        """A hand-built scan without the arguments, and no scheduler at all:
        the class default keeps the scan exactly as before."""
        tests = GeneralVariableScanMasterSatelliteTests()
        scan = tests._make_catch_error_scan(
            GeneralVariableScan_master_satellite_two_nodes)
        scan._selected_experiment_function = lambda experiment: None
        self.assertFalse(hasattr(scan, "scheduler"))
        scan._prepare_health_check()
        scan.run()
        written = {name for name, _, _ in scan.dataset_writes}
        self.assertFalse({n for n in written if "health_check" in n})

    def test_cadence_never_before_the_first_point(self):
        for every_n, checked_before in ((1, [1, 2, 3, 4, 5]), (3, [3])):
            with self.subTest(every_n=every_n):
                scan = self._prepared_scan(
                    every_n=every_n, sequence=(0.0, 1.0, 2.0, 3.0, 4.0, 5.0),
                    plan=["pass"] * 10)
                scan.run()
                order = self._points_and_checks(scan)
                checks = [index for index, entry in enumerate(order)
                          if entry[0] == "check"]
                points_before = [sum(1 for e in order[:i] if e[0] == "point")
                                 for i in checks]
                self.assertEqual(points_before, checked_before)
                self.assertEqual(
                    [e[1] for e in order if e[0] == "point"],
                    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])

    def test_cadence_counts_flat_points_of_paired_and_grid_scans(self):
        for paired, points in ((True, 3), (False, 6)):
            with self.subTest(paired=paired):
                scan = self._prepared_scan(every_n=2, sequence=(1.0, 2.0, 3.0),
                                           plan=["pass"] * 10)
                scan.scan_variable2 = "dummy_variable"
                scan.scan_sequence2 = [10.0, 20.0] if not paired else [10.0, 20.0, 30.0]
                scan.pair_scan_sequences = paired
                scan.run()
                order = self._points_and_checks(scan)
                self.assertEqual(sum(1 for e in order if e[0] == "point"), points)
                checks = [sum(1 for e in order[:i] if e[0] == "point")
                          for i, e in enumerate(order) if e[0] == "check"]
                self.assertEqual(checks, list(range(2, points, 2)))

    def test_not_repeated_on_an_underflow_retry(self):
        scan = self._prepared_scan(sequence=(1.0, 2.0), plan=["pass"] * 4)
        attempts = []

        def selected_function(experiment):
            attempts.append(experiment.f_FORT_Node2)
            experiment.log.append(("point", experiment.f_FORT_Node2))
            if attempts.count(2.0) == 1 and experiment.f_FORT_Node2 == 2.0:
                raise RTIOUnderflow("transient")

        scan._selected_experiment_function = selected_function
        scan.run()
        self.assertEqual(attempts, [1.0, 2.0, 2.0])
        self.assertEqual(len(self._submissions(scan)), 1)

    def test_one_at_a_time_in_node_then_transition_order(self):
        scan = self._prepared_scan(
            selection="both_nodes",
            flags=("Frequency_m11_Scan", "Frequency_00_Scan", "Frequency_11_Scan"),
            sequence=(1.0, 2.0), plan=["pass"] * 6)
        scan.run()
        checks = [(entry[2], entry[3][0]) for entry in self._submissions(scan)]
        self.assertEqual(checks, [
            ("Node1", "Frequency_00_Scan"), ("Node1", "Frequency_11_Scan"),
            ("Node1", "Frequency_m11_Scan"),
            ("Node2", "Frequency_00_Scan"), ("Node2", "Frequency_11_Scan"),
            ("Node2", "Frequency_m11_Scan")])
        # nothing of an earlier chain was still queued at any submission
        for entry in self._submissions(scan):
            self.assertEqual(entry[5], (), entry)
        for _, pipeline, expid, priority in scan.scheduler.submitted:
            self.assertEqual((pipeline, priority), ("main", 99))
            self.assertEqual(expid["arguments"]["parent_rid"], 500)

    def test_waits_for_the_scan_a_failed_check_submits(self):
        for at_once in (False, True):
            with self.subTest(scan_finished_before_the_gvs_resumed=at_once):
                scan = self._prepared_scan(
                    sequence=(1.0, 2.0), plan=["fail", "pass"],
                    new_center={("Node1", "Frequency_11_Scan"): 339.1e6},
                    prepare_scans_at_once=at_once)
                scan.run()
                ran = [entry for entry in scan.log if entry[0] in
                       ("optimizer", "point")]
                self.assertEqual(ran, [("point", 1.0), ("optimizer", 501, "fail"),
                                       ("optimizer", 502, "pass"), ("point", 2.0)])
                local = scan._HasEnvironment__dataset_mgr.local
                self.assertEqual(local["gvs_health_check_resubmissions"], [0])
                self.assertEqual(local["gvs_health_check_fidelity"], [0.9])
                self.assertEqual(local["gvs_health_check_chain_rids"],
                                 "scan point 1 Node1 freq11: 501,502")
                self.assertEqual(scan.f_microwaves_11_dds_Node1, 339.1e6)

    def test_core_released_before_every_pause(self):
        scan = self._prepared_scan(sequence=(1.0, 2.0, 3.0),
                                   plan=["fail", "fail", "pass", "pass"])
        scan.run()
        pauses = [i for i, entry in enumerate(scan.log) if entry[0] == "pause"]
        self.assertGreaterEqual(len(pauses), 4)
        for i in pauses:
            self.assertEqual(scan.log[i - 1], ("close",), scan.log[i - 3:i + 1])

    def test_resubmits_the_same_experiment_until_passed_and_counts(self):
        scan = self._prepared_scan(
            sequence=(1.0, 2.0, 3.0),
            plan=["fail", "fail", "fail", "crash", "pass",
                  "fail", "pass"])
        scan.run()
        checks = [s for s in scan.scheduler.submitted
                  if s[2]["arguments"]["run_health_check_and_optimize"]]
        # point 1: three check submissions (fail+fail = failed, fail+crash =
        # crashed: the scan left no result; pass); point 2: one (fail, then
        # the scan passes)
        self.assertEqual(len(checks), 4)
        first = checks[0][2]
        for _, _, expid, _ in checks[1:3]:
            self.assertEqual(expid, first)
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_resubmissions"], [2, 0])
        self.assertEqual(local["gvs_health_check_crashes"], [1, 0])
        self.assertEqual(local["gvs_health_check_scan_point"], [1, 2])
        self.assertEqual(scan._health_check_totals, {("Node1", "11"): 1})
        self.assertEqual(scan._health_check_crash_totals, {("Node1", "11"): 1})
        self.assertEqual(scan._health_check_resubmitted_points,
                         {("Node1", "11"): [1]})
        self.assertEqual(scan._health_check_crash_points,
                         {("Node1", "11"): [1]})
        status = scan._HasEnvironment__dataset_mgr._broadcaster.raw_view[
            "health_check_status"]
        self.assertEqual(status[0], False)
        self.assertIn("Node1 freq11: 1 failed (at scan point 1), "
                      "1 crashed (at scan point 1)", status[1])
        self.assertIn("total: 1 failed, 1 crashed", status[1])

    def test_a_crashed_chain_is_resubmitted_after_the_backoff(self):
        """No result from its RIDs -- here an older passed one is still in
        the master DB -- so the same check is submitted again, but only
        after health_check_crash_backoff_s, slept in steps of about a second
        with Terminate asked before each."""
        from unittest import mock
        scan = self._prepared_scan(sequence=(1.0, 2.0), plan=["crash", "pass"])
        scan.health_check_crash_backoff_s = 3.5
        scan.master.update({
            "health_check_result_rid": 42, "health_check_result_parent_rid": 41,
            "health_check_result_node": "Node1",
            "health_check_result_scan": "Frequency_11_Scan",
            "health_check_result_fidelity": 0.99,
            "health_check_result_passed": True,
            "health_check_result_for_parent_500":
                [42, "Node1", "Frequency_11_Scan", 0.99, 1]})
        asked = []
        original = scan.scheduler.check_termination

        def check_termination(rid=None):
            asked.append(len([s for s in sleeps if s > 0]))
            return original(rid)

        scan.scheduler.check_termination = check_termination
        sleeps = []
        with mock.patch.object(gvs_mixin.time, "sleep",
                               lambda s: sleeps.append(s)):
            scan.run()
        backoff = [s for s in sleeps if s > 0]
        self.assertEqual(backoff, [1.0, 1.0, 1.0, 0.5])
        # Terminate asked before the first sleep and after every one
        for n in range(len(backoff) + 1):
            self.assertIn(n, asked)
        # the resubmission came after the whole backoff
        submits = [i for i, entry in enumerate(scan.log) if entry[0] == "submit"]
        self.assertEqual(len(submits), 2)
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_resubmissions"], [1])
        self.assertEqual(local["gvs_health_check_crashes"], [1])
        self.assertEqual(scan._health_check_totals, {("Node1", "11"): 0})
        self.assertEqual(local["gvs_health_check_chain_rids"],
                         "scan point 1 Node1 freq11: 501,502")
        self.assertIn("ended with no result of its own (crashed); "
                      "resubmitting in 4 s", local["gvs_health_check_log"])

    def test_a_failed_check_is_resubmitted_at_once_at_99(self):
        """A genuine failure (the scan's fit check did not pass) waits for no
        backoff and goes in at 99, so it runs next."""
        from unittest import mock
        scan = self._prepared_scan(sequence=(1.0, 2.0),
                                   plan=["fail", "fail", "pass"])
        scan.health_check_crash_backoff_s = 60
        sleeps = []
        with mock.patch.object(gvs_mixin.time, "sleep",
                               lambda s: sleeps.append(s)):
            scan.run()
        self.assertEqual([s for s in sleeps if s > 0], [])
        self.assertEqual([priority for *_, priority in scan.scheduler.submitted],
                         [99, 99, 99])
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_resubmissions"], [1])
        self.assertEqual(local["gvs_health_check_crashes"], [0])

    def test_a_failed_check_whose_scan_crashed_is_a_crash(self):
        """A failed check always submits its scan, and the scan always
        publishes its fit check: a chain whose only result is the check's
        failure lost its scan (a FitError) and waits for the backoff."""
        scan = self._prepared_scan(sequence=(1.0, 2.0),
                                   plan=["fail", "crash", "pass"])
        outcomes = []
        original = scan._read_health_check_result

        def read(*args):
            outcome = original(*args)
            outcomes.append(outcome)
            return outcome

        scan._read_health_check_result = read
        scan.run()
        self.assertEqual([o for o, _ in outcomes], ["crashed", "passed"])
        self.assertEqual(outcomes[0][1], 0.3)      # the failed check's
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_crashes"], [1])
        self.assertEqual(local["gvs_health_check_chain_rids"],
                         "scan point 1 Node1 freq11: 501,502,503")

    def test_terminate_during_the_crash_backoff_ends_cleanly(self):
        from unittest import mock
        scan = self._prepared_scan(sequence=(1.0, 2.0, 3.0),
                                   plan=["crash", "pass"])
        scan.health_check_crash_backoff_s = 60
        sleeps = []

        def sleep(seconds):
            sleeps.append(seconds)
            if len([s for s in sleeps if s > 0]) == 5:
                scan.scheduler.termination_requested = True

        with mock.patch.object(gvs_mixin.time, "sleep", sleep):
            scan.run()     # ends cleanly, no exception
        self.assertEqual([s for s in sleeps if s > 0], [1.0] * 5)
        self.assertEqual(len(scan.scheduler.submitted), 1)
        self.assertEqual(scan.scheduler.deleted, [])
        self.assertEqual([e[1] for e in scan.log if e[0] == "point"], [1.0])
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_terminated_at_scan_point"], 1)
        self.assertEqual(scan.log[-1], ("save",))

    def test_an_optimizer_run_during_the_pause_does_not_change_the_outcome(self):
        """A run queued between the GVS's priority and 99 (here a user's
        failing Node1 freq11 optimizer, parent 777) runs after the chain while
        the GVS is still paused and overwrites health_check_result_*. The
        check passed and must count as passed."""
        scan = self._prepared_scan(sequence=(1.0, 2.0), plan=["pass", "fail"])
        original = scan.scheduler._run_optimizer

        def then_queue_an_unrelated_optimizer(rid, expid):
            original(rid, expid)
            if rid == 501:
                other = copy.deepcopy(expid)
                other["arguments"]["parent_rid"] = 777
                other["arguments"]["run_health_check_and_optimize"] = False
                scan.scheduler.runs[900] = {
                    "expid": other, "status": "prepare_done", "priority": 50,
                    "pipeline": "main"}

        scan.scheduler._run_optimizer = then_queue_an_unrelated_optimizer
        scan.run()
        ran = [entry for entry in scan.log if entry[0] in ("optimizer", "point")]
        self.assertEqual(ran, [("point", 1.0), ("optimizer", 501, "pass"),
                               ("optimizer", 900, "fail"), ("point", 2.0)])
        # the latest result is the other run's ...
        self.assertEqual(scan.master["health_check_result_rid"], 900)
        # ... and the check still counts as passed
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_resubmissions"], [0])
        self.assertEqual(local["gvs_health_check_fidelity"], [0.9])
        self.assertEqual(scan._health_check_totals, {("Node1", "11"): 0})
        self.assertEqual(len(scan.scheduler.submitted), 1)

    def test_rereads_only_that_center_from_the_master_db(self):
        scan = self._prepared_scan(
            selection="both_nodes", flags=("Frequency_01_Scan",),
            sequence=(1.0, 2.0), plan=["fail", "pass", "pass"],
            new_center={("Node1", "Frequency_01_Scan"): 336.97e6},
            overrides={"t_FORT_drop_Node1": 1e-5})
        before_node2 = scan.f_microwaves_01_dds_Node2
        reapplied = []
        scan._apply_run_wide_overrides = lambda: reapplied.append(True)
        scan.run()
        self.assertEqual(scan.f_microwaves_01_dds_Node1, 336.97e6)
        self.assertEqual(scan.f_microwaves_01_dds_Node2, before_node2)
        self.assertEqual(scan.base.reload_calls, 0)
        # mid-run, an override pass would re-write broadcast keys locally
        self.assertEqual(reapplied, [])
        self.assertFalse([name for name, _, _ in scan.dataset_writes
                          if name == "t_FORT_drop_Node1"])
        reads = [(name, archive) for name, archive in scan.dataset_reads
                 if not name.startswith("health_check_result_")]
        self.assertEqual(reads, [("f_microwaves_01_dds_Node1", False),
                                 ("f_microwaves_01_dds_Node2", False)])
        self.assertTrue(all(archive is False for _, archive in scan.dataset_reads))
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_center_after_Hz"],
                         [336.97e6, before_node2])
        self.assertEqual(local["gvs_health_check_node"], [1, 2])

    def test_single_node_other_node_check_leaves_this_run_alone(self):
        scan = self._prepared_scan(
            GeneralVariableScan_master_satellite_single_node, selection="Node1",
            flags=("Frequency_00_Scan",), sequence=(1.0, 2.0),
            plan=["fail", "pass"])
        scan.master["f_microwaves_00_dds_Node1"] = 334.6e6
        scan.run()
        # read from the master DB before and after, never set on this run
        self.assertEqual(
            [r for r in scan.dataset_reads if r[0].startswith("f_micro")],
            [("f_microwaves_00_dds_Node1", False)] * 2)
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_center_before_Hz"], [334.6e6])
        # and the other node is driven safe again: at the start and after
        self.assertEqual(scan.base.force_off_calls, 2)

    def test_republishes_this_runs_broadcasts_and_n_measurements_only(self):
        scan = self._prepared_scan(sequence=(1.0, 2.0), plan=["pass"])
        datasets = scan._HasEnvironment__dataset_mgr
        datasets.set("AllSPCMs_RO1", [0, 5], broadcast=True)
        datasets.set("scan_sequence1", [1.0, 2.0], broadcast=True)
        datasets.set("n_measurements", 12, broadcast=True, persist=True)
        datasets.set("p_AOM_A1_history_Node1", [-3.0], broadcast=True,
                     persist=True, archive=False)
        datasets.set("f_microwaves_11_dds_Node1", 1.0, broadcast=True,
                     persist=True, archive=False)
        mark = len(datasets._broadcaster.published)
        record = datasets.set

        def logged_set(key, value, **kwargs):
            scan.log.append(("dataset", key))
            record(key, value, **kwargs)

        datasets.set = logged_set
        scan.run()
        after = datasets._broadcaster.published[mark:]
        republished = {m[1]: m[2] for m in after if m[0] == "setitem"}
        # republished only after the event's last optimizer run
        last_optimizer = max(i for i, entry in enumerate(scan.log)
                             if entry[0] == "optimizer")
        for key in ("AllSPCMs_RO1", "scan_sequence1", "n_measurements"):
            writes = [i for i, entry in enumerate(scan.log)
                      if entry == ("dataset", key)]
            self.assertTrue(writes, key)
            self.assertTrue(all(i > last_optimizer for i in writes), key)
        self.assertEqual(republished["AllSPCMs_RO1"], False)
        self.assertEqual(republished["scan_sequence1"], False)
        self.assertEqual(republished["n_measurements"], True)
        self.assertNotIn("p_AOM_A1_history_Node1", republished)
        self.assertNotIn("f_microwaves_11_dds_Node1", republished)
        self.assertNotIn("delitem", [m[0] for m in after])
        self.assertEqual(datasets._broadcaster.raw_view["n_measurements"],
                         (True, 12))
        self.assertIn("AllSPCMs_RO1", datasets.local)
        self.assertNotIn("p_AOM_A1_history_Node1", datasets.local)

    def test_republish_run_datasets_keeps_flags(self):
        datasets = FakeDatasetManager()
        value = [1, 2]
        datasets.set("broadcast_only", value, broadcast=True, archive=False)
        datasets.set("both", [3], broadcast=True)
        datasets.set("local_only", [4])
        datasets.set("n_measurements", 7, broadcast=True, persist=True)
        datasets.set("calibration", 1.5, broadcast=True, persist=True)
        published = republish_run_datasets(datasets)
        self.assertEqual(set(published),
                         {"broadcast_only", "both", "n_measurements"})
        self.assertIs(datasets._broadcaster.raw_view["broadcast_only"][1], value)
        self.assertNotIn("broadcast_only", datasets.local)
        self.assertEqual(datasets._broadcaster.raw_view["both"],
                         (False, datasets.local["both"]))
        self.assertNotIn("local_only", datasets._broadcaster.raw_view)

    def test_records_are_local_and_unsuffixed(self):
        scan = self._prepared_scan(selection="both_nodes", sequence=(1.0, 2.0),
                                   plan=["pass", "fail", "pass"])
        scan.run()
        records = [(name, kwargs) for name, _, kwargs in scan.dataset_writes
                   if name.startswith("gvs_health_check")]
        self.assertTrue(records)
        for name, kwargs in records:
            self.assertFalse(name.endswith(("_Node1", "_Node2")), name)
            self.assertIs(kwargs.get("broadcast"), False, name)
            self.assertIs(kwargs.get("persist"), False, name)
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_transition_labels"],
                         "00,01,11,m10,m11")
        self.assertEqual(local["gvs_health_check_transition"], [2, 2])
        self.assertEqual(local["gvs_health_check_node"], [1, 2])
        self.assertEqual(len(local["gvs_health_check_time_end"]), 2)
        self.assertIn("forwarded override_ExperimentVariables_Node1",
                      local["gvs_health_check_log"])
        status_writes = [kwargs for name, _, kwargs in scan.dataset_writes
                         if name == "health_check_status"]
        self.assertTrue(all(k == {"broadcast": True, "persist": False}
                            for k in status_writes))
        self.assertGreaterEqual(scan.status_applets, 2)

    def test_prepare_refusals(self):
        cases = (
            ({"priority": 99}, "priority below 99"),
            ({"pipeline_name": "other"}, "'main' pipeline"),
            ({"flags": ()}, "no transition is ticked"),
            ({"every_n": 0}, "at least 1"),
            ({"overrides": {"f_microwaves_11_dds_Node1": 1.0}},
             "overrides or scans f_microwaves_11_dds_Node1"),
        )
        for kwargs, message in cases:
            with self.subTest(kwargs=kwargs):
                scan = self._make_scan(**kwargs)
                with self.assertRaisesRegex(ValueError, message):
                    scan._prepare_health_check()
        for check_n in (0, -1):
            with self.subTest(health_check_n_measurements=check_n):
                scan = self._make_scan()
                scan.health_check_n_measurements = check_n
                with self.assertRaisesRegex(
                        ValueError, "health_check_n_measurements must be at "
                                    "least 1"):
                    scan._prepare_health_check()
        scan = self._make_scan(selection="Node2")
        scan.scan_variable1 = "f_microwaves_11_dds_Node2"
        with self.assertRaisesRegex(ValueError, "f_microwaves_11_dds_Node2"):
            scan._prepare_health_check()

    def _forwarded_and_merged(self, scan, node):
        forwarded = scan._health_check_forwarded_overrides(node)
        merged = scan._validate_forwarded_overrides(node)
        return forwarded, merged

    def test_forwarded_overrides_two_nodes(self):
        scan = self._make_scan(overrides={
            "t_pumping_two_nodes": 2.0e-4,           # shared, overridden
            "f_FORT_Node1": 241.0e6,                 # Node1 per-node
            "target_780_QWP_Node2": 210.0,           # Node2 per-node
            "two_atom_threshold": 31000.0,           # plain global
            "n_measurements": 7,                     # never
            "t_MW_RF_pulse_Node1": 0.0,              # MW calibration: never
            "f_microwaves_dds_Node1": 1.0,           # the 00/01 scan variable
            "t_SPCM_first_shot_two_nodes": np.float64(0.02),
        })
        scan.scan_variable1 = "t_depumping_two_nodes"
        scan.t_pumping_two_nodes = 2.0e-4
        scan.t_SPCM_first_shot_two_nodes = np.float64(0.02)
        scan.t_depumping_two_nodes = 2.0e-6
        scan._prepare_health_check()
        scan._initialize_health_check_records()   # run(), before point 0
        scan.t_depumping_two_nodes = 7.0e-6       # the current scan point
        forwarded, merged = self._forwarded_and_merged(scan, "Node1")
        shared = set(scan.base.two_node_shared_variable_map)
        self.assertEqual(len(shared), 13)
        # all 13, the scanned one included
        self.assertEqual(set(forwarded) - {"f_FORT_Node1", "two_atom_threshold"},
                         shared)
        # ... at the run's value from before the scan, not the point's
        self.assertEqual(forwarded["t_depumping"], 2.0e-6)
        self.assertEqual(merged["t_depumping_Node1"], 2.0e-6)
        self.assertEqual(forwarded["t_pumping"], 2.0e-4)
        self.assertEqual(forwarded["t_OP_coil_relax_time"],
                         scan.t_OP_coil_relax_time_two_nodes)
        self.assertEqual(forwarded["f_FORT_Node1"], 241.0e6)
        self.assertEqual(forwarded["two_atom_threshold"], 31000.0)
        # resolved by the optimizer's merge to Node1's twins, nothing else
        self.assertEqual(merged["t_pumping_Node1"], 2.0e-4)
        self.assertEqual(merged["f_FORT_Node1"], 241.0e6)
        self.assertEqual(merged["two_atom_threshold"], 31000.0)
        self.assertTrue(all(not name.endswith("_Node2") for name in merged))
        self.assertEqual(len(merged), len(forwarded))
        text = format_override_dictionary(forwarded)
        self.assertIn("'t_SPCM_first_shot': 0.02", text)
        self.assertNotIn("float64", text)

        forwarded, merged = self._forwarded_and_merged(scan, "Node2")
        self.assertEqual(forwarded["target_780_QWP_Node2"], 210.0)
        self.assertNotIn("f_FORT_Node1", forwarded)
        self.assertEqual(merged["t_pumping_Node2"], 2.0e-4)
        self.assertTrue(all(not name.endswith("_Node1") for name in merged))

    def test_two_node_runs_cannot_override_a_shared_settings_twin(self):
        """Why the bare shared entries never meet a twin entry in the
        forwarded dictionary (the optimizer would refuse 'overridden twice'):
        a two-node run's own resolver refuses every per-node twin of the 13
        shared settings, in any override field and as a scan variable."""
        resolver = BaseExperimentMasterSatellite(
            experiment=None, experiment_mode="two_nodes")
        shared = resolver.two_node_shared_variable_map
        self.assertEqual(len(shared), 13)
        for bare in shared:
            for node in ("Node1", "Node2"):
                with self.subTest(twin=f"{bare}_{node}"):
                    with self.assertRaisesRegex(ValueError,
                                                "is not read in two_nodes"):
                        resolver.resolve_experiment_variable_target(
                            f"{bare}_{node}")

    def test_forwarded_overrides_single_node(self):
        scan = self._make_scan(
            GeneralVariableScan_master_satellite_single_node, overrides={
                "p_cooling_DP_MOT_Node2": -9.0,
                "advance_if_atom_in_first_shot": False,
                "f_FORT_Node2": 241.0e6,                 # scanned: never
                "t_microwave_01_pulse_Node2": 5e-6, "n_measurements": 3,
                "target_780_HWP_Node2": np.array([1.0, float("nan")]),
            })
        forwarded, merged = self._forwarded_and_merged(scan, "Node2")
        self.assertEqual(set(forwarded), {"p_cooling_DP_MOT_Node2",
                                          "advance_if_atom_in_first_shot",
                                          "target_780_HWP_Node2"})
        self.assertEqual(set(merged), set(forwarded))
        self.assertEqual(list(merged["target_780_HWP_Node2"][:1]), [1.0])
        # the other node gets the globals only
        forwarded, merged = self._forwarded_and_merged(scan, "Node1")
        self.assertEqual(forwarded, {"advance_if_atom_in_first_shot": False})
        self.assertEqual(merged, forwarded)

    def test_expid_carries_the_forwarded_dictionary_to_its_node_only(self):
        scan = self._prepared_scan(selection="Node2", flags=("Frequency_m10_Scan",))
        expid = scan._health_check_expid("Node2", "Frequency_m10_Scan", "{'a': 1}")
        arguments = expid["arguments"]
        self.assertEqual(expid["class_name"], _OPTIMIZER_CLASS)
        self.assertEqual(expid["file"],
                         "qn_artiq_routines\\MicrowaveScanOptimizer_master_satellite.py")
        self.assertEqual((expid["log_level"], expid["repo_rev"]), (30, "N/A"))
        self.assertEqual(arguments["override_ExperimentVariables_Node2"], "{'a': 1}")
        self.assertEqual(arguments["override_ExperimentVariables_Node1"], "{}")
        self.assertEqual([k for k, v in arguments.items() if v is True],
                         ["run_health_check_and_optimize", "enable_fitting",
                          "enable_geometric_frequency_scan",
                          "enable_node_specific_adaptive_scan",
                          "Frequency_m10_Scan"])
        self.assertIs(arguments["show_shared_applets"], False)
        self.assertEqual(arguments["n_measurements"], 50)
        self.assertEqual(arguments["n_measurements_optimization_scan"], 20)
        self.assertEqual(arguments["target_fidelity"], 0.8)
        self.assertEqual(arguments["adaptive_half_range_linewidths"], 1.5)
        self.assertEqual(arguments["adaptive_points_per_linewidth"], 5.0)

    def test_termination_deletes_nothing_and_logs_the_queue(self):
        """Terminate pressed as a chain run reaches prepare_done -- the check
        run (501), or the scan the failed check submitted (502). The GVS sees
        it BEFORE pause() and ends at once; the prepared run, and one of the
        chain already running elsewhere, stay in the queue (the user,
        2026-10-29: "don't delete the existing queue")."""
        for prepared in (501, 502):
            with self.subTest(terminated_as_prepared=prepared):
                scan = self._prepared_scan(sequence=(1.0, 2.0, 3.0),
                                           plan=["fail", "pass"],
                                           terminate_when_prepared=prepared)
                scan.scheduler.runs[499] = {
                    "expid": {"class_name": _OPTIMIZER_CLASS,
                              "arguments": {"parent_rid": 500.0}},
                    "status": "running", "priority": 99, "pipeline": "main"}
                scan.run()     # ends cleanly, no exception
                self.assertEqual(scan.scheduler.deleted, [])
                self.assertEqual(scan.scheduler.runs[prepared]["status"],
                                 "prepare_done")
                self.assertNotIn(prepared, [e[1] for e in scan.log
                                            if e[0] == "optimizer"])
                self.assertIn(499, scan.scheduler.runs)
                self.assertEqual([e[1] for e in scan.log if e[0] == "point"],
                                 [1.0])
                local = scan._HasEnvironment__dataset_mgr.local
                self.assertEqual(
                    local["gvs_health_check_terminated_at_scan_point"], 1)
                self.assertIn(
                    f"nothing deleted: optimizer run(s) still queued "
                    f"[{prepared}], running [499]",
                    local["gvs_health_check_log"])
                status = scan._HasEnvironment__dataset_mgr._broadcaster.raw_view[
                    "health_check_status"][1]
                self.assertIn(f"optimizer RIDs left in the queue: "
                              f"[{prepared}, 499]", status)
                self.assertEqual(scan.log[-1], ("save",))

    def test_termination_during_a_pause_still_ends_cleanly(self):
        """Terminate pressed while the GVS is already paused: ARTIQ runs the
        prepared chain run first, then pause() raises; the scan still ends
        cleanly."""
        scan = self._prepared_scan(sequence=(1.0, 2.0, 3.0),
                                   plan=["fail", "pass"])
        scan.scheduler.terminate_at_pause = 2
        scan.run()
        self.assertEqual([e for e in scan.log if e[0] == "optimizer"],
                         [("optimizer", 501, "fail"), ("optimizer", 502, "pass")])
        self.assertEqual(scan.scheduler.deleted, [])
        self.assertEqual([e[1] for e in scan.log if e[0] == "point"], [1.0])
        local = scan._HasEnvironment__dataset_mgr.local
        self.assertEqual(local["gvs_health_check_terminated_at_scan_point"], 1)
        self.assertEqual(scan.log[-1], ("save",))

    def test_terminate_pressed_before_a_submission(self):
        scan = self._prepared_scan(selection="both_nodes", sequence=(1.0, 2.0),
                                   plan=["pass"])
        original = scan.scheduler._run_optimizer

        def run_then_terminate(rid, expid):
            original(rid, expid)
            scan.scheduler.termination_requested = True

        scan.scheduler._run_optimizer = run_then_terminate
        scan.run()
        self.assertEqual(len(scan.scheduler.submitted), 1)
        self.assertEqual(scan.scheduler.deleted, [])

    def test_status_applet_request_is_never_fatal(self):
        """The real request, with a CCB that fails: a warning, no exception."""
        scan = self._make_scan()
        del scan._create_health_check_status_applet      # the class method
        scan.base = BaseExperimentMasterSatellite(
            experiment=None, experiment_mode="single_node", which_node="Node1")
        attempts = []

        class BrokenCCB:
            def issue(self, action, *args, **kwargs):
                attempts.append((action, args[0], kwargs.get("group")))
                raise ConnectionError("no dashboard")

        scan.ccb = BrokenCCB()
        with self.assertLogs(level="WARNING") as logs:
            scan._create_health_check_status_applet()
        self.assertEqual(attempts, [("create_applet", "MW health check status",
                                     ["Both nodes", "GVS health check"])])
        self.assertIn("could not create the health check applet",
                      "\n".join(logs.output))

    def test_dummy_scheduler_does_not_hang(self):
        """artiq_run / the compile sweep: RID 0, no status, no
        check_termination -- the check is skipped with a log line."""
        scan = self._prepared_scan(sequence=(1.0, 2.0, 3.0))
        scan.scheduler = types.SimpleNamespace(
            rid=0, priority=0, pipeline_name="main", get_status=lambda: {},
            check_pause=lambda rid=None: False, pause=lambda: None,
            submit=lambda *a, **k: self.fail("submitted"))
        scan._prepare_health_check()
        scan.run()
        self.assertEqual([e[1] for e in scan.log if e[0] == "point"],
                         [1.0, 2.0, 3.0])
        self.assertIn("no ARTIQ master (RID 0)",
                      scan._HasEnvironment__dataset_mgr.local["gvs_health_check_log"])


if __name__ == "__main__":
    unittest.main()
