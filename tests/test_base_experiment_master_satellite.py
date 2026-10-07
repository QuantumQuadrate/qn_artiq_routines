import sys
import types
import unittest
import json
from pathlib import Path


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
    for name, value in {
        "kernel": lambda function: function,
        "rpc": lambda function=None, **kwargs: (
            function if function is not None else (lambda decorated: decorated)
        ),
        "delay": lambda duration: None,
        "EnvExperiment": object,
        "NumberValue": lambda value, **kwargs: value,
        "BooleanValue": lambda value=False, **kwargs: value,
        "StringValue": lambda value, **kwargs: value,
        "MHz": 1e6,
        "kHz": 1e3,
        "ms": 1e-3,
        "us": 1e-6,
        "ns": 1e-9,
    }.items():
        if not hasattr(_experiment_stub, name):
            setattr(_experiment_stub, name, value)


from utilities.BaseExperiment_master_satellite import (  # noqa: E402
    BaseExperimentMasterSatellite,
)
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
        self.assertEqual(experiment.t_Node2_excitation_delay_mu, 189)
        self.assertEqual(experiment.t_Node2_rtio_offset_mu, 0)
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

    def test_n_measurements_is_broadcast_and_persisted_for_applets(self):
        """n_measurements has to reach applets as a DATASET, and stay persistent.

        Applets read broadcast datasets, never experiment attributes, and
        applets/plot_retention_and_loading.py subscripts
        data["n_measurements"] directly. A KeyError there lands in a bare
        `except:` that calls self.clear(), so the dataset going missing shows
        up as a permanently blank applet rather than as any error -- which is
        how it was found.

        persist must stay True as well. n_measurements is a DECLARED global
        that _load_experiment_variables requires to exist, so a
        broadcast-without-persist write would clear its persist flag and the
        master's next flush would drop it from dataset_db.pyon, breaking every
        master-satellite experiment at the following restart. That is exactly
        what the standalone Base does at utilities/BaseExperiment.py:786, and
        it is why the dataset vanished on 2026-10-05.
        """
        experiment, base = self.build_and_prepare("single_node", "Node1")
        # Stand in for the submitted GUI value, which every master-satellite
        # experiment re-asserts over the loaded global before this runs.
        experiment.n_measurements = 37
        experiment.dataset_writes.clear()

        base.initialize_result_state()

        writes = [
            write for write in experiment.dataset_writes
            if write[0] == "n_measurements"
        ]
        self.assertEqual(len(writes), 1)
        _, value, kwargs = writes[0]
        self.assertEqual(value, 37)
        self.assertTrue(kwargs.get("broadcast"))
        self.assertTrue(kwargs.get("persist"))

    def test_two_node_bindings(self):
        experiment, _ = self.build_and_prepare("two_nodes")
        self.assertEqual(experiment.dds_FORT_Node1.name, "urukul0_ch0")
        self.assertEqual(experiment.dds_FORT_Node2.name, "urukul4_ch3")
        self.assertEqual(experiment.sampler0_Node1.name, "sampler0")
        self.assertEqual(experiment.sampler0_Node2.name, "sampler3")
        self.assertEqual(experiment.SPCM_H1.name, "ttl0")
        # ttl_SPCM0 IS published in two-node mode, deliberately: all four
        # canonical detectors are master-local and already see both nodes
        # through the beamsplitter fan-out, so the legacy bare name is the one
        # real counter rather than anything node-specific. It was absent until
        # 2026-10-07, when two-node mode gained the broadcast aliases.
        self.assertIs(experiment.ttl_SPCM0, experiment.SPCM_H1)
        self.assertIs(
            experiment.ttl_SPCM1_OtherNode_counter,
            experiment.SPCM_V2_counter,
        )
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

    def test_two_node_mode_publishes_per_node_wiring_metadata(self):
        """Two-node code must say which crate it means, and the crates differ.

        The absent bare name is the assertion that matters: if coil_channels
        existed in two_nodes mode it would silently be one node's, and the
        other node's coils would be driven through the wrong DAC channels.
        """
        experiment, _ = self.build_and_prepare("two_nodes")

        self.assertEqual(experiment.coil_channels_Node1, [0, 1, 13, 14])
        self.assertEqual(experiment.coil_channels_Node2, [0, 1, 2, 3])
        self.assertEqual(experiment.AX_Zotino_channel_Node1, 13)
        self.assertEqual(experiment.AX_Zotino_channel_Node2, 2)
        self.assertEqual(experiment.UV_trig_channel_Node1, [8])

        self.assertFalse(hasattr(experiment, "coil_channels"))
        self.assertFalse(hasattr(experiment, "AX_Zotino_channel"))

        # Run-global dataset NAMES stay unsuffixed in both modes: they name run
        # state, not a node's wiring.
        self.assertEqual(
            experiment.measurements_progress, "measurements_progress"
        )
        self.assertEqual(experiment.scan_var_dataset, "scan_variables")

    def test_two_node_mode_computes_per_node_derived_amplitudes(self):
        """ampl_* are required in two_nodes, not merely convenient.

        The obvious alternative for a readout amplitude is the stabilizer's
        science setpoint, but prepare_laser_stabilizer refuses outside
        single_node AND FeedbackChannel.amplitudes[1] is 0.0 until feedback has
        run, so reaching for it open-loop would turn the FORT off.
        """
        experiment, _ = self.build_and_prepare("two_nodes")

        for node in ("Node1", "Node2"):
            loading = getattr(experiment, f"ampl_FORT_loading_{node}")
            self.assertAlmostEqual(
                loading,
                dB_to_V(getattr(experiment, f"p_FORT_loading_{node}")),
            )
            self.assertAlmostEqual(
                getattr(experiment, f"ampl_FORT_RO_{node}"),
                loading * getattr(experiment, f"p_FORT_RO_{node}"),
            )
            mot = getattr(experiment, f"ampl_cooling_DP_MOT_{node}")
            self.assertAlmostEqual(
                getattr(experiment, f"ampl_cooling_DP_RO_{node}"),
                mot * getattr(experiment, f"p_cooling_DP_RO_{node}"),
            )
            self.assertAlmostEqual(
                getattr(experiment, f"ampl_AOM_A1_{node}"),
                dB_to_V(getattr(experiment, f"p_AOM_A1_{node}")),
            )

        # The two nodes are genuinely calibrated differently; identical values
        # would mean the suffixing silently read one node twice.
        self.assertNotAlmostEqual(
            experiment.ampl_FORT_loading_Node1,
            experiment.ampl_FORT_loading_Node2,
        )
        self.assertFalse(hasattr(experiment, "ampl_FORT_loading"))
        self.assertFalse(hasattr(experiment, "ampl_cooling_DP_RO"))

    def test_each_stabilizer_binds_only_its_own_nodes_ambient_devices(self):
        """A per-node stabilizer must not touch the other node's hardware.

        AOMPowerStabilizer.run() switches devices that are not feedback
        channels: it blocks both repumps so their light does not contaminate
        the cooling PD reading, opens the cooling DP that the MOT channels
        need, and restores the six MOT AOMs at the end. Those were read as
        bare self.exp.<name> inside the kernel, which is right in standalone
        and single-node mode but wrong in two-node mode, where Base publishes
        the bare names as BROADCAST aliases -- so laser_stabilizer_Node2.run()
        drove Node1's repump, cooling DP and all six fiber AOMs, and vice
        versa.

        Found on hardware 2026-10-07: an RTIOUnderflow traceback inside
        laser_stabilizer_Node2.run() ran through _BroadcastTTLOut.on() ->
        self.node1.on() on channel 5, which is Node1's repump switch. The
        underflow itself was unrelated (slack erosion); it just happened to
        print the misbinding.

        It survived sequential feedback because the two runs leave the same
        end state, and because the broadcast DDS objects expose only .sw -- so
        no amplitude was ever written to the wrong crate. It would be actively
        destructive once the two nodes' feedback runs in parallel, which is
        planned.
        """
        import subroutines.aom_feedback as aom_feedback

        experiment, base = self.build_and_prepare("two_nodes")

        # aom_feedback captures cwd at import time and expects the artiq-master
        # directory; the suite runs from the repository root.
        artiq_master = Path(__file__).resolve().parents[3]
        self.assertTrue(
            (artiq_master / "repository" / "qn_artiq_routines" / "utilities"
             / "config" / "alice" / "feedback_channels.json").is_file(),
            "feedback config not found; adjust the artiq-master path",
        )
        original_cwd = aom_feedback.cwd
        aom_feedback.cwd = str(artiq_master) + "\\"
        try:
            for node, legacy_name in (("Node1", "alice"), ("Node2", "bob")):
                suffix = f"_{node}"
                other_suffix = "_Node2" if node == "Node1" else "_Node1"
                experiment.dds_defaults = base.node_resolvers[node].dds_defaults
                stabilizer = aom_feedback.AOMPowerStabilizer(
                    experiment=experiment,
                    dds_names=["dds_AOM_A1"],
                    iterations=1,
                    averages=1,
                    leave_AOMs_on=False,
                    leave_MOT_AOMs_on=True,
                    node_suffix=suffix,
                    config_node=legacy_name,
                )
                for attribute, device_name in (
                    ("ttl_repump_switch", "ttl_repump_switch"),
                    ("ttl_pumping_repump_switch", "ttl_pumping_repump_switch"),
                    ("dds_cooling_DP", "dds_cooling_DP"),
                    ("GRIN1and2_dds", "GRIN1and2_dds"),
                    ("ttl_exc0_switch", "ttl_exc0_switch"),
                    ("ttl_GRIN1_switch", "ttl_GRIN1_switch"),
                    ("mot_aom_1", "dds_AOM_A1"),
                    ("mot_aom_6", "dds_AOM_A6"),
                ):
                    bound = getattr(stabilizer, attribute)
                    self.assertIs(
                        bound,
                        getattr(experiment, device_name + suffix),
                        f"{node} stabilizer bound the wrong {attribute}",
                    )
                    self.assertIsNot(
                        bound,
                        getattr(experiment, device_name + other_suffix),
                        f"{node} stabilizer bound the OTHER node's "
                        f"{attribute}",
                    )
                    self.assertNotIn(
                        type(bound).__name__,
                        ("_BroadcastTTLOut", "_BroadcastDDS"),
                        f"{node} stabilizer bound a broadcast {attribute}, "
                        f"which drives both crates",
                    )
        finally:
            aom_feedback.cwd = original_cwd

    def test_two_node_retention_applet_thresholds_on_both_atoms(self):
        """The two-node retention applet must discriminate TWO atoms.

        plot_retention_and_loading computes cutoff = t_exposure * threshold
        and then calls RO1 above cutoff "loaded" and RO2 above cutoff
        "retained". In two-node mode loading is JOINT -- all four SPCMs see
        both traps through the beamsplitter fan-out -- so "loaded" has to mean
        an atom in BOTH traps, which is two_atom_threshold. Handing it
        single_atom_threshold would count one-atom events as loaded and make
        the retention denominator wrong.

        The exposure must also be the SUFFIXED name. Bare t_SPCM_first_shot
        still exists in dataset_db as a standalone leftover, so a two-node
        applet pointed at the bare name would scale a live threshold by a
        stale standalone value -- the same collision class as n_measurements.
        """
        from applets_master_satellite import (
            APPLET_SPECS,
            TWO_NODE_APPLET_SPECS,
            build_applet_command,
        )

        title = "retention and loading AllSPCMs"
        two_node_spec = [
            spec for spec in TWO_NODE_APPLET_SPECS if spec.title == title
        ]
        self.assertEqual(len(two_node_spec), 1, title)

        _, base = self.build_and_prepare("two_nodes")
        arguments = build_applet_command(two_node_spec[0], base).split()

        self.assertIn("two_atom_threshold", arguments)
        self.assertNotIn("single_atom_threshold", arguments)
        for argument in arguments:
            self.assertFalse(
                argument.startswith("single_atom_threshold"),
                "the joint criterion must not use a per-node single-atom "
                f"threshold: {argument}",
            )
        self.assertIn(
            "t_SPCM_first_shot_Node1", arguments,
            "the exposure must be the suffixed name; the bare one is a "
            "standalone leftover in dataset_db",
        )
        self.assertNotIn("t_SPCM_first_shot", arguments)

        # Single-node keeps the one-atom criterion. Pinned so a shared edit
        # cannot quietly change what single-node retention means.
        single_spec = [
            spec for spec in APPLET_SPECS if spec.title == title
        ][0]
        _, node1_base = self.build_and_prepare("single_node", "Node1")
        node1_arguments = build_applet_command(single_spec, node1_base).split()
        self.assertIn("single_atom_threshold_Node1", node1_arguments)
        self.assertNotIn("two_atom_threshold", node1_arguments)

    def test_each_mode_retires_only_the_conflicting_applet(self):
        """Starting one mode closes the other mode's retention applet -- and
        nothing else.

        "retention and loading AllSPCMs" exists in both spec sets under the
        SAME title, reading the same unsuffixed AllSPCMs_RO1/RO2, differing
        only in the threshold. Different top-level groups means the dashboard
        keeps both, because create_applet replaces by name only WITHIN a
        group. A leftover copy does not go blank -- it keeps plotting live
        data under the wrong criterion, which is the dangerous failure.

        But ONLY that one. The readout histograms and the atom loading time
        are just as meaningful during a two-node run, so retiring the whole
        idle group to fix one applet would throw away plots the owner wants
        to keep watching.
        """
        import applets_master_satellite as applets

        class RecordingCCB:
            def __init__(self):
                self.disabled_applets = []
                self.disabled_groups = []
                self.created = []

            def issue(self, action, *args, **kwargs):
                if action == "disable_applet":
                    self.disabled_applets.append((args[0], args[1]))
                elif action == "disable_applet_group":
                    self.disabled_groups.append(args[0])
                elif action == "create_applet":
                    self.created.append((args[0], kwargs.get("group")))

        conflicting = applets.conflicting_applet_titles()
        self.assertEqual(conflicting, {"retention and loading AllSPCMs"})

        # Applets that must survive a mode change, because they mean the same
        # thing either way.
        keep_open = {
            "All SPCMs RO1 histogram",
            "All SPCMs RO2 histogram",
            "Atom loading time (s)",
        }
        self.assertTrue(
            keep_open.isdisjoint(conflicting),
            "these must never be retired on a mode change",
        )

        experiment, base = self.build_and_prepare("single_node", "Node1")
        ccb = RecordingCCB()
        experiment.ccb = ccb
        applets.create_applets_for(
            experiment, base,
            specs=applets.APPLET_SPECS,
            shared_specs=applets.SHARED_APPLET_SPECS,
        )
        self.assertEqual(
            ccb.disabled_groups, [],
            "no whole group may be retired -- that takes the histograms with it",
        )
        self.assertIn(
            ("retention and loading AllSPCMs", applets.TWO_NODE_APPLET_GROUP),
            ccb.disabled_applets,
        )

        experiment, base = self.build_and_prepare("two_nodes")
        ccb = RecordingCCB()
        experiment.ccb = ccb
        applets.create_two_node_applets_for(
            experiment, base,
            specs=applets.TWO_NODE_APPLET_SPECS,
            shared_specs=applets.SHARED_APPLET_SPECS,
        )
        self.assertEqual(ccb.disabled_groups, [])
        for node in ("Node1", "Node2"):
            self.assertIn(
                ("retention and loading AllSPCMs", node), ccb.disabled_applets
            )
        # and nothing else was closed
        self.assertEqual(
            {title for title, _ in ccb.disabled_applets},
            {"retention and loading AllSPCMs"},
        )

        groups = dict(ccb.created)
        self.assertEqual(
            groups["retention and loading AllSPCMs"],
            [applets.TWO_NODE_APPLET_GROUP, "GVS and Cycler"],
        )

    def test_agreeing_scalars_are_rederived_after_overrides(self):
        """The bare shared scalars must follow the per-node values.

        They used to be published only from prepare(), which broke two things
        at once, because the sequence reads the BARE name:

          * overriding or scanning an agreeing scalar had NO EFFECT -- the bare
            name stayed frozen at its pre-override value;
          * an ASYMMETRIC override did not raise, because prepare() compared
            the two still-equal dataset values, and the later setattr of the
            suffixed names left the bare name holding Node1's value. Node2 then
            silently ran with Node1's exposure time, which is the single
            failure _publish_agreeing_scalars exists to prevent.

        refresh_variable_dependent_state is the one hook both the override path
        and the per-scan-point path run through, so the re-derivation belongs
        there and this test pins it.
        """
        experiment, base = self.build_and_prepare("two_nodes")

        # Agreeing override: the bare name must pick the new value up.
        experiment.t_SPCM_first_shot_Node1 = 0.042
        experiment.t_SPCM_first_shot_Node2 = 0.042
        base.refresh_variable_dependent_state()
        self.assertAlmostEqual(experiment.t_SPCM_first_shot, 0.042)

        # Asymmetric override of a JOINT scalar: one gate cannot have two
        # durations, so this must raise and name the offender.
        experiment.t_SPCM_first_shot_Node2 = 0.043
        with self.assertRaises(RuntimeError) as caught:
            base.refresh_variable_dependent_state()
        self.assertIn("t_SPCM_first_shot", str(caught.exception))

    def test_per_node_flags_are_not_forced_to_agree(self):
        """Genuinely per-node names must stay OFF the agreeing list.

        PGC_and_RO_with_on_chip_beams selects A5/A6, a per-node AOM pair;
        do_PGC_after_loading gates a per-node stage whose duration already
        differs between the nodes (t_PGC_after_loading is 0.6 vs 1.0 ms); and
        t_FORT_drop is a per-node duration and the knob retention is measured
        against. Sharing any of them would hand Node2 Node1's setting, so the
        two-node sequence addresses all three explicitly and they must not
        appear as bare names.
        """
        per_node_only = (
            "PGC_and_RO_with_on_chip_beams",
            "do_PGC_after_loading",
            "t_FORT_drop",
        )
        for name in per_node_only:
            self.assertNotIn(
                name,
                BaseExperimentMasterSatellite.TWO_NODE_AGREEING_SCALARS,
                f"{name} is per node; sharing it would make Node2 follow "
                f"Node1 silently.",
            )

        experiment, base = self.build_and_prepare("two_nodes")

        # Diverging them is legal and must not raise.
        experiment.PGC_and_RO_with_on_chip_beams_Node1 = True
        experiment.PGC_and_RO_with_on_chip_beams_Node2 = False
        experiment.do_PGC_after_loading_Node1 = True
        experiment.do_PGC_after_loading_Node2 = False
        experiment.t_FORT_drop_Node1 = 1e-6
        experiment.t_FORT_drop_Node2 = 0.0
        base.refresh_variable_dependent_state()

        for name in per_node_only:
            self.assertFalse(hasattr(experiment, name), name)
            for node in ("Node1", "Node2"):
                self.assertTrue(hasattr(experiment, f"{name}_{node}"))

    def test_single_node_wiring_and_amplitudes_are_unchanged(self):
        """Regression lock on the three extensions.

        _compute_derived_amplitudes and _install_wiring_metadata both became
        per-node by routing through _presentation_name, which returns the BARE
        name in single_node mode. If that ever stops being true, single-node
        physics changes silently, so pin the exact values here.
        """
        experiment, _ = self.build_and_prepare("single_node", "Node1")

        self.assertEqual(experiment.coil_channels, [0, 1, 13, 14])
        self.assertEqual(experiment.UV_trig_channel, [8])
        self.assertEqual(experiment.Magnetometer_X_ch, 1)
        self.assertAlmostEqual(
            experiment.ampl_FORT_loading, dB_to_V(experiment.p_FORT_loading)
        )
        self.assertAlmostEqual(
            experiment.ampl_FORT_RO,
            experiment.ampl_FORT_loading * experiment.p_FORT_RO,
        )

        # Nothing suffixed leaks into the single-node presentation.
        for name in ("coil_channels_Node1", "ampl_FORT_loading_Node1",
                     "ampl_FORT_RO_Node1"):
            self.assertFalse(hasattr(experiment, name), name)

    def test_two_node_result_state_seeds_the_shared_unsuffixed_surface(self):
        """A two-node run's results describe the JOINT measurement.

        One atom per trap, read out through one gate of the four master-local
        SPCMs, so there is no per-node quantity to suffix -- which is why this
        is the same method single_node uses rather than a parallel copy.
        """
        experiment, base = self.build_and_prepare("two_nodes")
        experiment.dataset_writes.clear()

        base.initialize_result_state()

        written = {name for name, _, _ in experiment.dataset_writes}

        # MEASUREMENT RESULTS are joint and unsuffixed -- one atom per trap
        # read out through one gate of the four master-local SPCMs, so there is
        # no per-node quantity to suffix.
        #
        # HARDWARE MONITORS are the exception, and must be suffixed: each node
        # has its own FORT, its own pickoff and its own photodiode, and the
        # applets read FORT_MM_monitor_Node1 / _Node2 separately
        # (applets_master_satellite.py:321). Seeding only the bare name left
        # both suffixed names absent and killed a two-node run on hardware with
        # "Cannot mutate nonexistent dataset 'FORT_MM_monitor_Node1'", after
        # the feedback had already run.
        per_node_monitors = {
            f"{monitor}_{node}"
            for monitor in ("FORT_MM_monitor", "FORT_APD_monitor")
            for node in ("Node1", "Node2")
        }
        self.assertEqual(
            {name for name in written if name.endswith(("_Node1", "_Node2"))},
            per_node_monitors,
            "two-node measurement results are joint and must not be "
            "node-suffixed; the only suffixed names here are the per-node "
            "hardware monitors",
        )
        # Seeded empty and broadcast, which is what makes the first
        # append_to_dataset legal rather than a KeyError.
        seeded = {name: value for name, value, _ in experiment.dataset_writes}
        for name in sorted(per_node_monitors):
            self.assertEqual(
                seeded.get(name), [],
                f"{name} must be seeded as an empty broadcast list, or the "
                f"first append_to_dataset raises",
            )
        for required in (
            "AllSPCMs_RO1", "AllSPCMs_RO2", "AllSPCMs_atom_check_in_loading",
            "Atom_loading_time", "time_without_atom", "atom_loading_wall_clock",
            "photocount_bins", "AllSPCMs_alternating_RO_alice",
            "AllSPCMs_alternating_RO_bob",
        ):
            self.assertIn(required, written, required)

        # The host scalars and per-measurement buffers the kernels index. These
        # must EXIST before compilation, not merely before running.
        for required in ("AllSPCMs_RO1", "SPCM0_RO1", "SPCM1_OtherNode_RO2",
                        "measurement", "atom_loading_time", "in_health_check"):
            self.assertTrue(hasattr(experiment, required), required)
        self.assertEqual(
            len(experiment.AllSPCMs_RO1_list), experiment.n_measurements
        )

    def test_progress_resets_in_both_modes_but_magnetometers_do_not(self):
        """measurements_progress is run-global; the magnetometers are not.

        In two_nodes mode resolve_result_dataset_name is the identity, so both
        nodes would resolve to one unsuffixed Magnetometer_MOT_X and the last
        writer would win. The two-node loader samples no magnetometers anyway.
        """
        for mode, node in (("single_node", "Node1"), ("two_nodes", None)):
            with self.subTest(mode=mode):
                experiment, base = self.build_and_prepare(mode, node)
                experiment.dataset_writes.clear()
                base.reset_result_state_for_scan_point()
                written = [name for name, _, _ in experiment.dataset_writes]

                self.assertIn("measurements_progress", written)
                magnetometer_writes = [
                    name for name in written if name.startswith("Magnetometer_")
                ]
                if mode == "single_node":
                    self.assertTrue(magnetometer_writes)
                else:
                    self.assertEqual(magnetometer_writes, [])

    def test_every_sampler_gets_a_deterministic_gain_latched(self):
        """Sampler.init() never writes the PGIA, so gain must be asserted.

        Node2's samplers used to get no gain write at all and inherited
        whatever the previous run latched. Gain code 0 is x1, which is what
        every channel already ran at -- the standalone idiom
        set_gain_mu(channel, 8) wrote 0b1000 into the NEXT channel's two-bit
        field and the next iteration's mask cleared it, so the word actually
        sent was 0x0000. So this changes no measured voltage; it removes
        inherited state.
        """
        for mode, node, expected_samplers in (
            ("single_node", "Node1", 3),
            ("single_node", "Node2", 3),
            ("two_nodes", None, 6),
        ):
            with self.subTest(mode=mode, node=node):
                experiment, base = self.build_and_prepare(mode, node)
                experiment.log.clear()
                base.initialize_hardware()

                gain_writes = [
                    name for name, operation in experiment.log
                    if operation == "set_gain_mu"
                ]
                # Eight channels on every sampler in the run, and no more.
                self.assertEqual(len(gain_writes), 8 * expected_samplers)
                self.assertEqual(
                    len(set(gain_writes)), expected_samplers,
                    "every sampler in the active node(s) must be written",
                )

    def test_missing_dataset_identifies_owner_and_name(self):
        experiment = FakeExperiment()
        del experiment.datasets["f_FORT_Node2"]
        base = BaseExperimentMasterSatellite(
            experiment, "single_node", "Node2"
        )

        with self.assertRaisesRegex(
            RuntimeError,
            r"ExperimentVariables_master_satellite_Node2\.py: f_FORT_Node2",
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
            r"ExperimentVariables_master_satellite_Node2\.py: f_FORT_Node2",
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
