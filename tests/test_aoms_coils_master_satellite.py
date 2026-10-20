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
        "kernel": _identity_decorator,
        "rpc": _identity_decorator,
        "delay": lambda duration: None,
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
        "TFloat": float,
        "TStr": str,
        "TInt32": int,
        "TInt64": int,
        "TBool": bool,
    }.items():
        if not hasattr(_experiment_stub, name):
            setattr(_experiment_stub, name, value)


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


def make_shared_hardware(log):
    devices = {"core": FakeCore("core", log)}
    with Path("utilities/config/master_satellite/device_aliases.json").open() as file:
        mapping = json.load(file)
    for node_mapping in mapping.values():
        for unified_name in node_mapping.values():
            devices.setdefault(unified_name, FakeDevice(unified_name, log))
    devices["core_dma"] = FakeDevice("core_dma", log)
    devices["scheduler"] = FakeDevice("scheduler", log)
    devices["k10cr1_ndsp"] = FakeRotator("k10cr1_ndsp", log)
    return devices


class ManualExperimentEnvironment:
    def initialize_environment(self, devices, datasets, log):
        self.devices = devices
        self.datasets = datasets
        self.log = log
        self.dataset_reads = []
        self.dataset_writes = []
        self.dataset_appends = []
        self.broadcast_sets = {}
    def setattr_argument(self, name, value, *args, **kwargs): setattr(self, name, value)
    def setattr_device(self, name): setattr(self, name, self.devices[name])
    def get_device(self, name): return self.devices[name]
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


class FakeStabilizerChannel:
    def __init__(self, dataset, dB_dataset):
        self.dataset = dataset
        self.dB_dataset = dB_dataset
        self.dB_history_dataset = dB_dataset + "_history"


class FakeStabilizer:
    def __init__(self, experiment, dds_names, iterations, averages, **kwargs):
        self.exp = experiment
        self.dds_names = dds_names
        self.iterations = iterations
        self.averages = averages
        self.kwargs = kwargs
        self.run_calls = 0
        self.monitor_calls = 0
        self.all_channels = [
            FakeStabilizerChannel("MOT1_monitor", "p_AOM_A1"),
            FakeStabilizerChannel("MOT2_monitor", "p_AOM_A2"),
        ]

    def run(self): self.run_calls += 1
    def monitor(self): self.monitor_calls += 1


# The harness sits after the experiment classes in the MRO so the dataset
# redirection layer resolves names before the harness records them.
class Node1ManualExperiment(AOMsCoils_master_satellite_Node1, ManualExperimentEnvironment):
    _stabilizer_factory = FakeStabilizer


class Node2ManualExperiment(AOMsCoils_master_satellite_Node2, ManualExperimentEnvironment):
    _stabilizer_factory = FakeStabilizer


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
        # Both nodes take their coil channels from Base's wiring table.
        self.assertEqual(node1.coil_channels, [0, 1, 13, 14])
        self.assertEqual(node2.coil_channels, [0, 1, 2, 3])

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
        for node, experiment_class, zotino_name, channels in (
            ("Node1", Node1ManualExperiment, "zotino0", (0, 1, 13, 14)),
            ("Node2", Node2ManualExperiment, "zotino1", (0, 1, 2, 3)),
        ):
            with self.subTest(node=node):
                experiment = self.make_experiment(
                    experiment_class, {"disable_coils": True}
                )
                experiment.apply_manual_state()
                zotino = self.devices[zotino_name]
                self.assertEqual(
                    [zotino.dac_values[channel] for channel in channels],
                    [0.0, 0.0, 0.0, 0.0],
                )

    def test_mot_coil_voltages_come_from_persistent_calibration(self):
        for node, experiment_class, zotino_name, channels in (
            ("Node1", Node1ManualExperiment, "zotino0", (0, 1, 13, 14)),
            ("Node2", Node2ManualExperiment, "zotino1", (0, 1, 2, 3)),
        ):
            with self.subTest(node=node):
                experiment = self.make_experiment(experiment_class)
                self.assertFalse(experiment.disable_coils)
                experiment.apply_manual_state()
                zotino = self.devices[zotino_name]
                self.assertEqual(
                    [zotino.dac_values[channel] for channel in channels],
                    [
                        self.datasets[f"AZ_bottom_volts_MOT_{node}"],
                        self.datasets[f"AZ_top_volts_MOT_{node}"],
                        self.datasets[f"AX_volts_MOT_{node}"],
                        self.datasets[f"AY_volts_MOT_{node}"],
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

    def test_each_dds_box_writes_its_authoritative_binding(self):
        """Each node's file names its DDS tones directly. Tie them back to
        Base's DDS binding table, the single source the old mixin read them
        from, so a wrong name in either file fails here.
        """
        from utilities.conversions import dB_to_V_kernel

        boxes = {
            "FORT_AOM_ON": "dds_FORT",
            "cooling_DP_ON": "dds_cooling_DP",
            "AOM_A1_ON": "dds_AOM_A1",
            "AOM_A2_ON": "dds_AOM_A2",
            "AOM_A3_ON": "dds_AOM_A3",
            "AOM_A4_ON": "dds_AOM_A4",
            "AOM_A5_ON": "dds_AOM_A5",
            "AOM_A6_ON": "dds_AOM_A6",
            "D1_pumping_DP_ON": "dds_D1_pumping_DP",
            "GRIN1and2_ON": "GRIN1and2_dds",
        }
        for node, experiment_class in (
            ("Node1", Node1ManualExperiment),
            ("Node2", Node2ManualExperiment),
        ):
            with self.subTest(node=node):
                experiment = self.make_experiment(
                    experiment_class, {box: True for box in boxes}
                )
                bindings = {
                    item["logical_alias"]: item
                    for item in experiment.base.node_resolvers[node].dds_bindings
                }
                # Not apply_manual_state(): initialize_hardware() would set
                # every DDS to these same values first.
                start = len(self.log)
                experiment.turn_on_AOMs()
                for box, alias in boxes.items():
                    if node == "Node1" and alias == "dds_D1_pumping_DP":
                        continue  # Node1's D1 runs from an external source
                    device = getattr(experiment, alias)
                    first_set = next(
                        value for name, operation, value in self.log[start:]
                        if name == device.name and operation == "set"
                    )
                    binding = bindings[alias]
                    self.assertEqual(first_set, {
                        "frequency": self.datasets[binding["frequency_attribute"]],
                        "amplitude": dB_to_V_kernel(
                            self.datasets[binding["power_attribute"]]
                        ),
                    }, box)

    def test_microwave_defaults_off_and_no_calibration_is_persisted(self):
        for experiment_class in (Node1ManualExperiment, Node2ManualExperiment):
            with self.subTest(node=experiment_class.NODE):
                experiment = self.make_experiment(experiment_class, {
                    "microwave_dds_ON": True,
                    "MW_RF_dds_ON": True,
                    "yes_Im_sure_I_want_MW_or_RF_dds_ON": False,
                })
                writes_after_prepare = len(experiment.dataset_writes)
                experiment.apply_manual_state()
                self.assertFalse(experiment.dds_microwaves.state)
                self.assertTrue(experiment.ttl_microwave_switch.state)
                self.assertFalse(experiment.dds_MW_RF.state)
                self.assertEqual(
                    len(experiment.dataset_writes), writes_after_prepare
                )

    def test_node2_combined_grin_modes_drive_dds_and_gate_together(self):
        from utilities.conversions import dB_to_V_kernel

        node2 = self.make_experiment(Node2ManualExperiment, {
            "Node2_GRIN1_AOM_ON": True, "Node2_GRIN2_AOM_ON": True,
        })
        node2.apply_manual_state()
        grin1and2_dds = self.devices["urukul5_ch2"]
        d1_dds = self.devices["urukul5_ch0"]
        self.assertTrue(grin1and2_dds.state)
        self.assertFalse(self.devices["ttl30"].state)  # GRIN1 gate open
        self.assertTrue(d1_dds.state)
        self.assertFalse(self.devices["ttl29"].state)  # GRIN2 gate open
        self.assertEqual(d1_dds.last_set, {
            "frequency": self.datasets["f_GRIN2_excitation_Node2"],
            "amplitude": dB_to_V_kernel(
                self.datasets["p_GRIN2_excitation_Node2"]
            ),
        })

    def test_node2_combined_modes_off_override_individual_grin_controls(self):
        node2 = self.make_experiment(Node2ManualExperiment, {
            "GRIN1and2_ON": True, "D1_pumping_DP_ON": True,
        })
        node2.apply_manual_state()
        # The unticked combined modes win, as in the standalone ordering:
        # both DDSes end OFF while the dedicated D1 gate stays open.
        self.assertFalse(self.devices["urukul5_ch2"].state)
        self.assertFalse(self.devices["urukul5_ch0"].state)
        self.assertFalse(self.devices["ttl28"].state)  # ttl_D1_pumping_switch open
        self.assertTrue(self.devices["ttl30"].state)   # GRIN1 gate closed
        self.assertTrue(self.devices["ttl29"].state)   # GRIN2 gate closed

    def test_node1_has_no_combined_grin_modes_and_keeps_grin_control(self):
        node1 = self.make_experiment(
            Node1ManualExperiment, {"GRIN1and2_ON": True}
        )
        self.assertFalse(hasattr(node1, "Node2_GRIN1_AOM_ON"))
        node1.apply_manual_state()
        self.assertTrue(self.devices["urukul0_ch3"].state)

    def test_node1_d1_is_switched_by_its_external_source_ttl(self):
        """Node1's D1 AOM runs from an external source since 2026-10-19.

        D1_pumping_DP_ON opens ttl_D1_pumping_switch (Node1 ttl15, active
        low: off() turns D1 on) and, as before, the GRIN1 gate D1 passes
        through. The Urukul channel that used to drive D1 (urukul2_ch2)
        stays off either way.
        """
        node1 = self.make_experiment(
            Node1ManualExperiment, {"D1_pumping_DP_ON": True}
        )
        node1.apply_manual_state()
        self.assertFalse(self.devices["ttl15"].state)        # D1 on
        self.assertFalse(self.devices["ttl14"].state)        # GRIN1 open
        self.assertFalse(self.devices["urukul2_ch2"].state)  # old D1 DDS off
        self.assertFalse(self.devices["urukul0_ch3"].state)  # GRIN tone off

        node1 = self.make_experiment(Node1ManualExperiment)
        node1.apply_manual_state()
        self.assertTrue(self.devices["ttl15"].state)         # D1 blocked
        self.assertTrue(self.devices["ttl14"].state)         # GRIN1 blocked
        self.assertFalse(self.devices["urukul2_ch2"].state)

        # The GRIN1 gate alone does not open D1.
        node1 = self.make_experiment(
            Node1ManualExperiment, {"GRIN1_AOM_switch_ON": True}
        )
        node1.apply_manual_state()
        self.assertTrue(self.devices["ttl15"].state)         # D1 blocked
        self.assertFalse(self.devices["ttl14"].state)        # GRIN1 open

    def test_laser_stabilizer_prepared_with_suffixed_persistence(self):
        node1 = self.make_experiment(Node1ManualExperiment)
        stabilizer = node1.laser_stabilizer
        self.assertIsInstance(stabilizer, FakeStabilizer)
        self.assertEqual(
            stabilizer.dds_names,
            eval(self.datasets["fast_feedback_dds_list_Node1"]),
        )
        self.assertEqual(
            node1.base.feedback_dataset_map["p_AOM_A1"], "p_AOM_A1_Node1"
        )
        self.assertEqual(
            node1.base.feedback_dataset_map["MOT1_monitor"],
            "MOT1_monitor_Node1",
        )
        self.assertEqual(
            node1.base.feedback_dataset_map["feedbackchannels"],
            "feedbackchannels_Node1",
        )
        self.assertIn("p_AOM_A1_Node1", node1.dataset_reads)
        self.assertEqual(
            node1.broadcast_sets["feedbackchannels_Node1"],
            ["p_AOM_A1_Node1", "p_AOM_A2_Node1"],
        )
        self.assertEqual(
            node1.broadcast_sets["p_AOM_A1_history_Node1"],
            [float(self.datasets["p_AOM_A1_Node1"])],
        )

        node2 = self.make_experiment(Node2ManualExperiment)
        self.assertEqual(node2.which_node, "bob")
        self.assertEqual(
            node2.base.feedback_dataset_map["p_AOM_A1"], "p_AOM_A1_Node2"
        )

    def test_run_feedback_gates_on_all_stabilized_aoms(self):
        all_on = {name: True for name in (
            "AOM_A1_ON", "AOM_A2_ON", "AOM_A3_ON",
            "AOM_A4_ON", "AOM_A5_ON", "AOM_A6_ON", "cooling_DP_ON",
        )}

        for experiment_class in (Node1ManualExperiment, Node2ManualExperiment):
            with self.subTest(node=experiment_class.NODE):
                feedback = self.make_experiment(
                    experiment_class, {**all_on, "run_laser_feedback": True}
                )
                feedback.run()
                self.assertEqual(feedback.laser_stabilizer.run_calls, 1)
                self.assertEqual(feedback.laser_stabilizer.monitor_calls, 0)

                monitor_only = self.make_experiment(experiment_class, all_on)
                monitor_only.run()
                self.assertEqual(monitor_only.laser_stabilizer.run_calls, 0)
                self.assertEqual(monitor_only.laser_stabilizer.monitor_calls, 1)

                partial = self.make_experiment(
                    experiment_class, {"run_laser_feedback": True}
                )
                partial.run()
                self.assertEqual(partial.laser_stabilizer.run_calls, 0)
                self.assertEqual(partial.laser_stabilizer.monitor_calls, 0)

                # Every one of the seven boxes is required.
                for missing in all_on:
                    almost = self.make_experiment(experiment_class, {
                        **all_on, missing: False, "run_laser_feedback": True,
                    })
                    almost.run()
                    self.assertEqual(
                        (almost.laser_stabilizer.run_calls,
                         almost.laser_stabilizer.monitor_calls),
                        (0, 0),
                        missing,
                    )

    def test_run_feedback_keeps_node2_combined_grin_state(self):
        """Feedback must end Node2 the way apply_manual_state left it.

        run_feedback re-applies the individual DDS and gate flags after the
        stabilizer; the Node2 combined GRIN modes have to follow, or the
        individual flags win: the GRIN DDSes come back on and GRIN1's gate
        opens, where the standalone utility ends them off.
        """
        all_on = {name: True for name in (
            "AOM_A1_ON", "AOM_A2_ON", "AOM_A3_ON",
            "AOM_A4_ON", "AOM_A5_ON", "AOM_A6_ON", "cooling_DP_ON",
        )}
        individual = {
            "GRIN1and2_ON": True, "D1_pumping_DP_ON": True,
            "GRIN1_AOM_switch_ON": True,
        }
        for feedback in (True, False):
            with self.subTest(run_laser_feedback=feedback, combined=False):
                node2 = self.make_experiment(Node2ManualExperiment, {
                    **all_on, **individual, "run_laser_feedback": feedback,
                })
                node2.run()
                self.assertFalse(self.devices["urukul5_ch2"].state)
                self.assertFalse(self.devices["urukul5_ch0"].state)
                self.assertTrue(self.devices["ttl30"].state)  # GRIN1 blocked
                self.assertTrue(self.devices["ttl29"].state)  # GRIN2 blocked
            with self.subTest(run_laser_feedback=feedback, combined=True):
                node2 = self.make_experiment(Node2ManualExperiment, {
                    **all_on, "run_laser_feedback": feedback,
                    "Node2_GRIN1_AOM_ON": True, "Node2_GRIN2_AOM_ON": True,
                })
                node2.run()
                self.assertTrue(self.devices["urukul5_ch2"].state)
                self.assertTrue(self.devices["urukul5_ch0"].state)
                self.assertFalse(self.devices["ttl30"].state)  # GRIN1 open
                self.assertFalse(self.devices["ttl29"].state)  # GRIN2 open

    def test_k10cr1_proxy_is_lazy_and_node_aware(self):
        from utilities.BaseExperiment_master_satellite import (
            _NodeAxisK10CR1Proxy,
        )

        node1 = self.make_experiment(Node1ManualExperiment)
        self.assertFalse(node1._k10cr1_requested)
        self.assertIsInstance(node1.k10cr1_ndsp, _NodeAxisK10CR1Proxy)
        # Lazy: publishing the proxy performs no device access.
        self.assertIsNone(node1.k10cr1_ndsp._client)

        node2 = self.make_experiment(
            Node2ManualExperiment, {"go_to_home_852QWP": True}
        )
        self.assertTrue(node2._k10cr1_requested)
        # Bare axis names resolve to the selected node; suffixed pass through.
        node2.k10cr1_ndsp.home("852_QWP")
        node2.k10cr1_ndsp.home("852_HWP_Node2")
        self.assertEqual(
            [entry for entry in self.log if entry[1] == "home"],
            [
                ("k10cr1_ndsp", "home", "852_QWP_Node2"),
                ("k10cr1_ndsp", "home", "852_HWP_Node2"),
            ],
        )

    def test_legacy_bare_axis_names_reach_the_selected_node(self):
        from subroutines.k10cr1_functions import move_to_target_deg

        node1 = self.make_experiment(Node1ManualExperiment)
        move_to_target_deg(node1, name="852_HWP", target_deg=5.0)
        rotator = self.devices["k10cr1_ndsp"]
        self.assertEqual(
            rotator.positions["852_HWP_Node1"],
            int(5.0 * self.datasets["deg_to_pos_Node1"]),
        )

    def test_k10cr1_missing_entry_fails_lazily_with_guidance(self):
        del self.devices["k10cr1_ndsp"]
        node1 = self.make_experiment(Node1ManualExperiment)  # prepare works
        with self.assertRaisesRegex(RuntimeError, "k10cr1_ndsp"):
            node1.k10cr1_ndsp.get_position("852_HWP")

    def test_852_target_moves_are_gated_on_852_booleans(self):
        # Each node's file is self-contained, so its own module holds the
        # K10CR1 calls, with that node's axis names written out.
        import AOMsCoils_master_satellite_Node1 as node1_module
        import AOMsCoils_master_satellite_Node2 as node2_module
        from unittest import mock

        for node, module, experiment_class in (
            ("Node1", node1_module, Node1ManualExperiment),
            ("Node2", node2_module, Node2ManualExperiment),
        ):
            with self.subTest(node=node):
                calls = []
                with mock.patch.object(
                    module, "go_to_home",
                    side_effect=lambda experiment, name: calls.append(
                        ("home", name)
                    ),
                ), mock.patch.object(
                    module, "move_to_target_deg",
                    side_effect=lambda experiment, name, target_deg: calls.append(
                        ("target", name, target_deg)
                    ),
                ), mock.patch.object(
                    module, "move_by_deg",
                    side_effect=lambda experiment, name, target_deg: calls.append(
                        ("by", name, target_deg)
                    ),
                ):
                    experiment = self.make_experiment(experiment_class, {
                        "go_to_target_780HWP": True,
                        "go_to_target_852QWP": True,
                    })
                    experiment.run()

                self.assertEqual(calls, [
                    ("target", f"780_HWP_{node}",
                     self.datasets[f"target_780_HWP_{node}"]),
                    ("target", f"852_QWP_{node}",
                     self.datasets[f"target_852_QWP_{node}"]),
                ])

    def test_every_k10cr1_box_moves_only_its_own_nodes_axis(self):
        """Each node's file writes its waveplate axis names out by hand, and
        the K10CR1 proxy passes an already-suffixed name through unchanged,
        so a stray suffix would rotate the other node's waveplate. Tick every
        K10CR1 box and check the axis of every call, in order.
        """
        import AOMsCoils_master_satellite_Node1 as node1_module
        import AOMsCoils_master_satellite_Node2 as node2_module
        from unittest import mock

        boxes = (
            "go_to_home_780HWP", "go_to_home_780QWP",
            "go_to_target_780HWP", "go_to_target_780QWP",
            "move_780HWP_by", "move_780QWP_by",
            "go_to_home_852HWP", "go_to_home_852QWP",
            "go_to_target_852HWP", "go_to_target_852QWP",
            "move_852HWP_by", "move_852QWP_by",
            "go_to_optimized_852_settings",
        )
        for node, module, experiment_class in (
            ("Node1", node1_module, Node1ManualExperiment),
            ("Node2", node2_module, Node2ManualExperiment),
        ):
            with self.subTest(node=node):
                axes = []
                with mock.patch.object(
                    module, "go_to_home",
                    side_effect=lambda experiment, name: axes.append(name),
                ), mock.patch.object(
                    module, "move_to_target_deg",
                    side_effect=lambda experiment, name, target_deg: axes.append(name),
                ), mock.patch.object(
                    module, "move_by_deg",
                    side_effect=lambda experiment, name, target_deg: axes.append(name),
                ):
                    experiment = self.make_experiment(
                        experiment_class, {box: True for box in boxes}
                    )
                    experiment.k10cr1_operations()

                # home, target, move-by for 780 then 852, then the optimized
                # 852 settings: HWP then QWP each time.
                self.assertEqual(axes, [
                    f"{beam}_{plate}_{node}"
                    for beam in ("780", "780", "780", "852", "852", "852", "852")
                    for plate in ("HWP", "QWP")
                ])

    def test_node_files_are_self_contained_experiments(self):
        """Each node's AOMsCoils is complete on its own (2026-10-19): one
        public experiment per file, and neither file imports the other or any
        shared AOMsCoils module, so a change for one node cannot reach the
        other. The mixin the two used to share is gone.
        """
        self.assertFalse(Path("AOMsCoils_master_satellite_mixin.py").exists())
        for node in ("Node1", "Node2"):
            filename = f"AOMsCoils_master_satellite_{node}.py"
            tree = ast.parse(Path(filename).read_text())
            public = {item.name for item in tree.body if isinstance(item, ast.ClassDef)
                      and not item.name.startswith("_")
                      and any(isinstance(base, ast.Name) and base.id == "EnvExperiment"
                              for base in item.bases)}
            self.assertEqual(public, {f"AOMsCoils_master_satellite_{node}"})
            imported = {item.module for item in ast.walk(tree)
                        if isinstance(item, ast.ImportFrom)}
            imported |= {alias.name for item in ast.walk(tree)
                         if isinstance(item, ast.Import) for alias in item.names}
            self.assertEqual(
                {name for name in imported if "AOMsCoils" in name}, set(),
                filename,
            )


if __name__ == "__main__":
    unittest.main()
