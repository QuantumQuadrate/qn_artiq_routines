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
    experiment_module.ms = 1e-3
    experiment_module.us = 1e-6
    artiq_module.experiment = experiment_module
    sys.modules["artiq"] = artiq_module
    sys.modules["artiq.experiment"] = experiment_module


from utilities.BaseExperiment_master_satellite import (  # noqa: E402
    BaseExperimentMasterSatellite,
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

        for node in ("Node1", "Node2"):
            config_path = (
                Path("utilities/config")
                / ("alice" if node == "Node1" else "bob")
                / "device_aliases.json"
            )
            with config_path.open() as config_file:
                logical = json.load(config_file)
            for defaults in logical["DDS_DEFAULTS"].values():
                setattr(self, f"{defaults['frequency']}_{node}", 80.0)
                setattr(self, f"{defaults['power']}_{node}", -10.0)

    def setattr_device(self, name):
        if name not in self.devices:
            raise KeyError(name)
        setattr(self, name, self.devices[name])


class BaseExperimentMasterSatelliteTests(unittest.TestCase):
    def build_and_prepare(self, mode, node=None):
        experiment = FakeExperiment()
        base = BaseExperimentMasterSatellite(experiment, mode, node)
        base.build()
        base.prepare()
        return experiment, base

    def test_single_node1_bindings(self):
        experiment, _ = self.build_and_prepare("single_node", "Node1")
        self.assertEqual(experiment.dds_FORT.name, "urukul0_ch0")
        self.assertEqual(experiment.sampler0.name, "sampler0")
        self.assertEqual(experiment.zotino0.name, "zotino0")
        self.assertEqual(experiment.ttl_SPCM0.name, "ttl0")

    def test_single_node2_uses_node2_devices_and_master_spcms(self):
        experiment, _ = self.build_and_prepare("single_node", "Node2")
        self.assertEqual(experiment.dds_FORT.name, "urukul4_ch3")
        self.assertEqual(experiment.sampler0.name, "sampler3")
        self.assertEqual(experiment.zotino0.name, "zotino1")
        self.assertEqual(experiment.ttl0.name, "ttl16")
        self.assertEqual(experiment.ttl_SPCM0.name, "ttl0")
        self.assertEqual(experiment.ttl_SPCM0_counter.name, "ttl0_counter")
        self.assertEqual(experiment.f_FORT, experiment.f_FORT_Node2)

    def test_two_node_bindings(self):
        experiment, _ = self.build_and_prepare("two_nodes")
        self.assertEqual(experiment.dds_FORT_Node1.name, "urukul0_ch0")
        self.assertEqual(experiment.dds_FORT_Node2.name, "urukul4_ch3")
        self.assertEqual(experiment.sampler0_Node1.name, "sampler0")
        self.assertEqual(experiment.sampler0_Node2.name, "sampler3")
        self.assertEqual(experiment.SPCM_H1.name, "ttl0")
        self.assertFalse(hasattr(experiment, "ttl_SPCM0"))

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
