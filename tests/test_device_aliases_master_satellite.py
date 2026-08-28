import unittest

from utilities.DeviceAliases_master_satellite import (
    DeviceAliasesMasterSatellite,
    DeviceResolutionError,
)


class FakeDevice:
    def __init__(self, name):
        self.name = name


class FakeCore:
    def __init__(self):
        self.reset_calls = 0

    def reset(self):
        self.reset_calls += 1


class FakeExperiment:
    def __init__(self):
        self.core = FakeCore()
        self.devices = {
            "urukul0_ch0": FakeDevice("urukul0_ch0"),
            "urukul4_ch3": FakeDevice("urukul4_ch3"),
        }
        self.f_FORT_Node1 = 80.0
        self.p_FORT_loading_Node1 = -10.0
        self.f_FORT_Node2 = 81.0
        self.p_FORT_loading_Node2 = -11.0

    def setattr_device(self, name):
        if name not in self.devices:
            raise KeyError(name)
        setattr(self, name, self.devices[name])


class DeviceAliasesMasterSatelliteTests(unittest.TestCase):
    def make_resolver(self, node):
        experiment = FakeExperiment()
        resolver = DeviceAliasesMasterSatellite(
            experiment,
            node,
            default_attribute_name_resolver=lambda name: f"{name}_{node}",
        )
        return experiment, resolver

    def test_node1_physical_resolution(self):
        _, resolver = self.make_resolver("Node1")
        self.assertEqual(
            resolver.resolve_physical_name("urukul0_ch0"), "urukul0_ch0"
        )

    def test_node2_physical_resolution(self):
        _, resolver = self.make_resolver("Node2")
        self.assertEqual(
            resolver.resolve_physical_name("urukul1_ch3"), "urukul4_ch3"
        )

    def test_node2_fort_two_stage_resolution(self):
        _, resolver = self.make_resolver("Node2")
        self.assertEqual(
            resolver.resolve_dds_physical_name("dds_FORT"),
            ("urukul1_ch3", "urukul4_ch3"),
        )

    def test_presentation_does_not_change_node2_fort_device(self):
        compatibility_experiment, compatibility_resolver = self.make_resolver(
            "Node2"
        )
        compatibility_device = compatibility_resolver.bind_dds(
            "dds_FORT", "dds_FORT"
        )

        two_node_experiment, two_node_resolver = self.make_resolver("Node2")
        two_node_device = two_node_resolver.bind_dds(
            "dds_FORT", "dds_FORT_Node2"
        )

        self.assertEqual(compatibility_device.name, "urukul4_ch3")
        self.assertEqual(two_node_device.name, "urukul4_ch3")
        self.assertIs(compatibility_experiment.dds_FORT, compatibility_device)
        self.assertIs(two_node_experiment.dds_FORT_Node2, two_node_device)
        self.assertEqual(
            compatibility_resolver.dds_bindings[0]["frequency_attribute"],
            "f_FORT_Node2",
        )

    def test_missing_mapping_fails_loudly(self):
        _, resolver = self.make_resolver("Node2")
        with self.assertRaisesRegex(
            DeviceResolutionError, "No master-satellite physical mapping"
        ):
            resolver.resolve_physical_name("missing_device")

    def test_resolver_never_resets_core(self):
        experiment, resolver = self.make_resolver("Node2")
        resolver.bind_dds("dds_FORT", "dds_FORT")
        self.assertEqual(experiment.core.reset_calls, 0)


if __name__ == "__main__":
    unittest.main()
