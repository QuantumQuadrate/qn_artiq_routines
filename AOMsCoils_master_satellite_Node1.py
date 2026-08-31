"""Node1-only manual AOM/DDS and coil control."""

from artiq.experiment import EnvExperiment

from AOMsCoils_master_satellite_mixin import _AOMsCoilsMasterSatelliteMixin


class AOMsCoils_master_satellite_Node1(
    _AOMsCoilsMasterSatelliteMixin, EnvExperiment
):
    """AOMsCoils_master_satellite_Node1

    Manually control only Node1 hardware without changing Node2 outputs.
    """

    NODE = "Node1"
    COIL_CHANNELS = (0, 1, 13, 14)

    def build(self):
        self._build_node_manual_controls()
