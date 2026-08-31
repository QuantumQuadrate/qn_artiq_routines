"""Node2-only manual AOM/DDS and coil control."""

from artiq.experiment import EnvExperiment

from AOMsCoils_master_satellite import _AOMsCoilsMasterSatelliteMixin


class AOMsCoils_master_satellite_Node2(
    _AOMsCoilsMasterSatelliteMixin, EnvExperiment
):
    """AOMsCoils_master_satellite_Node2

    Manually control only Node2 hardware without changing Node1 outputs.
    """

    NODE = "Node2"
    COIL_CHANNELS = (0, 1, 2, 3)

    def build(self):
        self._build_node_manual_controls()
