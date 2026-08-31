"""Native two-node GVS for master-satellite hardware."""

from artiq.experiment import EnvExperiment

from GeneralVariableScan_master_satellite_mixin import (
    _GeneralVariableScanMasterSatelliteMixin,
)


class GeneralVariableScan_master_satellite_two_nodes(
    _GeneralVariableScanMasterSatelliteMixin, EnvExperiment
):
    """GeneralVariableScan_master_satellite_two_nodes

    Run native two-node functions on both master-satellite nodes.
    """

    EXPERIMENT_MODE = "two_nodes"

    def build(self):
        self._build_master_satellite_scan()
