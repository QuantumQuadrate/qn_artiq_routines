"""Single-node compatibility GVS for master-satellite hardware."""

from artiq.experiment import EnumerationValue, EnvExperiment

from GeneralVariableScan_master_satellite import (
    _GeneralVariableScanMasterSatelliteMixin,
)


class GeneralVariableScan_master_satellite_single_node(
    _GeneralVariableScanMasterSatelliteMixin, EnvExperiment
):
    """GeneralVariableScan_master_satellite_single_node

    Run existing single-node experiment functions on selected Node1 or Node2
    hardware through the master-satellite stack.
    """

    EXPERIMENT_MODE = "single_node"

    def build(self):
        self.setattr_argument(
            "selected_node", EnumerationValue(self.VALID_NODES)
        )
        self._build_master_satellite_scan()
