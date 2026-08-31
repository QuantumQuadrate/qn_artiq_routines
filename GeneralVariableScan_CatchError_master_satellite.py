"""RTIO-underflow retry wrapper for the master-satellite GVS stack."""

import logging
import time

from artiq.experiment import BooleanValue, NumberValue
from artiq.coredevice.exceptions import RTIOUnderflow

from GeneralVariableScan_master_satellite import (
    GeneralVariableScanMasterSatellite,
)


class GeneralVariableScan_CatchError_master_satellite(
    GeneralVariableScanMasterSatellite
):
    """GeneralVariableScan_CatchError_master_satellite

    Retry an individual master-satellite scan point after an underflow.

    Function discovery, authoritative variable resolution, compatibility
    projection, device binding, result-state handling, and hardware
    initialization are inherited from GeneralVariableScanMasterSatellite.
    Only RTIOUnderflow is caught; DRTIO/link and device errors remain fatal so
    they are not hidden during master-satellite validation.
    """

    def build(self):
        super().build()
        self.setattr_argument(
            "underflow_max_retries",
            NumberValue(10, ndecimals=0, step=1, type="int"),
            "Catch Underflow",
        )
        self.setattr_argument(
            "underflow_backoff_ms",
            NumberValue(200.0, step=0.5),
            "Catch Underflow",
        )
        self.setattr_argument(
            "skip_only_that_iteration_if_exhausted",
            BooleanValue(True),
            "Catch Underflow",
        )

    def _report_underflow(self, message):
        logging.warning(message)
        print(message)

    def _underflow_backoff(self):
        """Pause on the host; the next point attempt resets core via Base."""
        time.sleep(max(0.0, float(self.underflow_backoff_ms)) / 1000.0)

    def run(self):
        self._initialize_run_state()

        iteration = 0
        self.set_dataset("iteration", iteration, broadcast=True)
        for variable1_value in self.scan_sequence1:
            for variable2_value in self.scan_sequence2:
                retries = 0
                while True:
                    try:
                        self._run_scan_point(
                            variable1_value, variable2_value, iteration
                        )
                        break
                    except RTIOUnderflow as error:
                        retries += 1
                        self._report_underflow(
                            f"RTIO underflow at iteration {iteration}, "
                            f"retry {retries}/{int(self.underflow_max_retries)}: "
                            f"{error}"
                        )
                        if retries >= int(self.underflow_max_retries):
                            message = (
                                f"RTIO underflow at iteration {iteration} "
                                f"exceeded max retries "
                                f"({int(self.underflow_max_retries)})."
                            )
                            self._report_underflow(message)
                            if not self.skip_only_that_iteration_if_exhausted:
                                raise
                            break
                        self._underflow_backoff()
                iteration += 1

        print(
            "**************** General Variable Scan master-satellite "
            "CatchError DONE ****************"
        )
