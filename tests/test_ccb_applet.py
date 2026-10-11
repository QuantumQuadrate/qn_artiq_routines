"""test_ccb_applet.py

Smallest possible check that an experiment can make an applet appear in the
dashboard on its own, via the CCB (the mechanism applets_master_satellite.py
is built on).

Touches NO hardware: it binds no core device, drives nothing, and only asks
the dashboard to display one number. Safe to submit at any time.

Expected result: an applet titled "CCB test - <dataset>" appears in the
dashboard under the group given below, showing that dataset's value.

If nothing appears, the usual causes are:
  * the experiment was run without a dashboard connected (CCB is a
    dashboard-side service -- artiq_run cannot create applets);
  * the applet dock's CCB policy is set to "Ignore requests" (it should be
    "Create applets"; the saved configuration had ccbp_global='create');
  * the dataset does not exist, in which case the applet appears but stays
    blank.
"""

from artiq.experiment import *


class test_ccb_applet(EnvExperiment):
    """test_ccb_applet"""

    def build(self):
        # "ccb" is a virtual device provided by the worker; no hardware.
        self.setattr_device("ccb")
        self.setattr_argument(
            "dataset_to_show",
            StringValue("n_measurements"),
            tooltip="Any existing scalar dataset.",
        )
        self.setattr_argument(
            "applet_group",
            StringValue("CCB test"),
            tooltip="Dashboard group to create the applet under.",
        )
        self.setattr_argument(
            "remove_instead",
            BooleanValue(False),
            tooltip="Tick to delete the whole group again instead of creating it.",
        )

    def run(self):
        group = str(self.applet_group)
        if self.remove_instead:
            self.ccb.issue("disable_applet_group", group)
            print(f"[ccb] requested removal of applet group {group!r}")
            return

        dataset = str(self.dataset_to_show)
        title = f"CCB test - {dataset}"
        command = "${artiq_applet}big_number " + dataset
        self.ccb.issue("create_applet", title, command, group=group)
        print(f"[ccb] requested applet {title!r} in group {group!r}")
        print(f"[ccb] command: {command}")
        print("[ccb] if no applet appeared, see this file's docstring.")
