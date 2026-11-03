"""
Show the master-satellite GVS microwave health check status as text.

The GVS (run_health_check on) writes it to the broadcast dataset
health_check_status: what the check is doing now, and per node and transition
how many times it had to resubmit the optimizer this run, after a failed
check and after a crashed chain counted apart. A count that keeps growing
means the lasers need relocking.

applet command:
python "C:\\...\\qn_artiq_routines\\applets\\health_check_status.py" health_check_status
"""

from PyQt5 import QtGui, QtWidgets

from artiq.applets.simple import SimpleApplet


class HealthCheckStatus(QtWidgets.QTextEdit):
    def __init__(self, args):
        QtWidgets.QTextEdit.__init__(self)
        self.dataset_name = args.status
        self.setReadOnly(True)
        self.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        # fixed width, so the per-transition counts line up
        self.setFont(QtGui.QFontDatabase.systemFont(
            QtGui.QFontDatabase.FixedFont))

    def data_changed(self, data, mods):
        try:
            text = str(data[self.dataset_name][1])
        except KeyError:
            text = "no health check status yet"
        self.setPlainText(text)


def main():
    applet = SimpleApplet(HealthCheckStatus)
    applet.add_dataset("status", "the GVS health check status text "
                                 "(health_check_status)")
    applet.run()


if __name__ == "__main__":
    main()
