"""
Example of how to make a "rolling" xy plot of two data channels
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph

from artiq.applets.simple import TitleApplet


class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.length_warning)
        self.mismatch = {'X values': False,
                         'Error bars': False,
                         'Fit values': False}

    def data_changed(self, data, mods, title):
        try:
            # should be a numpy array where each row is a different data channel
            y1 = data[self.args.y1][1]
            y2 = data[self.args.y2][1]

            # the data display will be rolling, only showing display_pts at a time
            pts = (data[self.args.pts][1])[0]
            print(pts)

        except KeyError:
            return
        x = data.get(self.args.x, (False, None))[1]
        if x is None:
            x = np.arange(len(y1))[-pts:]
        y1 = y1[-pts:]
        y2 = y2[-pts:]

        self.mismatch['y1'] = False
        self.mismatch['y2'] = False
        self.mismatch['x'] = False

        if len(y2) != len(y1):
            self.mismatch['y2'] = True
        if len(x) != len(y1):
            self.mismatch['x'] = True

        if any(self.mismatch.values()):
            if not self.timer.isActive():
                self.timer.start(1000)
            return
        else:
            self.timer.stop()

        self.clear()
        # self.plot(np.arange(len(y1)), y1, pen=(1,2), symbol="o")
        # self.plot(np.arange(len(y2)), y2, pen=(2,2), symbol="o")

        self.plot(x, y1, pen=(1, 2), symbol="o")
        self.plot(x, y2, pen=(2, 2), symbol="o")
        self.setTitle(title)

    def length_warning(self):
        self.clear()
        text = "⚠️ dataset lengths mismatch:\n"
        errors = ', '.join([k for k, v in self.mismatch.items() if v])
        text = ' '.join([errors, "should have the same length as Y values"])
        self.addItem(pyqtgraph.TextItem(text))


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("y1", "Y1 values")
    applet.add_dataset("y2", "Y2 values")
    applet.add_dataset("pts", "number of points to display")
    applet.add_dataset("x", "X values", required=False)
    # applet.add_dataset("error", "Error bars for each X value", required=False)
    # applet.add_dataset("fit", "Fit values for each X value", required=False)
    applet.run()

if __name__ == "__main__":
    main()
