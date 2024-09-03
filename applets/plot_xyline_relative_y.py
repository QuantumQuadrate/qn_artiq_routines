"""
same as plot_xy but there's a line and we have the option to divide by a specified y0

todo: properly dividing fit and error has not been dealt with
"""
#!/usr/bin/env python3

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
            y = np.array(data[self.args.y][1])
        except KeyError:
            return
        pts = data.get(self.args.pts, (False, None))[1]
        x = data.get(self.args.x, (False, None))[1]
        y0 = data.get(self.args.y0, (False, None))[1]

        if x is None:
            x = np.arange(len(y))
        if pts is not None:
            x = x[-pts:]
            y = y[-pts:]
        error = data.get(self.args.error, (False, None))[1]
        fit = data.get(self.args.fit, (False, None))[1]
        fitx = np.array(data.get(self.args.fitx, (False, None))[1])

        if y0 is not None:
            y /= y0
            # if error is not None:
            #     error /= y0
            # if fit is not None:
            #     fit /= y0

        marker_point = data.get(self.args.marker, (False, None))[1]

        if not len(y) or len(y) != len(x):
            self.mismatch['X values'] = True
        else:
            self.mismatch['X values'] = False
        if error is not None and hasattr(error, "__len__"):
            if not len(error):
                error = None
            elif len(error) != len(y):
                self.mismatch['Error bars'] = True
            else:
                self.mismatch['Error bars'] = False
        if fit is not None:
            if not len(fit):
                fit = None
            elif fitx is not None:
                if len(fit) != len(fitx):
                    self.mismatch['Fit values'] = True
            else:
                self.mismatch['Fit values'] = False
        if not any(self.mismatch.values()):
            self.timer.stop()
        else:
            if not self.timer.isActive():
                self.timer.start(1000)
            return

        self.clear()
        # self.plot(x, y, pen=None, symbol="x")
        self.plot(x, y, pen=(255, 255, 0),
                        symbol='x',
                        symbolBrush=(255, 255, 0),
                        symbolPen='w')
        self.setTitle(title)
        if error is not None:
            # See https://github.com/pyqtgraph/pyqtgraph/issues/211
            if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
                error = np.array(error)
            errbars = pyqtgraph.ErrorBarItem(
                x=np.array(x), y=np.array(y), height=error)
            self.addItem(errbars)
        if fit is not None:
            if fitx is not None:
                self.plot(fitx, fit, pen=(197, 5, 12))  # Badger Red
            else:
                self.plot(range(len(fit)), fit, pen=(197, 5, 12))  # Badger Red

        if marker_point is not None:
            print("whats wrong")
            markerx, markery = marker_point
            print(markerx, markery)
            self.plot([markerx], [markery], symbol='o', symbolPen=(0, 0, 255))

    def length_warning(self):
        self.clear()
        text = "⚠️ dataset lengths mismatch:\n"
        errors = ', '.join([k for k, v in self.mismatch.items() if v])
        text = ' '.join([errors, "should have the same length as Y values"])
        self.addItem(pyqtgraph.TextItem(text))


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("y", "Y values")
    applet.add_dataset("y0", "y0 value. plot will be normalized to this if given", required=False)
    applet.add_dataset("pts", "max number of pts to plot", required=False)
    applet.add_dataset("x", "X values", required=False)
    applet.add_dataset("marker", "a single plot point as a marker", required=False)
    applet.add_dataset("error", "Error bars for each X value", required=False)
    applet.add_dataset("fit", "Fit values", required=False)
    applet.add_dataset("fitx", "X values for the fit", required=False)

    applet.run()


if __name__ == "__main__":
    main()