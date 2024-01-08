#!/usr/bin/env python3

"""
better than the stock ARTIQ histogram which expects you to send it pre-binned data.
the dataset should not be pre-binned here. just send it in and add an optional
binning argument
"""

import PyQt5    # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph
import numpy as np

from artiq.applets.simple import TitleApplet


class HistogramPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args

    def data_changed(self, data, mods, title):
        try:
            y = data[self.args.y][1][1:]
            x = (data[self.args.x][1])[0]
            if self.args.x is None:
                bins = 100
            else:
                bins = x
            if self.args.t_exposure is None:
                t_exp = None
            else:
                t_exp = data[self.args.t_exposure][1]
            if self.args.color is None:
                color = 'b'
            else:
                color = data[self.args.color][1][0]
            if self.args.ignore_first_n_points is not None:
                n = data[self.args.ignore_first_n_points][1]
                y = y[n:]
        except KeyError:
            return

        hist, bin_edges = np.histogram(y, bins=bins)

        self.clear()
        self.plot(bin_edges, hist, fillLevel=0, stepMode=True,
                  brush=(0, 0, 255, 150), pen=color)

        self.setTitle(title)


def main():
    applet = TitleApplet(HistogramPlot)
    applet.add_dataset("y", "Y values")
    applet.add_dataset("x", "Bin boundaries", required=False)
    applet.add_dataset("t_exposure", "exposure time", required=False)
    applet.add_dataset("color", "Color", required=False)
    applet.add_dataset("ignore_first_n_points", "the number of leading data points to ignore", required=False)

    applet.run()

if __name__ == "__main__":
    main()
