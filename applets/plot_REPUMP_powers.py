"""
for a dataset 'counts' which is updated every experiment measurement compute the iteration-wise mean
 and scale it by 't_exposure'. 'measurements' is the number of measurements per iteration.


applet command:

This is from plot_iteration_wise_variable.py

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
        # self.timer.timeout.connect(self.length_warning)
        self.mismatch = {'X values': False,
                         'Error bars': False,
                         'Fit values': False}

    def data_changed(self, data, mods, title):

        try:  # not all of these are persistent
            counts1 = np.array(data[self.args.counts1][1][1:])
            counts2 = np.array(data[self.args.counts2][1][1:])
            measurements = data[self.args.measurements][1]

            iteration = len(counts1)//measurements

            zoom = 100
            offset = 0.02

            if iteration > 0:

                mean1_by_iteration = [
                    np.mean((counts1[i * measurements:(i + 1) * measurements]+offset)*zoom)
                    for i in range(iteration)]

                mean2_by_iteration = [
                    np.mean((counts2[i * measurements:(i + 1) * measurements])*zoom)
                    for i in range(iteration)]

                self.clear()

                self.plot(range(iteration), mean1_by_iteration,
                          pen=(0, 0, 255),
                          symbol='o',
                          symbolBrush=(0, 0, 255),
                          symbolPen='w',
                          name='REPUMP1_PD')

                self.plot(range(iteration), mean2_by_iteration,
                          pen=(255, 0, 0),
                          symbol='o',
                          symbolBrush=(255, 0, 0),
                          symbolPen='w',
                          name='REPUMP2_PD')

                self.addLegend()
            else:
                self.clear()
        except:
            self.clear()


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts1", "repump1 power")
    applet.add_dataset("counts2", "repump2 power")
    applet.add_dataset("measurements", "how many measurements per data point")

    applet.run()


if __name__ == "__main__":
    main()
