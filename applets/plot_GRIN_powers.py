"""
for a dataset 'counts' which is updated every experiment measurement compute the iteration-wise mean
 and scale it by 't_exposure'. 'measurements' is the number of measurements per iteration.


applet command:

This is from plot_iteration_wise_variable.py

"""
# todo: make this code universal for REPUMP, GRIN1, PR...?

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

            norm_factor1 = 0.03    # D1
            norm_factor2 = 0.1  # Exc

            # norm_factor1 = 1.0  # D1
            # norm_factor2 = 1.0  # Exc

            offset1 = 0.0
            offset2 = 0.0

            if iteration > 0:

                mean1_by_iteration = [
                    np.mean((counts1[i * measurements:(i + 1) * measurements]+offset1)/norm_factor1)
                    for i in range(iteration)]

                mean2_by_iteration = [
                    np.mean((counts2[i * measurements:(i + 1) * measurements]+offset2)/norm_factor2)
                    for i in range(iteration)]

                self.clear()

                self.plot(range(iteration), mean1_by_iteration,
                          pen='tomato',
                          symbol='o',
                          symbolBrush='tomato',
                          symbolPen='w',
                          name='GRIN1_D1')

                self.plot(range(iteration), mean2_by_iteration,
                          pen='dodgerblue',
                          symbol='o',
                          symbolBrush='dodgerblue',
                          symbolPen='w',
                          name='GRIN1_EXC')

                self.addLegend()
            else:
                self.clear()
        except:
            self.clear()


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts1", "GRIN1_D1 power")
    applet.add_dataset("counts2", "GRIN1_EXC power")
    applet.add_dataset("measurements", "how many measurements per data point")

    applet.run()


if __name__ == "__main__":
    main()