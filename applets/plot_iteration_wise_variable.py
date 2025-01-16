"""
for a dataset 'counts' which is updated every experiment measurement compute the iteration-wise mean
 and scale it by 't_exposure'. 'measurements' is the number of measurements per iteration.

Example: if you record the amount of counts on the SPCM from the FORT on each measurement in a GeneralVariableScan,
then after each iteration (one scan step, say, 100 measurements), this will update with the computed mean for that
iteration.

applet command:
ppython "C:\...\artiq-master\repository\qn_artiq_routines\applets\plot_loading_background.py"
SPCM0_RO1 SPCM0_RO2 n_measurements single_atom_counts_per_s t_SPCM_first_shot
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
            counts = np.array(data[self.args.counts][1][1:])
            measurements = data[self.args.measurements][1]
            t_exposure = data[self.args.t_exposure][1]
            color = data.get(self.args.color, (False, None))[1]
            if color is None:
                color = (0, 0, 255)

            iteration = len(counts)//measurements
            if iteration > 0:
                mean_by_iteration = [
                    np.mean(counts[i * measurements:(i + 1) * measurements] / t_exposure)
                    for i in range(iteration)]

                self.clear()
                self.plot(range(iteration), mean_by_iteration,
                          pen=color,
                          symbol='o',
                          symbolBrush=color,
                          symbolPen='w',
                          name='shot 1 background')
                self.addLegend()
            else:
                self.clear()
        except:
            self.clear()


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts", "first readout photocounts")
    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("t_exposure", "the atom readout exposure time")
    applet.add_dataset("color", "color", required=False)

    applet.run()


if __name__ == "__main__":
    main()