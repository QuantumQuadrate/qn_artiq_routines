"""
arguments are counts, counts2, and a defined threshold to discriminate between atom
and background. todo: make the threshold parameter optional and instead determine it
with a fit

applet command:
ppython "C:\...\artiq-master\repository\qn_artiq_routines\applets\plot_loading_background.py"
SPCM0_RO1 SPCM0_RO2 n_measurements single_atom_threshold t_SPCM_first_shot
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

        try: # not all of these are persistent
            counts_shot1 = np.array(data[self.args.counts_shot1][1][1:])
            counts_shot2 = np.array(data[self.args.counts_shot2][1][1:])
            measurements = data[self.args.measurements][1]
            threshold_cts_per_s = data[self.args.threshold_cts_per_s][1]
            t_exposure = data[self.args.t_exposure][1]
            cutoff = t_exposure*threshold_cts_per_s

            iteration = len(counts_shot1)//measurements
            if iteration > 0:
                if len(counts_shot1) == len(counts_shot2):
                    mean1_by_iteration = [
                        np.mean(counts_shot1[i * measurements:(i + 1) * measurements][
                                counts_shot1[i * measurements:(i + 1) * measurements] < cutoff] / t_exposure)
                        for i in range(iteration)]

                    mean2_by_iteration = [
                        np.mean(counts_shot2[i * measurements:(i + 1) * measurements][
                                counts_shot2[i * measurements:(i + 1) * measurements] < cutoff] / t_exposure)
                        for i in range(iteration)]

                    self.clear()
                    self.plot(range(iteration), mean1_by_iteration,
                              pen=(0, 0, 255),
                              symbol='o',
                              symbolBrush=(0, 0, 255),
                              symbolPen='w',
                              name='shot 1 background')
                    self.plot(range(iteration), mean2_by_iteration,
                              pen=(255, 0, 0),
                              symbol='o',
                              symbolBrush=(255, 0, 0),
                              symbolPen='w',
                              name='shot 2 background')
                    
                # todo: add std error
                # if error is not None:
                #     # See https://github.com/pyqtgraph/pyqtgraph/issues/211
                #     if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
                #         error = np.array(error)
                # errbars = pyqtgraph.ErrorBarItem(
                #     x=x, y=retention_array, height=2*error, pen=(255, 0, 0)) # error should be +/- the std, hence 2*
                # self.addItem(errbars)
                self.addLegend()
            else:
                self.clear()
        except:
            self.clear()

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts_shot1", "first readout photocounts")
    applet.add_dataset("counts_shot2", "second readout photocounts")
    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("threshold_cts_per_s", "the atom signal counts/s threshold")
    applet.add_dataset("t_exposure", "the atom readout exposure time")

    applet.run()


if __name__ == "__main__":
    main()