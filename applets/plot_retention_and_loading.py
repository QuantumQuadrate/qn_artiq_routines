"""
plot the atom retention and loading rate

retention: the fraction of atoms detected in the first readout which
survive the science phase and are hence detected in the second readout

loading rate: the fraction of successful attempts to load an atom which

applet command:
python "C:\..\qn_artiq_routines\applets\plot_retention_and_loading.py"
photocounts photocounts2 n_measurements iteration single_atom_counts_per_s t_SPCM_first_shot
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

        counts_shot1 = data[self.args.counts_shot1][1][1:]
        counts_shot2 = data[self.args.counts_shot2][1][1:]
        measurements = data[self.args.measurements][1]
        iteration = int(data[self.args.iteration][1])

        if (self.args.threshold_cts_per_s is not None and
                self.args.t_exposure is not None):
            threshold_cts_per_s = data[self.args.threshold_cts_per_s][1]
            t_exposure = data[self.args.t_exposure][1]
            cutoff = int(t_exposure*threshold_cts_per_s)
        else:
            cutoff = 220 # for testing purposes # todo: compute this

        retention_array = np.zeros(iteration)
        loading_rate_array = np.zeros(iteration)
        n_atoms_loaded_array = np.zeros(iteration)

        if iteration > 0:
            for i in range(iteration):
                shot1 = counts_shot1[i * measurements:(i + 1) * measurements]
                shot2 = counts_shot2[i * measurements:(i + 1) * measurements]
                atoms_loaded = [x > cutoff for x in shot1]
                n_atoms_loaded = sum(atoms_loaded)
                n_atoms_loaded_array[i] = n_atoms_loaded
                atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
                loading_fraction = n_atoms_loaded / measurements
                loading_rate_array[i] = loading_fraction
                retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / n_atoms_loaded
                retention_array[i] = retention_fraction
            #
            error = np.array([1/np.sqrt(n) if n > 0 else 0 for n in n_atoms_loaded_array])

            x = np.arange(iteration)
            self.plot(x, retention_array, pen=(255, 0, 0),
                      symbol='o',
                      symbolBrush=(255, 0, 0),
                      symbolPen='w',
                      name='retention')
            self.plot(x, loading_rate_array, pen=(0, 100, 100),
                      symbol='o',
                      symbolBrush=(0, 100, 100),
                      symbolPen='w',
                      name='loading')
            self.setYRange(-0.05, 1.05, padding=0)
            self.setTitle(title)
            # self.addLegend()

            # the text that shows up is huge. # todo
            # xaxis = pyqtgraph.AxisItem('bottom')#,text='my x axis')
            # xaxis.setLabel(text='test')
            # self.addItem(xaxis)

            #
            if error is not None:
                # See https://github.com/pyqtgraph/pyqtgraph/issues/211
                if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
                    error = np.array(error)
            errbars = pyqtgraph.ErrorBarItem(
                x=x, y=retention_array, height=error, pen=(255, 0, 0))
            self.addItem(errbars)
            # if fit is not None:
            #     xi = np.argsort(x)
            #     self.plot(x[xi], fit[xi])

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts_shot1", "first readout photocounts")
    applet.add_dataset("counts_shot2", "second readout photocounts")
    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("iteration", "the current experiment iteration")
    applet.add_dataset("threshold_cts_per_s", "the atom signal counts/s threshold")
    applet.add_dataset("t_exposure", "the atom readout exposure time")
    applet.run()

if __name__ == "__main__":
    main()

    # self.plot(x, y, pen=(255, 138, 0),
    #                 symbol='x',
    #                 symbolBrush=(255, 138, 0),
    #                 symbolPen='w')