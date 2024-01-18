"""
plot the atom retention, i.e., the fraction of atoms detected in the first readout which
survive the science phase and are hence detected in the second readout

arguments are counts, counts2, and a defined threshold to discriminate between atom
and background. todo: make the threshold parameter optional and instead determine it
with a fit
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
        counts_shot1 = data[self.args.counts_shot1][1][1:]
        counts_shot2 = data[self.args.counts_shot2][1][1:]
        measurements = 100 #data[self.args.measurements][1]
        iteration = data[self.args.iteration][1]
        cutoff = 180 #data[self.args.cutoff][1]

        retention_array = np.zeros(iteration)
        loading_rate_array = np.zeros(iteration)
        n_atoms_loaded_array = np.zeros(iteration)

        for i in range(iteration):
            shot1 = counts_shot1[i * measurements:(i + 1) * measurements]
            shot2 = counts_shot2[i * measurements:(i + 1) * measurements]
            atoms_loaded = [x > cutoff for x in shot1]
            n_atoms_loaded = sum(atoms_loaded)
            atoms_retained = [x > cutoff and y for x, y in zip(shot2, atoms_loaded)]
            retention_fraction = 0 if not n_atoms_loaded > 0 else sum(atoms_retained) / sum(atoms_loaded)
            retention_array[i] = retention_fraction
        #
        error = [1/np.sqrt(n) if n > 0 else 0 for n in n_atoms_loaded_array]

        # fit = data.get(self.args.fit, (False, None))[1]

        # if not len(y) or len(y) != len(x):
        #     self.mismatch['X values'] = True
        # else:
        #     self.mismatch['X values'] = False
        # if error is not None and hasattr(error, "__len__"):
        #     if not len(error):
        #         error = None
        #     elif len(error) != len(y):
        #         self.mismatch['Error bars'] = True
        #     else:
        #         self.mismatch['Error bars'] = False
        # if fit is not None:
        #     if not len(fit):
        #         fit = None
        #     elif len(fit) != len(y):
        #         self.mismatch['Fit values'] = True
        #     else:
        #         self.mismatch['Fit values'] = False
        # if not any(self.mismatch.values()):
        #     self.timer.stop()
        # else:
        #     if not self.timer.isActive():
        #         self.timer.start(1000)
        #     return
        #
        # self.clear()
        # # self.plot(x, y, pen=None, symbol="x")
        x = range(iteration)
        self.plot(x, retention_array, pen=(255, 0, 0),
                        symbol='o',
                        symbolBrush=(255, 0, 0),
                        symbolPen='w')
        self.setTitle(title)
        # if error is not None:
        #     # See https://github.com/pyqtgraph/pyqtgraph/issues/211
        #     if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
        #         error = np.array(error)
        # errbars = pyqtgraph.ErrorBarItem(
        #     x=x, y=retention_array, height=error)
        # self.addItem(errbars)
        # if fit is not None:
        #     xi = np.argsort(x)
        #     self.plot(x[xi], fit[xi])

    def length_warning(self):
        self.clear()
        text = "⚠️ dataset lengths mismatch:\n"
        errors = ', '.join([k for k, v in self.mismatch.items() if v])
        text = ' '.join([errors, "should have the same length as Y values"])
        self.addItem(pyqtgraph.TextItem(text))

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts_shot1", "first readout photocounts")
    applet.add_dataset("counts_shot2", "second readout photocounts")
    # applet.add_dataset("threshold", "the atom signal threshold", required=False)
    # applet.add_dataset("threshold", "the atom signal threshold", required=False)
    applet.run()

if __name__ == "__main__":
    main()

    # self.plot(x, y, pen=(255, 138, 0),
    #                 symbol='x',
    #                 symbolBrush=(255, 138, 0),
    #                 symbolPen='w')