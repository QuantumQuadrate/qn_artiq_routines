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
        try:
            counts0 = np.array(data[self.args.exc0_counts][1][1:])
            counts1 = np.array(data[self.args.exc0_counts][1][1:])
            readout_counts = np.array(data[self.args.readout_counts][1][1:])


            measurements = data[self.args.measurements][1]
            n_exc_cycles = data[self.args.n_exc_cycles][1]

            threshold_cts_per_s = data[self.args.threshold_cts_per_s][1]
            t_exposure = data[self.args.t_exposure][1]
            cutoff = int(t_exposure*threshold_cts_per_s)

            count_len = measurements * n_exc_cycles

            iteration = len(counts0) // count_len

            if iteration > 0:
                counts_no_atom = np.zeros(iteration)
                counts_with_atom = np.zeros(iteration)

                x = np.arange(iteration)

                try:
                    nsteps = len(data.get(self.args.scan_sequence1, (False, None))[1])
                    scan_sequence1 = data[self.args.scan_sequence1][1]
                    if nsteps > 1 or scan_sequence1 != [0.0]:
                        x = np.array(scan_sequence1[:iteration])
                except:
                    pass


                for i in range(iteration):
                    counts0_i = counts0[i * count_len : (i+1) * count_len]
                    counts1_i = counts0[i * count_len: (i + 1) * count_len]
                    readout_counts_i = readout_counts[i * count_len: (i + 1) * count_len]


                    counting_no_atom = 0
                    counting_with_atom = 0


                    for a, b, c in zip(counts0_i, counts1_i, readout_counts_i):
                        if (a > 0 or b > 0) and c > cutoff:
                            counting_with_atom += 1
                        elif (a > 0 or b > 0) and c <= cutoff:
                            counting_no_atom += 1


                    # total counts from SPCM0 & SPCM1

                    counts_no_atom[i] = counting_no_atom
                    counts_with_atom[i] = counting_with_atom

                self.clear()


                self.plot(x, counts_no_atom,
                          pen=None,
                          symbol='o',
                          symbolBrush='royalblue',
                          symbolPen='royalblue',
                          name='no atom')
                self.plot(x, counts_with_atom,
                          pen=None,
                          symbol='o',
                          symbolBrush='crimson',
                          symbolPen='crimson',
                          name='with atom')

                # self.setYRange(-0.0, 1.0, padding=0)

                title = str(data[self.args.scan_vars][1])

                self.setTitle(title)
                self.addLegend()
            else:
                self.clear()

        except:
            self.clear()





def main():
    applet = TitleApplet(XYPlot)

    applet.add_dataset("exc0_counts", "excitation counts from SPCM0")
    applet.add_dataset("exc1_counts", "excitation counts from SPCM1")
    applet.add_dataset("readout_counts", "readout counts from SPCM0")
    applet.add_dataset("n_exc_cycles", "how many excitation cycles per measurement")

    applet.add_dataset("measurements", "how many measurements per data point")
    applet.add_dataset("threshold_cts_per_s", "the atom signal counts/s threshold")
    applet.add_dataset("t_exposure", "the atom readout exposure time")
    applet.add_dataset("scan_vars", "the names of the scan variable(s)", required=False)
    applet.add_dataset("scan_sequence1", "the scan steps", required=False)

    applet.run()


if __name__ == "__main__":
    main()