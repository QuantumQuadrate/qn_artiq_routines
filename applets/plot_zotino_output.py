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

            setpoint1 = data[self.args.setpoint1][1]
            setpoint2 = data[self.args.setpoint2][1]
            measurements = data[self.args.measurements][1]

            iteration = len(counts1)//measurements


            offset1 = 0.0 # PD1 background @ 0 power
            offset2 = 0.0 # PD2 background @ 0 power

            norm_factor1 = setpoint1
            norm_factor2 = setpoint2


            # norm_factor1 = 0.0003 # Repump 1 normalizing factor
            # norm_factor2 = 0.00036 # Repump 2 normalizing factor

            if iteration > 0:

                mean1_by_iteration = [
                    np.mean((counts1[i * measurements:(i + 1) * measurements]+offset1)/norm_factor1)
                    for i in range(iteration)]
                mean2_by_iteration = [
                    np.mean((counts2[i * measurements:(i + 1) * measurements]+offset2)/norm_factor2)
                    for i in range(iteration)]
                # mean3_by_iteration = [
                #     np.mean((counts3[i * measurements:(i + 1) * measurements]+offset3)/norm_factor3)
                #     for i in range(iteration)]
                # mean4_by_iteration = [
                #     np.mean((counts4[i * measurements:(i + 1) * measurements]+offset4)/norm_factor4)
                #     for i in range(iteration)]
                # mean5_by_iteration = [
                #     np.mean((counts5[i * measurements:(i + 1) * measurements]+offset5)/norm_factor5)
                #     for i in range(iteration)]
                # mean6_by_iteration = [
                #     np.mean((counts6[i * measurements:(i + 1) * measurements]+offset6)/norm_factor6)
                #     for i in range(iteration)]

                self.clear()

                self.plot(range(iteration), mean1_by_iteration,
                          pen='aquamarine',
                          symbol='o',
                          symbolBrush='aquamarine',
                          symbolPen='w',
                          name='test1')

                self.plot(range(iteration), mean2_by_iteration,
                          pen='salmon',
                          symbol='o',
                          symbolBrush='salmon',
                          symbolPen='w',
                          name='test2')

                # self.plot(range(iteration), mean3_by_iteration,
                #           pen='tomato',
                #           symbol='o',
                #           symbolBrush='tomato',
                #           symbolPen='w',
                #           name='REPUMP_AO3')
                #
                # self.plot(range(iteration), mean4_by_iteration,
                #           pen='dodgerblue',
                #           symbol='o',
                #           symbolBrush='dodgerblue',
                #           symbolPen='w',
                #           name='REPUMP_AO4')
                #
                # self.plot(range(iteration), mean5_by_iteration,
                #           pen='g',
                #           symbol='o',
                #           symbolBrush='g',
                #           symbolPen='w',
                #           name='REPUMP_AO5')
                #
                # self.plot(range(iteration), mean6_by_iteration,
                #           pen='orange',
                #           symbol='o',
                #           symbolBrush='orange',
                #           symbolPen='w',
                #           name='REPUMP_AO6')

                self.addLegend()
            else:
                self.clear()
        except:
            self.clear()


def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("counts1", "AZ bottom")
    applet.add_dataset("counts2", "AZ top")
    applet.add_dataset("setpoint1", "AZ bottom setpoint")
    applet.add_dataset("setpoint2", "AZ top setpoint")
    # applet.add_dataset("counts3", "repump3 power")
    # applet.add_dataset("counts4", "repump4 power")
    # applet.add_dataset("counts5", "repump5 power")
    # applet.add_dataset("counts6", "repump6 power")
    applet.add_dataset("measurements", "how many measurements per data point")

    applet.run()


if __name__ == "__main__":
    main()
