"""
Plot the signals representing the laser beam powers. units are determined by the datasets

applet command:
python "C:\...\qn_artiq_routines\applets\bar_plot_MOT_powers.py"
MOT1_monitor MOT2_monitor MOT3_monitor MOT4_monitor MOT5_monitor MOT6_monitor
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph

from artiq.applets.simple import TitleApplet


class XYPlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        # self.labels = [f'MOT{i + 1}' for i in range(6)]

        self.labels = ["health_check_uw_freq00", "health_check_uw_freq01", "health_check_uw_freq11"]
        self.labels_to_show = ["Freq 00\nFidelity", "Freq 01\nFidelity", "Freq 11\nFidelity"]
        # red-green color-blind friendly RGB colors from ChatGPT
        self.colors = [
                        (230, 97, 0),  # orange - yellow
                        (255, 194, 10), # yellow
                        (112, 173, 71), # green
                        (64, 176, 166), # oval
                        (12, 123, 220),  # blue
                        (93, 58, 155),  # violet

                        (255, 138, 0),  # orange
                        (149, 0, 210),  # purple
                        (0, 170, 150),  # teal
                        (128, 128, 128),  # gray
                        (139, 69, 19),  # brown
                        (255, 105, 180)]  # pink


    def data_changed(self, data, mods, title):
        try:
            # the data display will be rolling, only showing display_pts at a time
            # should be a numpy array where each row is a different data channel
            uw_fidelity_data = []

            for i in range(3):
                uw_fidelity_data.append(data[getattr(self.args, self.labels[i])][1])

            # MOT_switchyard_input = data.get(self.args.MOT_switchyard_input, (False, None))[1]
            #
            # if MOT_switchyard_input is not None:
            #     MOT_data.append(data[self.args.MOT_switchyard_input][1][-1])
            #     self.labels += ['MOT_switchyard_input']

        except KeyError:
            return

        x = range(len(uw_fidelity_data))

        self.clear()
        w = 0.8
        xpts = range(-1,len(uw_fidelity_data)+1)
        bar_graph = pyqtgraph.BarGraphItem(x = x, height = uw_fidelity_data, width = w, brushes=self.colors)
        self.setYRange(0.0, 1.1, padding=0)
        self.setXRange(-w, len(uw_fidelity_data)-1+w, padding=0)

        self.plot(xpts, np.full(len(xpts), 0.8), pen='red', style=PyQt5.QtCore.Qt.DashLine, width=0.2)
        self.addItem(bar_graph)
        self.setTitle(title)

        for i,datum in enumerate(uw_fidelity_data):
            if datum > 1.0:
                self.text = pyqtgraph.TextItem(str(round(datum, 3)) + "\n" + self.labels_to_show[i], color=(0, 0, 0))
                self.addItem(self.text)
                self.text.setPos(i - w / 3, 0.95)
            elif datum > 0.5:
                self.text = pyqtgraph.TextItem(str(round(datum,3))+"\n"+self.labels_to_show[i],color=(0,0,0))
                self.addItem(self.text)
                self.text.setPos(i- w/3, round(datum,3)-0.05)
            else:
                self.text = pyqtgraph.TextItem(str(round(datum, 3))+"\n"+self.labels_to_show[i], color=(255, 255, 255))
                self.addItem(self.text)
                self.text.setPos(i - w / 3, round(datum, 3) + 0.4)

def main():
    applet = TitleApplet(XYPlot)
    applet.add_dataset("health_check_uw_freq00", "health_check_uw_freq00")
    applet.add_dataset("health_check_uw_freq01", "health_check_uw_freq01")
    applet.add_dataset("health_check_uw_freq11", "health_check_uw_freq11")
    # applet.add_dataset("MOT_switchyard_input", "MOT PD0 voltage", required=False)
    applet.run()

if __name__ == "__main__":
    main()
