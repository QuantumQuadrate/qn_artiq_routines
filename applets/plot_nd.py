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
        try:
            # should be a numpy array where each row is a different data channel
            y_data = data[self.args.y][1]
            print(y_data)
        except KeyError:
            return
        x = data.get(self.args.x, (False, None))[1]
        if x is None:
            x = np.arange(y_data.shape[1])

        self.clear()
        nrows = y_data.shape[0]
        for i,y in enumerate(y_data):
            self.plot(x, y, pen=(i,nrows), symbol="o")
        self.setTitle(title)
        # if error is not None:
        #     # See https://github.com/pyqtgraph/pyqtgraph/issues/211
        #     if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
        #         error = np.array(error)
        #     errbars = pyqtgraph.ErrorBarItem(
        #         x=np.array(x), y=np.array(y), height=error)
        #     self.addItem(errbars)
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
    applet.add_dataset("y", "Y values")
    applet.add_dataset("x", "X values", required=False)
    # applet.add_dataset("error", "Error bars for each X value", required=False)
    # applet.add_dataset("fit", "Fit values for each X value", required=False)
    applet.run()

if __name__ == "__main__":
    main()
