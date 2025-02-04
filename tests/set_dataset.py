from artiq.experiment import *

class SetDataset(EnvExperiment):
    def build(self):
        pass

    def run(self):
        self.set_dataset('t_MOT_loading', 1000*ms, persist=True)