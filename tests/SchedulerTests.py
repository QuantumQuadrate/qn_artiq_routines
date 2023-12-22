from artiq.experiment import *

class SchedulerTest(EnvExperiment):

    def build(self):
        self.setattr_device("scheduler")

    def prepare(self):
        pass

    def run(self):
        print(self.scheduler.rid)
        print(self.scheduler.expid)
        # print(dir(self.scheduler))
        print(self.scheduler.pipeline_name)