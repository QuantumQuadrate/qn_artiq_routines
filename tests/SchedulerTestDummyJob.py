"""
try submitting this run with scheduler test
"""

from artiq.experiment import *
from time import sleep

class SchedulerTestDummyJob(EnvExperiment):

    def build(self):
        self.setattr_device("scheduler")
        self.setattr_argument("test_variable",NumberValue(42))
        self.setattr_device("core")
        self.setattr_device("led0")

    def prepare(self):
        pass

    @rpc(flags={'async'})
    def print_job_attrs(self):
        print("test job submitted has the following:")
        print(self.scheduler.rid)
        print(self.scheduler.expid)
        print(self.scheduler.pipeline_name)

    @kernel
    def run(self):
        self.core.reset()
        self.print_job_attrs()

        for i in range(5):
            self.led0.pulse(10 * us)
            delay(0.5*s)
        print("higher priority experiment finished")

"""
print(dir(self.scheduler))

print:['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', 
       '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', 
       '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__',
       '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', 
       '__str__', '__subclasshook__', '__weakref__', '_check_pause', 
       '_check_termination', '_submit', 'check_pause', 'check_termination', 
       'delete', 'expid', 'get_status', 'pause', 'pause_noexc', 'pipeline_name', 
       'priority', 'request_termination', 'rid', 'set_run_info', 'submit']"""