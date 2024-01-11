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