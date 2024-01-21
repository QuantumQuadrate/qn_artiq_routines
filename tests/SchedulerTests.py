from artiq.experiment import *
from time import sleep

class SchedulerTest(EnvExperiment):

    def build(self):
        self.setattr_device("scheduler")
        self.setattr_device("core")
        self.setattr_device("led0")

    def prepare(self):
        self.new_job_expid = {'log_level': 30,
                             'file': 'qn_artiq_routines\\tests\\SchedulerTestDummyJob.py',
                             'class_name': 'SchedulerTestDummyJob',
                             # 'arguments': {}, # this will just use the defaults
                             'arguments': {'test_variable': 3.14},
                             'repo_rev': 'N/A'}

        print(self.scheduler.rid)
        print(self.scheduler.priority)
        print(self.scheduler.expid)
        print(self.scheduler.pipeline_name)

    @rpc(flags={'async'})
    def submit_higher_priority_job(self):
        print("submitting another experiment")
        self.scheduler.submit(pipeline_name=self.scheduler.pipeline_name,
                              expid=self.new_job_expid,
                              priority=self.scheduler.priority + 1,  # higher priority than the current run
                              due_date=None,
                              flush=False)

    @kernel
    def my_function(self):
        self.core.reset()
        self.led0.pulse(10*us)

    def run(self):

        for i in range(10):

            print("parent experiment, i=",i)
            self.my_function()
            sleep(1)

            # check some condition, e.g. i
            if i == 5:
                self.submit_higher_priority_job()

            if self.scheduler.check_pause():
                self.core.comm.close()  # put the hardware in a safe state before checking pause
                self.scheduler.pause()  # check if we need to run a new experiment*

            # *technically we already know we will need to pause, so we could pull it
            # inside the if statement in anticipation of that, but this simulates a
            # more a general use case

        print("finish this run")
        sleep(2)
        print("scheduler test done")



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