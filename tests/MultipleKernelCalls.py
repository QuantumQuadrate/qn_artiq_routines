from artiq.experiment import *

class MultipleKernelCalls(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("led0")

    def prepare(self):
        pass

    @kernel
    def setup(self):
        self.core.break_realtime()

    # trivial example of function run on kernel with a return
    @kernel
    def my_function(self) -> TFloat:
        self.core.reset()
        self.led0.pulse(10*us)
        return 0.5

    @kernel
    def my_function2(self):
        self.core.reset()
        self.led0.pulse(10 * us)

    def caller(self, func):
        res = func()
        print(res)

    def run(self):
        # self.setup()
        res = self.my_function()
        res = self.my_function()
        self.my_function2()

        print(res)
        print("experiment finished")

    # similar to previous example but use a wrapper function for the kernels
    def run(self):
        for i in range(10):
            self.caller(self.my_function)

        print("experiment finished")