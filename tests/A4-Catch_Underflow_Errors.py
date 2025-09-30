"""
Catching underflow errors or managing slack dynamically to avoid underflow errors completely. Enable only one of the parts in the code below.

Akbar 2025-03-12
"""
from artiq.experiment import *

class Catch_Underflow_errors(EnvExperiment):

    def build(self):
       self.setattr_device("core")
       self.setattr_device("ttl7")


    @kernel
    def run(self):
        self.core.reset()
        # self.ttl7.output()

        delay(1 * us)
        # self.ttl7.off()



        #### increasing the delay dynamically only when the slack is smaller than a threshold.
        #### This works well and adds a delay when necessary to avoid Underflow error.
        slack_threshold = 1e4 ### acceptable slack in mu. 1e4 = 10us.
        for x in range(100):
            # self.ttl7.pulse(90 * ns)
            delay(10*ns)
            slack = now_mu() - self.core.get_rtio_counter_mu()
            if slack<slack_threshold:
                delay(1*ms)




        # ### Catching the underflow using try/except. This works. However, note that it catches the error after it has occured.
        # ### Thus, the ttl signal sometimes stays high during the 1ms delay.
        # for x in range(100):
        #     try:
        #         self.ttl7.pulse(90 * ns)
        #         delay(10*ns)
        #     except RTIOUnderflow:
        #         delay(1*ms)






        # #### prints out the slack in mu. You will see that the slack decreases gradually until we get Underflow.
        # for x in range(100):
        #     self.ttl7.pulse(90 * ns)
        #     delay(10 * ns)
        #     slack = now_mu() - self.core.get_rtio_counter_mu()
        #     self.print_async(slack)






        # #### gives Underflow unless the delay is 1000ns:
        # for x in range(100):
        #     self.ttl7.pulse(20 * ns)
        #     delay(100 * ns)


        print("*************   Experiment finished   *************")

    @rpc(flags={"async"})
    def print_async(self, x):
        """print asynchronously so we don't block the RTIO counter.
        useful for debugging"""
        print(x)
