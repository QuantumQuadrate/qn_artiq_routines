### Test test2

# from artiq.experiment import *
# # import pandas as pd
# import time
# from artiq.language.core import *
# from artiq.language.types import *
#
#
# St_t = time.time()
# file0 = open('C:/Users/QC/OneDrive - UW-Madison/Desktop/Akbar/AkSTimeList0.txt', 'a')
# file1 = open('C:/Users/QC/OneDrive - UW-Madison/Desktop/Akbar/AkSTimeList1.txt', 'a')
# class TimeStamp_ttl(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl0")
#         self.setattr_device("ttl1")
#
#     @rpc(flags={"async"})
#     def saveTime0(self, x):
#         file0.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl0.input()
#         self.ttl1.input()
#         timestamps = np.range(100)
#
#         delay(1 * us)
#
#         for x in range(10):
#             delay(10 * ms)
#             with parallel:
#                 tend0 = self.ttl0.gate_rising(10 * us)
#                 tend1 = self.ttl1.gate_rising(10 * us)
#
#             # ### this saves the timestamp of the first event only:
#             # timePuls0 = self.ttl0.timestamp_mu(tend0)
#             # timePuls1 = self.ttl1.timestamp_mu(tend1)
#             # self.saveTime0(self.core.mu_to_seconds(timePuls0))
#             # self.saveTime1(self.core.mu_to_seconds(timePuls1))
#             # delay(1 * ms)
#
#
#
#         self.saveTime0(timestamps)
#         file0.close()
#         # file1.close()
#
#
#     def analyze(self):
#         End_t = time.time()
#         # self.saveTime1(TimeList)
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))
#
#


from artiq.experiment import *

class TimeStamp_ttl(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl0")

    @kernel
    def run(self):
        self.core.reset()
        self.ttl0.input()

        # Preallocate an array for timestamps
        num_iterations = 10
        max_events = 5  # Maximum number of events per iteration
        timestamps = [[-1.0] * max_events for _ in range(num_iterations)]  # Nested arrays

        for i in range(num_iterations):
            delay(1*ms)
            tend = self.ttl0.gate_rising(100 * us)

            # Read and store timestamps for each detected event
            event_counter = 0
            while event_counter < max_events:
                event_time = self.ttl0.timestamp_mu(tend)
                if event_time == -1:  # No more events detected
                    break
                timestamps[i][event_counter] = self.core.mu_to_seconds(event_time)
                event_counter += 1

        self.save_data(timestamps)
        delay(10 * ms)
        print("************  ALL DONE  ************")

    @rpc(flags={"async"})
    def save_data(self, data):
        with open("AkStimestamps.txt", "a") as file:
            for iteration in data:
                file.write(", ".join(map(str, iteration)) + "\n")

    @rpc(flags={"async"})
    def saveTime1(self, x):
        file1.write(str(x) + '\n')
        # print(x)
