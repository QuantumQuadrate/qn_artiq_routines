"""
An independent code to monitor ttl0 for signals timestamp the signals and count the total number of time stamps.
Saves the timestamps as well as the total counts in txt file next to the HDF5 file. Also counts the signals using
ttl.count and edge counter for comparison.

When timestamp after openning the gate, artiq stamps only one event (I think the first event). To time stamp all the
events received within  the window I have used a while loop to tag up to max_clicks counts in each window.

timestamping and ttl.count both consume the window, thus, we cannot do ttl.count and timestamp with the same window.

"""
### time tagging and counting the tags:
from artiq.experiment import *
import numpy as np  # Use numpy for random number generation

class TimeStamp_ttl(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl0")
        self.setattr_device("ttl0_counter")

        self.setattr_argument("num_iterations", NumberValue(100, ndecimals=0, step=1))
        self.setattr_argument("t_SPCM_exposure", NumberValue(100 * us, unit='us'))

    @kernel
    def run(self):
        self.core.reset()
        self.ttl0.input()

        max_clicks = 10 ### maximum number of clicks that will be time tagged in each gate window
        timestamps = [[-1.0] * max_clicks for _ in range(self.num_iterations)]
        ttl_counts = [0] * self.num_iterations ### counts we measure using ttl.count
        edge_counts = [0] * self.num_iterations ### counts we measure using edge counter
        stamp_counts = [0] * self.num_iterations ### counts we measure by adding up the timestamps

        for i in range(self.num_iterations):
            delay(10*ms)
            t_end_SPCM0 = self.ttl0.gate_rising(self.t_SPCM_exposure)

            click_counter = 0
            while click_counter < max_clicks:
                click_time = self.ttl0.timestamp_mu(t_end_SPCM0)
                if click_time == -1:
                    break
                timestamps[i][click_counter] = self.core.mu_to_seconds(click_time)
                click_counter += 1

            ### counting by adding up the timestamps
            stamp_counts[i] = click_counter

            ### counting using ttl.count:
            delay(10 * us)
            t_end_SPCM0 = self.ttl0.gate_rising(self.t_SPCM_exposure)
            ttl_counts[i] = self.ttl0.count(t_end_SPCM0)

            ### counting using edge counter:
            delay(10 * us)
            self.ttl0_counter.gate_rising(self.t_SPCM_exposure)
            edge_counts[i] = self.ttl0_counter.fetch_count()


        self.save_data(timestamps)
        self.save_counts("ttl_counts", ttl_counts)
        self.save_counts("edge_counts", edge_counts)
        self.save_counts("stamp_counts", stamp_counts)

        print("************  ALL DONE  ************")


    @rpc(flags={"async"})
    def save_data(self, data):
        with open("Click_timestamps.txt", "a") as file:
            for iteration in data:
                file.write(", ".join(map(str, iteration)) + "\n")

    @rpc(flags={"async"})
    def save_counts(self, title, data):
        with open("Click_Counts.txt", "a") as file:
            file.write(title + ": ")
            file.write(", ".join(map(str, data)) + "\n")

