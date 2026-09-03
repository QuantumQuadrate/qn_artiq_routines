#########  This is the practice code that I used to implement different things in artiq mainly focusing on
####    real-time input/out features with TTLs. Examples include time tagging and counting TTL pulses recieved
####    at two ttls in parallel, generating a TTL pulse from one channel, and detecting on the other channel with
####    just a few ns delay, etc.

from artiq.experiment import *
# import pandas as pd
import time
from artiq.language.core import *
from artiq.language.types import *



######  Generating a single TTL pulse from channel 0 and reading it from channel 5:
# class Tutorial(EnvExperiment):
#
#     def build(self):
#
#         self.setattr_device("core")
#         self.setattr_device("ttl0")
#         self.setattr_device("ttl5")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl0.input()
#         self.ttl5.output()
#
#         self.ttl5.on()
#         delay(10 * us)
#
#         self.ttl0.sample_input()  # read the value of the TTL input at the position of the time cursor
#         TTLValue = self.ttl0.sample_get()  # Returns the value of a sample previously obtained with sample_input()
#
#         delay(5 * ms)
#         self.ttl5.off()
#
#         delay(5 * ms)
#
#         print(TTLValue)
#         print("core done!")






##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################                  My time tagging attemp starts here                ############################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################



# #### A TTL pulse was generated with a AWG and sent to TTL2. It reads the signal and sends a triggering
# ####    TTL signal to TTL6. It worked.
# #### The time tracking below helps to understand how the gating and finding the edge interact with the
# ####    timeline.
# #### Basically, when openning the gate for say 60ms, it puts the time cursor at 60ms later than "now".
# ####    When a rising edge is detected during this 60ms time, it records the time of the edge and "unblocks"
# ####    the core so it can do other things. However, it does not move the cursor. If we want to generate a
# ####    TTL pulse at a specific time after the rising edge of the input TTL, we have to use at_mu(t_edge).
# ####    Otherwise, the pulse will be generated at a random time (it should be generated at now+60ms, I guess.
# ####    But for some reason it shows up at a random point when I checked that on the scope).
# ####
# class Tutorial(EnvExperiment):
#
#     def build(self):
#
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#         self.setattr_device("ttl6")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#         self.ttl6.output()
#
#         delay(1 * ms)
#
#         time1=now_mu()
#         t_end = self.ttl2.gate_rising(60 * ms)  # opens gate for rising edges to be detected on TTL2 for 60ms
#         # sets variable t_end as time(in MUs) at which detection stops
#         time2=now_mu()
#
#         t_edge = self.ttl2.timestamp_mu(t_end)  # sets variable t_edge as time(in MUs) at which first edge is detected
#         # if no edge is detected, sets t_edge to -1
#         time3 = now_mu()
#
#         if t_edge > 0:  # runs if an edge has been detected
#             time4 = now_mu()
#             at_mu(t_edge)  # set time cursor to position of edge
#             time5 = now_mu()
#
#             delay(5 * us)  # 5us delay, to prevent underflow
#             time6 = now_mu()
#             self.ttl6.pulse(20 * us)  # outputs 5ms pulse on TTL6
#             time7 = now_mu()
#
#             print('t2-t1=',10**6 * self.core.mu_to_seconds(time2 - time1))
#             print('t3-t2=',10**6 * self.core.mu_to_seconds(time3 - time2))
#             print('t_end-t2=', 10**6 * self.core.mu_to_seconds(t_end - time2))
#             print('t_edge-t1=',10**6 * self.core.mu_to_seconds(t_edge - time1))
#             print('t4-t3=',10**6 * self.core.mu_to_seconds(time4 - time3))
#             print('t5-t4=', 10 ** 6 * self.core.mu_to_seconds(time5 - time4))
#             print('t6-t5=', 10 ** 6 * self.core.mu_to_seconds(time6 - time5))
#             print('t7-t6=', 10 ** 6 * self.core.mu_to_seconds(time7 - time6))
#
#         # self.ttl2.count(t_end)  #### discard remaining edges and close gate (I am not sure what this does in terms
#         #### of closing the gate and if it is necessary. Below, I realize that it makes a big difference when I want
#         #### to use a loop and time tag many events in a loop).
#         #### Since the "print"s above consume time, you need to disable the "print"s before enabeling this line.
#         # print(t_edge)
#         print("core done!")







######  I used a AWG to generate TTL pulses sent to TTL2. This code counts or time tag the recieved pulses. It worked.
class Tutorial(EnvExperiment):

    def build(self):

        self.setattr_device("core")
        self.setattr_device("ttl1")
        self.setattr_device("ttl2")

    @kernel
    def run(self):
        self.core.reset()
        self.ttl1.input()
        self.ttl2.input()

        delay(1 * us)

        #### This works. But a shorter version is a few lines below
        # t_end = self.ttl2.gate_rising(100 * ms)
        # numbPuls = self.ttl2.count(t_end)

        #### Instead of the two lines above, the manual says you can use the following.
        #### But it does not count when I tried:
        # numbPuls = self.ttl2.count(now_mu()+10000000000) ### This does not work!

        #### Instead of the two lines above, the manual says you can use the following.
        #### It works:
        # numbPuls1 = self.ttl1.count(self.ttl1.gate_rising(10 * ms))

        #### However, i cannot read the counts from two ttls at the same time in parallel. This does not work:
        # with parallel:
        #     numbPuls1 = self.ttl1.count(self.ttl1.gate_rising(100 * ms))
        #     numbPuls2 = self.ttl2.count(self.ttl2.gate_rising(100 * ms))


        #### The following seems to work. However, I need to check the precision and if it misses any signals:
        with parallel:
            tend1 = self.ttl1.gate_rising(50 * ms)
            tend2 = self.ttl2.gate_rising(50 * ms)

        # numbPuls1 = self.ttl1.count(tend1)
        # numbPuls2 = self.ttl2.count(tend2)

        timePuls1 = self.ttl1.timestamp_mu(tend1)
        timePuls2 = self.ttl2.timestamp_mu(tend2)
        dt = self.core.mu_to_seconds(timePuls2 - timePuls1)

        # delay(10 * us)
        # timePuls = self.ttl2.timestamp_mu(self.ttl2.gate_rising(2 * ms))


        # print(numbPuls1)
        # print(numbPuls2)
        print(timePuls1)
        print(timePuls2)
        print(10**6 * dt)
        print("core done!")







######  I used a AWG to generate TTL pulses sent to TTL2. This code time tags a few of the received TTL pulses
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         # numbPuls = self.ttl2.count(self.ttl2.gate_rising(10 * ms))
#         # delay(10 * us)
#         timePuls1 = self.ttl2.timestamp_mu(self.ttl2.gate_rising(2 * ms))
#         ### It is weird that adding a delay here does not do anything!
#         delay(10 * s)
#         ### It is also weird that printing the values here gives an error unless you add a long delay
#         ### (which does not affect the performance at all!!)
#         print(timePuls1)
#         timePuls2 = self.ttl2.timestamp_mu(self.ttl2.gate_rising(2 * ms))
#
#         delay(1 * ms)
#
#         # print(timePuls1)
#         print(timePuls2)
#         dt = self.core.mu_to_seconds(timePuls2 - timePuls1)
#         print(dt)
#
#         # print(numbPuls)
#         # print(timePuls)
#         print("core done!")







####   The code below helps to understand the behaviour of the timeline cursor in a loop. I am using an AWG to
####   generate TTL pulses of 100us period sent to ttl2. In the loop, I add a delay(1000*ms), then open the gate for
####   10ms to detect a rising edge, time tag the edge, and generate a 20ms pulse at ttl6. The time is recoreded
####   at different points and saved on a file on desktop. This is what happens to the time cursor:
####   In the first run of the loop, the program waits for 1000ms (the first delay) because the cursor is
####   placed at 1000ms later than "now". Then the time cursor is placed at 10ms later. However, the program
####   does not wait for extra 10ms (you can see that by increasing this time to 3s, for example); it unblocks
####   the core as soon as it detects a rising edge. It time stamps the edge, and generate the ttl6 pulse that
####   I can see on the scope. Now, goes to the 2nd loop. However, remember that the time cursor is at 1000ms later
####   which has not come yet. The code moves the time cursor to another 1000ms later and another 10ms later. However,
####   it does not generate any other ttl6 pulses because it is going to generate the ttl6 pulse at the position
####   of the time cursor which is a long time later. The code does not wait until that time arrives. It goes through
####   all the 10 loops very quickly (it moves the time cursor in each loop by an extra 1010ms) and finishes the code
####   without geenrating any more ttl6 pulses. So, ttl6 is generated only once (in the first loop).
####     A few other interesting points:
####       * It time tags all the pulses of ttl2 in all 10 loops, although the time cursor is several seconds
####         in future. This is because we never close the gate. So, it is open for 10ms which is long
####         enough to recieve ten ttl2 pulses.
####       * If the gate time is reduced from 10ms to 0.5ms, then it recieves 5 ttl2 pulses during the 0.5ms
####         gate openning time (in the first 5 loops), then in the 6th loop, since the gate is closed, it waits
####         for 1000ms (the requested delay in the loop), it opens the gate again and detects 5 more ttl2
####         pulses. In this case, two ttl6 pulses are generated; one in the first loop, and another in the
####         6th loop. So, I conclude here that, when the gate is open, it goes through the loops and ignores
####         the delay(). If the gate is closed, it waits when it sees a delay. In either case, it never waits
####         to get to the position of the time cursor.
####       * If the gate openning time is 10ms, a ttl6 pulse gets generated in the first loop. However, if the
####         gate openning time is 100ms, a ttl6 never gets negerated! So, it is not patient enough to wait more
####         than 50-70ms!! If a gate openning time of 60ms is chosen, the ttl6 remains on high level forever
####         (until the next time you generate a ttl6 pulse); it is not even patient enough to turn it off.
####         It turns the ttl6 ON and exits the code!!
####       * When a print is requested in the loop, it issues an overflow error after a couple of loops
####         (unless a short gate time of 1ms is used). I am not sure why! I guess because it is because
####         of the way print interacts with the time cursor; maybe it has to happen at the position of the time
####         cursor or something like that.
# TimeList = []
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#         self.setattr_device("ttl6")
#
#
#     @rpc(flags={"async"})
#     def saveTime(self, x):
#         TimeList.append(x)
#         TimeData = pd.DataFrame(TimeList)
#         TimeData.to_csv('~/Desktop/AkSTimeList.csv')
#         # print(TimeList)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#         self.ttl6.output()
#
#         delay(1 * us)
#         for x in range(10):
#             time1 = now_mu()
#             self.saveTime(self.core.mu_to_seconds(time1))
#
#             delay(1000 * ms)
#
#             self.saveTime(self.core.mu_to_seconds(now_mu()))
#             t_end = self.ttl2.gate_rising(10 * ms)
#             time2 = now_mu()
#
#             timePuls = self.ttl2.timestamp_mu(t_end)
#             # self.saveTime(self.core.mu_to_seconds(timePuls))
#             self.ttl6.pulse(20 * ms)
#             # print('t2=', time2)
#
#         print("core done!")







####   To further understand the timming in a loop, here I am using the same code as above, but with core.reset()
####   and core.count() to see how they affect the timming. Using either after detecting a rising edge, closes the
####   gate. So, the loop runs as expected with 1000ms delay in each (so it takes about 10s to finish the code).
####   Some details:
####     reset() closes the gate and I think empties the FIFO list. After detecting the rising edge, reset() moves
####             back the time cursor to about 5ms after the rising edge. This time is not only large, but also
####             uncertain. Sometimes it moves the cursor to 4.99ms after the rising edge, sometimes to 5.34ms
####             after the rising edge.
####     count() closes the gate too. But moves the cursor to about 9.935ms after the detected rising edge.
####             Sometimes it puts the cursor at 9.925ms after the edge. So, not completely reliable.
####     Changing the gate openning time from 10ms to 1ms, reduces the lag time of count() from ~10ms to ~1ms.
####             But it has no effect for the lag time of the reset() method. So, reset() moves the cursor always
####             to ~5ms after the rising edge. However, when using count(), the code waits until the end of the
####             gate openning time arrives (because it wants to count the pulses). So, if a delay of 1000ms is
####             used in the loop, and the gate is open for another 1000ms, the code will take about 20s to run
####             10 loops.

# TimeList = []
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#         self.setattr_device("ttl6")
#
#     ##### This is super slow. So, I will use open file (later below), instead of pandas.
#     @rpc(flags={"async"})
#     def saveTime(self, x):
#         TimeList.append(x)
#         TimeData = pd.DataFrame(TimeList)
#         TimeData.to_csv('~/Desktop/AkSTimeList.csv')
#         # print(TimeList)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#         self.ttl6.output()
#         self.saveTime(self.core.mu_to_seconds(now_mu()))
#
#         delay(1 * us)
#         for x in range(10):
#             self.saveTime(self.core.mu_to_seconds(now_mu()))  #1
#             delay(1000 * ms)
#             self.saveTime(self.core.mu_to_seconds(now_mu()))  #2
#             t_end = self.ttl2.gate_rising(10 * ms)
#             timePuls = self.ttl2.timestamp_mu(t_end)
#             self.saveTime(self.core.mu_to_seconds(now_mu()))  #3
#             self.saveTime(self.core.mu_to_seconds(timePuls))  #4
#
#             # self.core.reset()
#             self.ttl2.count(t_end)
#             self.saveTime(self.core.mu_to_seconds(now_mu()))  #5
#
#         print("core done!")









####  The same as above, but in a for loop and save the times in a list on the desktop: I used a AWG to generate
####  TTL pulses sent to ttl2. This code time tags a few of the received TTL pulses in a "for loop". It works, but
####  for a limited number of pulses. The AWG generates ttl at 10kHz, and the gate_rising is 100ms.
####  If I  increase the number of loops to more than 1024, I get overflow error. If I decrease the ttl rate to 1kHz,
####  or decrease the gate_rising to 10ms, then I can run until 1034 loops without overflow error.
####  With ttl generated at 100kHz rate, I can run only 54 loops.
# St_t = time.time()
# TimeList = []
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def saveTime(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         for x in range(1024):
#             t_end = self.ttl2.gate_rising(100 * ms)
#             timePuls = self.ttl2.timestamp_mu(t_end)
#             self.saveTime(self.core.mu_to_seconds(timePuls)*10**3)
#
#         file1.close()
#
#     def analyze(self):
#         End_t = time.time()
#         print("$$$ %s seconds $$$" % round((End_t - St_t),3))







#### A stupidly easy fix to the problem above: Instead of openning the gate in the loop (which moves the cursor further
####    back each time), open the gate before the loop. Then, I am able to open the gate for as long as I want and
####    time tag as many signal as I want and save them on a file very fast.
#### If the ttl signal is generated at 5kHz rate, I can monitor the channel, time tag, and save on desktop for as long
####    as I want (I tested and recorded for 100s). With a ttl signal at 10kHz rate, I can time tag for about 7s, i.e.
####    70,000 signals before it raises overflow error. So, if I want to use this to time tag very fast signals, I
####    have to reset() the memory every few seconds (which is fine for the single-photon generation experiment).
####    The reason for overflow error is that the CPU reads the timetags from the cache at a certain rate at maximum.
####    When the input ttl rate is faster that the CPU reading rate, data starts to accumulate in the cache. Then,
####    after some time the cache (which can hold about 64, maybe 50 time tags) is full and overflow occurs.
# St_t = time.time()
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def saveTime(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         t_end = self.ttl2.gate_rising(10 * s)
#         for x in range(1000):
#             timePuls = self.ttl2.timestamp_mu(t_end)
#             self.saveTime(timePuls)
#
#
#         file1.close()
#
#     def analyze(self):
#         End_t = time.time()
#         print("$$$ %s seconds $$$" % round((End_t - St_t),3))





#### Obviously, to time tag for a longer time, we can use two loops and reset() the buffer in the outer loop which
####    runs, say, every 1s. So, we time tag the photons every 1s, for example. This is absolutely fine, as we
####    want to have a refresh rate for the coincidence measurement anyway.
# St_t = time.time()
# TimeList = []
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList.txt', 'a')
#
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def saveTime(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         for x in range(10):
#             t_end = self.ttl2.gate_rising(1 * s)
#             for x in range(5000):
#                 timePuls = self.ttl2.timestamp_mu(t_end)
#                 self.saveTime(timePuls)
#             self.core.reset()
#
#         file1.close()
#
#     def analyze(self):
#         End_t = time.time()
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))







#### The following code shows that when two "for loops" are used in parallel, artiq finishes the first
#### loop completely and then goes to run the 2nd loop. So, we should not do it this way to time tag
#### the two ttl signals. An AWG was used to generate the signals.
# TimeList1 = []
# TimeList2 = []
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         TimeList1.append(x)
#         TimeData1 = pd.DataFrame(TimeList1)
#         TimeData1.to_csv('~/Desktop/AkSTimeList1.csv')
#         # print(TimeList1)
#
#     def saveTime2(self, x):
#         TimeList2.append(x)
#         TimeData2 = pd.DataFrame(TimeList2)
#         TimeData2.to_csv('~/Desktop/AkSTimeList2.csv')
#         # print(TimeList2)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         with parallel:
#             for i in range(5):
#                 with sequential:
#                     tend1 = self.ttl1.gate_rising(100 * ms)
#                     timePuls1 = self.ttl1.timestamp_mu(tend1)
#                     self.saveTime1(timePuls1)
#                     self.core.reset()
#                     delay(500 * ms)
#                     print('i=', i)
#             for j in range(5):
#                 with sequential:
#                     self.core.reset()
#                     tend2 = self.ttl2.gate_rising(100 * ms)
#                     timePuls2 = self.ttl2.timestamp_mu(tend2)
#                     # self.saveTime2(timePuls2)
#                     delay(500 * ms)
#                     print('j=', j)
#
#         print("core DONE")







#### Sending 1 Hz at ttl1 and x2 Hz at ttl2.
#### The part without loop shows that:
####    Artiq waits to recieves both signals and then read them both from the cache at the same time. In the meanwhile,
####    it saves all the recieved timetags in the cache (which apparently can keep only 64 timetags). So, if x2 is
####    smaller than 64 Hz, it works fine. However, if x2 is larger than 64, depending when it recieves the two
####    signals relative to each other, it may issue overflow error; if it has just missed the previous ttl1 signal,
####    it has to wait for almost 1s to get the next ttl1 signal and in the meantime tag all the ttl2 signals
####    and keep them in cache.
#### The part with loop:
####    If the loop runs 4 times, and x2 is 10Hz, then, this saves 4 signals from ttl1 and
####    4 signals from ttl2. Since ttl2 receives signals much faster, the core keeps the time tags in cache
####    until it receives a signal from ttl1 and then save all of them together. In other words, it doesn't
####    miss ttl4 signals while waiting for ttl1; it saves the first 4 signals from ttl1 and the first 4
####    signals from ttl2 that it recieves after openning the gate. If we want to run the loop 10 times, for example,
####    it means it has to keep about 100 timetags from ttl2 while waiting to get 10 ttl1 signals, and thus,
####    runs out of cache and give overflow error.
# St_t = time.time()
# TimeList = []
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList1.txt', 'a')
# file2 = open('/home/fiber-cavity/Desktop/AkSTimeList2.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime2(self, x):
#         file2.write(str(x) + '\n')
#         # print(x)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#         time1 = now_mu()
#
#         with parallel:
#             tend1 = self.ttl1.gate_rising(10 * s)
#             tend2 = self.ttl2.gate_rising(10 * s)
#
#
#         timePuls1 = self.ttl1.timestamp_mu(tend1)
#         timePuls2 = self.ttl2.timestamp_mu(tend2)
#         self.saveTime1(10 ** 3 * self.core.mu_to_seconds(time1))
#         self.saveTime1(10 ** 3 * self.core.mu_to_seconds(timePuls1))
#         self.saveTime2(10 ** 3 * self.core.mu_to_seconds(time1))
#         self.saveTime2(10 ** 3 * self.core.mu_to_seconds(timePuls2))
#
#         # for x in range(4):
#         #     timePuls1 = self.ttl1.timestamp_mu(tend1)
#         #     timePuls2 = self.ttl2.timestamp_mu(tend2)
#         #
#         #     self.saveTime1(10**3 * self.core.mu_to_seconds(timePuls1))
#         #     self.saveTime2(10**3 * self.core.mu_to_seconds(timePuls2))
#
#         file1.close()
#         # file2.close()
#
#     def analyze(self):
#         End_t = time.time()
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))









#### If the signals sent to ttl1 and ttl2 are different by "delta", then the cache gets full in about 64/delta,
####    because pulses are kept in the cache until they are read by the CPU, which apparently, reads the two ttl
####    signals one by one (one from ttl1 and one from tt2). If there are, for example, 50 ttl2 signals but no ttl1
####    signals, the core keep all these 50 signals in the cache until it recieves a signal from ttl1. Then,
####    timestamps and releases the first ttl2 signal that is in the cache. Apparently, the cache can keep only
####    64 timestamps.
#### Again an AWG sends 10kHz and 9kHz to ttl1 and ttl2, respectively. With this code I learned that when the
####    gate-open-time ends, i.e. when the gate closes, the core resets the cache. To see this, I am using
####    two loops. The first loop runs for 60ms and saves 601 timestamps. If the gate-open-time is 70ms instead
####    of 60ms, an overflow error occurs because the cache gets full. So, the first loop almost fill up the
####    cache, but not completely. Moreover, the gate time is chosen precisely such that the gate gets closed
####    after the loop. Then, I start another loop but without reset() in between. If the cache remains full,
####    the second loop should not run and should give overflow error. However, it works fine. But I need to
####    use a delay of more than 5ms. This is another indication that closing the gate uses reset() automatically
####    as it requires a 5ms delay, at least.
# St_t = time.time()
# TimeList = []
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList1.txt', 'a')
# file2 = open('/home/fiber-cavity/Desktop/AkSTimeList2.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime2(self, x):
#         file2.write(str(x) + '\n')
#         # print(x)
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#         time1 = now_mu()
#
#         with parallel:
#             tend1 = self.ttl1.gate_rising(60 * ms)
#             tend2 = self.ttl2.gate_rising(60 * ms)
#
#         tt=0
#         for x in range(601):
#             timePuls1 = self.ttl1.timestamp_mu(tend1)
#             timePuls2 = self.ttl2.timestamp_mu(tend2)
#
#             self.saveTime1(10**3 * self.core.mu_to_seconds(timePuls1))
#             self.saveTime2(10**3 * self.core.mu_to_seconds(timePuls2))
#
#             if timePuls1 > 0:
#                 tt = timePuls1
#
#         at_mu(tt)
#         delay(6*ms)
#         with parallel:
#             tend1 = self.ttl1.gate_rising(60 * ms)
#             tend2 = self.ttl2.gate_rising(60 * ms)
#
#
#         for x in range(600):
#             timePuls1 = self.ttl1.timestamp_mu(tend1)
#             timePuls2 = self.ttl2.timestamp_mu(tend2)
#
#             self.saveTime1(10**3 * self.core.mu_to_seconds(timePuls1))
#             self.saveTime2(10**3 * self.core.mu_to_seconds(timePuls2))
#
#         file1.close()
#         # file2.close()
#
#     def analyze(self):
#         End_t = time.time()
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))







#### A simple code just to count the signals in a time interval from two ttls simultaneouslly. It does not work
####    to count the signals in parallel!!
####    So, I am counting the signals in series (see the last loop):
# St_t = time.time()
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#
# # #####   One of the channels has to have a short gate time! otherwise overflows:
# #         t1 = now_mu()
# #         with parallel:
# #             tend1 = self.ttl1.gate_rising(100 * ms)
# #             tend2 = self.ttl2.gate_rising(7 * ms)
# #
# #         t2 = now_mu()
# #         count1 = self.ttl1.count(tend1)
# #         count2 = self.ttl2.count(tend2)
# #
# #         print(count1, count2)
# #         print(self.core.mu_to_seconds(t2 - t1))
# #         print(self.core.mu_to_seconds(tend1 - t1))
# #         print(self.core.mu_to_seconds(tend2 - t1))
#
#
# # #####   not working (underflow error):
# #         with parallel:
# #             with sequential:
# #                 tend1 = self.ttl1.gate_rising(5 * ms)
# #                 count1 = self.ttl1.count(tend1)
# #                 delay(10*ms)
# #             with sequential:
# #                 tend2 = self.ttl2.gate_rising(50 * ms)
# #                 count2 = self.ttl2.count(tend2)
#
#
# ##### Instead of counting the two channels in parallel, I am counting in series:
#         for x in range(100):
#             tend1 = self.ttl1.gate_rising(1000 * ms)
#             count1 = self.ttl1.count(tend1)
#
#             delay(5 * ms)
#
#             tend2 = self.ttl2.gate_rising(10 * ms)
#             count2 = self.ttl2.count(tend2)
#
#             print(count1,count2)
#             delay(10*ms)
#
#
#
#     def analyze(self):
#         End_t = time.time()
#         # self.saveTime2(TimeList)
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))











#### Here I am using two loops. The inner loop runs 500 times and timetags 500 pulses from each channel. The outer
####    loop runs 1000 times. I don't have to use the outer loop if the two rates are similar and we don't
####    run out of cache. However, if the two rates for ttl1 and tt2 are different, we have to reset()
####    every while (when the cache is full). That is why I am using two loops.
#### I tested which one is faster: saving every single timetag in the inner loop on a file, or saving them in a
####    list and save the list later on in the outer loop. I realized that there is no a big difference. When no
####    RPC is called, i.e. the timetags are not saved in any way, the reset() function is very reliable and moves
####    the time cursor always less than 5ms (between 4.6 to 4.9 ms). However, when I use one RPC in the inner loop
####    to either save the tags in a file or in a list, the reset() moves the cursor often by 5ms, but sometimes by
####    200ms. If I use two RPC calls to save both channels in either a file or a list, the reset() moves the cursor
####    by 200ms more often.
# St_t = time.time()
# TimeList1 = []
# TimeList2 = []
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList1.txt', 'a')
# file2 = open('/home/fiber-cavity/Desktop/AkSTimeList2.txt', 'a')
# file3 = open('/home/fiber-cavity/Desktop/AkSTimeList3.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def AppendList1(self, x):
#         TimeList1.append(x)
#
#     @rpc(flags={"async"})
#     def AppendList2(self, x):
#         TimeList2.append(x)
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime2(self, x):
#         file2.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime3(self, x):
#         file3.write(str(x) + '\n')
#         # print(x)
#
#
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         tt=0
#         t1 = 0
#         t2 = 0
#         for x in range(100):
#             with parallel:
#                 tend1 = self.ttl1.gate_rising(600 * ms)
#                 tend2 = self.ttl2.gate_rising(600 * ms)
#             for y in range(500):
#                 timePuls1 = self.ttl1.timestamp_mu(tend1)
#                 timePuls2 = self.ttl2.timestamp_mu(tend2)
#
#                 self.saveTime1(10 ** 3 * self.core.mu_to_seconds(timePuls1))
#                 self.saveTime2(10 ** 3 * self.core.mu_to_seconds(timePuls2))
#
#                 # self.AppendList1(10 ** 3 * self.core.mu_to_seconds(timePuls1))
#                 # self.AppendList2(10 ** 3 * self.core.mu_to_seconds(timePuls2))
#                 if timePuls1 > 0:
#                     tt = timePuls1
#             # t1 = now_mu()
#             self.core.reset()
#             # at_mu(tt)
#             # delay(20 * ms)
#             t2 = now_mu()
#             self.saveTime3(10 ** 3 * self.core.mu_to_seconds(t2 - tt))
#
#
#         file1.close()
#         # file2.close()
#
#         # tt = 0
#         # with parallel:
#         #     tend1 = self.ttl1.gate_rising(60 * ms)
#         #     tend2 = self.ttl2.gate_rising(60 * ms)
#         # for y in range(500):
#         #     timePuls1 = self.ttl1.timestamp_mu(tend1)
#         #     timePuls2 = self.ttl2.timestamp_mu(tend2)
#         #
#         #     # self.saveTime1(10 ** 3 * self.core.mu_to_seconds(timePuls1))
#         #     # self.saveTime2(10 ** 3 * self.core.mu_to_seconds(timePuls2))
#         #     if timePuls1 > 0:
#         #         tt = timePuls1
#         # self.core.reset()
#         # at_mu(tt)
#         # delay(10 * ms)
#         #
#         # #
#         # tend1 = self.ttl1.gate_rising(60 * ms)
#         # # with parallel:
#         # #     tend1 = self.ttl1.gate_rising(60 * ms)
#         # #     tend2 = self.ttl2.gate_rising(60 * ms)
#         # # for y in range(500):
#         # #     timePuls1 = self.ttl1.timestamp_mu(tend1)
#         # #     timePuls2 = self.ttl2.timestamp_mu(tend2)
#
#
#
#
#
#     def analyze(self):
#         End_t = time.time()
#         # self.saveTime2(TimeList)
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))












####    I am using AWG to generate TTLs at about 1 MHz rate sent to ttl1 and ttl2. Then, I am openning the gates
####    in parallel for 100ns and timestamp the two signals, save them in two files. It works well, even if I
####    run the loop 1 million times. I use a delay of 1ms representing the atom tomography. If the delay is
####    smaller than 300 us, I get overflow error. The code works well even if I use two different rates
####    for the two channels, say 1 MHz at ttl1 and 0.1 MHz for ttl2.
#### I think this is closer to the experimental condition of atom-photon entanglemen: we want to timestamp the two
####    ttls every 1-5 ms (limited mainly by the atom tomography and cooling cycle). In each cycle, the gate is open
####    for about 100ns to collect the photons after the excitation.
# St_t = time.time()
# file1 = open('/home/fiber-cavity/Desktop/AkSTimeList1.txt', 'a')
# file2 = open('/home/fiber-cavity/Desktop/AkSTimeList2.txt', 'a')
# class Tutorial(EnvExperiment):
#
#     def build(self):
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#
#     @rpc(flags={"async"})
#     def saveTime1(self, x):
#         file1.write(str(x) + '\n')
#         # print(x)
#
#     @rpc(flags={"async"})
#     def saveTime2(self, x):
#         file2.write(str(x) + '\n')
#         # print(x)
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#
#         delay(1 * us)
#
#         for x in range(10000):
#             with parallel:
#                 tend1 = self.ttl1.gate_rising(100 * ns)
#                 tend2 = self.ttl2.gate_rising(100 * ns)
#
#             timePuls1 = self.ttl1.timestamp_mu(tend1)
#             timePuls2 = self.ttl2.timestamp_mu(tend2)
#
#             self.saveTime1(10 ** 3 * self.core.mu_to_seconds(timePuls1))
#             self.saveTime2(10 ** 3 * self.core.mu_to_seconds(timePuls2))
#             delay(1 * ms)
#
#
#         file1.close()
#         # file2.close()
#
#
#     def analyze(self):
#         End_t = time.time()
#         # self.saveTime2(TimeList)
#         print("$$$ %s seconds $$$" % round((End_t - St_t), 3))












#### Using one TTL as output and another as input in parallel to generate and detect the
#### TTLs all with artiq, i.e. no AWG. It works.
# class Tutorial(EnvExperiment):
#
#     def build(self):
#
#         self.setattr_device("core")
#         self.setattr_device("ttl1")
#         self.setattr_device("ttl2")
#         self.setattr_device("ttl5")
#         self.setattr_device("ttl6")
#
#     @kernel
#     def run(self):
#         self.core.reset()
#         self.ttl1.input()
#         self.ttl2.input()
#         self.ttl5.output()
#         self.ttl6.output()
#
#
#         # #### The following code works. It generates a signal on ttl 5 and measures the rising time on ttl1.
#         # #### However, there are some delicate delays required:
#         # #### 1st, there must be a delay of minimum 8ns before doing anything (at the very begining before parallel).
#         # ####    Otherwise, the signal wont get detected.
#         # #### 2nd, no other delays are required (in the sequences for example). But if using any delays, the delay in
#         # ####    the end sequence has to be no more than 70ns longer than the first delay. Otherwise, the signal wont
#         # ####    get detected. The 2nd delay could be zero.
#         # #### 3rd, the gating sequence has to come after the pulse generation sequence even when there is no delay in
#         # ####    those sequences! If the two sequentials are swapped, an undeflow error will be raised.
#         # #### In any case, there is always a constant time difference of either 198ns or 199ns between the generation
#         # ####    and detection of the ttl pulse as seen in the results below.
#         # delay(1 * us)
#         # with parallel:
#         #     with sequential:
#         #         # delay(100 * ns)
#         #         Time1 = now_mu()
#         #         self.ttl5.pulse(2 * us)
#         #     with sequential:
#         #         # delay(170 * ns)
#         #         TTLtime = self.ttl1.timestamp_mu(self.ttl1.gate_rising(10 * us))
#         # print(TTLtime)
#         # print('time delay =', 10**9 * self.core.mu_to_seconds(TTLtime - Time1), 'ns')
#
#         # #### The following will do the same as above (in case you have to separate the gate openning from the timestamp)
#         # delay(1 * us)
#         # with parallel:
#         #     with sequential:
#         #         delay(100 * ns)
#         #         Time1 = now_mu()
#         #         self.ttl5.pulse(2 * us)
#         #     with sequential:
#         #         delay(170 * ns)
#         #         tend = self.ttl1.gate_rising(5000 * ms)
#         # TTLtime = self.ttl1.timestamp_mu(tend)
#         # print(TTLtime)
#         # print('time delay =', 10**9 * self.core.mu_to_seconds(TTLtime - Time1), 'ns')
#         #
#         # print("core done!")
#
#
#         #### Generating two signals from the two TTLs and detecting the timing. If two TTL signals are generated
#         #### at the same time, it detects the 2nd one about 3ns later than the 1st TTL. This delay changes from 2ns to
#         #### 4ns (on average 3ns). It works.
#         delay(1 * us)
#         with parallel:
#             with parallel:
#                 self.ttl5.pulse(50 * ns)
#                 with sequential:
#                     delay(5 * ns)
#                     self.ttl6.pulse(50 * ns)
#
#             tend1 = self.ttl1.gate_rising(50 * ms)
#             tend2 = self.ttl2.gate_rising(50 * ms)
#
#         TTLtime1 = self.ttl1.timestamp_mu(tend1)
#         TTLtime2 = self.ttl2.timestamp_mu(tend2)
#         dt = self.core.mu_to_seconds(TTLtime2 - TTLtime1)
#         print('time delay =', 10 ** 9 * dt, 'ns')
#
#
#         print("CORE DONE")