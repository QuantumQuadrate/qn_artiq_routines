"""
This code allows for tuning the coils with the homemade potentiometer box
by reading its output into the Sampler and outputting a voltage
from a corresponding Zotino channel. In this way the MOT can be optimized
manually and then the ARTIQ variables for the coil voltages can be updated
automatically.
"""

from artiq.experiment import *
import numpy as np

from utilities.BaseExperiment import BaseExperiment

class SamplerMOTCoilTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("FORT_AOM_on", BooleanValue(False))
        self.setattr_argument("MOT_AOMs_on", BooleanValue(True))
        self.setattr_argument("set_current_coil_volts_at_finish", BooleanValue(False))
        self.setattr_argument("set_best_coil_volts_at_finish", BooleanValue(False)) # this will override the previous value
        self.setattr_argument("leave_coils_on_at_finish", BooleanValue(True))
        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("coil_volts_multiplier",
                              NumberValue(3.3)) # scales the value read by the Sampler
        self.setattr_argument("differential_mode",
                              BooleanValue(False),"differential mode (tune voltage wrt current coil settings)") # scan the coils with respect to the current settings
        self.setattr_argument("differential_multiplier",
                              NumberValue(1),"differential mode (tune voltage wrt current coil settings)") # scales the value read by the Sampler

        group = "SPCM settings"
         # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(300 * ms, unit='ms'), group)
        # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1 * 10 ** 5), group)  # saturation limit in counts/dt.
        self.setattr_argument("print_count_rate", BooleanValue(True), group)

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(200), "Laser feedback")
        self.setattr_argument("enable_laser_feedback", BooleanValue(True), "Laser feedback")

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        print(self.dt_exposure)
        self.n_steps = int(60*self.run_time_minutes/self.dt_exposure+0.5)
        print(self.n_steps)
        self.sampler_buffer = np.zeros(8)
        self.control_volts_channels = [0,1,2,3] # the sampler channels to read
        self.default_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]

        self.count_rate_dataset = 'SPCM_count_rate_Hz'
        self.set_dataset(self.count_rate_dataset,
                             [0.0],
                             broadcast=True)

        self.maxcount_dataset = 'SPCM_max_counts_and_volts'
        self.set_dataset(self.count_rate_dataset,
                         [0.0, *self.default_volts],
                         broadcast=True)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        if self.MOT_AOMs_on:
            # Turn on AOMs to load the MOT.
            self.dds_cooling_DP.sw.on()
            self.dds_cooling_SP.sw.on()
            self.dds_MOT_RP.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
        else:
            if self.MOT_AOMs_on:
                # Turn on AOMs to load the MOT.
                self.dds_cooling_DP.sw.off()
                self.dds_cooling_SP.sw.off()
                self.dds_MOT_RP.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A6.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
        if self.FORT_AOM_on:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()

            # wait for AOMs to thermalize
            delay(3000 * ms)

        control_volts = [0.0]*4
        best_volts = [0.0]*4 # the volts that we had when the best MOT was found
        max_count_rate = 0.0

        delay(10 * ms)

        # for outputting a voltage proportional to counts read from the SPCM on ttl0
        ch = 6
        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        Satdt = self.sat1s * self.dt_exposure  # saturation limit in counts/dt.
        delay(1 * ms)

        if self.enable_laser_feedback:
            # takes several iterations for process value to stabilize
            for i in range(20):
                print("running feedback step", i)
                self.laser_stabilizer.run()  # must come after relevant DDS's have been set
                delay(500 * ms)

        print("ready!")

        for i in range(self.n_steps):

            if self.enable_laser_feedback:
                if (i % self.AOM_feedback_period_cycles) == 0:
                    print("running feedback")
                    self.core.break_realtime()
                    self.laser_stabilizer.run()
                    delay(100*ms)

            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            count_rate_Hz = count1 / self.dt_exposure
            if self.print_count_rate:
                print(round(count_rate_Hz))
            delay(10 * ms)
            volt1 = count1 * 5 / Satdt  # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
            self.zotino0.write_dac(ch, volt1)
            self.zotino0.load()
            delay(1 * ms)

            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, count_rate_Hz)
            if count_rate_Hz > max_count_rate:
                max_count_rate = count_rate_Hz
                best_volts = control_volts
                self.set_dataset(self.maxcount_dataset,[count_rate_Hz]+best_volts)
            delay(1 * ms)

            self.sampler1.sample(self.sampler_buffer)
            if self.differential_mode:
                control_volts = [self.sampler_buffer[ch]*self.differential_multiplier + self.default_volts[ch]
                                 for ch in self.control_volts_channels]
            else:
                control_volts = [self.sampler_buffer[ch] * self.coil_volts_multiplier
                                 for ch in self.control_volts_channels]

            delay(1*ms)

            # set coils based on the sampler values we read
            self.zotino0.set_dac(
                control_volts,
                channels=self.coil_channels)

        if not self.leave_coils_on_at_finish:
            self.zotino0.write_dac(ch, 0.0)
            self.zotino0.load()
            delay(1 * ms)

        volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
        if self.set_best_coil_volts_at_finish:
            for i in range(4):
                self.set_dataset(volt_datasets[i], best_volts[i], broadcast=True)
        elif self.set_current_coil_volts_at_finish:
            for i in range(4):
                self.set_dataset(volt_datasets[i], control_volts[i], broadcast=True)

        print("Best volts [VZ_bottom,VZ_top,Vx,Vy]:")
        print(best_volts)
        print("Best count rate (kHz):")
        print(max_count_rate/1e3)


        print("Experiment finished.")
