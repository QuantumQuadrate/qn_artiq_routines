"""
This code allows for tuning the coils with the homemade potentiometer box
by reading its output into the Sampler and outputting a voltage
from a corresponding Zotino channel. In this way the MOT can be optimized
manually and then the ARTIQ variables for the coil voltages can be updated
automatically. 
"""

from artiq.experiment import *

from utilities.BaseExperiment import BaseExperiment

class SamplerMOTCoilTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("FORT_AOM_on", BooleanValue(True))
        self.setattr_argument("MOT_AOMs_on", BooleanValue(True))
        self.setattr_argument("update_coil_volts_at_finish", BooleanValue(False))
        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("coil_volts_multiplier", NumberValue(1)) # scales the value read by the Sampler

        group = "SPCM settings"
         # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(300 * ms, unit='ms'), group)
        # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1 * 10 ** 5), group)  # saturation limit in counts/dt.
        self.setattr_argument("print_count_rate", BooleanValue(True), group)

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
        self.sampler_buffer = [0.0]*8
        self.control_volts_channels = [0,1,2,3] # the sampler channels to read

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

        delay(10 * ms)

        # for outputting a voltage proportional to counts read from the SPCM on ttl0
        ch = 6
        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        Satdt = self.sat1s * self.dt_exposure  # saturation limit in counts/dt.
        delay(1000 * ms)

        print("ready!")

        for i in range(self.n_steps):

            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            if self.print_count_rate:
                print(round(count1 / self.dt_exposure))
            delay(10 * ms)
            volt1 = count1 * 5 / Satdt  # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
            self.zotino0.write_dac(ch, volt1)
            self.zotino0.load()
            delay(1 * ms)

            self.sampler1.sample(self.sampler_buffer)
            control_volts = [self.sampler_buffer[ch]*self.coil_volts_multiplier for ch in self.control_volts_channels]

            delay(1*ms)

            # set coils based on the sampler values we read
            self.zotino0.set_dac(
                control_volts,
                channels=self.coil_channels)

            delay(1 * ms)

        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        if self.update_coil_volts_at_finish:
            volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
            for i in range(4):
                self.set_dataset(volt_datasets[i], control_volts[i], broadcast=True)

        print("Experiment finished.")
