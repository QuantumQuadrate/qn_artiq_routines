"""
This code allows for tuning the coils with the homemade potentiometer box
by reading its output into the Sampler and outputting a voltage
from a corresponding Zotino channel. In this way the MOT can be optimized
manually and then the ARTIQ variables for the coil voltages can be updated
automatically.
"""

from artiq.experiment import *
import numpy as np

import sys
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class SamplerCoilTuneFORTLoading(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("continues_loading", BooleanValue(True))
        self.setattr_argument("set_current_coil_volts_at_finish", BooleanValue(False))
        self.setattr_argument("leave_coils_on_at_finish", BooleanValue(True))
        self.setattr_argument("run_time_minutes", NumberValue(2))
        self.setattr_argument("coil_volts_multiplier",
                              NumberValue(3)) # scales the value read by the Sampler
        self.setattr_argument("differential_mode",
                              BooleanValue(True),"differential mode (tune voltage wrt current coil settings)") # scan the coils with respect to the current settings
        self.setattr_argument("differential_multiplier",
                              NumberValue(0.5),"differential mode (tune voltage wrt current coil settings)") # scales the value read by the Sampler

        group = "SPCM settings"
         # exposure time of the SPCM

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(20), "Laser feedback")
        self.setattr_argument("enable_laser_feedback", BooleanValue(True), "Laser feedback")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        self.n_steps = int(60*self.run_time_minutes/(
                self.t_MOT_loading+self.t_FORT_loading+0.5/self.AOM_feedback_period_cycles)+0.5)
        self.sampler_buffer = np.zeros(8)
        self.control_volts_channels = [0,1,2,3] # the sampler channels to read
        self.default_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT]

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # todo: these are going to be regularly used, so put these in the base experiment
        self.set_dataset("photocounts_per_s", [0], broadcast=True)

        self.set_dataset("photocount_bins", [50], broadcast=True)

        # turn off AOMs we aren't using, in case they were on previously
        self.dds_D1_pumping_SP.sw.off()
        self.dds_excitation.sw.off()

        # turn on cooling MOT AOMs
        self.dds_cooling_DP.sw.on()  # cooling double pass
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A5.sw.on()

        delay(2000 * ms)  # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1 * ms)

        counts = 0

        saturated_coils = [False]*4

        control_volts = [0.0]*4

        print("ready!")

        delay(10*ms)

        for i in range(self.n_steps):

            self.sampler1.sample(self.sampler_buffer)
            if self.differential_mode:
                control_volts = [self.sampler_buffer[ch]*self.differential_multiplier + self.default_volts[ch]
                                 for ch in self.control_volts_channels]
            else:
                control_volts = [self.sampler_buffer[ch] * self.coil_volts_multiplier
                                 for ch in self.control_volts_channels]

            delay(1*ms)

            self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

            ############################
            # load the MOT
            ############################
            self.zotino0.set_dac(
                control_volts,
                channels=self.coil_channels)
            self.dds_cooling_DP.sw.on()

            if self.continues_loading:
                if self.enable_laser_feedback:
                    if (i % (self.AOM_feedback_period_cycles*10)) == 0:
                        print("running feedback")
                        self.core.break_realtime()
                        self.laser_stabilizer.run()
                        delay(1 * ms)
                        self.dds_FORT.sw.on()
                else:
                    delay(1 * ms)
                delay(50 * ms)
                t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
                counts = self.ttl0.count(t_gate_end)

                delay(3 * ms)

                # update the datasets
                self.append_to_dataset('photocounts_per_s', counts / self.t_SPCM_first_shot)


                for j in range(4):
                    try:
                        self.zotino0.write_dac(self.coil_channels[j], control_volts[j])
                        if saturated_coils[j]:  # shouldn't get here unless we're no longer saturated
                            self.print_async("no longer saturated")
                            saturated_coils[j] = False
                    except ValueError:
                        if not saturated_coils[j]:
                            self.print_async("warning: voltage saturated for", self.coil_names[j])
                            saturated_coils[j] = True
                        sign = control_volts[j] / (control_volts[j] ** 2) ** (1 / 2)
                        self.zotino0.write_dac(self.coil_channels[j], sign * 9.9)
                    self.zotino0.load()

            else:
                if self.enable_laser_feedback:
                    if (i % self.AOM_feedback_period_cycles) == 0:
                        print("running feedback")
                        self.core.break_realtime()
                        self.laser_stabilizer.run()
                        delay(1 * ms)
                        self.dds_FORT.sw.on()
                else:
                    delay(1 * ms)
                # wait for the MOT to load
                delay_mu(self.t_MOT_loading_mu)

                # turn on the dipole trap and wait to load atoms
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
                delay_mu(self.t_FORT_loading_mu)

                # turn off the coils
                self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)

                delay(3 * ms)  # should wait for the MOT to dissipate

                # set the cooling DP AOM to the readout settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

                ############################
                # take the first shot
                ############################
                self.dds_cooling_DP.sw.on()
                t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
                counts = self.ttl0.count(t_gate_end)
                delay(1 * ms)
                self.dds_cooling_DP.sw.off()

                # delay to mimic a plausible real experiment
                delay(10 * ms)

                # effectively turn the FORT AOM off
                self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
                # set the cooling DP AOM to the MOT settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

                delay(2 * ms)

                # update the datasets
                self.append_to_dataset('photocounts_per_s', counts/self.t_SPCM_first_shot)

                delay(1*ms)

                for j in range(4):
                    try:
                        self.zotino0.write_dac(self.coil_channels[j], control_volts[j])
                        if saturated_coils[j]: # shouldn't get here unless we're no longer saturated
                            self.print_async("no longer saturated")
                            saturated_coils[j] = False
                    except ValueError:
                        if not saturated_coils[j]:
                            self.print_async("warning: voltage saturated for", self.coil_names[j])
                            saturated_coils[j] = True
                        sign = control_volts[j]/(control_volts[j]**2)**(1/2)
                        self.zotino0.write_dac(self.coil_channels[j], sign*9.9)
                    self.zotino0.load()

        if not self.leave_coils_on_at_finish:
            for ch in self.coil_channels:
                self.zotino0.write_dac(ch, 0.0)
                self.zotino0.load()
                delay(1 * ms)

        volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
        if self.set_current_coil_volts_at_finish:
            for i in range(4):
                self.set_dataset(volt_datasets[i], control_volts[i], broadcast=True, persist=True)

        print("Experiment finished.")
