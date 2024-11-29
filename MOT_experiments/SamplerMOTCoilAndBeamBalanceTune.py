"""
This code is a combination of SamplerMOTCoilTune and SamplerMOTBeamBalance,
as the name implies. With only one potentiometer box we can have both
functionalities by using an external switch to toggle between the two
during the same run of this experiment. At the time of writing this, this is
done using ttl14 set to high, and a switch to toggle the connection of ttl14
to ttl3, the state of which we read to determine the which variables to tune.
"""

from artiq.experiment import *
import numpy as np
import logging

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class SamplerMOTCoilAndBeamBalanceTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("FORT_AOM_on", BooleanValue(False))

        # this can be used to modulate the FORT while atom loading in nothing mode, to find the
        # trap frequencies by watching the atom loading go away
        self.setattr_argument("FORT_modulation_switch_on", BooleanValue(False))

        self.setattr_argument("set_current_parameters_at_finish", BooleanValue(False))
        self.setattr_argument("leave_coils_on_at_finish", BooleanValue(True))

        self.both_mode = "both"
        self.beam_mode = "beams only"
        self.coil_mode = "coils only"
        self.nothing_mode = "nothing (for debugging)"
        self.setattr_argument("what_to_tune", EnumerationValue([self.both_mode, self.coil_mode,
                                                                self.beam_mode, self.nothing_mode]))

        self.setattr_argument("run_time_minutes", NumberValue(10))

        # if False, the beams are balanced pairwise, e.g. MOT1 goes up the same amount MOT3 goes down
        # if True, all four potentiometer knobs are used to adjust MOT1-4 (hence disable_z_beam is
        # ignored)
        self.setattr_argument("individual_beam_mode",
                              BooleanValue(False), "beam balance tune settings")
        self.setattr_argument("max_setpoint_percent_deviation", # nominal
                              NumberValue(0.1),"beam balance tune settings")

        # we can balance the z beams with confidence by measuring the powers outside the chamber,
        # so unless we are fine tuning loading, we may want trust our initial manual balancing.
        # this has no effect if individual_beam_mode is True, which disallows z beam tuning anyway
        self.setattr_argument("disable_z_beam_tuning",
                              BooleanValue(True), "beam balance tune settings")
        self.setattr_argument("coil_volts_multiplier",
                              NumberValue(3.0),"coil tune settings")  # scales the value read by the Sampler
        self.setattr_argument("differential_mode",
                              BooleanValue(True),
                              "coil tune settings")  # scan the coils with respect to the current settings

        self.setattr_argument("differential_multiplier",
                              NumberValue(0.5),
                              "coil tune settings")  # scales the value read by the Sampler

        self.setattr_argument("change_z_offset_and_grad_B_mode",
                              BooleanValue(True),
                              "coil tune settings")  # scan the coils with respect to the current settings

        group = "SPCM settings"
         # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(15 * ms, unit='ms'), group)

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(500), "Laser feedback")
        self.setattr_argument("monitor_only", BooleanValue(False), "Laser feedback")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        self.beam_tuning_disabled = not (self.what_to_tune == self.beam_mode or self.what_to_tune == self.both_mode)

        self.n_steps = int(60*self.run_time_minutes/self.dt_exposure+0.5)
        self.sampler_buffer = np.zeros(8)

        self.setpoint_datasets = ["set_point_PD1_AOM_A1","set_point_PD2_AOM_A2","set_point_PD3_AOM_A3",
                                  "set_point_PD4_AOM_A4","set_point_PD5_AOM_A5","set_point_PD6_AOM_A6"]
        self.default_setpoints = [getattr(self,dataset) for dataset in self.setpoint_datasets]

        self.sampler_channels = [0, 1, 2, 3]  # the sampler channels to read

        self.default_volts = [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT,
                              self.AX_volts_MOT, self.AY_volts_MOT]
        self.volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]

        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)
        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # Turn on AOMs to load the MOT.
        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        delay(1*ms)

        if self.FORT_AOM_on:
            self.dds_FORT.sw.on()

            # wait for AOMs to thermalize
            delay(3000 * ms)

        delay(1 * ms)

        # warm up to get make sure we get to the setpoints
        for i in range(10):
            self.laser_stabilizer.run(monitor_only=self.monitor_only)
            if self.FORT_AOM_on:
                self.dds_FORT.sw.on()

        delay(1*ms)

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        saturated_coils = [False] * 4
        control_volts = [0.0] * 4
        volts = 0.0

        ampl1_factor = 1.0
        ampl2_factor = 1.0
        ampl3_factor = 1.0
        ampl4_factor = 1.0
        ampl5_factor = 1.0
        ampl6_factor = 1.0

        last_state = False

        delay(1 * ms)
        if self.FORT_modulation_switch_on:
            self.FORT_mod_switch.on()  # transmit the modulation to the VCA
            delay(1 * ms)

        print("ready!")

        for i in range(self.n_steps):

            if (i % self.AOM_feedback_period_cycles) == 0:
                print("running feedback")
                self.core.break_realtime()
                self.laser_stabilizer.run(monitor_only=self.monitor_only)
                delay(100*ms)
                if self.FORT_AOM_on:
                    self.dds_FORT.sw.on()
                if self.FORT_modulation_switch_on:
                    self.FORT_mod_switch.on()  # transmit the modulation to the VCA
                    delay(1 * ms)

            else:
                delay(10*ms)

            # todo: probably want to some error handling here so we don't set the amplitude too high
            self.dds_AOM_A1.set(amplitude=self.stabilizer_AOM_A1.amplitude * ampl1_factor,
                                frequency=self.AOM_A1_freq)
            self.dds_AOM_A2.set(amplitude=self.stabilizer_AOM_A2.amplitude * ampl2_factor,
                                frequency=self.AOM_A2_freq)
            self.dds_AOM_A3.set(amplitude=self.stabilizer_AOM_A3.amplitude * ampl3_factor,
                                frequency=self.AOM_A3_freq)
            self.dds_AOM_A4.set(amplitude=self.stabilizer_AOM_A4.amplitude * ampl4_factor,
                                frequency=self.AOM_A4_freq)
            self.dds_AOM_A5.set(amplitude=self.stabilizer_AOM_A5.amplitude * ampl5_factor,
                                frequency=self.AOM_A5_freq)
            self.dds_AOM_A6.set(amplitude=self.stabilizer_AOM_A6.amplitude * ampl6_factor,
                                frequency=self.AOM_A6_freq)

            t_end = self.ttl0.gate_rising(self.dt_exposure)
            counts = self.ttl0.count(t_end)
            count_rate_per_s = counts / self.dt_exposure

            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, count_rate_per_s)

            # check the current state of 14 read by ttl3
            self.ttl3.sample_input()
            delay(0.1 * ms)
            tune_beam_balance = bool(self.ttl3.sample_get())

            if self.what_to_tune == self.beam_mode:
                if last_state != tune_beam_balance:
                    if tune_beam_balance:
                        self.print_async("beam balance mode")
                    else:
                        self.print_async("coil tune mode")
                    last_state = tune_beam_balance

            self.sampler1.sample(self.sampler_buffer)

            if self.what_to_tune != self.nothing_mode:
                # only defer to tune_beam_balance if the mode is 'both'

                if (tune_beam_balance and self.what_to_tune == self.both_mode) or self.what_to_tune == self.beam_mode:
                    # print("updating amplitudes") # todo remove
                    if self.individual_beam_mode:
                        ampl1_factor = 1.0 + self.sampler_buffer[self.sampler_channels[0]] * self.max_setpoint_percent_deviation/3.5
                        ampl2_factor = 1.0 + self.sampler_buffer[self.sampler_channels[1]] * self.max_setpoint_percent_deviation/3.5
                        ampl3_factor = 1.0 + self.sampler_buffer[self.sampler_channels[2]] * self.max_setpoint_percent_deviation/3.5
                        ampl4_factor = 1.0 + self.sampler_buffer[self.sampler_channels[3]] * self.max_setpoint_percent_deviation/3.5

                    else:
                        ampl1_factor = 1.0 + self.sampler_buffer[self.sampler_channels[0]] * self.max_setpoint_percent_deviation/3.5
                        ampl3_factor = 1.0 - self.sampler_buffer[self.sampler_channels[0]] * self.max_setpoint_percent_deviation/3.5
                        ampl2_factor = 1.0 + self.sampler_buffer[self.sampler_channels[1]] * self.max_setpoint_percent_deviation/3.5
                        ampl4_factor = 1.0 - self.sampler_buffer[self.sampler_channels[1]] * self.max_setpoint_percent_deviation/3.5

                        if not self.disable_z_beam_tuning:
                            ampl5_factor = 1.0 + self.sampler_buffer[self.sampler_channels[2]] * self.max_setpoint_percent_deviation/3.5
                            ampl6_factor = 1.0 - self.sampler_buffer[self.sampler_channels[2]] * self.max_setpoint_percent_deviation/3.5
                else:

                    if self.differential_mode:
                        if self.change_z_offset_and_grad_B_mode:
                            # the anti-Helmholtz contribution
                            control_volts[0] = self.sampler_buffer[0] * self.differential_multiplier + self.default_volts[0]
                            control_volts[1] = self.sampler_buffer[0] * self.differential_multiplier + self.default_volts[1]

                            # an offset to the Z coils. the -1 sign is because of how the Z coils are
                            # wired. this is the Helmholtz contribution to shim the field
                            control_volts[0] += self.sampler_buffer[1] * self.differential_multiplier
                            control_volts[1] += -1 * self.sampler_buffer[1] * self.differential_multiplier

                            control_volts[2:] = [self.sampler_buffer[ch] * self.differential_multiplier + self.default_volts[ch]
                                                 for ch in self.sampler_channels[2:]]

                        else:
                            control_volts = [self.sampler_buffer[ch] * self.differential_multiplier + self.default_volts[ch]
                                         for ch in self.sampler_channels]

                    else:
                        if self.change_z_offset_and_grad_B_mode:

                            # the anti-Helmholtz contribution
                            control_volts[0] = self.sampler_buffer[0] * self.coil_volts_multiplier
                            control_volts[1] = self.sampler_buffer[0] * self.coil_volts_multiplier

                            # an offset to the Z coils. the -1 sign is because of how the Z coils are
                            # wired. this is the Helmholtz contribution to shim the field
                            control_volts[0] += self.sampler_buffer[1] * self.coil_volts_multiplier
                            control_volts[1] += -1 * self.sampler_buffer[1] * self.coil_volts_multiplier


                            control_volts[2:] = [self.sampler_buffer[ch] * self.coil_volts_multiplier
                                             for ch in self.sampler_channels[2:]]

                        else:
                            control_volts = [self.sampler_buffer[ch] * self.coil_volts_multiplier
                                             for ch in self.sampler_channels]


                    delay(1 * ms)

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

                delay(1*ms)

        if not self.leave_coils_on_at_finish:
            for ch in self.coil_channels:
                self.zotino0.write_dac(ch, 0.0)
                self.zotino0.load()
                delay(1 * ms)

        if self.what_to_tune != self.nothing_mode:
            if self.set_current_parameters_at_finish:

                if self.what_to_tune == self.both_mode or self.what_to_tune == self.beam_mode:
                    # run with monitor only will update the values it read from the detectors
                    # without doing feedback
                    self.laser_stabilizer.run(monitor_only=True, defaults_at_start=False)

                    current_power_setpoints = [
                        self.stabilizer_AOM_A1.value,
                        self.stabilizer_AOM_A2.value,
                        self.stabilizer_AOM_A3.value,
                        self.stabilizer_AOM_A4.value,
                        self.stabilizer_AOM_A5.value,
                        self.stabilizer_AOM_A6.value
                    ]

                    print("setting current laser setpoints")
                    print(current_power_setpoints)

                    if self.disable_z_beam_tuning:
                        n_beams = 4
                    else:
                        n_beams = 6

                    ampl_factors = [ampl1_factor,ampl2_factor,ampl3_factor,ampl4_factor,ampl5_factor,ampl6_factor]
                    for i in range(n_beams):
                        print("amplitude factor", i, ampl_factors[i])

                    for i in range(n_beams):
                        print("MOT ", i+1, "setpoint % change:", current_power_setpoints[i]/self.default_setpoints[i])
                        self.set_dataset(self.setpoint_datasets[i], current_power_setpoints[i], broadcast=True, persist=True)

                if self.what_to_tune == self.both_mode or self.what_to_tune == self.coil_mode:
                    print("setting current coil values")
                    for i in range(4):
                        volts = control_volts[i]
                        if (volts**2)**(1/2) > 10.0:
                            sign = volts / (volts ** 2) ** (1 / 2)
                            volts = sign * 9.9
                        self.set_dataset(self.volt_datasets[i], volts, broadcast=True, persist=True)

        print("Experiment finished.")
