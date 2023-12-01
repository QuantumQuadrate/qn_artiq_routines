"""
Two-shot experiment with a varying time between the shots over which the atom is held in the trap
"""

from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

from utilities.BaseExperiment import BaseExperiment


class SingleAtomTrapLifetime(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # this is an argument for using a scan package, maybe
        self.scan_datasets = ["t_delay_between_shots"]
        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value))
        except KeyError as e:
            print(e)
            self.setattr_argument("t_delay_between_shots", StringValue(
                'np.array([0.0, 1.0, 10.0, 100.0,1000.])*ms'))

        self.setattr_argument("n_measurements", NumberValue(50, ndecimals=0, step=1))
        self.setattr_argument("atom_counts_threshold", NumberValue(260, ndecimals=0, step=1))
        self.setattr_argument("no_first_shot", BooleanValue(default=False))
        self.setattr_argument("MOT_AOMs_always_on", BooleanValue(default=False)) # good for diagnosing readout heating
        self.setattr_argument("do_PGC_in_MOT", BooleanValue(False))
        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), "Histogram setup (set bins=0 for auto)")
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True),"Laser power stabilization")

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.t_exp_trigger = 1*ms

        self.sampler_buffer = np.full(8, 0.0)
        self.cooling_volts_ch = 7

        self.t_delay_between_shots_list = eval(self.t_delay_between_shots)
        self.n_iterations = len(self.t_delay_between_shots_list)

        self.atom_loaded = False
        self.atoms_loaded = 0
        self.atoms_retained = 0
        self.atom_retention = [0.0]*self.n_iterations

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.expt()
        print("Experiment finished.")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        # todo: these are going to be regularly used, so put these in the base experiment
        self.set_dataset("photocounts", [0], broadcast=True)
        self.set_dataset("photocounts2", [0], broadcast=True)

        self.set_dataset("photocount_bins", [self.bins], broadcast=True)
        self.set_dataset("atom_retention", [0.0], broadcast=True)

        # turn on cooling MOT AOMs
        self.dds_cooling_DP.sw.on() # cooling double pass
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        delay(1 * ms)

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1*ms)

        counts = 0
        counts2 = 0

        iteration = 0
        for t_delay_between_shots in self.t_delay_between_shots_list:

            # for computing atom loading and retention statistics, though I usually just do this in post
            self.atom_loaded = False
            self.atoms_loaded = 0
            self.atoms_retained = 0

            # loop the experiment sequence
            for measurement in range(self.n_measurements):

                if self.enable_laser_feedback:
                    if measurement % 10 == 0:
                        self.laser_stabilizer.run()
                        delay(1 * ms)
                    self.dds_FORT.sw.on()
                    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)

                self.ttl7.pulse(self.t_exp_trigger) # in case we want to look at signals on an oscilloscope

                # Turn on the MOT coils and cooling light
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
                # delay(2 * ms)
                self.dds_cooling_DP.sw.on()

                # wait for the MOT to load
                delay_mu(self.t_MOT_loading_mu)

                # todo: try loading from a PGC phase

                # turn on the dipole trap and wait to load atoms
                self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_loading)
                delay_mu(self.t_FORT_loading_mu)

                # turn off the coils
                self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                                     channels=self.coil_channels)
                if self.do_PGC_in_MOT:
                    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
                    delay(self.t_PGC_in_MOT)

                if not self.MOT_AOMs_always_on:
                    self.dds_cooling_DP.sw.off()

                delay(3*ms) # should wait several ms for the MOT to dissipate

                # set the cooling DP AOM to the readout settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

                if not self.no_first_shot:
                    # take the first shot
                    self.dds_cooling_DP.sw.on()
                    t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
                    counts = self.ttl0.count(t_gate_end)
                    delay(1*ms)
                    self.dds_cooling_DP.sw.off()

                delay(t_delay_between_shots)

                # take the second shot
                self.dds_cooling_DP.sw.on()
                t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
                counts2 = self.ttl0.count(t_gate_end)
                delay(1*ms)

                if not self.MOT_AOMs_always_on:
                    self.dds_cooling_DP.sw.off()

                # todo: check the FORT extinction ratio here
                # effectively turn the FORT AOM off
                self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.ampl_FORT_loading)
                # set the cooling DP AOM to the MOT settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

                # analysis
                if counts > self.atom_counts_threshold:
                    self.atoms_loaded += 1
                    self.atom_loaded = True
                else:
                    self.atom_loaded = False

                if counts2 > self.atom_counts_threshold and self.atom_loaded:
                    self.atoms_retained += 1

                delay(2*ms)

                # update the datasets
                if not self.no_first_shot:
                    self.append_to_dataset('photocounts', counts)
                self.append_to_dataset('photocounts2', counts2)


            # compute the retention fraction based and update the dataset
            if self.atoms_loaded > 0:
                retention_fraction = self.atoms_retained/self.atoms_loaded
            else:
                retention_fraction = 0.0
            self.atom_retention[iteration] = retention_fraction
            self.set_dataset('atom_retention',self.atom_retention[:iteration+1],broadcast=True)
            iteration += 1

        delay(1*ms)
        # leave MOT on at end of experiment, but turn off the FORT
        self.dds_cooling_DP.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)