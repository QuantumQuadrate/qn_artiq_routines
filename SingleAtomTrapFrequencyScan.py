"""
Two-shot experiment in which the FORT is modulated between the two shots.
By scanning over the modulation frequency we can measure the trap vibrational frequencies
in order to extract the trap depth and Gaussian waist.
"""

from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import load_MOT_and_FORT


class SingleAtomTrapFrequencyScan(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1))
        # self.setattr_argument("require_loading", BooleanValue(True))

        # the voltage steps. the voltage controls the frequency of a Rigol DG1022Z
        rigol_group = "Rigol DG1022Z settings"
        self.scan_datasets = ["f_modulation_sequence"]
        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value),rigol_group)
        except KeyError as e:
            print(e)
            self.setattr_argument("f_modulation_sequence", StringValue(
                'np.linspace(5.0, 20.0, 10)*kHz'),rigol_group)

        # these should correspond to the settings on the Rigol DG1022Z in order to correctly
        # set the Zotino voltage below
        self.setattr_argument("FM_deviation", NumberValue(40.0*kHz, ndecimals=0, unit='kHz',step=1),rigol_group)
        self.setattr_argument("carrier_frequency", NumberValue(50.0*kHz, ndecimals=0, unit='kHz',step=1),rigol_group)

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        self.sampler_buffer = np.full(8, 0.0)
        self.cooling_volts_ch = 7

        self.f_modulation_list = eval(self.f_modulation_sequence)

        # make sure we aren't going to generate a voltage outside of the -5 to 5V range of the Rigol input
        assert min(self.f_modulation_list) >= self.carrier_frequency - self.FM_deviation, \
            "f_modulation lower bound should be > carrier_frequency - FM_deviation!"
        assert max(self.f_modulation_list) <= self.carrier_frequency + self.FM_deviation, \
            "f_modulation upper bound should be < carrier_frequency + FM_deviation!"

        # conversion to voltage assuming a linear response. this is not exactly true.
        self.V_modulation_list = (self.f_modulation_list - self.carrier_frequency)* 5/self.FM_deviation
        self.n_iterations = len(self.f_modulation_list)

        # this should never be triggered because the asserts above should catch this,
        # but it serves as redundancy in case someone changes the frequency modulation definitions
        assert min(self.V_modulation_list) >= -5, f"{min(self.V_modulation_list)} is invalid voltage for Rigol"
        assert max(self.V_modulation_list) <= 5, f"{max(self.V_modulation_list)} is invalid voltage for Rigol"


        print(self.V_modulation_list)

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
        self.set_dataset("photocount_bins", [50], broadcast=True)
        self.set_dataset("photocounts", [0], broadcast=True)
        self.set_dataset("photocounts2", [0], broadcast=True)


        # turn off AOMs we aren't using, in case they were on previously
        self.dds_D1_pumping_SP.sw.off()
        self.dds_excitation.sw.off()

        # turn on cooling MOT AOMs
        self.dds_cooling_DP.sw.on() # cooling double pass
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1*ms)

        counts = 0
        counts2 = 0

        # Zotino channel for biasing the VCA when we aren't modulating the FORT
        self.zotino0.write_dac(5, 0.62)
        self.zotino0.load()

        # ttl channel for toggling between what we send to the VCA: the zotino voltage or the Rigol signal
        self.ttl12.output()
        self.ttl12.off()

        iteration = 0

        # this is the iteration loop, i.e. the variable steps for the experiment
        for V_modulation in self.V_modulation_list:

            # for computing atom loading and retention statistics
            self.atom_loaded = False
            self.atoms_loaded = 0
            self.atoms_retained = 0

            # these are the datasets for plotting only, an we restart them each iteration
            self.set_dataset("photocounts_current_iteration", [0], broadcast=True)
            self.set_dataset("photocounts2_current_iteration", [0], broadcast=True)
            delay(1*ms)

            # loop the experiment sequence
            for measurement in range(self.n_measurements):

                # make sure any beams that can mess up the fW measurement (e.g. OP, excitation) are turned off
                if self.enable_laser_feedback:
                    if measurement % 10 == 0:
                        self.laser_stabilizer.run()
                        delay(1 * ms)
                    self.dds_FORT.sw.on()
                    self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

                ############################
                # load the MOT and FORT
                ############################

                load_MOT_and_FORT(self)

                ############################
                # take the first shot
                ############################
                self.dds_cooling_DP.sw.on()
                t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
                counts = self.ttl0.count(t_gate_end)
                delay(1*ms)
                self.dds_cooling_DP.sw.off()

                delay(3*ms)

                ############################
                # FORT amplitude modulation phase
                ############################

                self.zotino0.write_dac(4, V_modulation) # todo: set to V_modulation
                self.zotino0.load()
                self.ttl12.on() # toggle the modulation to the VCA
                delay(5*ms)
                self.ttl12.off()
                delay(1*ms)

                ############################
                # take the second shot
                ############################

                self.dds_cooling_DP.sw.on()
                t_gate_end = self.ttl0.gate_rising(self.t_SPCM_second_shot)
                counts2 = self.ttl0.count(t_gate_end)
                delay(1*ms)
                self.dds_cooling_DP.sw.off()

                # set the cooling DP AOM to the MOT settings
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

                delay(2*ms)

                # update the datasets
                self.append_to_dataset('photocounts', counts)
                self.append_to_dataset('photocounts2', counts2)
                self.append_to_dataset('photocounts_current_iteration', counts)
                self.append_to_dataset('photocounts2_current_iteration', counts2)
                self.set_dataset("iteration", iteration, broadcast=True)

            iteration += 1

        delay(1*ms)
        # leave MOT on at end of experiment
        self.dds_cooling_DP.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        # effectively turn the FORT AOM off
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

        # finally, in case the worker refuses to die
        self.write_results()
