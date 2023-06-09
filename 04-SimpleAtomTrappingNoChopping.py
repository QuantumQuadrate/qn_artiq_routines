"""
A simple experiment to look for a trapped single atom signal

1. Turn on cooling and RP AOMs
Experiment cycle (repeats n times)
2. Turn on magnetic fields (Zotino)
3. Turn on fiber AOMs (Urukul), wait some time to load the MOT
4. Turn on the dipole trap AOM
5. Turn off the quadrupole fields (PGC phase)
6. Turn off the fiber AOMs
7. Turn on the fiber AOMs and read from a single SPCM (TTL) for a certain exposure time
8. Store the number of counts registered by the SPCM in an array
End of experiment
9. Save the array of counts to a file


* Tested that all the delays and timings are functioning as expected by pulsing TTLs at different points and monitoring
    with an oscilloscope.

"""

from artiq.experiment import *
import csv
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

from BaseExperiment import BaseExperiment


class SimpleAtomTrapNoChop(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1))
        self.setattr_argument("datadir",
                              StringValue('C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\'),"File to save data")
        self.setattr_argument("datafile", StringValue('atom_loading_counts.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), "Histogram setup")
        self.setattr_argument("counts_per_bin", NumberValue(10, ndecimals=0, step=1), "Histogram setup")
        self.setattr_argument("print_counts", BooleanValue(True))
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True),"Laser power stabilization")

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        if self.prepend_date_to_datafile:
            self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        else:
            self.datafile = self.datadir + self.datafile

        # experiment trigger pulse width
        self.t_exp_trigger = 1*ms

        self.sampler_buffer = [0.0] * 8
        self.cooling_volts_ch = 7

        self.hist_bins = np.zeros(self.bins, dtype=int)
        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.expt()
        print("Experiment finished.")

    @rpc(flags={"async"}) # means this code runs asynchronously; won't block the rtio counter
    def file_setup(self, rowheaders=[]):
        with open(self.datafile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(rowheaders)
            f.close()

    @rpc(flags={"async"})
    def file_write(self, data):
        with open(self.datafile, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)
            f.close()

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        self.set_dataset("photocounts", self.hist_bins, broadcast=True)
        # self.set_dataset("mot_photocounts", self.hist_bins, broadcast=True)

        self.file_setup(rowheaders=['counts'])

        # turn on cooling/RP AOMs
        self.dds_cooling_DP.sw.on() # cooling double pass
        self.dds_cooling_SP.sw.on()  # cooling single pass
        self.dds_MOT_RP.sw.on()  # MOT repump

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        self.AOMservo.run()

        # loop the experiment sequence
        for measurement in range(self.n_measurements):
            self.ttl7.pulse(self.t_exp_trigger) # in case we want to look at signals on an oscilloscope

            # Set magnetic fields for MOT loading
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # Set and turn on fiber AOMs to load the MOT. The MOT AOMs upstream are assumed to be on.
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            delay(1 * ms)

            # wait for the MOT to load
            delay_mu(self.t_MOT_loading_mu)

            # change double pass power and frequency to PGC settings
            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

            # turn on the dipole trap and wait to load atoms
            self.dds_FORT.sw.on()
            delay_mu(self.t_FORT_loading_mu)

            ### change the magnetic fields for imaging
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_RO, self.AZ_top_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                channels=self.coil_channels)
            delay(10 * ms)

            # change AOMs to "imaging" settings
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_RO)
            # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
            delay(10 * ms)

            # # take the shot
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            counts = self.ttl0.count(t_gate_end)
            delay(1*ms)
            if self.print_counts:
                print(counts)
            delay(10 * ms)

            # reset parameters
            # self.dds_AOM_A2.sw.off()  # fiber AOMs off
            # self.dds_AOM_A3.sw.off()
            # self.dds_AOM_A1.sw.off()
            # self.dds_AOM_A6.sw.off()
            # self.dds_AOM_A4.sw.off()
            # self.dds_AOM_A5.sw.off()
            self.dds_FORT.sw.off()  # FORT AOM off
            # self.dds_AOM_A3.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
            # self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
            #                      channels=self.coil_channels)

            bin_idx = int(counts / self.counts_per_bin)
            if bin_idx < self.bins:
                self.hist_bins[bin_idx] += 1
                self.mutate_dataset("photocounts", bin_idx, self.hist_bins[bin_idx])

            if self.print_measurement_number:
                print("measurement", measurement)
            if self.save_data:
                self.file_write([counts])
            delay(10*ms)

        delay(1*ms)
        # leave MOT on at end of experiment, but turn off the FORT
        self.dds_FORT.sw.off()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)