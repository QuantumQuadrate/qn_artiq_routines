"""
This code scans the coils to find the optimum coil parameters that put the MOT at the focus of the
parabolic mirror. It saves the photon counts in a file.

"""
import os

from artiq.experiment import *
import csv
from datetime import datetime as dt
import numpy as np
import sys
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment


class CoilScanSPCMCount(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()


        self.scan_datasets = ["Vz_bottom_array", "Vz_top_array", "Vx_array", "Vy_array"]
        group = "Coil steps"

        self.setattr_argument("print_time_estimate_only", BooleanValue(True))

        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value), group)
                # print("retrieved dataset", dataset, "=", value)
        except KeyError as e:
            print(e)
            self.setattr_argument("Vz_bottom_array", StringValue(
                '[6 - l*(6 - 8)/20 for l in range(20)]'), "Coil steps")
            self.setattr_argument("Vz_top_array", StringValue(
                '[6.75 - i*(6.75 - 8.75)/20 for i in range(20)]'), "Coil steps")
            self.setattr_argument("Vx_array", StringValue(
                '[-0.5 - j*(-0.5 - 0.5)/20 for j in range(20)]'), "Coil steps")
            self.setattr_argument("Vy_array", StringValue(
                '[-0.5 - j*(-0.5 - 0.5)/20 for k in range(20)]'), "Coil steps")

        self.setattr_argument("differential_scan", BooleanValue(True))
        self.setattr_argument("FORT_on", BooleanValue(False))
        self.setattr_argument("coils_enabled", BooleanValue(True))
        self.setattr_argument("datafile", StringValue('coil_scan.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(30), "Laser feedback")
        self.setattr_argument("enable_laser_feedback", BooleanValue(True), "Laser feedback")

        # dev ops
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

        self.base.set_datasets_from_gui_args()
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
        self.datadir = os.getcwd()
        if self.prepend_date_to_datafile:
            self.datafile = os.path.join(self.datadir, self.t_experiment_run + '_' + self.datafile)
        else:
            self.datafile = os.path.join(self.datadir + self.datafile)

        # save a copy of the strings defining the scan variables
        self.V_array_strings = [
            self.Vz_bottom_array,self.Vz_top_array,self.Vx_array, self.Vy_array]

        # evaluate the strings we used to define the coil steps in the GUI.
        self.Vz_bottom_array = eval(self.Vz_bottom_array) #.replace('zbottom_steps','self.zbottom_steps'))
        self.Vz_top_array = eval(self.Vz_top_array)
        self.Vx_array = eval(self.Vx_array)
        self.Vy_array = eval(self.Vy_array)

        if self.differential_scan:
            self.Vz_bottom_array += self.AZ_bottom_volts_MOT
            self.Vz_top_array += self.AZ_top_volts_MOT
            self.Vx_array += self.AX_volts_MOT
            self.Vy_array += self.AY_volts_MOT

        assert max(abs(self.Vz_bottom_array)) <= 10, "Zotino can not output > 10V. Decrease VZ_bottom scan range!"
        assert max(abs(self.Vz_top_array)) <= 10, "Zotino can not output > 10V. Decrease VZ_top scan range!"
        assert max(abs(self.Vx_array)) <= 10, "Zotino can not output > 10V. Decrease VX scan range!"
        assert max(abs(self.Vy_array)) <= 10, "Zotino can not output > 10V. Decrease VY scan range!"

        self.zbottom_steps = len(self.Vz_bottom_array)
        self.ztop_steps = len(self.Vz_top_array)
        self.xsteps = len(self.Vx_array)
        self.ysteps = len(self.Vy_array)

        self.sampler_buffer = np.full(8, 0.0)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        print("estimated duration of scan:",
              self.xsteps * self.ysteps * self.ztop_steps * self.zbottom_steps * (
                      self.t_MOT_loading + self.t_SPCM_exposure) / 3600,
              "hours")
        if not self.print_time_estimate_only:

            for i in range(4):
                # takes the value from the GUI and updates the dataset so next time we recompute the arguments
                # in the GUI, these values will be the defaults. you can find the values in the hdf file
                # for each experiment.
                self.set_dataset(self.scan_datasets[i], self.V_array_strings[i], broadcast=True, persist=True)

            self.file_setup(rowheaders=['counts','AZ_bottom V','AZ_top V','AY V','AX V','cooling PD V'])

            # Turn on the magnetic fields
            if self.coils_enabled:
                delay(10 * ms)
                # the starting point
                self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                                     channels=self.coil_channels)
            else:
                self.zotino0.set_dac(
                    [0.0, 0.0, 0.0, 0.0],
                    channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # Turn on AOMs to load the MOT.
            self.dds_cooling_DP.sw.on()


            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()

            # wait for AOMs to thermalize
            delay(4000 * ms)

            step = 0
            i_vz_step = 1

            rtio_log("zotino0")

            delay(100*ms)
            if self.enable_laser_feedback:
                self.laser_stabilizer.run()
                delay(100 * ms)

            if self.FORT_on:
                self.dds_FORT.sw.on()

            for Vz_top in self.Vz_top_array:
                print(i_vz_step, "out of ", len(self.Vz_top_array), "outer loop steps")

                # set MOT coils to setting we know generates a MOT, as a check
                coil_volts = [self.AZ_bottom_volts_MOT,
                              self.AZ_top_volts_MOT,
                              self.AX_volts_MOT,
                              self.AY_volts_MOT]

                self.zotino0.set_dac(
                    coil_volts,
                    channels=self.coil_channels)

                delay(1500 * ms)

                # trigger Luca to save an image
                self.zotino0.write_dac(6, 4.0)
                self.zotino0.load()
                delay(5 * ms)
                self.zotino0.write_dac(6, 0.0)
                self.zotino0.load()

                delay(200*ms)

                for Vz_bottom in self.Vz_bottom_array:
                    for Vx in self.Vx_array:
                        for Vy in self.Vy_array:

                            if self.enable_laser_feedback:
                                if (step % self.AOM_feedback_period_cycles) == 0:
                                    self.core.break_realtime()
                                    self.laser_stabilizer.run()
                                    delay(10*ms)

                                    # the feedback turns the FORT off at the end
                                    if self.FORT_on:
                                        self.dds_FORT.sw.on()

                            # do the experiment sequence
                            delay(10 * ms)

                            # update coil values
                            coil_volts = [Vz_bottom,
                                          Vz_top,
                                          Vx,
                                          Vy]

                            delay(1*ms)
                            # print(coil_volts)
                            self.zotino0.set_dac(
                                coil_volts,
                                channels=self.coil_channels)

                            # wait for the MOT to load
                            delay_mu(self.t_MOT_loading_mu)

                            # take the shot
                            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
                            counts = self.ttl0.count(t_gate_end)

                            if self.print_meas_result:
                                print("counts", counts)
                            delay(1 * ms)
                            if self.save_data:
                                ### all items in the list must be the same type that is why 0.0 is added to counts
                                self.sampler0.sample(self.sampler_buffer) # check cooling laser power
                                self.file_write([counts+0.0,
                                                 coil_volts[0],
                                                 coil_volts[1],
                                                 coil_volts[2],
                                                 coil_volts[3],
                                                 self.sampler_buffer[self.cooling_volts_ch]])
                            delay(10 * ms)
                            step += 1
                i_vz_step += 1 # the outer loop counter
            ### reset parameters
            delay(10*ms)

            self.dds_FORT.sw.off()

            self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                                 channels=self.coil_channels)

            self.cleanup()  # asynchronous stuff like closing files
        print("Experiment finished.")


    @rpc(flags={"async"})  # means this code runs asynchronously; won't block the rtio counter
    def cleanup(self):
        self.file_obj.close()

    @rpc(flags={"async"})
    def file_setup(self, rowheaders=[]):
        # open file once, then close it at end of the experiment
        self.file_obj = open(self.datafile, 'a', newline='')
        self.csvwriter = csv.writer(self.file_obj)
        if rowheaders != []:
            self.csvwriter.writerow(rowheaders)

    @rpc(flags={"async"})
    def file_write(self, data):
        self.csvwriter.writerow(data)

