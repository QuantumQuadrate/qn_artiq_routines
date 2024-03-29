"""
This code scans the coils to find the MOT. It sends triggering signals to save TIFFs by Andor Luca.

"""

from artiq.experiment import *
import csv
from datetime import datetime as dt

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class CoilScanFindMOT(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()


        self.scan_datasets = ["Vz_bottom_array", "Vz_top_array", "Vx_array", "Vy_array"]
        group = "Coil steps"
        try:
            for dataset in self.scan_datasets:
                value = self.get_dataset(dataset)
                self.setattr_argument(dataset, StringValue(value), group)
                print("retrieved dataset", dataset, "=", value)
        except KeyError as e:
            print(e)
            self.setattr_argument("Vz_bottom_array", StringValue(
                '[0.6 - l*(0.6 - 1)/20 for l in range(20)]'), "Coil steps")
            self.setattr_argument("Vz_top_array", StringValue(
                '[-1.5 - i*(3.3 - 1.5)/20 for i in range(20)]'), "Coil steps")
            self.setattr_argument("Vx_array", StringValue(
                '[0.15 - j*(0.8 + 0.15)/20 for j in range(20)]'), "Coil steps")
            self.setattr_argument("Vy_array", StringValue(
                '[0.025 - k*(0.9 + 0.025)/20 for k in range(20)]'), "Coil steps")

        # self.setattr_argument("coils_enabled", BooleanValue(True))
        self.setattr_argument("datadir", StringValue('C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\'),"File to save data")
        self.setattr_argument("datafile", StringValue('coil_scan.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")

        ### when to run the AOM feedback (after how many iterations in the for loops)
        # self.setattr_argument("AOM_feedback_period_cycles", NumberValue(30), "Laser feedback")
        # self.setattr_argument("enable_laser_feedback", BooleanValue(True), "Laser feedback")


        ### dev ops
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        # self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
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

        ### where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        if self.prepend_date_to_datafile:
            self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        else:
            self.datafile = self.datadir + self.datafile

        # save a copy of the strings defining the scan variables
        self.V_array_strings = [
            self.Vz_bottom_array,self.Vz_top_array,self.Vx_array, self.Vy_array]

        # evaluate the strings we used to define the coil steps in the GUI.
        self.Vz_bottom_array = eval(self.Vz_bottom_array) #.replace('zbottom_steps','self.zbottom_steps'))
        self.Vz_top_array = eval(self.Vz_top_array)
        self.Vx_array = eval(self.Vx_array)
        self.Vy_array = eval(self.Vy_array)

        self.zbottom_steps = len(self.Vz_bottom_array)
        self.ztop_steps = len(self.Vz_top_array)
        self.xsteps = len(self.Vx_array)
        self.ysteps = len(self.Vy_array)

        # self.sampler_buffer = [0.0]*8
        # self.cooling_volts_ch = 7 # we'll read this channel later and save it to the file

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()

        for i in range(4):
            # takes the value from the GUI and updates the dataset so next time we recompute the arguments
            # in the GUI, these values will be the defaults. you can find the values in the hdf file
            # for each experiment.
            self.set_dataset(self.scan_datasets[i], self.V_array_strings[i], broadcast=True, persist=True)

        self.file_setup(rowheaders=['AZ_bottom V','AZ_top V','AY V','AX V'])

        # Turn on the magnetic fields
        delay(10 * ms)
        # the starting point
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
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

        print("estimated duration of scan:",
              self.xsteps*self.ysteps*self.ztop_steps*self.zbottom_steps*(self.t_MOT_loading+self.t_SPCM_exposure)/3600,
              "hours")
        step = 0
        i_vz_step = 1

        rtio_log("zotino0")

        for Vz_top in self.Vz_top_array:
            print(i_vz_step, "out of ", len(self.Vz_top_array), "outer loop steps")

            for Vz_bottom in self.Vz_bottom_array:
                for Vx in self.Vx_array:
                    for Vy in self.Vy_array:

                        # self.core.break_realtime()

                        # if self.enable_laser_feedback:
                        #     if (step % self.AOM_feedback_period_cycles) == 0:
                        #         print("running feedback")
                        #         self.core.break_realtime()
                        #         self.laser_stabilizer.run()
                        #         delay(10 * ms)


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

                        ### take the shot by triggering Luca to save an image
                        self.zotino0.write_dac(6, 4.0)
                        self.zotino0.load()
                        delay(5 * ms)
                        self.zotino0.write_dac(6, 0.0)
                        self.zotino0.load()
                        delay(200 * ms)


                        if self.save_data:
                            self.file_write([coil_volts[0],
                                             coil_volts[1],
                                             coil_volts[2],
                                             coil_volts[3]])
                        delay(10 * ms)
                        step += 1
            i_vz_step += 1 # the outer loop counter
        ### reset parameters
        delay(10*ms)

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

