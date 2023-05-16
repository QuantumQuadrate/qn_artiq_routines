"""
This code scans the coils to find the optimum coil parameters that put the MOT at the focus of the
parabolic mirror. It saves the photon counts in a file.

"""

from artiq.experiment import *
import csv
from artiq.coredevice import ad53xx # for converting volts to mu for the zotino
from artiq.coredevice.exceptions import RTIOUnderflow
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

import importlib
from subroutines.stabilizer import *


class CoilScanSPCMCount(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        # declare the hardwire devices/channels we will use
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")

        self.setattr_device("sampler0")  # for reading photodetector for laser power monitoring
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("ttl0")  # input for counting SPCM clicks
        self.setattr_device("ttl7")  # output for experiment trigger

        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(0.6, unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")

        self.setattr_argument("AZ_top_volts_MOT", NumberValue(-2.9, unit="V", ndecimals=3, step=0.025),
                              "A-Z shim/quad top coil")
        self.setattr_argument("AX_volts_MOT", NumberValue(-0.35, unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        self.setattr_argument("AY_volts_MOT", NumberValue(-0.43, unit="V", ndecimals=3, step=0.025),"A-Y shim coils")
        self.setattr_argument("f_780DP_MOT", NumberValue(111.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_780DP_MOT", NumberValue(-4, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("AOM3_freq", NumberValue(130.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM3, MOT cooling single pass")
        self.setattr_argument("AOM3_power", NumberValue(1, unit="dBm", scale=1, ndecimals=1),
                              "AOM3, MOT cooling single pass")
        self.setattr_argument("AOM4_freq", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("AOM4_power", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")

        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency, empirically
        self.setattr_argument("AOM_A1_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A1")
        self.setattr_argument("AOM_A1_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A1")
        # self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "AOM A1")

        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        # self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        # self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A4_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A4")
        self.setattr_argument("AOM_A4_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A4")
        # self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "AOM A4")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.47 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        # self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.52 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        # self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(0.6, unit="V", ndecimals=3, step=0.01),
                              "A-Z shim/quad bottom coils")
        self.setattr_argument("AZ_top_volts_MOT", NumberValue(-2.5, unit="V", ndecimals=3, step=0.01),
                              "A-Z shim/quad top coil")
        self.setattr_argument("AX_volts_MOT", NumberValue(-0.19, unit="V", ndecimals=3, step=0.01), "A-X shim coils")
        self.setattr_argument("AY_volts_MOT", NumberValue(-0.25, unit="V", ndecimals=3, step=0.01), "A-Y shim coils")
        self.setattr_argument("disable_coils", BooleanValue(default=False))

        self.setattr_argument("xsteps", NumberValue(20, type='int', ndecimals=0, step=1, scale=1), "Coil steps")
        self.setattr_argument("ysteps", NumberValue(20, type='int', ndecimals=0, step=1, scale=1), "Coil steps")
        self.setattr_argument("ztop_steps", NumberValue(20, type='int', ndecimals=0, step=1, scale=1), "Coil steps")
        self.setattr_argument("zbottom_steps", NumberValue(20, type='int', ndecimals=0, step=1, scale=1), "Coil steps")

        self.setattr_argument("coils_enabled", BooleanValue(True))
        self.setattr_argument("t_MOT_loading", NumberValue(500 * ms, unit="ms", ndecimals=0, step=10 * ms))
        self.setattr_argument("t_SPCM_exposure", NumberValue(50 * ms, unit="ms", ndecimals=1, step=5 * ms))
        self.setattr_argument("datadir", StringValue('C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\'),"File to save data")
        self.setattr_argument("datafile", StringValue('coil_scan.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")

        # when to run the AOM feedback (after how many iterations in the for loops)
        self.setattr_argument("AOM_feedback_period_cycles", NumberValue(30))

        # dev ops
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        if self.prepend_date_to_datafile:
            self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        else:
            self.datafile = self.datadir + self.datafile

        # experiment trigger pulse width
        self.t_exp_trigger = 1*ms

        # convert times to machine units
        self.t_MOT_loading_mu = self.core.seconds_to_mu(self.t_MOT_loading)
        self.t_SPCM_exposure_mu = self.core.seconds_to_mu(self.t_SPCM_exposure)

        # converts RF power in dBm to amplitudes in V
        self.ampl_780DP_MOT = math.sqrt(2 * 50 * 10 ** (self.p_780DP_MOT / 10 - 3))
        # self.ampl_780DP_PGC = math.sqrt(2 * 50 * 10 ** (self.p_780DP_PGC / 10 - 3))
        # self.ampl_780DP_imaging = math.sqrt(2 * 50 * 10 ** (self.p_780DP_imaging / 10 - 3))
        self.AOM3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM3_power / 10 - 3))
        self.AOM4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM4_power / 10 - 3))

        self.AOM_A1_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A1_power / 10 - 3))
        self.AOM_A2_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A2_power / 10 - 3))
        self.AOM_A3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A4_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A5_power / 10 - 3))
        self.AOM_A6_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A6_power / 10 - 3))

        self.coil_channels = [0,1,2,3]
        # convert frequencies to frequency-tuning-words... when is this necessary?

        # setup stabilization for the cooling laser power

        # todo: eventually read conversion functions such as this from a config file
        def volts_to_optical_mW(x: TFloat) -> TFloat:
            """
            the conversion of PD voltage to cooling light power at the switchyard MOT 1 path
            """
            x += 0.011  # this accounts for a mismatch between what the Sampler reads and what
            # the multimeter that I used for the fit reads
            return -0.195395 + 17.9214 * x

        self.sampler_buffer = [0.0]*8
        self.cooling_volts_ch = 7
        self.AOMservo = AOMPowerStabilizer(experiment=self,
                                           dds_names=["urukul0_ch1"],
                                           sampler_name="sampler0",
                                           sampler_channels=[self.cooling_volts_ch],
                                           transfer_functions=[volts_to_optical_mW],
                                           setpoints=[0.7],  # in mW
                                           proportionals=[0.04],
                                           iters=5,  # keep iters/t_meas_delay small or rtio underflow
                                           t_meas_delay=20 * ms)

        print("prepare - done")

    @kernel
    def run(self):
        self.init_hardware()

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

        # Set and turn on AOMs to load the MOT.
        self.urukul0_ch1.sw.on()
        self.urukul0_ch2.sw.on()
        self.urukul0_ch3.sw.on()
        self.urukul1_ch0.sw.on()
        self.urukul1_ch1.sw.on()
        self.urukul1_ch2.sw.on()
        self.urukul1_ch3.sw.on()
        self.urukul2_ch0.sw.on()
        self.urukul2_ch1.sw.on()

        # wait for AOMs to thermalize
        delay(4000 * ms)

        print("estimated duration of scan:",
              self.xsteps*self.ysteps*self.ztop_steps*self.zbottom_steps*(self.t_MOT_loading+self.t_SPCM_exposure)/3600,
              "hours")
        step = 0
        for i in range(self.ztop_steps):
            print(i, "out of ", self.ztop_steps, "outer loop steps")
            for l in range(self.zbottom_steps):
                for j in range(self.xsteps):
                    for k in range(self.ysteps):

                        self.core.break_realtime()

                        if (step % self.AOM_feedback_period_cycles) == 0:
                            print("running feedback")
                            self.core.break_realtime()
                            self.AOMservo.run()
                            delay(10*ms)

                        # do the experiment sequence
                        self.ttl7.pulse(self.t_exp_trigger)

                        # update coil values
                        delay(1 * ms)

                        Vz_top = -1.5 - i * (3.3 - 1.5)/self.ztop_steps
                        Vx = 0.15 - j * (0.8 + 0.15) / self.xsteps
                        Vy = 0.025 - k * (0.9 + 0.025) / self.ysteps
                        Vz_bottom = 0.6 - l * (0.6 - 1) /self.zbottom_steps

                        coil_volts = [Vz_bottom,
                                      Vz_top,
                                      Vx,
                                      Vy]

                        delay(1*ms)
                        # print(coil_volts)
                        self.zotino0.set_dac(
                            coil_volts,
                            channels=self.coil_channels)
                        delay(1*ms)

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
                        delay(1 * ms)
                        step += 1

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

    @kernel
    def init_hardware(self):
        """
        Sets amplitudes and frequencies for the Urukul0 channels
        :return:
        """

        self.core.reset()
        self.ttl0.input()  # for reading pulses from SPCM
        self.ttl7.output()  # for outputting a trigger each cycle

        self.urukul0_cpld.init()
        self.urukul1_cpld.init()
        self.urukul2_cpld.init()

        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()
        self.urukul2_ch0.init()
        self.urukul2_ch1.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))
        self.urukul2_ch0.set_att(float(0))
        self.urukul2_ch1.set_att(float(0))
        self.zotino0.init()

        self.core.break_realtime()

        # URUKUL 0 - FORT, MOT and D2 state prep AOMs:

        delay(1 * ms)
        self.urukul0_ch1.set(frequency=self.f_780DP_MOT, amplitude=self.ampl_780DP_MOT)

        delay(1 * ms)
        self.urukul0_ch2.set(frequency=self.AOM3_freq, amplitude=self.AOM3_ampl)

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.AOM4_freq, amplitude=self.AOM4_ampl)

        # URUKUL 1 and 2 - MOT arm fiber AOMs:
        delay(1 * ms)
        self.urukul1_ch0.set(frequency=self.AOM_A2_freq, amplitude=self.AOM_A2_ampl)
        self.urukul1_ch1.set(frequency=self.AOM_A3_freq, amplitude=self.AOM_A3_ampl)
        self.urukul1_ch2.set(frequency=self.AOM_A1_freq, amplitude=self.AOM_A1_ampl)
        self.urukul1_ch3.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)
        self.urukul2_ch0.set(frequency=self.AOM_A4_freq, amplitude=self.AOM_A4_ampl)
        self.urukul2_ch1.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)

        self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set

        print("init_hardware - done")