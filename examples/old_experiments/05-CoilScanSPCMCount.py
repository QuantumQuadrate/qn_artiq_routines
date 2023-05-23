"""
This code scans the coils to find the optimum coil parameters that put the MOT at the focus of the
parabolic mirror. It saves the photon counts on a file.

"""

from artiq.experiment import *
import csv
from artiq.coredevice import ad53xx # for converting volts to mu for the zotino
from artiq.coredevice.exceptions import RTIOUnderflow
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

class CoilScan_SPCMCount(EnvExperiment):


    def build(self):
        """
        declare hardware and user-configurable independent variables
        """

        # declare the hardwire devices/channels we will use
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")
        self.setattr_device("zotino0") # for controlling coils
        self.setattr_device("ttl0")  # input for counting SPCM clicks
        self.setattr_device("ttl7")  # output for experiment trigger

        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(0.6, unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        # self.setattr_argument("AZ_bottom_volts_PGC", NumberValue(0*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        # self.setattr_argument("AZ_bottom_volts_RO", NumberValue(0*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        self.setattr_argument("AZ_top_volts_MOT", NumberValue(-2.9, unit="V", ndecimals=3, step=0.025),
                              "A-Z shim/quad top coil")
        # self.setattr_argument("AZ_top_volts_PGC", NumberValue(0 * (-1.64 / 2.5), unit="V", ndecimals=3, step=0.025),
        #                       "A-Z shim/quad top coil")
        # self.setattr_argument("AZ_top_volts_RO", NumberValue(0 * (-1.64 / 2.5), unit="V", ndecimals=3, step=0.025),
        #                       "A-Z shim/quad top coil")
        self.setattr_argument("AX_volts_MOT", NumberValue(-0.35, unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        # self.setattr_argument("AX_volts_PGC", NumberValue(0.5*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        # self.setattr_argument("AX_volts_RO", NumberValue(0.5*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        self.setattr_argument("AY_volts_MOT", NumberValue(-0.43, unit="V", ndecimals=3, step=0.025),"A-Y shim coils")
        # self.setattr_argument("AY_volts_PGC", NumberValue(0.4*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-Y shim coils")
        # self.setattr_argument("AY_volts_RO", NumberValue(0.4*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-Y shim coils")

        # define user-configurable independent variables in human-readable units
        # these will show up in the GUI
        # self.setattr_argument("f_FORT", NumberValue(100.0 * MHz, unit="MHz", ndecimals=1),
        #                       "AOM1, FORT switching AOM")
        # self.setattr_argument("p_FORT_loading", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
        #                       "AOM1, FORT switching AOM")
        # self.setattr_argument("p_FORT_RO", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
        #                       "AOM1, FORT switching AOM")
        # self.setattr_argument("p_FORT_PGC", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
        #                       "AOM1, FORT switching AOM")
        self.setattr_argument("f_cooling_DP_MOT", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_cooling_DP_PGC", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_cooling_DP_RO", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_MOT", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_PGC", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_RO", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        # self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")

        self.setattr_argument("f_cooling_SP", NumberValue(130.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM3, MOT cooling single pass")
        self.setattr_argument("p_cooling_SP", NumberValue(1, unit="dBm", scale=1, ndecimals=1),
                              "AOM3, MOT cooling single pass")
        # self.setattr_argument("Cooling_SP_AOM_ON", BooleanValue(default=False), "AOM3, MOT cooling single pass")

        self.setattr_argument("f_cooling_RP", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("p_cooling_RP", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")
        # self.setattr_argument("MOT_RP_AOM_ON", BooleanValue(default=False), "AOM4, MOT RP/Exc")

        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency, empirically
        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(-5, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        # self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        # self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        # self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        # self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")


        self.setattr_argument("coils_enabled", BooleanValue(True))
        self.setattr_argument("t_MOT_loading", NumberValue(500 * ms, unit="ms", ndecimals=0, step=10 * ms))
        self.setattr_argument("t_SPCM_exposure", NumberValue(50 * ms, unit="ms", ndecimals=1, step=5 * ms))
        self.setattr_argument("datadir", StringValue('/home/fiber-cavity/NetworkExperimentData/CoilScanSPCMCount/'),"File to save data")
        self.setattr_argument("datafile", StringValue('coil_scan_zmax_---mV_xmax_ymax_---mV.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

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
        self.ampl_cooling_DP_MOT = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_MOT / 10 - 3))
        self.ampl_cooling_DP_PGC = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_PGC / 10 - 3))
        self.ampl_cooling_DP_RO = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_RO / 10 - 3))
        self.AOM3_ampl = math.sqrt(2 * 50 * 10 ** (self.p_cooling_SP / 10 - 3))
        self.AOM4_ampl = math.sqrt(2 * 50 * 10 ** (self.p_cooling_RP / 10 - 3))

        self.AOM_A2_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A2_power / 10 - 3))
        self.AOM_A3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A5_power / 10 - 3))
        self.AOM_A6_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A6_power / 10 - 3))

        self.coil_channels = [0,1,2,3]
        # convert frequencies to frequency-tuning-words... when is this necessary?

        # convert voltage for coils to machine units?
        # self.zero_volts = NumberValue(0, unit="V") #ad53xx.voltage_to_mu(0)
        # print(self.zero_volts)

    @kernel
    def run(self):
        self.init_hardware()

        self.file_setup(rowheaders=['counts','AZ_bottom V','AZ_top V','AY V','AX V'])

        # Turn on the magnetic fields
        if self.coils_enabled:
            delay(10 * ms)
            # the starting point
            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                                 channels=self.coil_channels)
        else:
            self.zotino0.set_dac(
                [0.0,0.0,0.0,0.0],
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
        delay(1 * ms)


        #### this was PH's method for scanning, taking into account the impedance mismatch
        # xsteps = ysteps = 25
        # zsteps = 25

        # xmax = ymax = 0.3 # Volts, amplitude of the scan, i.e. we scan +/-xmax V
        # zmax = 0.25 # Volts, """

        # dARTIQdMonitorX = (0.25/0.176) # impedance mismatch coefficients
        # dARTIQdMonitorY = (0.18/0.232)
        # dARTIQdMonitorZtop = (2.9/1.8)
        # dx_Volts =  dARTIQdMonitorX*(2*xmax/(xsteps-1)) # 10 mV steps
        # dy_Volts =  dARTIQdMonitorY*(2*ymax/(ysteps-1)) # 10 mV steps
        # dz_top_Volts = dARTIQdMonitorZtop*(2*zmax/zsteps)  # 20 mV steps
        # print("scanning Z top +/-",1/dARTIQdMonitorZtop*10**3*dz_top_Volts*(zsteps-1)/2,"mV in",zsteps,"steps")
        # print("scanning X +/-",1/dARTIQdMonitorX*10**3*dx_Volts*(xsteps-1)/2,"mV in",xsteps,"steps")
        # print("scanning Y +/-",1/dARTIQdMonitorY*10**3*dy_Volts*(ysteps-1)/2,"mV in",ysteps,"steps")
        # print("estimated duration of scan:",xsteps*ysteps*zsteps*(self.t_MOT_loading+self.t_SPCM_exposure)/3600,"hours")
        # delay(1*ms)

        # for i in range(1): #zsteps):
        #     print(i, "out of ", zsteps)
        #     for j in range(1):#xsteps):
        #         for k in range(1):#ysteps):
        ####

        # xmax = ymax = 0.3 # Volts, amplitude of the scan, i.e. we scan +/-xmax V
        # zmax = 0.25 # Volts, """

        zsteps = 45
        ysteps = 30
        xsteps = 30

        # Vx_steps, Vy_steps, Vz_steps = self.get_xyz_steps(xsteps, ysteps, zsteps)

        print("estimated duration of scan:",xsteps*ysteps*zsteps*(self.t_MOT_loading+self.t_SPCM_exposure)/3600,"hours")

        for i in range(zsteps):
            print(i, "out of ", zsteps)
            for j in range(xsteps):
                for k in range(ysteps):

                    self.core.break_realtime()

                    # do the experiment sequence
                    self.ttl7.pulse(self.t_exp_trigger)

                    # update coil values
                    delay(1 * ms)
                    # coil_volts = [self.AZ_bottom_volts_MOT,
                    #     self.AZ_top_volts_MOT + dz_top_Volts * (-(zsteps - 1) / 2 + i),
                    #     self.AX_volts_MOT + dx_Volts * (-(xsteps - 1) / 2 + j),
                    #     self.AY_volts_MOT + dy_Volts * (-(ysteps - 1) / 2 + k)]

                    Vz = -1.5 - i * (3.3-1.5)/zsteps
                    Vx = 0.15 - j * (0.8 + 0.15) / xsteps
                    Vy = 0.025 - k * (0.9 + 0.025) / ysteps

                    coil_volts = [self.AZ_bottom_volts_MOT,
                                  Vz,
                                  Vx,
                                  Vy]

                    # coil_volts = [self.AZ_bottom_volts_MOT,
                    #               Vz_steps[i],
                    #               Vx_steps[j],
                    #               Vy_steps[k]]

                    # print(coil_volts)

                    delay(1*ms)
                    # print(coil_volts)
                    self.zotino0.set_dac(
                        coil_volts,
                        channels=self.coil_channels)
                    delay(1*ms)

                    # wait for the MOT to load
                    delay_mu(self.t_MOT_loading_mu)
                    # print("i,j,k=",i,j,k)

                    # take the shot
                    t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
                    counts = self.ttl0.count(t_gate_end)
                    # print(counts)  # To-do print value to file instead.
                    # delay(10 * ms)

                    if self.print_meas_result:
                        print("counts", counts)
                    delay(1 * ms)
                    if self.save_data:
                        # all items in the list must be the same type
                        self.file_write([counts+0.0, coil_volts[0], coil_volts[1], coil_volts[2], coil_volts[3]])
                    delay(1 * ms)

                    self.urukul0_ch1.sw.off() # cooling/RP AOMs off
                    self.urukul0_ch2.sw.off()
                    self.urukul0_ch3.sw.off()
                    self.urukul1_ch0.sw.off()  # fiber AOMs off
                    self.urukul1_ch1.sw.off()
                    self.urukul1_ch2.sw.off()
                    self.urukul1_ch3.sw.off()
                    delay(10*ms)
                    self.urukul0_ch1.sw.on()
                    self.urukul0_ch2.sw.on()
                    self.urukul0_ch3.sw.on()
                    self.urukul1_ch0.sw.on()
                    self.urukul1_ch1.sw.on()
                    self.urukul1_ch2.sw.on()
                    self.urukul1_ch3.sw.on()
                    delay(1*ms)

        ### reset parameters
        delay(10*ms)
        self.urukul0_ch1.sw.off() # cooling/RP AOMs off
        self.urukul0_ch2.sw.off()
        self.urukul0_ch3.sw.off()
        self.urukul1_ch0.sw.off() # fiber AOMs off
        self.urukul1_ch1.sw.off()
        self.urukul1_ch2.sw.off()
        self.urukul1_ch3.sw.off()
        self.urukul0_ch0.sw.off()  # FORT AOM off
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=[0, 1, 2, 3])

        self.cleanup() # asynchronous stuff like closing files
        print("Experiment finished.")


    @rpc(flags={"async"})  # means this code runs asynchronously; won't block the rtio counter
    def cleanup(self):
        self.file_obj.close()

    @rpc(flags={"async"})
    def get_xyz_steps(self, xsteps, ysteps, zsteps):
        """the voltage steps to output with the zotino"""
        Vz_steps = np.linspace(-1.5, -3.3, zsteps)
        Vx_steps = np.linspace(-0.15, -0.8, xsteps)
        Vy_steps = np.linspace(0.025, -0.9, ysteps)
        return Vx_steps, Vy_steps, Vz_steps

    @rpc(flags={"async"})
    def file_setup(self,rowheaders=[]):
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

        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))
        self.zotino0.init()
        self.urukul0_cpld.init()
        self.urukul1_cpld.init()

        self.core.break_realtime()

        self.urukul0_cpld.set_profile(0)
        self.urukul1_cpld.set_profile(0)

        self.core.break_realtime()

        # URUKUL 0 - FORT, MOT and D2 state prep AOMs:

        delay(1 * ms)
        self.urukul0_ch1.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(1 * ms)
        self.urukul0_ch2.set(frequency=self.f_cooling_SP, amplitude=self.AOM3_ampl)

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.f_cooling_RP, amplitude=self.AOM4_ampl)

        # URUKUL 1 - MOT arm fiber AOMs:
        delay(1 * ms)
        self.urukul1_ch0.set(frequency=self.AOM_A2_freq, amplitude=self.AOM_A2_ampl)

        delay(1 * ms)
        self.urukul1_ch1.set(frequency=self.AOM_A3_freq, amplitude=self.AOM_A3_ampl)

        delay(1 * ms)
        self.urukul1_ch2.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)

        delay(1 * ms)
        self.urukul1_ch3.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)