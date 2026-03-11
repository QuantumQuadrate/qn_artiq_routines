"""
The laser stabilizer seems to have some issues.
Here I am trying to write a new stabilizer and testing it in an independent experiment.
Have not finished it yet.

Akbar 2026-01-28

"""


from artiq.experiment import *
import numpy as np
import sys, os

# repo path boilerplate (same style you showed)
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd + "\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from utilities.conversions import dB_to_V_kernel


class SimpleFiberAOMStabilizer(EnvExperiment):
    """
    Sequential (series) P-ish stabilizer for 6 fiber AOMs:
      - only one AOM on at a time
      - Sampler averaging (n_avg samples with dt between)
      - up to max_iters updates, then give up
      - RF power clamped to [p_min_dBm, p_max_dBm]
      - updates in ~step_dB increments
      - logs:
          MOTi_monitor (PD history)
          p_AOM_Ai_history (RF dBm history)
    """

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # timing / averaging
        self.setattr_argument("n_avg", NumberValue(10, type="int", ndecimals=0, step=1))

        # loop controls
        self.setattr_argument("max_iters", NumberValue(100, type="int", ndecimals=0, step=1))
        self.setattr_argument("tol_frac", NumberValue(0.02, ndecimals=3, step=0.001))  # 2%


        # proportional-ish gain in dB per fractional error (quantized to step_dB)
        # e.g. err_frac=0.10 and kp_dB=1.0 -> delta=0.10 dB
        # self.setattr_argument("kp_dB", NumberValue(1.0, unit="dB", ndecimals=2, step=0.1))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        # Device + channel mapping (from your existing feedback_channels.json)
        # A1 ch7, A2 ch5, A3 ch3, A4 ch4, A5 ch1, A6 ch2
        self.fiber_dds = [
            self.dds_AOM_A1,
            self.dds_AOM_A2,
            self.dds_AOM_A3,
            self.dds_AOM_A4,
            self.dds_AOM_A5,
            self.dds_AOM_A6,
        ]
        self.sampler_ch = [7, 5, 3, 4, 6, 0]

        # Dataset names
        self.mon_ds = ["MOT1_monitor", "MOT2_monitor", "MOT3_monitor",
                       "MOT4_monitor", "MOT5_monitor", "MOT6_monitor"]
        self.p_hist_ds = ["p_AOM_A1_history", "p_AOM_A2_history", "p_AOM_A3_history",
                          "p_AOM_A4_history", "p_AOM_A5_history", "p_AOM_A6_history"]

        # Setpoint datasets (specified in ExperimentVariables)
        self.setpoint_ds = ["set_point_PD1_AOM_A1", "set_point_PD2_AOM_A2", "set_point_PD3_AOM_A3",
                            "set_point_PD4_AOM_A4", "set_point_PD5_AOM_A5", "set_point_PD6_AOM_A6"]

        # Power datasets used as starting points (in dBm)
        self.p_start_ds = ["p_AOM_A1", "p_AOM_A2", "p_AOM_A3", "p_AOM_A4", "p_AOM_A5", "p_AOM_A6"]

        # AOM frequency attributes from ExperimentVariables (BaseExperiment exposes them on self)
        self.aom_freqs = [
            float(self.AOM_A1_freq),
            float(self.AOM_A2_freq),
            float(self.AOM_A3_freq),
            float(self.AOM_A4_freq),
            float(self.AOM_A5_freq),
            float(self.AOM_A6_freq),
        ]

        # Make kernel-safe
        self.kernel_invariants |= {"aom_freqs"}

        # Pre-create datasets (lists) so you can plot after the run
        for ds in self.mon_ds + self.p_hist_ds:
            self.set_dataset(ds, [], broadcast=True)

        # Make these usable inside kernels
        self.kernel_invariants = getattr(self, "kernel_invariants", set())
        self.kernel_invariants |= {"fiber_dds", "sampler_ch"}



        # Resolve setpoints on host (NO strings / get_dataset in kernel)
        self.setpoints = []
        for ds in self.setpoint_ds:
            v = self.get_dataset(ds, default=None)
            if v is None:
                raise RuntimeError(f"Missing dataset '{ds}'. Set it in GUI or create it before running.")
            self.setpoints.append(float(v))

        # Resolve start powers on host too (in dBm)
        self.p_start = []
        for ds in self.p_start_ds:
            v = self.get_dataset(ds, default=None)
            if v is None:
                raise RuntimeError(f"Missing dataset '{ds}'. Set it in GUI or create it before running.")
            self.p_start.append(float(v))

        # Make kernel-safe
        self.kernel_invariants |= {"setpoints", "p_start"}


        self.dt_between_avg_us = 100
        self.t_settle_us = 1000

        self.p_min_dBm = -30.0
        self.p_max_dBm = -3.0
        self.step_dB = 0.1
        self.kp_dB = 1.0


        print("prepare - done")

    # ---------------- Kernel helpers ----------------

    @kernel
    def _all_fiber_aoms_off(self):
        for d in self.fiber_dds:
            d.sw.off()

    @kernel
    def _set_fiber_aom_dbm(self, idx: TInt32, freq_hz: TFloat, p_dBm: TFloat):
        # Convert dBm to volts-amplitude for DDS amplitude
        amp_v = dB_to_V_kernel(p_dBm)
        dds = self.fiber_dds[idx]
        dds.set(frequency=freq_hz, amplitude=amp_v)

    @kernel
    def _read_sampler0_avg(self, ch: TInt32, n_avg: TInt32, dt_us: TInt32) -> TFloat:
        acc = 0.0
        buf = np.full(8, 0.0)
        for _ in range(n_avg):
            self.sampler0.sample(buf)
            acc += buf[ch]
            delay(dt_us * us)
        return acc / n_avg


    @kernel
    def run(self):
        self.base.initialize_hardware()

        # Ensure starting clean
        self.core.reset()
        self._all_fiber_aoms_off()
        delay(self.t_settle_us * us)

        for i in range(6):
            # Reset per-AOM histories
            self.set_dataset(self.mon_ds[i], [0.0], broadcast=True)
            self.set_dataset(self.p_hist_ds[i], [0.0], broadcast=True)

            setpoint = self.setpoints[i]
            freq = self.aom_freqs[i]


            # Start from the usual GUI dataset p_AOM_Ai (dBm), clamped
            p_dBm = self.p_start[i]

            if p_dBm < self.p_min_dBm:
                p_dBm = self.p_min_dBm
            elif p_dBm > self.p_max_dBm:
                p_dBm = self.p_max_dBm

            # Only this AOM on
            self._all_fiber_aoms_off()
            self._set_fiber_aom_dbm(i, freq, p_dBm)
            self.fiber_dds[i].sw.on()
            delay(self.t_settle_us * us)

            # Feedback iterations (<= max_iters)
            for k in range(int(self.max_iters)):
                pv = float(self._read_sampler0_avg(self.sampler_ch[i], int(self.n_avg), int(self.dt_between_avg_us)))

                # log history (so you can plot after)
                self.append_to_dataset(self.mon_ds[i], pv)
                self.append_to_dataset(self.p_hist_ds[i], p_dBm)

                # stop if within tolerance
                err_frac = (setpoint - pv) / setpoint
                if abs(err_frac) <= float(self.tol_frac):
                    break

                # proportional-ish step in dB, quantized to step_dB
                delta_dB = float(self.kp_dB) * err_frac
                step = float(self.step_dB)

                # quantize; ensure at least one step if still out of tol
                q = round(delta_dB / step) * step
                if q == 0.0:
                    if delta_dB > 0.0:
                        q = step
                    elif delta_dB < 0.0:
                        q = -step
                    else:
                        q = 0.0

                if p_dBm < self.p_min_dBm:
                    p_dBm = self.p_min_dBm
                elif p_dBm > self.p_max_dBm:
                    p_dBm = self.p_max_dBm

                # apply updated RF power
                self._set_fiber_aom_dbm(i, freq, p_dBm)
                delay(self.t_settle_us * us)

            # turn off before moving to the next AOM (series isolation)
            self.fiber_dds[i].sw.off()

        # done
        self._all_fiber_aoms_off()
