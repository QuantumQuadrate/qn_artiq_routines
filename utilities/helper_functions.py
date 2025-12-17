from artiq.experiment import *

@rpc(flags={"async"})
def print_async(*x):
    """print asynchronously so we don't block the RTIO counter.
    useful for debugging"""
    print(*x)

@kernel
def _rand_shift(self, step: TFloat, salt: TInt32) -> TFloat:
    """
    Kernel-safe pseudo-random dither for grid shifting, used in tune_shims_for_atom_loading.

    Returns a random offset in [-step/2, +step/2), using:
      1) a time-based seed (now_mu()) mixed with a user-provided salt
      2) a small xorshift32 scrambler to improve bit mixing
      3) conversion of a 31-bit integer to a float in [0, 1)
      4) mapping to [-0.5, +0.5) * step

    Notes:
    - 'salt' lets you generate independent random shifts for X/Y (use different salts).
    """
    # --- (A) Build a 32-bit seed from time and salt ---
    # now_mu() changes with RTIO time, so later calls naturally get different seeds.
    # Multiplying salt by an LCG constant helps spread different salts across bits.
    x = (now_mu() ^ (salt * 1103515245) ^ 12345) & 0xffffffff

    # --- (B) Scramble the seed with xorshift32 (good enough for dithering) ---
    x ^= ((x << 13) & 0xffffffff)
    x ^= (x >> 17)
    x ^= ((x << 5) & 0xffffffff)
    x &= 0xffffffff  # keep in 32-bit space

    # --- (C) Convert to uniform float u in [0, 1) ---
    # Use only 31 bits to avoid sign issues; multiply by float constant to force float math.
    x = x & 0x7fffffff
    u = x * (1.0 / 2147483648.0)

    # --- (D) Map u to [-step/2, +step/2) ---
    return (u - 0.5) * step

@kernel
def _count_threshold_crossings_in_loading_window(self,
                                                 threshold_per_s,
                                                 gate_time,
                                                 n_samples) -> TInt32:
    """
    Helper function for automatic shim tuning during atom loading.
    Counts the number of above-threshold *intervals* in a time series (your "crossing" metric).
    Equivalent to the AtomLoadingOptimizer logic: count falling edges + final high state.
    """
    atoms_loaded = 0

    # First sample
    delay(100 * us)
    with parallel:
        self.ttl_SPCM0_counter.gate_rising(gate_time)
        self.ttl_SPCM1_counter.gate_rising(gate_time)
    c0 = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)
    q_last = (c0 / gate_time) > threshold_per_s

    # Remaining samples
    for _ in range(n_samples - 1):
        delay(100 * us)
        with parallel:
            self.ttl_SPCM0_counter.gate_rising(gate_time)
            self.ttl_SPCM1_counter.gate_rising(gate_time)
        c = int((self.ttl_SPCM0_counter.fetch_count() + self.ttl_SPCM1_counter.fetch_count()) / 2)
        q = (c / gate_time) > threshold_per_s

        # Count a "loaded atom interval" when we go from above->below (falling edge)
        if q_last and (not q):
            atoms_loaded += 1
        q_last = q

    # If we end above threshold, count that final interval
    if q_last:
        atoms_loaded += 1

    return atoms_loaded