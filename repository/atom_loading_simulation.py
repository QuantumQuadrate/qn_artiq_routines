from artiq.experiment import *
import numpy as np
import time
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.special import factorial
from scipy import stats
from scipy.stats import poisson
from repository.run_modeling import *
import matplotlib.pyplot as plt
# a no-hardware simulation that we can use to test plotting

class SingleAtomLoading(EnvExperiment):


    def build(self):
        self.setattr_argument("measurements", NumberValue(100, ndecimals=0, step=1))
        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1))
        self.setattr_argument("counts_per_bin", NumberValue(10, ndecimals=0, step=1))
        
    def sample_photocounts(self):
        """generator function to sample counts from background + single atom distribution"""
    
        n = self.measurements
    
        mu_bg = 100 # the mean background
        mu_sig = 100 # the signal mean 
        poisson = lambda x, m: (m**x/np.math.factorial(x))*np.exp(-m)
        count_dist = lambda counts: poisson(counts, mu_bg) + poisson(counts, mu_sig+mu_bg)
        
        domain = [0,1000] # assume we don't measure fewer than 0 or more than 500 counts
        x1,x2 = domain

        fmax = poisson(mu_bg, mu_bg) # the maximum
        
        # these arrays exist so we could plot the points (x,f) and (x,y) to verify the distribution sampling
        y_dist = np.empty(n) 
        f_dist = np.empty(n)
        x_dist = np.empty(n) # this is the distribution we want
        j = 0 # dist index
        
        while j < n:
        
            x = int((x2-x1)*np.random.rand()+0.5) # rand val on domain of f(x)
            f = count_dist(x)
            y = np.random.rand()*fmax # rand val on range of f(x)
            if y <= f:
                y_dist[j]=y
                f_dist[j]=f
                x_dist[j]=x # x vals with approximate gaussian pdf
                j+=1
                yield x # the number of counts

    def sample_cts_fit(self):
        """generator function to sample counts from background + single atom distribution"""

        n = self.measurements

        mu_bg = 100  # the mean background
        mu_sig = 100  # the signal mean

        def count_dist(x_data, bg, sig):
            return poisson.pmf(x_data, bg) + poisson.pmf(x_data-sig,bg)

        domain = [0, 500]  # assume we don't measure fewer than 0 or more than 500 counts
        x1, x2 = domain

        fmax = poisson.pmf(mu_bg, mu_bg)  # the maximum

        # these arrays exist so we could plot the points (x,f) and (x,y) to verify the distribution sampling
        y_dist = np.empty(n)
        f_dist = np.empty(n)
        x_dist = np.empty(n)  # this is the distribution we want
        j = 0  # dist index

        while j < n:

            x = int((x2 - x1) * np.random.rand() + 0.5)  # rand val on domain of f(x)
            f = count_dist(x, mu_bg, mu_sig)
            y = np.random.rand() * fmax  # rand val on range of f(x)
            if y <= f:
                y_dist[j] = y
                f_dist[j] = f
                x_dist[j] = x  # x vals with approximate gaussian pdf
                j += 1



        def otsu_intraclass_variance(image, threshold):
            """
            Otsu’s intra-class variance.
            If all pixels are above or below the threshold, this will throw a warning that can safely be ignored.
            """
            return np.nansum([
                np.mean(cls) * np.var(image, where=cls)
                #   weight   ·  intra-class variance
                for cls in [image >= threshold, image < threshold]
            ])

        otsu_threshold = min(
            range(int(np.min(x_dist)) + int(self.counts_per_bin), int(np.max(x_dist))),
            key=lambda th: otsu_intraclass_variance(x_dist, th)
        )
        error = []
        #load_error = np.array([1 / np.sqrt(n) if n > 0 else 0 for n in y_dist])
        #print(len(load_error))
        atoms_loaded = np.sum(x_dist[:] >= otsu_threshold)

        #args = x_dist, y_dist
        #print(start_modeling("count_dist", args))
        self.set_dataset("otsu_threshold", otsu_threshold, broadcast=True)
        self.set_dataset("trapped_ct", (np.sum(x_dist[:] >= otsu_threshold)), broadcast=True)
        # self.set_dataset("atoms_loaded", atoms_loaded, broadcast = True)
        self.set_dataset("x_dist", x_dist, broadcast=True)
        self.set_dataset("y_dist", y_dist, broadcast=True)
        self.set_dataset("f_dist", f_dist, broadcast=True)
        #self.set_dataset("err_bars", load_error, broadcast = True)

        # result is a scipy optimize result object, the fit parameters
        # are stored in result.x

    def run(self):
        
        hist_bins = np.zeros(self.bins,dtype=int)
        self.set_dataset("photocounts", hist_bins, broadcast=True)

        #hdbins = np.array(range(0,np.max(),self.counts_per_bin))

        counts = self.sample_photocounts()
        domain = [0,300]
        counts_pruned = np.array([x for x in counts if x < domain[1]])
        ypts, bins, _ = plt.hist(counts_pruned, bins=40)
        xpts = np.linspace(min(counts_pruned), max(counts_pruned), len(bins) - 1)
        self.set_dataset("x1", xpts, broadcast=True)
        self.set_dataset("y1", ypts, broadcast=True)






        print("starting")
        self.sample_cts_fit()
        max_val = -1
        for counts in self.sample_photocounts():
            if counts > max_val:
                max_val = counts
            bin = int(counts/self.counts_per_bin)
            if bin < self.bins:
                hist_bins[bin] += 1
                self.mutate_dataset("photocounts", bin, hist_bins[bin], )

            #print("counts=",counts,"bin=",bin)
            time.sleep(0.1)
        print("finished")
        hdbins = np.array(range(0, max_val, int(max_val/self.counts_per_bin)+1))
        self.set_dataset("hdbins", hdbins, broadcast=True)
        print("creating Applet")




        # NaNs only arise if the class is empty, in which case the contribution should be zero, which `nansum` accomplishes.
class SingleAtomTest(EnvExperiment):
    """
    Count Fiting
    """

    def build(self):
        self.name = "Atom Loading Fit Testing"
        self.setattr_argument("measurements", NumberValue(100, ndecimals=0, step=1))
        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1))
        self.setattr_argument("counts_per_bin", NumberValue(10, ndecimals=0, step=1))
        self.counts = self.sample_photocounts()

    def run(self):
        print(f"START: {self.name}")
        domain = [0,500]
        counts_pruned = np.array([x for x in self.counts if (x < domain[1])] )
        #print("counts pruned")

        p0 = [5, 10, 100, 200, 10, 1]
        bin_count = self.bins
        args = [counts_pruned, p0, bin_count]

        params = start_modeling("count_dist", args=args)
        atoms_loaded = params['atoms_loaded']
        otsu_threshold = params['otsu_threshold']
        modeled_func = params['modeled_func']
        opt_params = params['opt_params']

        #print(opt_params)
        self.set_dataset("otsu_threshold", otsu_threshold, broadcast=True)
        self.set_dataset("atoms_loaded", atoms_loaded, broadcast = True)
        self.set_dataset("real_dat", counts_pruned, broadcast = True)
        self.set_dataset("height_dat", self.y_dist, broadcast=True)
        self.set_dataset("count_dat", modeled_func[0], broadcast=True)
        self.set_dataset("bin_dat", modeled_func[1], broadcast=True)


        print(f"END: {self.name}")

    def sample_photocounts(self):
        """generator function to sample counts from background + single atom distribution"""

        n = self.measurements

        mu_bg = 100  # the mean background
        mu_sig = 100  # the signal mean
        poisson = lambda x, m: (m ** x / np.math.factorial(x)) * np.exp(-m)
        count_dist = lambda counts: poisson(counts, mu_bg) + poisson(counts, mu_sig + mu_bg)

        domain = [0, 1000]  # assume we don't measure fewer than 0 or more than 500 counts
        x1, x2 = domain

        fmax = poisson(mu_bg, mu_bg)  # the maximum

        # these arrays exist so we could plot the points (x,f) and (x,y) to verify the distribution sampling
        self.y_dist = np.empty(n)
        f_dist = np.empty(n)
        x_dist = np.empty(n)  # this is the distribution we want
        j = 0  # dist index

        while j < n:

            x = int((x2 - x1) * np.random.rand() + 0.5)  # rand val on domain of f(x)
            f = count_dist(x)
            y = np.random.rand() * fmax  # rand val on range of f(x)
            if y <= f:
                self.y_dist[j] = y
                f_dist[j] = f
                x_dist[j] = x  # x vals with approximate gaussian pdf
                j += 1
                yield x  # the number of counts



