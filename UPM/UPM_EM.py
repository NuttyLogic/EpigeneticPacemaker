#!/usr/bin/env python3

import numpy as np
import scipy.stats as stats
from Utilities.UPMHelpers import PMEM


class MethylationEM:

    def __init__(self, methylation_array=None, sample_list=None, site_list=None, states=None, iter_limit=100,
                 error_tolerance=0.0001):
        self.meth_matrix = np.asarray(methylation_array)
        self.number_sites, self.number_times = self.meth_matrix.shape
        self.sample_list = sample_list
        self.site_list = site_list
        self.states = states
        self.UPM_EC_EM_results = self.get_em_results(iter_limit=iter_limit, error_tolerance=error_tolerance)
        self.UPM_EC_EM_results['MC_times'] = list(self.states)
        self.run_statistics = self.get_run_statistics()

    def get_em_results(self, iter_limit=None, error_tolerance=None):
        results = PMEM(states=self.states,
                       meth_matrix=self.meth_matrix,
                       iter_limit=iter_limit,
                       error_tolerance=error_tolerance)
        return results

    def get_run_statistics(self):
        run_stats = {'MC_rss': (self.UPM_EC_EM_results['MC_err'] ** 2) * self.number_sites * self.number_times,
                     'mhtn_rss': (self.UPM_EC_EM_results['PM_err'] ** 2) * self.number_sites * self.number_times}
        run_stats['chi2'] = np.log(run_stats['MC_rss'] / run_stats['mhtn_rss']) * self.number_times * self.number_sites
        run_stats['pval'] = stats.chi2.cdf(x=run_stats['chi2'],  df=self.number_times)
        return run_stats
