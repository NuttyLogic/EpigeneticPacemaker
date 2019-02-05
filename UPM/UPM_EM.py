
# ! usr/bin/env python3

import numpy as np
import scipy.stats as stats
from Utilities.UPMHelpers import PMEM


class MethylationEM:

    def __init__(self, methylation_table=None, sample_list=None, site_list=None, times=None, iter_limit=100,
                 err_tolerance=0.0001):
        self.methylation_table = methylation_table
        self.sample_list = sample_list
        self.site_list = site_list
        self.times = times
        self.number_sites = len(methylation_table)
        self.number_times = len(times)
        self.UPM_EC_EM_results = self.get_em_results(iter_limit=iter_limit, err_tolerance=err_tolerance)
        self.UPM_EC_EM_results['MC_times'] = list(self.times)
        self.run_statistics = self.get_run_statistics()

    def get_em_results(self, iter_limit=None, err_tolerance=None):
        results = PMEM(times=self.times,
                       table=self.methylation_table,
                       itr_limit=iter_limit,
                       err_tolerance=err_tolerance)
        return results

    def get_run_statistics(self):
        run_stats = {'MC_rss': (self.UPM_EC_EM_results['MC_err'] ** 2) * self.number_sites * self.number_times,
                     'mhtn_rss': (self.UPM_EC_EM_results['PM_err'] ** 2) * self.number_sites * self.number_times}
        run_stats['chi2'] = np.log(run_stats['MC_rss'] / run_stats['mhtn_rss']) * self.number_times * self.number_sites
        run_stats['pval'] = stats.chi2.cdf(x=run_stats['chi2'],  df=self.number_times)
        return run_stats

    def output_results(self):
        times_output = open(f'{self.output}-MC_PM-times.tsv', 'w')
        times_output.write('Sample\tMC-age\tPM-age\tMC-age/PM-age\n')
        for sample, mc_age, pm_age in zip(self.sample_list,
                                          self.UPM_EC_EM_results['MC_times'],
                                          self.UPM_EC_EM_results['PM_times']):
            mc_pm_ratio = mc_age / pm_age
            times_output.write(f'{sample}\t{mc_age}\t{pm_age}\t{mc_pm_ratio:0.6f}\n')
        times_output.close()
        rate_output = open(f'{self.output}-MC_PM-rates.tsv', 'w')
        rate_output.write('Sample\tMC-rate\tPM-rate\tMC-rate/PM-rate\tMC-d\tPM-d\tMC-d/PM-d\n')
        for site, mc_rate, pm_rate, mc_d, pm_d in zip(self.site_list,
                                                      self.UPM_EC_EM_results['MC_rates'],
                                                      self.UPM_EC_EM_results['PM_rates'],
                                                      self.UPM_EC_EM_results['MC_d'],
                                                      self.UPM_EC_EM_results['PM_d']):
            mc_pm_ratio = mc_rate / pm_rate
            mc_pm_d_ratio = mc_d / pm_d
            rate_output.write(f'{site}\t{mc_rate:0.6f}\t{pm_rate:0.6f}\t{mc_pm_ratio:0.6f}'
                              f'\t{mc_d:0.6f}\t{pm_d:0.6f}\t{mc_pm_d_ratio:0.6f}\n')
        rate_output.close()
