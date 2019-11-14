from typing import Dict, Union
import numpy as np
import scipy.stats as stats
from EpigeneticPacemaker.EPMBase import EPMBase
from EpigeneticPacemaker.EPMCompute import epm_expectation_maximization


class EpigeneticPacemaker(EPMBase):

    def __init__(self, iter_limit: int = 100, error_tolerance: float = 0.00001):
        EPMBase.__init__(self)
        self.iter_limit = iter_limit
        self.error_tolerance = error_tolerance

    def fit(self, meth_array: np.ndarray = None, states: np.ndarray = None):
        assert isinstance(meth_array, (np.ndarray, np.generic)), 'Pass numpy array'
        assert isinstance(states, (np.ndarray, np.generic)), 'Pass numpy array'
        self.EPM = epm_expectation_maximization(meth_array=meth_array,
                                                states=states,
                                                iter_limit=self.iter_limit,
                                                error_tolerance=self.error_tolerance)
        self.EPM.update(self.get_epm_statistics(meth_array))

    def get_epm_statistics(self, meth_array: np.ndarray) -> Union[Dict[str, float], None]:
        number_sites, number_samples = meth_array.shape
        if not self.EPM:
            print('Fit model before calculating model statistics')
            return None
        run_stats = {'MC_rss': (self.EPM['MC_error'] ** 2) * number_sites * number_samples,
                     'EPM_rss': (self.EPM['EPM_error'] ** 2) * number_sites * number_samples}
        run_stats['chi2'] = np.log(run_stats['MC_rss'] / run_stats['EPM_rss']) * number_sites * number_samples
        run_stats['pval'] = stats.chi2.cdf(x=run_stats['chi2'], df=number_samples)
        return run_stats
