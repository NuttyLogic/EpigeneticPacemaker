import os
import unittest
import numpy as np
from EpigeneticPacemaker.Helpers import pearson_correlation

data_dir = os.path.dirname(os.path.realpath(__file__))

# import simulated control data
phenos = np.genfromtxt(f'{data_dir}/test_data/val_phenos.tsv', delimiter='\t')
meth = np.genfromtxt(f'{data_dir}/test_data/val_meth.tsv', delimiter='\t')
site_info = np.genfromtxt(f'{data_dir}/test_data/val_site_info.tsv', delimiter='\t')

# unpack data
ages, expected_states = phenos[:, 0], phenos[:, 1]
m_nots, rates = site_info[:, 0], site_info[:, 1]






class TestEPMCV(unittest.TestCase):

    def setUp(self):
        pass

    def test_fold_prediction(self):
        """Check that cv gives same age prediction"""
        for cv_age, epm_age in zip(epm_cv.predictions[:, 0], fold_predictions):
            self.assertAlmostEqual(cv_age, epm_age)

    def test_ran_fold_prediction(self):
        """Check that cv gives same age prediction"""
        for cv_age, ran_cv_age in zip(epm_cv.predictions[:, 0],
                                      epm_cv_ran.predictions[:, 0]):
            self.assertAlmostEqual(cv_age, ran_cv_age, 2)

    def test_site_rates(self):
        for epm_rate, sim_rate in zip(epm_cv._coefs[:, 0], rates):
            self.assertAlmostEqual(epm_rate, sim_rate, 1)

    def test_site_intercepts(self):
        for epm_intercept, m_not in zip(epm_cv._intercepts, m_nots):
            self.assertAlmostEqual(epm_intercept, m_not, 1)


if __name__ == '__main__':
    unittest.main()