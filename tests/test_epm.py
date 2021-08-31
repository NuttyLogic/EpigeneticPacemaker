import os
import unittest
import numpy as np
from EpigeneticPacemaker import EpigeneticPacemaker

data_dir = os.path.dirname(os.path.realpath(__file__))

# import simulated control data
phenos = np.genfromtxt(f'{data_dir}/test_data/val_phenos.tsv', delimiter='\t')
meth = np.genfromtxt(f'{data_dir}/test_data/val_meth.tsv', delimiter='\t')
site_info = np.genfromtxt(f'{data_dir}/test_data/val_site_info.tsv', delimiter='\t')

# unpack data
ages, expected_states = phenos[:, 0], phenos[:, 1]
m_nots, rates = site_info[:, 0], site_info[:, 1]

# fit epm
epm = EpigeneticPacemaker(learning_rate=0.1, scale_X=True, verbose=True)
epm.fit(ages, meth)

states = epm.predict(meth)[:, 0]


class TestEPM(unittest.TestCase):

    def setUp(self):
        pass

    def test_epm_states(self):
        for state, expected_state in zip(states, expected_states):
            self.assertAlmostEqual(state, expected_state, 1)

    def test_site_rates(self):
        for epm_rate, sim_rate in zip(epm._coefs[:, 0], rates):
            self.assertAlmostEqual(epm_rate, sim_rate, 1)

    def test_site_intercepts(self):
        for epm_intercept, m_not in zip(epm._intercepts, m_nots):
            self.assertAlmostEqual(epm_intercept, m_not, 1)


if __name__ == '__main__':
    unittest.main()
