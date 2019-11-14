import os
from typing import List
import unittest
import numpy as np
from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker
from EpigeneticPacemaker.EpigeneticPacemakerCV import EpigeneticPacemakerCV


def import_upm_table(upm_file: str = None):
    sample_names = None
    sample_ages = None
    methylation_sites = []
    methylation_table = []
    with open(upm_file, 'r') as upm:
        for line in upm:
            line_split = line.replace('\n', '').replace(' ', '').split('\t')
            if line_split[0] == 'ID':
                sample_names = line_split[1:]
            elif line_split[0] == 'Age':
                sample_ages = [float(x) for x in line_split[1:]]
            else:
                methylation_sites.append(line_split[0])
                methylation_table.append([float(x) for x in line_split[1:]])
    return sample_names, sample_ages, methylation_sites, methylation_table


def import_control_csv(file_path: str) -> List:
    con_data = []
    with open(file_path, 'r') as file:
        for line in file:
            line_split = line.strip().split(',')
            control_data.append([float(value) for value in line_split])
    return con_data


test_data_folder = os.path.dirname(os.path.realpath(__file__)) + '/test_data/'
test_input_file_path = f'{test_data_folder}meth-tab-n100-m200.txt'

test_samples, test_ages, test_meth_sites, test_meth_table = import_upm_table(upm_file=test_input_file_path)

# fit epm training training model
test_em = EpigeneticPacemaker(iter_limit=100, error_tolerance=0.00001)
test_em.fit(np.array(test_meth_table), np.array(test_ages))
test_ages = test_em.predict(np.array(test_meth_table))

control_data = import_control_csv(f'{test_data_folder}epm_control_data.csv')

control_ages, control_rates, control_intercepts = control_data[0], control_data[1], control_data[2]

# fit epm_cv model
epm_cv = EpigeneticPacemakerCV(cv_folds=4, iter_limit=100, error_tolerance=0.00001, verbose=True)
epm_cv.fit(np.array(test_meth_table), np.array(test_ages))
cv_ages = epm_cv.predicted_states

test_indices = epm_cv.models['iter_0']['test_indices']
train_indices = [count for count in range(len(test_samples)) if count not in test_indices]

# fit a separate model for the first fold to evaluate it against the cv model
epm = EpigeneticPacemaker()
epm.fit(np.array(test_meth_table)[:, train_indices], np.array(test_ages)[train_indices])
epm_test_ages = epm.predict(np.array(test_meth_table)[:, test_indices])

cv_model_ages = epm_cv.predict(np.array(test_meth_table))


class TestEPM(unittest.TestCase):

    def setUp(self):
        pass

    def test_epm_ages(self):
        for control_age, test_age in zip(control_ages, test_ages):
            assert np.isclose(control_age, test_age, rtol=0.001)

    def test_epm_rates(self):
        for control_rate, test_rate in zip(control_rates, test_em.EPM['EPM_rates']):
            assert np.isclose(control_rate, test_rate, rtol=0.001)

    def test_epm_intercepts(self):
        for control_intercept, test_intercept in zip(control_intercepts, test_em.EPM['EPM_intercepts']):
            assert np.isclose(control_intercept, test_intercept, rtol=0.001)

    def test_age_len(self):
        self.assertEqual(len(test_ages), len(test_samples))

    def test_fold_prediction(self):
        """Check that cv gives same age prediction"""
        for cv_age, epm_age in zip(cv_ages, epm_test_ages):
            self.assertAlmostEqual(cv_age, epm_age)

    def test_age_prediction(self):
        """Check that cv age prediction within .1 of control age"""
        for cv_age, control_age in zip(cv_model_ages, control_ages):
            age_diff = abs(control_age - cv_age)
            self.assertLess(age_diff, .1)


if __name__ == '__main__':
    unittest.main()
