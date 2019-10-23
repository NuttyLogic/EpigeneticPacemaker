import os
from typing import Dict
import unittest
import numpy as np
from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker


def import_upm_table(upm_file=None):
    assert(isinstance(upm_file, str))
    upm_file = open(upm_file, 'r')
    sample_names = None
    sample_ages = None
    methylation_sites = []
    methylation_table = []
    with upm_file as upm:
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


def convert_value(value: str) -> float or str:
    try:
        return float(value)
    except ValueError:
        return value


def import_control_csv(file_path: str) -> Dict[str, np.array]:
    categories = []
    control_data = None
    with open(file_path, 'r') as file:
        for count, line in enumerate(file):
            line_split = line.strip().split(',')
            if count == 0:
                categories = line_split
                control_data = {cat: [] for cat in categories}
            else:
                for category, value in zip(categories, line_split):
                    control_data[category].append(convert_value(value))
    return control_data


test_data_folder = os.path.dirname(os.path.realpath(__file__))
test_input_file_path = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt'
control_rates = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-rates-MCvsPM.csv'
control_ages = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-times-MCvsPM.csv'

test_samples, test_ages, test_meth_sites, test_meth_table = import_upm_table(upm_file=test_input_file_path)


test_em = EpigeneticPacemaker(iter_limit=100, error_tolerance=0.00001)
test_em.fit(np.array(test_meth_table), np.array(test_ages))
test_ages = test_em.predict(np.array(test_meth_table))

control_epm_ages = import_control_csv(control_ages)
control_epm_rates = import_control_csv(control_rates)


class TestUPM(unittest.TestCase):

    def setUp(self):
        pass

    def test_pm_ages(self):
        for control_age, test_age in zip(control_epm_ages['PM-age'], test_ages):
            assert np.isclose(control_age, test_age, rtol=0.001)

    def test_pm_rates(self):
        for control_rate, test_rate in zip(control_epm_rates['PM-rate'], test_em.EPM['EPM_rates']):
            assert np.isclose(control_rate, test_rate, rtol=0.001)

    def test_rate_len(self):
        self.assertEqual(len(test_ages), len(test_samples))


if __name__ == '__main__':
    unittest.main()
