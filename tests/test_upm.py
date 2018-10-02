import unittest
import os
import pandas as pd
import numpy as np
from UPMRun import MethylationEM
from UPMHelpers import import_upm_table

test_data_folder = dir_path = os.path.dirname(os.path.realpath(__file__))
test_input_file_path = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt'
control_rates = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-rates-MCvsPM.csv'
control_ages = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-times-MCvsPM.csv'

test_samples, test_ages, test_meth_sits, test_meth_table = import_upm_table(upm_file=test_input_file_path)


test_em = MethylationEM(methylation_table=test_meth_table,
                        sample_list=test_samples,
                        site_list=test_meth_sits,
                        times=test_ages,
                        iter_limit=100,
                        err_tolerance=0.00001,
                        output_path=f'{test_data_folder}/test_data/test_output')

control_age_df = pd.read_csv(control_ages, sep=',', index_col=0)
control_rate_df = pd.read_csv(control_rates, sep=',', index_col=0)

class TestUPM(unittest.TestCase):

    def setUp(self):
        pass

    def test_pm_ages(self):
        for control_age, test_age in zip(list(control_age_df['PM-age']), test_em.UPM_EC_EM_results['PM_times']):
            assert np.isclose(control_age, test_age)

    def test_pm_rates(self):
        for control_rate, test_rate in zip(list(control_rate_df['PM-rate']), test_em.UPM_EC_EM_results['PM_rates']):
            print(control_rate, test_rate)
            assert np.isclose(control_rate, test_rate)

    def test_rate_len(self):
        self.assertEqual(len(list(control_age_df['PM-age'])), len(test_em.sample_list))


if __name__ == '__main__':
    unittest.main()