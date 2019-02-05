import os
import unittest
import numpy as np
from UPM.UPM_EM import MethylationEM
from Utilities.UPMHelpers import import_upm_table
from UPM.UPM_EM_CV import UPM_CV
from UPM.PredictAges import predict_age

test_data_folder = dir_path = os.path.dirname(os.path.realpath(__file__))
test_input_file_path = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt'
control_rates = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-rates-MCvsPM.csv'
control_ages = f'{test_data_folder}/test_data/meth-tab-n100-m200.txt.CEM-times-MCvsPM.csv'

test_samples, test_ages, test_meth_sites, test_meth_table = import_upm_table(upm_file=test_input_file_path)


all_samples = MethylationEM(methylation_table=test_meth_table,
                            sample_list=test_samples,
                            site_list=test_meth_sites,
                            times=test_ages,
                            iter_limit=100,
                            err_tolerance=0.00001)

# test age prediction
sample_0_age = predict_age(np.asarray(test_meth_table)[:, 0:1],
                           r_rates=all_samples.UPM_EC_EM_results['PM_rates'],
                           r_d=all_samples.UPM_EC_EM_results['PM_d'])

# test leave one out
lov_cv = UPM_CV(methylation_array=test_meth_table, upm_signal=test_ages, cv_size=1, sample_labels=test_samples,
                collect_stats=True, verbose=True)
lov_cv.cv_upm_run()
lov_cv.get_cv_rates()


class TestUPM(unittest.TestCase):

    def setUp(self):
        pass

    def test_age_prediction(self):
        self.assertAlmostEqual(sample_0_age['PM_times'][0], all_samples.UPM_EC_EM_results['PM_times'][0])

    def test_lov(self):
        for lov_value, all_value in zip(all_samples.UPM_EC_EM_results['PM_times'], lov_cv.predicted_ages.values()):
            pm_difference = abs(lov_value - all_value)
            self.assertLess(pm_difference, 1)


if __name__ == '__main__':
    unittest.main()
