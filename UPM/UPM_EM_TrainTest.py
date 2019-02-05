from UPM.PredictAges import predict_age
from UPM.UPM_EM import MethylationEM


class MethylationEM_TrainTest:

    def __init__(self, training_ages=None, training_methylation=None, training_samples=None, testing_ages=None,
                 testing_methylation=None, testing_samples=None, cpg_site_labels=None):
        self.training_ages = training_ages
        self.training_methylation = training_methylation
        self.training_samples = training_samples
        self.testing_ages = testing_ages
        self.testing_methylation = testing_methylation
        self.testing_samples = testing_samples
        self.cpg_site_labels = cpg_site_labels
        self.PM_model = None
        self.test_ages = None

    def run(self):
        self.fit_model()
        self.predict_test()

    def fit_model(self):
        self.PM_model = MethylationEM(methylation_table=self.training_methylation,
                                      sample_list=self.training_samples,
                                      site_list=self.cpg_site_labels,
                                      times=self.training_ages,
                                      iter_limit=100,
                                      err_tolerance=0.00001)

    def predict_test(self):
        self.test_ages = predict_age(table=self.testing_methylation,
                                     r_rates=self.PM_model.UPM_EC_EM_results['PM_rates'],
                                     r_d=self.PM_model.UPM_EC_EM_results['PM_d'])
