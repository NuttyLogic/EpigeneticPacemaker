import numpy as np
from tqdm import tqdm
from UPM.UPM_EM_TrainTest import MethylationEM_TrainTest


class UPM_CV:
    """Cross validated Universal Pacemaker model fitting. Returns age prediction for samples outside training
    set and average site rates and starting methylation values
    Keywords:
        methylation_array (np.array / list): array of methylation values with rows as CpG sites and samples as columns,
                                             can accept lists or np.array
        upm_signal (np.array / list): UPM signal (ie. age) to fit model
        cv_size (int): cross validation testing sample size
        sample_labels (list of str): ordered sample label labels to return predicted UPM age
        methylation_sites (np.array or list): if passed subset of methylation array rows are used to fit the model,
                                              else all sites in the array are used
        collect_stats (bool): default=True
        verbose (bool): default=False
    Attributes:
        self.methylation_array (np.array / list): array of methylation values with rows as CpG sites
                                                  and samples as columns, can accept lists or np.array
        self.upm_signal (np.array / list): UPM signal (ie. age) to fit model
        self.cv_size (int): cross validation testing sample size
        self.sample_labels (list of str): ordered sample label labels to return predicted UPM age
        self.methylation_sites (np.array or list): if passed subset of methylation array rows are used to fit the model,
                                                   else all sites in the array are used
        self.collect_stats (bool): default=True
        self.tqdm_disable (bool): disable TQDM output
        self.upm_run_stats (dict): pm_rates and pm_d for each CV fold
        self.cv_rates (dict): average pm_rates and pm_d across CV folds
        self.predicted_agas (dict): dict of predicted UPM signals
        """

    def __init__(self, methylation_array=None, upm_signal=None, cv_size=1,
                 sample_labels=None, methylation_sites=None, collect_stats=True, verbose=False):
        self.methylation_array = np.asarray(methylation_array)
        self.upm_signal = upm_signal
        self.cv_size = cv_size
        self.sample_indicies = [x for x in range(len(methylation_array[0]))]
        self.sample_labels = sample_labels
        self.methylation_sites = methylation_sites
        if not self.methylation_sites:
            self.methylation_sites = [x for x in range(len(self.methylation_array[:, 0]))]
        self.collect_stats = collect_stats
        self.tqdm_disable = True if not verbose else False
        self.upm_run_stats = {}
        self.cv_results = {}
        self.predicted_ages = {}

    def cv_upm_run(self):
        steps = int(len(self.sample_indicies) / self.cv_size) + 1
        if not len(self.sample_indicies) % self.cv_size:
            steps = int(len(self.sample_indicies) / self.cv_size)
        iter_count = 0
        for step in tqdm(range(steps), desc='Processing Folds', disable=self.tqdm_disable):
            start = step * self.cv_size
            end = start + self.cv_size

            # select test indicies
            test_samples = self.sample_indicies[start:end]
            test_methylation_array = self.methylation_array[:, start:end]
            test_signal = self.upm_signal[start:end]

            # select training indicies
            training_samples = [x for x in self.sample_indicies if x not in test_samples]
            training_methylation_array = self.methylation_array[:, training_samples]
            training_signal = [self.upm_signal[x] for x in training_samples]

            # subest methylation array
            test_methylation_array = test_methylation_array[self.methylation_sites, :]
            training_methylation_array = training_methylation_array[self.methylation_sites, :]
            # get upm object
            upm = MethylationEM_TrainTest(training_ages=training_signal,
                                          training_methylation=training_methylation_array,
                                          training_samples=training_samples,
                                          testing_ages=test_signal,
                                          testing_methylation=test_methylation_array,
                                          testing_samples=test_samples)
            upm.run()
            if self.collect_stats:
                self.upm_run_stats[f'iter_{iter_count}'] = dict(test_samples=test_samples,
                                                                PM_rates=upm.PM_model.UPM_EC_EM_results['PM_rates'],
                                                                PM_d=upm.PM_model.UPM_EC_EM_results['PM_d'])
            for sample_index, sample_age in zip(test_samples, upm.test_ages['PM_times']):
                self.predicted_ages[self.sample_labels[sample_index]] = sample_age
            iter_count += 1

    def get_cv_rates(self):
        pm_rates = []
        pm_d = []
        for count, values in enumerate(self.upm_run_stats.values()):
            pm_rates.append(values['PM_rates'])
            pm_d.append(values['PM_d'])
        cv_pm = []
        cv_d = []
        for site_pm in zip(*pm_rates):
            cv_pm.append(np.mean(site_pm))
        for site_d in zip(*pm_d):
            cv_d.append(np.mean(site_d))
        self.cv_results['PM_rates'] = cv_pm
        self.cv_results['PM_d'] = cv_d
