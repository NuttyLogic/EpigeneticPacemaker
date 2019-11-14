import random
from typing import List
import numpy as np
from tqdm import tqdm
from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker
from EpigeneticPacemaker.EPMBase import EPMBase


class EpigeneticPacemakerCV(EPMBase):
    """
        """

    def __init__(self, verbose: bool = False,
                 cv_folds: int = 3, randomize_order: bool = False, iter_limit: int = 100,
                 error_tolerance: float = 0.00001):
        EPMBase.__init__(self)
        self.tqdm_disable = True if not verbose else False
        self.cv_folds = cv_folds
        self._epm = EpigeneticPacemaker(iter_limit=iter_limit, error_tolerance=error_tolerance)
        self.randomize = randomize_order
        self.models = {}
        self.predicted_states = {}

    def fit(self, meth_array: np.ndarray, states: np.ndarray):
        assert isinstance(meth_array, (np.ndarray, np.generic)), 'Pass numpy array'
        assert isinstance(states, (np.ndarray, np.generic)), 'Pass numpy array'
        cv_groups = self.get_cv_folds(len(states))
        fold_count = 0
        for test_indices in tqdm(cv_groups, desc='Processing Folds', disable=self.tqdm_disable):
            train_indices = [index for index in range(len(states)) if index not in test_indices]
            train_array = meth_array[:, train_indices]
            train_states = states[train_indices]

            test_array = meth_array[:, test_indices]

            self._epm.fit(meth_array=train_array, states=train_states)
            test_states = self._epm.predict(meth_array=test_array)
            for index, state in zip(test_indices, test_states):
                self.predicted_states[index] = state

            self.models[f'iter_{fold_count}'] = dict(test_indices=test_indices,
                                                     EPM_rates=np.copy(self._epm.EPM['EPM_rates']),
                                                     EPM_intercepts=np.copy(self._epm.EPM['EPM_intercepts']))
            fold_count += 1
        self.set_cv_model()

    def get_cv_folds(self, sample_number: int) -> List[List[int]]:
        if self.cv_folds < 0:
            self.cv_folds = sample_number
        sample_indices = [count for count in range(sample_number)]
        if self.randomize:
            random.shuffle(sample_indices)
        step_size = int(sample_number / self.cv_folds)
        sample_remainder = sample_number % self.cv_folds
        if sample_remainder:
            step_size += int(sample_remainder / (self.cv_folds - 1))
        test_indices = []
        for fold in range(self.cv_folds):
            if fold + 1 == self.cv_folds:
                test_indices.append(sample_indices[fold * step_size:])
            else:
                test_indices.append(sample_indices[fold * step_size: fold * step_size + step_size])
        return test_indices

    def set_cv_model(self):
        rates, intercepts = [], []
        for fold, model in self.models.items():
            rates.append(model['EPM_rates'])
            intercepts.append(model['EPM_intercepts'])
        cv_rates: np.ndarray = np.array([np.mean(rate) for rate in zip(*rates)])
        cv_intercepts: np.ndarray = np.array([np.mean(intercept) for intercept in zip(*intercepts)])
        self.EPM = dict(EPM_rates=cv_rates, EPM_intercepts=cv_intercepts)
        predicted_states = []
        for index in range(len(self.predicted_states)):
            predicted_states.append(self.predicted_states[index])
        self.predicted_states = np.array(predicted_states)
