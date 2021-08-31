import random
import numpy as np
from tqdm import tqdm
from EpigeneticPacemaker.EPMBase import EPMBase


class EpigeneticPacemakerCV(EPMBase):
    """
        """

    def __init__(self,
                 cv_folds: int = 3, randomize_sample_order: bool = False,
                 iter_limit=100, n_jobs=1,
                 error_tolerance=0.001, learning_rate=0.01,
                 scale_X=True, verbose=False, cv_predictions=False
                 ):
        EPMBase.__init__(self)
        self.cv_folds = cv_folds
        self.randomize = randomize_sample_order
        self.iter_limit = iter_limit
        self.n_jobs = n_jobs
        self.error_tolerance = error_tolerance
        self.learning_rate = learning_rate
        self.scale_X = scale_X
        self.verbose = verbose
        self.cv_predictions = cv_predictions
        self.predictions = {}

    def fit(self, X, Y, sample_weights=None):
        cv_groups = self.get_cv_folds(X.shape[0])
        fold_count = 0
        # reshape X if one dimensional
        X_fit = X if len(X.shape) > 1 else X.reshape(-1, 1)
        coefs, intercepts, errors = np.zeros((Y.shape[0], X_fit.shape[1])), np.zeros(Y.shape[0]), 0.0
        training_sample_count = 0
        for test_indices in tqdm(cv_groups, disable=True if not self.verbose else False, desc='CV Folds'):
            train_indices = [index for index in range(X.shape[0]) if index not in test_indices]
            training_sample_count += len(train_indices)
            train_Y = Y[:, train_indices]
            train_X = X_fit[train_indices, :]

            test_Y = Y[:, test_indices]

            self.fit_epm(train_X, train_Y, sample_weights=sample_weights)
            test_states = self.predict(test_Y)
            if self.cv_predictions:
                for index, state in zip(test_indices, test_states):
                    self.predictions[index] = state

            # weight the contribution of each fold by the number of samples in the fold
            coefs += self._coefs * len(train_indices)
            intercepts += self._intercepts * len(train_indices)
            errors += self._error * len(train_indices)
            fold_count += 1
        self._coefs = coefs / training_sample_count
        self._intercepts = intercepts / training_sample_count
        self._error = errors / training_sample_count
        self.unpack_out_of_fold_predictions()

    def get_cv_folds(self, sample_number):
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

    def unpack_out_of_fold_predictions(self):
        if self.cv_predictions:
            self.predictions = np.array([self.predictions[index] for index in range(len(self.predictions))])

