import numpy as np
from tqdm import tqdm
from EpigeneticPacemaker.EPMCompute import one_epm_step, predict_epm_states
from EpigeneticPacemaker.EPMScaler import EPMScaler
from EpigeneticPacemaker.Helpers import pearson_correlation


class EPMBase:

    def __init__(self, iter_limit=100, n_jobs=1,
                 error_tolerance=0.001, learning_rate=0.01,
                 scale_X=True, verbose=False):
        self._coefs = None
        self._intercepts = None
        self._error = None
        self.iter_limit = iter_limit
        self.n_jobs = n_jobs
        self.error_tolerance = error_tolerance
        self.learning_rate = learning_rate
        self.scale_X = scale_X
        self.verbose = verbose

    def predict(self, Y: np.ndarray):
        if self._coefs is None:
            print("EPM model not trained\nRun .fit method to train model")
            return 1
        return predict_epm_states(self._coefs,
                                  self._intercepts,
                                  Y)

    def fit_epm(self, X, Y, sample_weights=None, verbose=False):
        fit_X, fit_Y = np.copy(X, order='k'), np.copy(Y, order='k')
        if len(fit_X.shape) == 1:
            fit_X = fit_X.reshape(-1, 1)
        error = None
        scaler = EPMScaler()
        scaler.fit(fit_X)
        for _ in tqdm(range(self.iter_limit), disable=True if not verbose else False):
            _iter_sys, _iter_update = one_epm_step(fit_X, fit_Y, n_jobs=self.n_jobs)
            _iter_error = sum(_iter_sys[2])
            if not error:
                error = _iter_error
            else:
                if error - _iter_error < self.error_tolerance:
                    break
                else:
                    error = _iter_error
            fit_X += self.learning_rate * _iter_update[1].T
            if self.scale_X:
                fit_X = scaler.transform(fit_X)
        self._coefs = _iter_sys[0]
        self._intercepts = _iter_sys[1]
        self._error = error

    def score(self, X, Y):
        predictions = self.predict(Y)
        scores = pearson_correlation(predictions, X.T)
        return np.array([scores[i][i] for i in range(scores.shape[0])])