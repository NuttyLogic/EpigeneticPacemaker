from typing import Tuple, Union
import numpy as np
import scipy.stats as stats
from EpigeneticPacemaker.EPMCompute import predict_epm_states


class EPMBase:

    def __init__(self):
        self.EPM = None

    def predict(self, meth_array: np.ndarray) -> Union[np.array, None]:
        if not self.EPM:
            print('Fit model before predicting epigenetic state')
            return None
        return predict_epm_states(meth_array=meth_array,
                                  rates=self.EPM['EPM_rates'],
                                  intercepts=self.EPM['EPM_intercepts'])['EPM_states']

    def score(self, meth_array: np.ndarray, states: np.ndarray) -> Union[Tuple[float, float], None]:
        if not self.EPM:
            print('Fit model before predicting epigenetic state')
            return None
        predicted_states: np.ndarray = predict_epm_states(meth_array=meth_array,
                                                          rates=self.EPM['EPM_rates'],
                                                          intercepts=self.EPM['EPM_intercepts'])['EPM_states']
        return stats.pearsonr(predicted_states, states)
