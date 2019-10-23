import scipy.stats as stats
from EpigeneticPacemaker.EPMCompute import predict_epm_states


class EPMBase:

    def __init__(self):
        self.EPM = None

    def predict(self, meth_array):
        if not self.EPM:
            print('Fit model before predicting epigenetic state')
            return None
        return predict_epm_states(meth_array=meth_array,
                                  r_rates=self.EPM['EPM_rates'],
                                  r_d=self.EPM['EPM_intercepts'])['EPM_states']

    def score(self, meth_array, states):
        if not self.EPM:
            print('Fit model before predicting epigenetic state')
            return None
        predicted_states = predict_epm_states(meth_array=meth_array,
                                              r_rates=self.EPM['EPM_rates'],
                                              r_d=self.EPM['EPM_intercepts'])['EPM_states']
        return stats.pearsonr(predicted_states, states)
