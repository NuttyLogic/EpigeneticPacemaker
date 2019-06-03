import numpy as np
from Utilities.UPMHelpers import calc_error


def predict_age(meth_matrix=None, r_rates=None, r_d=None):
    assert (len(r_rates) == meth_matrix.shape[0]), f'Not a rate for every site, {r_rates} ' \
        f'!= number of sites, {meth_matrix.shape[0]})'
    assert (len(r_d) == meth_matrix.shape[0]), f'Not a d for every site, {r_d} ' \
        f'!= number of sites, {meth_matrix.shape[0]})'

    sum_r_sq1 = np.sum(r_rates ** 2)
    sum1 = np.sum(r_rates * r_d)

    sum2 = np.dot(r_rates, meth_matrix)
    states_updated = (sum2 - sum1) / sum_r_sq1

    new_err = calc_error(meth_matrix=meth_matrix, states=states_updated, rates=r_rates, d=r_d)

    results_dict = {'PM_err': new_err,
                    'PM_times': states_updated}

    return results_dict
