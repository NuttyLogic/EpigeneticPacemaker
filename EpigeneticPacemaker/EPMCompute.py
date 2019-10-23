
import numpy as np


def find_solution_direct(meth_array, states):
    sum_states = np.sum(states)
    sum_states_sqr = np.sum(states ** 2)

    b1_1 = 1 / sum_states_sqr - sum_states ** 2 / (
                (sum_states_sqr ** 2) * ((sum_states ** 2) / sum_states_sqr - len(states)))
    b1_2 = sum_states / (sum_states_sqr * (sum_states ** 2 / sum_states_sqr - len(states)))

    top = states * b1_1 + b1_2

    b2_1 = (sum_states_sqr * (sum_states ** 2 / sum_states_sqr - len(states)))
    b2_2 = 1 / ((sum_states ** 2) / sum_states_sqr - len(states))

    bottom = states * sum_states / b2_1 - b2_2

    r_rates = np.dot(meth_array, top)

    r_d = np.dot(meth_array, bottom)

    return r_rates, r_d


def EPM_expectation_maximization(states=None, meth_array=None, iter_limit=100, error_tolerance=.00001):
    states = np.copy(np.asarray(states))
    meth_array = np.asarray(meth_array)
    em_iter: int = 0

    init_err = None
    init_rates = None
    init_d = None
    new_err = None

    while True:

        r_rates, r_d = find_solution_direct(meth_array=meth_array, states=states)

        prev_err = calc_error(meth_array=meth_array, states=states, rates=r_rates, d=r_d)

        if not em_iter:
            init_err = prev_err
            init_rates = r_rates
            init_d = r_d

        em_iter += 1

        assert (len(r_rates) == meth_array.shape[0]), f'Not a rate for every site, {r_rates} ' \
            f'!= number of sites, {meth_array.shape[0]})'
        assert (len(r_d) == meth_array.shape[0]), f'Not a d for every site, {r_d} ' \
            f'!= number of sites, {meth_array.shape[0]})'

        sum_r_sq1 = np.sum(r_rates ** 2)
        sum1 = np.sum(r_rates * r_d)

        sum2 = np.dot(r_rates, meth_array)
        states_updated = (sum2 - sum1) / sum_r_sq1

        new_err = calc_error(meth_array, states=states_updated, rates=r_rates, d=r_d)

        imp = prev_err - new_err

        assert (new_err < prev_err), f'new_err > prev_err: {new_err} vs {prev_err}'

        states = states_updated

        if em_iter == iter_limit:
            break
        elif imp < error_tolerance:
            break

    model_params = {'MC_error': init_err,
                    'MC_rates': init_rates,
                    'MC_intercepts': init_d,
                    'EPM_error': new_err,
                    'EPM_rates': r_rates,
                    'EPM_intercepts': r_d,
                    'EPM_iter': em_iter}

    return model_params


def calc_error(meth_array=None, states=None, rates=None, d=None):
    total_error = 0.0
    number_sites, number_states = meth_array.shape
    for count, site in enumerate(meth_array):
        total_error += sum((site - states * rates[count] - d[count]) ** 2)
    return np.sqrt(total_error / (number_sites * number_states))


def predict_epm_states(meth_array=None, r_rates=None, r_d=None):
    assert (len(r_rates) == meth_array.shape[0]), f'Not a rate for every site, {r_rates} ' \
        f'!= number of sites, {meth_array.shape[0]})'
    assert (len(r_d) == meth_array.shape[0]), f'Not a d for every site, {r_d} ' \
        f'!= number of sites, {meth_array.shape[0]})'

    sum_r_sq1 = np.sum(r_rates ** 2)
    sum1 = np.sum(r_rates * r_d)

    sum2 = np.dot(r_rates, meth_array)
    states_updated = (sum2 - sum1) / sum_r_sq1

    new_err = calc_error(meth_array=meth_array, states=states_updated, rates=r_rates, d=r_d)

    results_dict = {'EPM_error': new_err,
                    'EPM_states': states_updated}

    return results_dict

