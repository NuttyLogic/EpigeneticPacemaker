from typing import Dict, Tuple
import numpy as np


def find_solution_direct(meth_array: np.ndarray, states: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fit linear model at each site of interest"""
    rates: np.ndarray = (np.sum((meth_array - np.mean(meth_array, axis=0)) * (states - np.mean(states)), axis=1)
                         / np.sum((states - np.mean(states))**2))
    intercepts: np.ndarray = np.mean(meth_array, axis=1) - rates * np.mean(states)
    return rates, intercepts


def update_states(meth_array: np.ndarray = None,
                  rates: np.ndarray = None, intercepts: np.ndarray = None) -> np.ndarray:
    """Update state given by S = (m_hat_{i,j} - m_0) / r_i. For numerical simplicity actually calculate
    S = (r_i*\\hat{m}_{i,j} - m^{0}r_i)/(r_i^2))"""
    assert (rates.shape[0] == meth_array.shape[0]), f'Not a rate for every site, {rates.shape[0]} ' \
                                                    f'!= number of sites, {meth_array.shape[0]})'
    assert (intercepts.shape[0] == meth_array.shape[0]), f'Not an intercept for every site, {intercepts.shape[0]} ' \
                                                         f'!= number of sites, {meth_array.shape[0]})'
    return (np.dot(rates, meth_array) - np.sum(rates * intercepts)) / np.sum(rates ** 2)


def epm_expectation_maximization(meth_array: np.ndarray = None, states: np.ndarray = None,
                                 iter_limit: int = 100, error_tolerance: float = .00001) -> Dict[str, np.ndarray]:
    states = np.copy(np.asarray(states))
    em_iter = 0

    init_err, init_rates, init_intercepts = None, None, None

    while True:

        rates, intercepts = find_solution_direct(meth_array=meth_array, states=states)

        prev_err: float = calc_error(meth_array=meth_array, states=states, rates=rates, intercepts=intercepts)

        if not em_iter:
            init_err, init_rates, init_intercepts = prev_err, rates, intercepts

        em_iter += 1

        states_updated = update_states(meth_array=meth_array, rates=rates, intercepts=intercepts)

        new_err = calc_error(meth_array, states=states_updated, rates=rates, intercepts=intercepts)

        # calculate model improvement
        imp = prev_err - new_err

        assert (new_err < prev_err), f'new_err > prev_err: {new_err} vs {prev_err}'

        states = states_updated

        if em_iter == iter_limit:
            break
        elif imp < error_tolerance:
            break

    model_params = {'MC_error': init_err,
                    'MC_rates': init_rates,
                    'MC_intercepts': init_intercepts,
                    'EPM_error': new_err,
                    'EPM_rates': rates,
                    'EPM_intercepts': intercepts,
                    'EPM_iter': em_iter}

    return model_params


def calc_error(meth_array: np.array = None, states: np.array = None,
               rates: np.array = None, intercepts: np.array = None) -> float:
    total_error = 0.0
    number_sites, number_states = meth_array.shape
    for count, site in enumerate(meth_array):
        total_error += sum((site - states * rates[count] - intercepts[count]) ** 2)
    return np.sqrt(total_error / (number_sites * number_states))


def predict_epm_states(meth_array: np.ndarray = None,
                       rates: np.ndarray = None, intercepts: np.ndarray = None) -> Dict[str, np.ndarray]:

    states_updated = update_states(meth_array=meth_array, rates=rates, intercepts=intercepts)

    new_err = calc_error(meth_array=meth_array, states=states_updated, rates=rates, intercepts=intercepts)

    return {'EPM_error': new_err, 'EPM_states': states_updated}

