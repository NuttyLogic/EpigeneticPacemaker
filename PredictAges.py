import numpy as np
from UPMHelpers import calc_error


def predict_age(table=None, r_rates=None, r_d=None):
    n_s = len(table)
    n_t = len(table[0])

    assert (len(r_rates) == n_s), f'len(r_rates1), {r_rates} != number of sites, {n_s})'
    assert (len(r_d) == n_s), f'len(r_d1), {r_d} != number of sites, {n_s})'

    sum_r_sq1 = sum([x ** 2 for x in r_rates])
    sum1 = sum([r_rates[j] * r_d[j] for j in range(n_s)])

    sum2 = np.zeros(n_t)
    times_n = np.zeros(n_t)

    for i in range(n_t):
        sum2[i] = sum([r_rates[j] * table[j][i] for j in range(n_s)])
        times_n[i] = (sum2[i] - sum1) / sum_r_sq1

    new_err = calc_error(table=table, times=times_n, rates=r_rates, d=r_d, n_s=n_s, n_t=n_t)

    results_dict = {'PM_err': new_err,
                    'PM_times': times_n}

    return results_dict
