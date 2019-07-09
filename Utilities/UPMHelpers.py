
import numpy as np


def find_solution_direct(meth_matrix, states):
    sum_states = np.sum(states)
    sum_states_sqr = np.sum(states ** 2)

    b1_1 = 1 / sum_states_sqr - sum_states ** 2 / (
                (sum_states_sqr ** 2) * ((sum_states ** 2) / sum_states_sqr - len(states)))
    b1_2 = sum_states / (sum_states_sqr * (sum_states ** 2 / sum_states_sqr - len(states)))

    top = states * b1_1 + b1_2

    b2_1 = (sum_states_sqr * (sum_states ** 2 / sum_states_sqr - len(states)))
    b2_2 = 1 / ((sum_states ** 2) / sum_states_sqr - len(states))

    bottom = states * sum_states / b2_1 - b2_2

    r_rates = np.dot(meth_matrix, top)

    r_d = np.dot(meth_matrix, bottom)

    return r_rates, r_d


def PMEM(states=None, meth_matrix=None, iter_limit=100, error_tolerance=.00001):
    states = np.asarray(states)
    meth_matrix = np.asarray(meth_matrix)
    em_iter: int = 0

    init_err = None
    init_rates = None
    init_d = None

    while True:

        r_rates, r_d = find_solution_direct(meth_matrix=meth_matrix, states=states)

        prev_err = calc_error(meth_matrix=meth_matrix, states=states, rates=r_rates, d=r_d)

        if not em_iter:
            init_err = prev_err
            init_rates = r_rates
            init_d = r_d

        em_iter += 1

        assert (len(r_rates) == meth_matrix.shape[0]), f'Not a rate for every site, {r_rates} ' \
            f'!= number of sites, {meth_matrix.shape[0]})'
        assert (len(r_d) == meth_matrix.shape[0]), f'Not a d for every site, {r_d} ' \
            f'!= number of sites, {meth_matrix.shape[0]})'

        sum_r_sq1 = np.sum(r_rates ** 2)
        sum1 = np.sum(r_rates * r_d)

        sum2 = np.dot(r_rates, meth_matrix)
        states_updated = (sum2 - sum1) / sum_r_sq1

        new_err = calc_error(meth_matrix, states=states_updated, rates=r_rates, d=r_d)

        imp = prev_err - new_err

        assert (new_err < prev_err), f'new_err > prev_err: {new_err} vs {prev_err}'

        states = states_updated
        results_dict = {'MC_err': init_err,
                        'MC_rates': init_rates,
                        'MC_d': init_d,
                        'PM_err': new_err,
                        'PM_times': states,
                        'PM_rates': r_rates,
                        'PM_d': r_d,
                        'EM_iter': em_iter}

        if em_iter == iter_limit:
            break
        elif imp < error_tolerance:
            break

    return results_dict


def calc_error(meth_matrix=None, states=None, rates=None, d=None):
    total_error = 0.0
    try:
        number_sites, number_states = meth_matrix.shape
    except ValueError:
        number_sites, number_states = len(meth_matrix), 1
    for count, site in enumerate(meth_matrix):
        total_error += sum((site - states * rates[count] - d[count]) ** 2)
    return np.sqrt(total_error / (number_sites * number_states))


def import_upm_table(upm_file=None):
    assert(isinstance(upm_file, str))
    upm_file = open(upm_file, 'r')
    sample_names = None
    sample_ages = None
    methylation_sites = []
    methylation_table = []
    with upm_file as upm:
        for line in upm:
            line_split = line.replace('\n', '').replace(' ', '').split('\t')
            if line_split[0] == 'ID':
                sample_names = line_split[1:]
            elif line_split[0] == 'Age':
                sample_ages = [float(x) for x in line_split[1:]]
            else:
                methylation_sites.append(line_split[0])
                methylation_table.append([float(x) for x in line_split[1:]])
    return sample_names, sample_ages, methylation_sites, methylation_table
