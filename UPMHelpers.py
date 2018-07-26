
import numpy as np


def find_solution_direct(y=None, times=None, n_t=None, n_s=None):

    sum_time = sum(times)
    sum_time_sqr = sum([time ** 2 for time in times])

    time_matrix = []

    b1_1 = 1 / sum_time_sqr - sum_time ** 2 / ((sum_time_sqr ** 2) * ((sum_time ** 2) / sum_time_sqr - n_t))
    b1_2 = sum_time / (sum_time_sqr * (sum_time ** 2 / sum_time_sqr - n_t))

    for i in range(n_s):
        row0 = (n_s * n_t) * [0.]
        for j in range(n_t):
            row0[i * n_t + j] = times[j] * b1_1 + b1_2
        time_matrix.append(row0)

    b2_1 = (sum_time_sqr * (sum_time ** 2 / sum_time_sqr - n_t))
    b2_2 = 1 / ((sum_time ** 2) / sum_time_sqr - n_t)

    for i in range(n_s):
        row0 = (n_s * n_t) * [0.]
        for j in range(n_t):
            row0[i * n_t + j] = sum_time * times[j] / b2_1 - b2_2
        time_matrix.append(row0)

    r_rates = []
    for i in range(n_s):
        temp = 0.
        for j in range(n_t):
            temp += time_matrix[i][i * n_t + j] * y[i * n_t + j]
        r_rates.append(temp)

    r_d = []
    for i in range(n_s):
        temp = 0.
        for j in range(n_t):
            temp += time_matrix[n_s + i][i * n_t + j] * y[i * n_t + j]
        r_d.append(temp)
    return r_rates, r_d


def PMEM(times=None, table=None, itr_limit=100, err_tolerance=.00001):

    times = list(times)
    n_s = len(table)
    n_t = len(table[0])
    y = [table[i][j] for i in range(n_s) for j in range(n_t)]
    imp = 1
    EM_itr = 0

    init_err = None
    init_rates = None
    init_d = None

    while True:

        EM_itr += 1

        r_rates1, r_d1 = find_solution_direct(y=y, times=times, n_t=n_t, n_s=n_s)

        prev_err = calc_error(table=table, times=times, rates=r_rates1, d=r_d1, n_s=n_s, n_t=n_t)

        if EM_itr == 1:
            init_err = prev_err
            init_rates = r_rates1
            init_d = r_d1

        assert (len(r_rates1) == n_s), f'len(r_rates1), {r_rates1} != number of sites, {n_s})'
        assert (len(r_d1) == n_s), f'len(r_d1), {r_d1} != number of sites, {n_s})'

        sum_r_sq1 = sum([x ** 2 for x in r_rates1])
        sum1 = sum([r_rates1[j] * r_d1[j] for j in range(n_s)])

        sum2 = n_t * [0.]
        times_n = n_t * [0.]

        for i in range(n_t):
            sum2[i] = sum([r_rates1[j] * table[j][i] for j in range(n_s)])
            times_n[i] = (sum2[i] - sum1) / sum_r_sq1

        new_err = calc_error(table=table, times=times_n, rates=r_rates1, d=r_d1, n_s=n_s, n_t=n_t)

        imp = prev_err - new_err

        assert(new_err < prev_err), f'new_err > prev_err: {new_err} vs {prev_err}'

        times = times_n

        results_dict = {'MC_err': init_err,
                        'MC_rates': init_rates,
                        'MC_d': init_d,
                        'PM_err': new_err,
                        'PM_times': times_n,
                        'PM_rates': r_rates1,
                        'PM_d': r_d1,
                        'EM_iter': EM_itr}

        if EM_itr == itr_limit:
            break
        elif imp < err_tolerance:
            break

    return results_dict


def calc_error(table=None, times=None, rates=None, d=None, n_s=None, n_t=None):
    tot_err = 0.

    for i in range(n_s):
        for j in range(n_t):
            t_j = times[j]
            err = (table[i][j] - t_j * rates[i] - d[i]) ** 2
            tot_err += err

    return np.sqrt(tot_err / (n_s * n_t))


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
