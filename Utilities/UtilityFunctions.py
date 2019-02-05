from UPM.UPM_EM import MethylationEM


def output_results(upm_object, output_path):
    assert isinstance(upm_object, MethylationEM)
    with open(f'{output_path}-MC_PM-times.tsv', 'w') as times_output:
        times_output.write('Sample\tMC-age\tPM-age\tMC-age/PM-age\n')
        for sample, mc_age, pm_age in zip(upm_object.sample_list,
                                          upm_object.UPM_EC_EM_results['MC_times'],
                                          upm_object.UPM_EC_EM_results['PM_times']):
            mc_pm_ratio = mc_age / pm_age
            times_output.write(f'{sample}\t{mc_age}\t{pm_age}\t{mc_pm_ratio:0.6f}\n')
    with open(f'{output_path}-MC_PM-rates.tsv', 'w') as rate_output:
        rate_output.write('Sample\tMC-rate\tPM-rate\tMC-rate/PM-rate\tMC-d\tPM-d\tMC-d/PM-d\n')
        for site, mc_rate, pm_rate, mc_d, pm_d in zip(upm_object.site_list,
                                                      upm_object.UPM_EC_EM_results['MC_rates'],
                                                      upm_object.UPM_EC_EM_results['PM_rates'],
                                                      upm_object.UPM_EC_EM_results['MC_d'],
                                                      upm_object.UPM_EC_EM_results['PM_d']):
            mc_pm_ratio = mc_rate / pm_rate
            mc_pm_d_ratio = mc_d / pm_d
            rate_output.write(f'{site}\t{mc_rate:0.6f}\t{pm_rate:0.6f}\t{mc_pm_ratio:0.6f}'
                              f'\t{mc_d:0.6f}\t{pm_d:0.6f}\t{mc_pm_d_ratio:0.6f}\n')
