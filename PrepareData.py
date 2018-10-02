#! usr/bin/env python3

import sys
import numpy as np
from UPMTableIterator import UPMTableIterator

class UPMTable:
    """Reads a Structure Table for the Universal Pace Maker. The structured table is a tsv with identifiers for the
    ID and AGE information.
     table structure:
        Sample_ID
        AGE
        xxx_site_labels"""

    def __init__(self, beta_matrix_path='None', sample_ids=None, age_data=None,
                 number_of_sites_to_keep=None, pearson_age_correlation_threshold=None):
        self.beta_matrix_path = beta_matrix_path
        self.number_of_sites_to_keep = number_of_sites_to_keep
        self.pearson_age_correlation_threshold = pearson_age_correlation_threshold
        if isinstance(sample_ids, list):
            self.sample_ids = sample_ids
        else:
            self.get_sample_ids()
        if isinstance(age_data, list):
            self.age_data = [float(age) for age in age_data]
            self.age_average = np.mean(self.age_data)
            self.age_std = np.std(self.age_data)
            self.x_var = sum([(self.age_data[i] - self.age_average) ** 2 for i in range(len(self.age_data))])
        else:
            self.age_data = age_data
            self.age_average = None
            self.age_std = None
            self.x_var = None
            self.get_age_data()
        self.covariance_sites = []
        self.methylation_array = []
        self.methylation_sites = []
        self.get_beta_age_covariance()
        self.get_variant_sites()
        self.get_beta_matrix()

    def get_age_data(self):
        for count, line in enumerate(UPMTableIterator(self.beta_matrix_path)):
            if line[0] == 'AGE':
                self.age_data = [float(age) for age in line[1:]]
                self.age_average = np.mean(self.age_data)
                self.age_std = np.std(self.age_data)
                self.x_var = sum([(self.age_data[i] - self.age_average) ** 2 for i in range(len(self.age_data))])
                break
            elif count > 10000:
                sys.exit('AGE descriptor not found')

    def get_sample_ids(self):
        for count, line in enumerate(UPMTableIterator(self.beta_matrix_path)):
            if line[0] == 'Sample_ID':
                self.sample_ids = line[1:]
                break
            elif count > 10000:
                sys.exit('Sample_ID descriptor not found')

    def get_beta_age_covariance(self):
        for line in UPMTableIterator(self.beta_matrix_path):
            beta_values = self.convert_beta_values(line[1:])
            beta_average = np.mean(beta_values)
            sum_xy = sum([(self.age_data[i] - self.age_average) * (beta_values[i] - beta_average)
                          for i in range(len(self.age_data))])
            y_var = sum([(beta_values[i] - beta_average) ** 2 for i in range(len(self.age_data))])
            bottom = np.sqrt(self.x_var) * np.sqrt(y_var)
            r = sum_xy / bottom
            self.covariance_sites.append([line[0], abs(r)])

    def get_variant_sites(self):
        self.covariance_sites.sort(key=lambda x: x[1], reverse=True)
        if self.number_of_sites_to_keep:
            self.covariance_sites = self.covariance_sites[0:self.number_of_sites_to_keep]
        elif self.pearson_age_correlation_threshold:
            self.covariance_sites = [site for site in self.covariance_sites if site[1]
                                     >= self.pearson_age_correlation_threshold]
        self.covariance_sites = {site[0]:site[1] for site in self.covariance_sites}

    def get_beta_matrix(self):
        for line in UPMTableIterator(self.beta_matrix_path):
            try:
                self.covariance_sites[line[0]]
            except KeyError:
                continue
            else:
                self.methylation_sites.append(line[0])
                beta_values = self.convert_beta_values(line[1:])
                self.methylation_array.append(beta_values)
        self.methylation_array = np.asarray(self.methylation_array)

    @staticmethod
    def convert_beta_values(values):
        value_list = []
        for value in values:
            try:
                value_list.append(float(value))
            except ValueError:
                value_list.append(np.nan)
        return value_list
