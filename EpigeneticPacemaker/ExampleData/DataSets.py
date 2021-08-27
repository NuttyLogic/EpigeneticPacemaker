import gzip
import io
import os
import numpy as np

example_dir = os.path.dirname(os.path.realpath(__file__)) +'/'


def convert_to_float(line):
    converted_line = []
    for value in line:
        try:
            converted_line.append(float(value))
        except ValueError:
            converted_line.append(value)
    return converted_line


def load_data_set(file_path):
    formatted_data = []
    with io.BufferedReader(gzip.open(file_path, 'rb')) as data:
        for line in data:
            formatted_data.append(convert_to_float(line.decode('utf-8').strip().split('\t')))
    sample_names = formatted_data[0]
    cpg_sites = [line[0] for line in formatted_data[1:-1]]
    phenotypes = np.array(formatted_data[-1][1:])
    methylation_values = np.array([line[1:] for line in formatted_data[1:-1]])
    return sample_names, cpg_sites, phenotypes, methylation_values


def get_example_data():
    test_data = load_data_set(f'{example_dir}GSE74193_test.tsv.gz')
    train_data = load_data_set(f'{example_dir}GSE74193_train.tsv.gz')
    return test_data, train_data
