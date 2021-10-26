import os
import sys
import time
from datetime import timedelta
import numpy as np
import pandas as pd
from argument_deconv import parse_args
from utils_deconv import get_tissue_index, load_np_array, get_tissues_by_col_idx
from single_save_of_probabilistic_data import DISTRIBUTION_SUFFIX_HARD, P_T_GIVEN_R_SUFFIX, PAT_GZ_SUFFIX

MODELS = ['Zero', 'First', 'Second', 'Third', 'Fourth', 'Fifth', 'Sixth', 'Seventh', 'Eighth', 'Ninth']


def validate_tissue_name(args):
    assert (args.tissue_name is not None)


def get_predicted_rates(full_names_lst, tissue_idx, suffix):
    """
    :param full_names_lst: the names of the files
    :param tissue_idx: the index of the relevant (target) tissue
    :param suffix: the suffix of the given files
    :return: a list of the probabilities of the target tissue in each file , the file names (in the same order)
    """
    predicted_rates = []
    file_names = []
    for full_name in full_names_lst:
        cur_dist_vec = load_np_array(full_name)
        predicted_rates.append(cur_dist_vec[tissue_idx])
        filename = os.path.basename(full_name)
        file_names.append(filename.split(suffix)[0] + PAT_GZ_SUFFIX)
    return np.array(predicted_rates)*100, np.array(file_names, dtype=str)


def aggregate_files(dir_name):
    """
    :param dir_name: the name of the directory that contains the relevant files
    :return: list of hard-distribution files, list of posterior files
    """
    files_lst = list(sorted(os.listdir(dir_name)))
    hard_lst, prob_matrix_lst = [], []
    for filename in files_lst:
        full_name = os.path.join(dir_name, filename)
        if filename.endswith(DISTRIBUTION_SUFFIX_HARD):
            hard_lst.append(full_name)
        elif filename.endswith(P_T_GIVEN_R_SUFFIX):
            prob_matrix_lst.append(full_name)
    return hard_lst, prob_matrix_lst


def save_summary_file(args):
    """
    save a csv file of 2 columns: file name , P(target) according to hard method
    """
    relevant_tissues =  get_tissues_by_col_idx(args.groups_file, args.tissues)
    tissue_idx = get_tissue_index(args.tissue_name, relevant_tissues, args.groups_file, args.column_index)
    hard_lst, prob_matrix_lst = aggregate_files(args.directory_name)
    predicted_list_hard, file_names_1 = get_predicted_rates(hard_lst, tissue_idx, DISTRIBUTION_SUFFIX_HARD)
    data_to_save = np.hstack((file_names_1.reshape(-1,1), predicted_list_hard.reshape(-1,1)))
    df_to_save = pd.DataFrame(data_to_save, columns=['file name', 'predicted {} (%)'.format(args.tissue_name)])
    df_to_save.to_csv(os.path.join(args.directory_name, "deconv_summary.csv"), sep=',', encoding='utf-8', index=False)
    return


if __name__ == '__main__':
    start_time = time.time()
    args = parse_args(add_target=True)
    validate_tissue_name(args)
    save_summary_file(args)
    #print('time:', timedelta(seconds=time.time() - start_time), file=sys.stderr)
