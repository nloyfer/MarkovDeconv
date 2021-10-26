import gc
import pandas as pd
import os
import sys
import time
from datetime import timedelta
import numpy as np
from markov_model_analyzer import get_a_matrix_of_p_of_tissue_given_read
from argument_deconv import parse_args, get_transitions_path, parse_tissues_distribution, validate_prior_input
from utils_deconv import save_np_array, create_directory_if_does_not_exist, PAT_GZ_SUFFIX, get_tissues_by_col_idx

DISTRIBUTION_SUFFIX_HARD = "_tissue_distribution_by_hard_assignment.binary"
P_T_GIVEN_R_SUFFIX = "_p_of_tissue_given_read_matrix.binary"
READS_INFO_SUFFIX = "_reads_info.binary"


def get_tissues_counter_vector(src_vec, rep_vec, number_of_tissues):
    """
    :param src_vec: np array of size: number of reads.
                    the i'th element is the index of the tissue that the i'th read came from
    :param rep_vec: np array of size: number of reads.
                    the i'th element is the number of the times the i'th read appears
    :param number_of_tissues: number of tissues
    :return: np array of size: number of tissue.
                the i'th element is the number of reads that came from the i'th tissue (including read repetition)
    """
    counter_vec = np.zeros(number_of_tissues, dtype=np.int)
    tissue_idx_list = np.unique(src_vec)
    for idx in tissue_idx_list:
        reads_from_target = src_vec == idx
        tissue_counter = np.sum(reads_from_target.reshape(-1, 1) * rep_vec)  # numerator
        counter_vec[idx] = tissue_counter
    return counter_vec


def get_tissue_distribution_by_hard_assignment(prob_data, reads_repetition):
    """
    :param prob_data: numpy array of size : (reads_number, tissue_number).
                The [j,i] element is the probability of the j'th tissue given the i'th read
    :param reads_repetition: a vector of size number of unique reads.
                                the i'th element is the number of times the i'th read appears
    :return: np array of size: number of tissue.
                the distribution over the tissues according to hard assignment logic.
    """
    total_reads = np.sum(reads_repetition)  # denominator
    read_src_vec = np.argmax(prob_data, axis=1).reshape(-1, 1)
    tissues_counter_vec = get_tissues_counter_vector(read_src_vec, reads_repetition, prob_data.shape[1])  # numerator
    assert (int(np.sum(tissues_counter_vec)) == int(total_reads))
    return tissues_counter_vec / total_reads


def validate_prob_matrix(probabilities_matrix):
    """
    validate the given matrix
    :param probabilities_matrix: matrix of P(tissues|read)
    """
    assert ((0.99999 < np.sum(probabilities_matrix, axis=1)).all())
    assert ((np.sum(probabilities_matrix, axis=1) < 1.00001).all())
    assert 1.00001 >= np.max(probabilities_matrix) >= np.min(probabilities_matrix) >= 0.0


def save_distribution_vectors(args, tissues_names_lst, transitions_path, t_distribution):
    """
    This function save the matrix of P(tissues|read) and the tissue distribution (hard assignment).
    :param args: ArgumentParser object
    :param tissues_names_lst: the names of the tissues that we want to check the probabilities for them.
    :param transitions_path: the file that contains the transition probabilities matrix. (The i'th element is
                                    a transition matrix of size (2**model_orderX#blocks), related to the i'th tissue.
    :param t_distribution: a prior vector, the i'th element is the probability of the i'th tissue.
    """
    log_probabilities_matrix, reads_repetition = \
        get_a_matrix_of_p_of_tissue_given_read(args, args.groups_file, tissues_names_lst, transitions_path,
                                               t_distribution, args.read_length, args.pat_path)
    probabilities_matrix = np.exp(log_probabilities_matrix)  # the sum of each row is 1 # the sum of the all matrix is the number of read
    log_probabilities_matrix = None
    gc.collect()
    validate_prob_matrix(probabilities_matrix)
    # save P(tissues|read)
    file_name = os.path.basename(args.pat_path).split(PAT_GZ_SUFFIX)[0]
    output_matrix = os.path.join(args.directory_name, str(file_name) + P_T_GIVEN_R_SUFFIX)
    save_np_array(probabilities_matrix, output_matrix)

    # tissue distribution vector
    distribution = get_tissue_distribution_by_hard_assignment(probabilities_matrix, reads_repetition.reshape(-1,1))
    output_distribution = os.path.join(args.directory_name, str(file_name) + DISTRIBUTION_SUFFIX_HARD)
    save_np_array(distribution, output_distribution)
    # TODO- save it to csv file (see example below ) and not in binary file
    # df_to_save = pd.DataFrame(data_to_save, columns=['file name', 'predicted {} (%)'.format(args.tissue_name)])
    # df_to_save.to_csv(os.path.join(args.directory_name, "deconv_summary.csv"), sep=',', encoding='utf-8', index=False)
    return probabilities_matrix


if __name__ == '__main__':
    start_time = time.time()
    args = parse_args(add_pat_path=True, add_markers_file=True, add_delay=True, add_plabels=True)
    assert (not os.path.isdir(args.pat_path))
    create_directory_if_does_not_exist(args.directory_name)
    selected_tissues_lst = get_tissues_by_col_idx(args.groups_file, args.tissues)
    trans_matrix_file = get_transitions_path(args.transitions_dir, args.model_order, args.column_index)
    prior_pai = parse_tissues_distribution(args.prior_tissues_distribution, selected_tissues_lst)
    validate_prior_input(prior_pai, selected_tissues_lst)
    probs = save_distribution_vectors(args, selected_tissues_lst, trans_matrix_file, prior_pai)
    df = pd.DataFrame(columns=selected_tissues_lst, data=probs)
    df['label'] = np.array(selected_tissues_lst)[np.argmax(probs, axis=1)]
    df['label'].to_csv(args.labels_path, sep='\t', header=None, index=None)
    print('time:', timedelta(seconds=time.time() - start_time), file=sys.stderr)
