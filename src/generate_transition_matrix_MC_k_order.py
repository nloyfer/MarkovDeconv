from argparse import ArgumentParser
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from utils_deconv import TRANSITIONS_FILE_PATTERN

BASE = 2
GZ_FILE_SUFFIX = ".m{}.gz"
PSEUDO_COUNT_ANY_GROUP = 1
DEFAULT_GROUP = 1
SAVE_FACTOR = 50000
OUTPUT_TRANSITIONS_DIR = None
INPUT_TISSUES_FILE = None
INPUT_SUFF_STATS_DIR = None 


def parse_args():
    p = ArgumentParser()
    g = p.add_argument_group('I/O')

    g.add_argument('--output_transitions_dir', '-otrans', type=str,
                   help="the directory in which the transitions matrix files will written.")  # out_dir , o , '.'

    g.add_argument('--input_tissues_file', '-in_tissue', type=str,
                   help="a path to a file in which the samples files are mapping to the tissues")  # group_file, g  , group.csv

    g.add_argument('--column_index', '-g', type=int, default=DEFAULT_GROUP,
                   help='the index of the column that represent the relevant tissues division')  # remove - always 1

    g.add_argument('--input_directory', '-i', type=str,
                   help='a path to the directory which contains files of the sufficient statistics of the model.')  #  remove == out_dir

    g.add_argument('--pseudo_count', '-pseudo', type=int, default=PSEUDO_COUNT_ANY_GROUP,
                   help='a pseudo count for any possible group of k+1 nucleotides')  # pseudo_count'

    g.add_argument('--file_suffix', '-suffix', type=str, default=GZ_FILE_SUFFIX,
                   help='the suffix of the files in the input_directory which contains'
                        'the sufficient statistics of the model. default: .m{k}.gz')  # remove

    # markers_file , L , requiered = True, no default
    # input_dir, -i,   requiered = True, no default

    g.add_argument('--markov_order', '-k', type=int, help='the order of the Markov Chain model')  # default 4, create all 0...k
    args = p.parse_args()
    return args


def get_groups_of_tissues(tissues_groups_path, group_idx):
    """
    :param tissues_groups_path - a path to a csv file (separated by ',') in which each sample from the first column
                        map to some tissue (according to the group_idx column)
    :param group_idx - the index of the column (in tissues_groups_path) that represent the relevant tissues division.
    :return: groups_df - a pandas data frame (size #tissuesX2) in which the indexes are the tissues (unique).
                        The value in the first column for index i is a list of samples that related to the i'th tissue.
    """
    map_df = pd.read_csv(tissues_groups_path, header=0, sep=',').iloc[:, [0, group_idx]]
    groups_df = map_df.groupby(map_df.columns[1]).agg(lambda x: x.tolist())
    groups_df.sort_index(inplace=True)
    return groups_df


def generate_transitions_matrix(dir_input, k, pseudo, file_suffix, tissues_groups_path, group_idx):
    """
    create a matrix of the Markov model parameters, using multi-processing (process for each tissue)
    :param dir_input: a path to a directory that contains files with the sufficient statistics of the model.
                        In each file the order of the columns is a binary order and each row related to some block.
    :param k: the order of the markov model
    :param file_suffix: the suffix of the files that related to the given model_order ('k')
                            and contain the sufficient statistics of the model.
    :param pseudo: a pseudo count to add for each group of k+1 nucleotides
    :param tissues_groups_path - a path to a csv file (separated by ',') in which each sample from the first column
                        map to some tissue (according to the group_idx column)
    :param group_idx - the index of the column (in tissues_groups_path) that represent the relevant tissues division.
    :return: a matrix of the Markov model parameters, type: float32, size: (#tissues)X(2**k)X(#blocks)
                and an array of tissues in the same order of the first dimension of the returned matrix
    """
    groups_df = get_groups_of_tissues(tissues_groups_path, group_idx)
    processes = []
    n_processes = 7 if k == 5 else 16  # TODO - verify with Netanel. Do we need to add parameters for this?
    with Pool(n_processes) as p:
        for i, group in enumerate(np.array(groups_df.iloc[:, 0])):
            processes.append(p.apply_async(thread_function, (i, group, k, pseudo, dir_input, file_suffix)))
        p.close()
        p.join()

    prob_matrix = [pr.get() for pr in processes]
    return np.array(prob_matrix).astype(np.float32) if prob_matrix != [] else None


def get_full_file_names(dir_input, file_suffix, samples_prefix):
    """
    :param dir_input: a path to a directory that contains files with the sufficient statistics of the model.
                        In each file the order of the columns is a binary order and each row related to some block.
    :param file_suffix: the suffix of the files that related to the given model_order ('k')
                            and contains the sufficient statistics of the model.
    :param samples_prefix: a list of samples that related to the same tissue
    :return: a list of (sufficient statistics) files that exist in the given directory (dir_input),
                and have the given file_suffix, and have a prefix that appears in the given samples_prefix list.
    """
    all_files_lst = list(sorted(os.listdir(dir_input)))
    relevant_file_list = []
    for prefix in samples_prefix:
        for file_name in all_files_lst:
            if file_name.endswith(file_suffix) and file_name.startswith(prefix):
                relevant_file_list.append(file_name)
                break
    return relevant_file_list


def thread_function(process_id, samples_lst, k, pseudo, dir_input, file_suffix):
    """
    The process task
    :param process_id: the id of the process
    :param samples_lst: a list of samples that related to the same tissue
    :param k: the order of the markov model
    :param pseudo:  a pseudo count to  add for each group of k+1 nucleotides
    :param dir_input: a path to a directory that contains files with the sufficient statistics of the model.
                        In each file the order of the columns is a binary order and each row related to some block.
    :param file_suffix: the suffix of the files that related to the given model_order ('k')
                            and contains the sufficient statistics of the model.
    :return: a matrix of the Markov model transition parameters for a specific tissue,
                type: float32, size: (2**k)X(#blocks)
    """
    # print('in process', process_id)
    ss_files_lst = get_full_file_names(dir_input, file_suffix, samples_lst)

    col_numbers = np.arange(2**(k+1)) + BASE
    num_of_blocks = get_num_of_blocks(dir_input, ss_files_lst[0], col_numbers)
    current_counts = np.zeros((num_of_blocks, 2**(k+1)), dtype=np.float32)
    pseudo_counts_vec = np.full(2 ** (k+1), pseudo)
    # sum up the all sufficient statistics from the files that related to the same tissue, and to the k model
    for filename in ss_files_lst:
        filename = os.path.join(dir_input, filename)
        df = pd.read_csv(filename, header=None, sep='\t', usecols=col_numbers, dtype=np.float32)
        current_counts += np.array(df)

    total_tabel = current_counts[:, 0::2] + current_counts[:, 1::2] + (pseudo * 2)  # size: num_of_blocks, 2**k
    denominator = np.zeros((num_of_blocks, 2**(k+1)), dtype=np.float32)
    denominator[:, 0::2], denominator[:, 1::2] = total_tabel, total_tabel

    numerators = (current_counts + pseudo_counts_vec).T
    cur_probs = (numerators / denominator.T).astype(np.float32)  # size ( 2**(k+1))X(#blocks)
    return cur_probs[0::2, :]  #  size (2**k)X(#blocks)  # only P(T|...)


def get_num_of_blocks(dir_input, filename, col_numbers):
    """
    This function assumes that in each sufficient statistics file from the given dir_input directory
    there is the same number of blocks (rows).
    :param dir_input: a path to a directory that contains files with the sufficient statistics of the model.
                        In each file the order of the columns is a binary order and each row related to some block.
    :param filename: a path to some sufficient statistics file.
    :param col_numbers: the columns in the given 'filename' file that contains the sufficient statistics values
    :return: the number of blocks (rows) in the sufficient statistics files that exist in the given dir_input directory
    """
    filename = os.path.join(dir_input, filename)
    df = pd.read_csv(filename, header=None, sep='\t', usecols=col_numbers, dtype=np.float32)
    return df.shape[0]


def save_transitions_matrix_file(prob_output_file, prob_matrix):
    """
    This function saves the given prob_matrix to given prob_output_file in binary format (dtype=np.uint16)
    :param prob_output_file - the path of the output file
    :param prob_matrix - the probabilities matrix to be saved in the given output file
    :return: None
    """
    with open(prob_output_file, mode='wb') as file_obj:
        (prob_matrix*SAVE_FACTOR).astype(dtype=np.uint16).tofile(file_obj)


def validate_matrix(transitions_mat):
    """
    verify that the values in the given matrix are between 0 to 1 (including)
    """
    assert (np.all(transitions_mat <= 1))
    assert (np.all(transitions_mat >= 0))


if __name__ == '__main__':
    args = parse_args()
    suffix_of_files = args.file_suffix.format(args.markov_order)
    transitions_matrix = generate_transitions_matrix(args.input_directory.format(args.markov_order), args.markov_order,
                                                     args.pseudo_count, suffix_of_files, args.input_tissues_file, args.column_index)
    validate_matrix(transitions_matrix)
    output_file = os.path.join(args.output_transitions_dir, TRANSITIONS_FILE_PATTERN.format(args.markov_order, args.column_index))
    save_transitions_matrix_file(output_file, transitions_matrix)
