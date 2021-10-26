import gc
import math
from k_order_MC_counting_strict_only import k_order_MC_counting_strict_only as k_order_MC
import numpy as np
import pandas as pd
from multiprocessing import Pool
from time import sleep
from utils_deconv import get_tissues_by_col_idx, IllegalArgumentError, is_file_exists, BED_FILE_COLUMNS

NUM_OF_CHROMOSOMES = 22
EXTRA_CHR = ['M', 'X', 'Y']
CHROMOSOME_PREFIX = "chr"
DELAY_FACTOR = 1.3


def get_delay_per_model(model_order, add_delay):
    """
    :param model_order: 1 for first order markov chain
    :param add_delay: boolean flag- if true add delay, else define the delay as 0
    :return: the delay for this model in seconds
    """
    return 0 if not add_delay else math.floor(model_order * DELAY_FACTOR) if model_order < 5 else 6.5


def is_empty_region(markers_path, region):
    """
    :param markers_path:
    :param region:
    :return: True if there are no markers in the given markers_path that related to the given region
    """
    is_file_exists(markers_path)
    markers_df = pd.read_csv(markers_path, header=None, sep='\t', comment='#', names=BED_FILE_COLUMNS,
                             usecols=list(range(len(BED_FILE_COLUMNS))))
    region_df = markers_df[markers_df['chr'] == region]
    return region_df.empty


def get_a_matrix_of_p_of_tissue_given_read(args, tissue_file, relevant_tissues, transitions_matrix_file,
                                           t_distribution, min_len, pat_path):
    """
    :param args: ArgumentParser object
    :param relevant_tissues: the names of the tissues that we want to check the probabilities for them.
    :param pat_path: a path to a pat file that contains the reads that we want to check the probabilities for them.
    :param min_len: the minimum length of read for the analysis
    :param tissue_file: a path to the file that contains the names of the tissues (sorted), according to
                                their order in probs_matrix_file and trans_matrix_file.
    :param transitions_matrix_file: the file that contains the transition probabilities matrix. (The i'th element is
                                    a transition matrix of size (2**model_orderX#blocks), related to the i'th tissue.
    :param t_distribution: a prior vector, the i'th element is the probability of the i'th tissue.
    :return: numpy array of size : (reads_number, tissue_number).
                if is_p_t_given_r==True, The [j,i] element is the log probability of the j'th tissue given the i'th read
                else, The [j,i] element is the log probability of the j'th read given the i'th cell type (tissue)
    """
    tissues_list = get_tissues_by_col_idx(tissue_file, args.column_index)
    processes = []
    with Pool(10) as p:  # number of processes running at time is 10 maximum
        all_chr = list(range(NUM_OF_CHROMOSOMES)) + EXTRA_CHR
        for i, val in enumerate(all_chr):
            region = CHROMOSOME_PREFIX + str(i + 1) if i < NUM_OF_CHROMOSOMES else CHROMOSOME_PREFIX + val
            if is_empty_region(args.markers_bed_path, region):
                continue
            processes.append(p.apply_async(thread_function, (args.model_order, i + 1, pat_path, tissues_list,
                                                             relevant_tissues, region, min_len, transitions_matrix_file,
                                                             t_distribution, args.markers_bed_path, args.genome)))
            # delay in order to prevent a situation of loading the transition file  in all processes simultaneously:
            sleep(get_delay_per_model(args.model_order, args.delay))
        p.close()
        p.join()

    gc.collect()
    prob_matrix, counts = [], []
    for j, pr in enumerate(processes):
        return_val = pr.get()
        if return_val is not None:
            prob_matrix.append(return_val[0])
            counts.append(return_val[1])

    if not prob_matrix:
        raise IllegalArgumentError('Probabilities matrix is empty!')
    else:
        return np.concatenate(prob_matrix, axis=0).astype(np.float32), np.concatenate(counts).astype(np.int16)


def thread_function(model_order, process_id, reads_pat_file, tissues_list, relevant_tissues, region,
                    min_len, transitions_matrix_file, t_distribution, markers_file, genome):
    """
    A function that each process perform.
    :param model_order: the order of the markov chain model
    :param relevant_tissues: the names of the tissues that we want to check the probabilities for them.
    :param process_id: the number that represents the current process
    :param min_len: the minimum length of read for the analysis
    :param reads_pat_file: a path to a pat file that contains the reads that we want to check the probabilities for them.
    :param tissues_list: the names of the tissues (sorted), according to their order in probs_matrix_file and transitions_matrix_file.
    :param region: the name of the chromosome for this thread (process)
    :param transitions_matrix_file: the file that contains the transition probabilities matrix.
                            (The i'th element is a matrix of size (4X#blocks), related to the i'th tissue.
                            The [k,j] element in this matrix is the the probability to see the k'th pair in the
                            i'th tissue and in the j'th block. There are 4 pairs: TT, TC, CT, CC)
    :param t_distribution: a prior vector, the i'th element is the probability of the i'th tissue.
    :param markers_file: a bed file that contains 3 columns: name of chromosome, start block in bp, end block in bp.
    :return: the process results -
                return None if there are no reads in the 'reads_pat_file' for the given region
                else - return a tuple of (probabilities_matrix , reads_repetition_col)
                    probabilities matrix - numpy array of size : (reads_number, tissue_number) that related to the
                    reads which are relevant to this process only.
                    if is_p_t_given_r==True, The [j,i] element is the probability of the j'th tissue given the i'th read
                    else, The [j,i] element is the probability of the j'th read given the i'th cell type (tissue)
                    reads_repetition_col- a vector of size number of unique reads (which are relevant to this process
                    only). The i'th element is the number of times the i'th read appears
    """
    MC_model = k_order_MC(reads_pat_file, tissues_list, relevant_tissues, transitions_matrix_file,
                          region, min_len, t_distribution, markers_file, genome, model_order)
    return MC_model.get_matrix_of_log_p_tissue_given_read()

