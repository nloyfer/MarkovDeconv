import re
import os
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from utils_deconv import read_shell, is_file_exists, PAT_FILE_COLUMNS, BED_FILE_COLUMNS


def check_executable(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(os.path.join(p, cmd), os.X_OK):
            return True
    eprint(f'executable {cmd} not found in PATH')
    return False

def get_wgbs_tools_path(wtpath=None):
    # user did not specify wgbstools path:
    if not wtpath:
        wtpath = 'wgbstools'

    if wtpath == 'wgbstools':
        if not check_executable(wtpath):
            exit(1)
    elif op.isfile(wtpath):
        if not os.access(wtpath, os.X_OK):
            eprint('Could not run wgbstools:', wtpath)
            exit(1)
    else:
        eprint('Could not run wgbstools:', wtpath)
        exit(1)

    return wtpath

class Markov_Model_strict_only:  # TODO- change the name of the class
    STRICT = '--strict'
    C_BASE, C_AS_STR = 'C', '1'
    T_BASE, T_AS_STR = 'T', '0'
    UNKNOWN_BASE = '.'  # CpG site but maybe methylated and maybe not methylated.

    def __init__(self, pat_file, tissues_list, relevant_tissues, region, read_min_len,
                 tissues_distribution, markers_bed_file, genome):
        """
        :param pat_file: a path to a pat file that contains the reads that we want to check the probabilities for them.
        :param tissues_list: the names of the tissues (sorted), according to their order in probs_matrix_file and trans_matrix_file.
        :param relevant_tissues: the names of the tissues that we want to check the probabilities for them.
        :param region: a name of chromosome.
        :param read_min_len: the minimum length of read for the analysis
        :param tissues_distribution: a vector of size number of tissues. the i'th element is the probability od the i'th tissue.
        :param markers_bed_file: a file that contains 3 columns: name of chromosome, start block in bp, end block in bp.
        """
        self.wt = get_wgbs_tools_path()
        self.min_read_len = read_min_len
        self.region = region
        self.genome = genome
        self.tissues_list = tissues_list
        self.num_of_tissues = len(self.tissues_list)
        self.markers_of_region_df = self.load_markers_file_by_region(markers_bed_file)
        # create a list of the indices of the relevant tissues (self.tissue_indexes_list)
        if len(relevant_tissues) != self.num_of_tissues:
            self.tissue_indexes_list = [self.tissues_list.index(t_name.lower()) for t_name in relevant_tissues]
            self.tissue_indexes_list.sort()
        else:
            self.tissue_indexes_list = None

        self.first_marker_index, self.last_marker_index = self.get_start_and_end_markers_idx()
        self.cur_block = None  # the index of the block regarding to 6M blocks (all of the blocks)
        self.start_cur_block, self.end_cur_block = None, None  # in CpG sites
        self.pat_data_generator = self.load_pat_file_by_chr_and_markers(pat_file, genome)
        self.pat_data = []

        if len(relevant_tissues) != self.num_of_tissues:
            tissue_idx_list = [self.tissues_list.index(t_name.lower()) for t_name in relevant_tissues]
            self.cell_types_distribution = [prob for _, prob in sorted(zip(tissue_idx_list, tissues_distribution))]
        else:
            self.cell_types_distribution = tissues_distribution

        self.num_of_relevant_tissues = len(relevant_tissues)
        self.reads_repetition_col = np.zeros(0, dtype=np.int16)


    ############################################## tissue_given_read ###################################################

    def get_matrix_of_log_p_tissue_given_read(self):
        """
        This function calculates for each pair of (read, tissue) the probability of the tissue given the read
        :return: None if there are no reads in the given pat_file that are located in the blocks from the
                        markers_bed_file and in self.region and the length of them is at least self.min_read_len.
                else- return a tuple of size 2:
                    the first element is a  numpy array of size : (reads_number, tissue_number).
                    The [j,i] element is the probability of the j'th tissue given the i'th read
                    The second element is a numpy array that contains the repetition column of the pat_file,
                    filtered according to self.region, self.min_read_len, and to the markers_bed_file.
        """
        all_probabilities_data = []
        for df in self.pat_data_generator:
            data = df[df[PAT_FILE_COLUMNS.index('pat')].str.len() >= self.min_read_len]
            self.pat_data = np.array(data)
            self.reads_repetition_col = np.hstack((self.reads_repetition_col, self.pat_data[:, PAT_FILE_COLUMNS.index('count')]))
            reads_number = self.pat_data.shape[0]
            probabilities_data = np.zeros((reads_number, self.num_of_relevant_tissues), dtype=float)
            for read_idx in range(reads_number):
                read_row_array = self.pat_data[read_idx, :len(PAT_FILE_COLUMNS)]
                assert (self.start_cur_block <= read_row_array[1] < self.end_cur_block)
                assert((read_row_array[1] + len(read_row_array[2]) - 1) < self.end_cur_block)
                probabilities_data[read_idx, :] = self.log_p_cell_type_given_read(read_row_array[2])
            all_probabilities_data.append(probabilities_data)

        if all_probabilities_data:
            return np.concatenate(all_probabilities_data, axis=0), self.reads_repetition_col
        return None

    def log_p_cell_type_given_read(self, read_pattern):
        """
        :param read_pattern: the string of the read pattern (in CpG sites)
        :return: np-array of size number of tissues.
                    the i'th element is the log probability of the i'th tissue given the given read.
        """
        # Numerators - the i'th element is the numerator of the i'th tissue
        log_p_cell_type_vec = np.log(self.cell_types_distribution)
        log_p_read_given_cell_type_vec = self.log_p_read_given_cell_type(read_pattern, self.cur_block)
        log_numerator = log_p_cell_type_vec + log_p_read_given_cell_type_vec
        # Denominator
        log_denominator = logsumexp(log_numerator)
        return log_numerator - log_denominator

    ############################################## read_given_tissue ###################################################

    def log_p_read_given_cell_type(self, read, block_idx):
        """
        assuming that the given read contained only in one block - the given block.
        :param read: the string of the read
        :param block_idx: the index of the block that contains this read
        :return: np-array of size num of tissues.
                    The i'th element is the the log of the probability of the read given the i'th tissue
        """
        raise NotImplementedError

    def load_pat_file_by_chr_and_markers(self, pat_file, genome):
        """
        :param pat_file: a path to a pat file that contains reads.
        :return: generator of pandas data frames, each value is a data frame that contains reads that related to
                    a certain marker (from the given bed file) and related to the given chromosome (region param)
        """
        is_file_exists(pat_file)
        if self.markers_of_region_df.empty:
            return iter([])
        first_index = self.markers_of_region_df.index[0]
        for index, row in self.markers_of_region_df.iterrows():
            cur_marker = '{}:{}-{}'.format(row['chr'], row['start_bp'], row['end_bp'])
            cur_cmd = f'{self.wt} cview --strip {pat_file} --genome {genome} -r {cur_marker} {self.STRICT}' # TODO cmd of mix!!! rate!
            df_of_markers = read_shell(cur_cmd, header=None, index_col=False, sep='\t')
            if df_of_markers is None:
                continue

            self.cur_block = index - first_index  # in order to start the counting from zero
            self.start_cur_block, self.end_cur_block = self.update_start_end_block_site(row['chr'], row['start_bp'], row['end_bp'], genome)  # TODO - use also cols 3,4 (ask Netanel if can we assume this)
            yield df_of_markers

    def update_start_end_block_site(self, chrom, start_bp, end_bp, genome):
        """
        :param chrom: the chromosome of the marker
        :param start_bp: start block in bp
        :param end_bp: end block in bp
        :return: the start of the block (in site), the end of the block (in site)
        """
        marker_as_region = f'{chrom}:{start_bp}-{end_bp}'
        cmd = f'{self.wt} convert --genome {genome} --no_anno -r {marker_as_region}'
        str_output = read_shell(cmd, header=None, sep='\t')[0][0]
        start, end = re.split('-',re.findall(r'(\s\d+-\d+)',str_output)[0])
        return int(start), int(end)

    def get_start_and_end_markers_idx(self):
        """
        :return: 2 values.
                (1) The index of the first marker that related to self.region,
                (2) The index of the last marker that related to self.region (including).
                The indexes regards to the order in the original markers file.
                We don't need to save the all markers' indexes that related to the self.region because the markers'
                indexes are continuous (the markers file is sorted by chromosome and then by the location in the chromosome)
                If there are no markers that related to the self.region chromosome return 0, 0.
        """
        assert (True not in np.array(self.markers_of_region_df.duplicated(subset=None, keep='first')))
        return (self.markers_of_region_df.index[0], self.markers_of_region_df.index[-1]) if not self.markers_of_region_df.empty else (0, 0)

    def load_markers_file_by_region(self, markers_path):
        """
        :param markers_path: a path to the file that contains the markers. Each row related to
                    a marker-block, format: chromosome, start location, end location, block_index.
        :return: pandas data frame of size . len(BED_FILE_COLUMNS) X #markers that came from self.region
                    The columns are BED_FILE_COLUMNS, and the rows are the markers that related to self.region
        """
        is_file_exists(markers_path)
        markers_df = pd.read_csv(markers_path, header=None, sep='\t', comment='#', names=BED_FILE_COLUMNS, usecols=list(range(len(BED_FILE_COLUMNS))))
        return markers_df[markers_df['chr'] == self.region]
