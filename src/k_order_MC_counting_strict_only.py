import gc
import re
import pandas as pd
import numpy as np
from utils_deconv import is_file_exists, BED_FILE_COLUMNS
from Markov_Model_strict_only import Markov_Model_strict_only


class k_order_MC_counting_strict_only(Markov_Model_strict_only):  # TODO- change the name of the class
    LOAD_FACTOR = 50000
    ZERO_CORRECT = 0.00001
    TRANSITION_FILE_SEPARATOR = "MC"

    def __init__(self, pat_file, tissues_list, relevant_tissues, transitions_file, region,
                 read_min_len, tissues_distribution, markers_bed_file, genome, k):
        """
        :param pat_file: a path to a pat file that contains the reads that we want to check the probabilities for them.
        :param tissues_list: the names of the tissues (sorted), according to their order in probs_matrix_file and trans_matrix_file.
        :param relevant_tissues: the names of the tissues that we want to check the probabilities for them.
        :param transitions_file: the file that contains the probabilities matrix
                            (The i'th element is a matrix of size (4X#blocks), related to the i'th tissue.
                            The [k,j] element in this matrix is the the probability to see the k'th pair in the
                            i'th tissue and in the j'th block. There are 4 pairs: TT, TC, CT, CC)
        :param region: a name of chromosome.
        :param read_min_len: the minimum length of read for the analysis
        :param tissues_distribution: a vector of size number of tissues. the i'th element is the probability od the i'th tissue.
        :param markers_bed_file: a file that contains 3 columns: name of chromosome, start block in bp, end block in bp.
        :param k: the order of the markov chain
        """
        self.k = k
        # invoking the __init__ of the parent class
        super().__init__(pat_file, tissues_list, relevant_tissues, region, read_min_len, tissues_distribution, markers_bed_file, genome)
        # initialize the probabilities matrices
        self.transitions_matrix_array = self.generate_transitions_matrix_array(transitions_file)
        # validate that the blocks dimension in transitions matrix is equal to the number of markers in markers_bed_file
        blocks_num = (self.transitions_matrix_array[0]).shape[2]
        markers_df = pd.read_csv(markers_bed_file, header=None, sep='\t', comment='#', names=BED_FILE_COLUMNS, usecols=list(range(len(BED_FILE_COLUMNS))))
        region_df = markers_df[markers_df['chr'] == self.region]
        markers_num = np.array(region_df).shape[0]
        assert (markers_num == blocks_num)

    def get_log_of_all_transition_vec(self, k, block_idx):
        """
        :param k: the order of the model
        :param block_idx: the index of the block that contains the read
        :return: a vector of size (2**(k+1), #tissues) that contains the transitions probabilities to move from
                some k sites to some specific site in the given block and in all tissues
        """
        trans_prob_vec = self.transitions_matrix_array[k][:, :, block_idx].T  # shape: (2**k, #tissues)
        log_trans_prob_vec_all = np.zeros((trans_prob_vec.shape[0] * 2, trans_prob_vec.shape[1]), dtype=np.float32)
        log_trans_prob_vec_all[0::2, :] = np.log(trans_prob_vec)  # T given
        log_trans_prob_vec_all[1::2, :] = np.log(1 - trans_prob_vec)  # C given
        return log_trans_prob_vec_all

    def get_log_prob_of_first_site(self, block_idx, first_nuc):
        """
        :param block_idx: the index of the block that contains the read
        :param first_nuc: the first site of the read "C" or "T"
        :return: the log probability to see the first site in the given block and in all tissues
        """
        cur_prob_k_0 = self.transitions_matrix_array[0][:, 0, block_idx]  # vector of size num of tissues
        return np.log(1 - cur_prob_k_0) if first_nuc == self.C_BASE else np.log(cur_prob_k_0)

    def get_log_prob_of_1_to_k_sites(self, binary_read, read_len, block_idx):
        """
        :param binary_read: the string of the read where "T" replaced by "0" and "C" replaced by "1"
        :param read_len: the length of the read
        :param block_idx: the index of the block that contains the read
        :return: the log probability of the first 1,..,k sites (includes k):
                    log(P(read[1]|read[0])*P(read[2]|read[0:2])*...)
        """
        cur_log_prob = np.zeros(self.transitions_matrix_array[self.k].shape[0])
        i = 1
        while i < self.k and i < read_len:
            loc = int(binary_read[:i + 1], 2)
            log_trans_prob_vec_all = self.get_log_of_all_transition_vec(i, block_idx)
            cur_log_prob += log_trans_prob_vec_all[loc, :]
            i += 1
        return cur_log_prob

    def log_p_read_given_cell_type(self, read, block_idx):
        """
        assuming that the given read contained only in one block - the given block.
        :param read: the string of the read
        :param block_idx: the index of the block that contains this read
        :return: np-array of size num of tissues.
                    The i'th element is the the log of the probability of the read given the i'th tissue
        """
        read_len = len(read)
        assert (read_len != 0)
        binary_read = read.replace(self.C_BASE, self.C_AS_STR).replace(self.T_BASE, self.T_AS_STR)
        cur_log_prob = np.zeros(self.transitions_matrix_array[self.k].shape[0])  # shape: number of tissues
        # adding the log probability of the first site
        if read[0] != self.UNKNOWN_BASE:
            cur_log_prob += self.get_log_prob_of_first_site(block_idx, read[0])
        # adding the log probability of the first 1,..,k sites (includes k): P(read[1]|read[0])*P(read[2]|read[0:2])*...
        read_prefix = read[:self.k]  # if (read_len >= self.k) else read
        if self.UNKNOWN_BASE not in read_prefix:
            cur_log_prob += self.get_log_prob_of_1_to_k_sites(binary_read, read_len, block_idx)
        if read_len <= self.k:
            return cur_log_prob
        # create the counts vector of sequences of length 2**(self.k+1)
        log_trans_prob_vec_all = self.get_log_of_all_transition_vec(self.k, block_idx)
        read_known_parts = re.split(r"[.]+", binary_read)
        counts = np.zeros(log_trans_prob_vec_all.shape[0])
        for cur_read in read_known_parts:
            if len(cur_read) <= self.k:
                continue
            for i in range(self.k, len(cur_read)):
                loc = int(cur_read[i-self.k:i+1], 2)
                counts[loc] += 1
        # calculate the log probability for i in 2**(self.k+1): cur_log_prob += (counts[i]*log_trans_prob_vec[i, :])
        cur_log_prob += np.sum(np.multiply(counts.reshape(-1, 1), log_trans_prob_vec_all), axis=0)
        return cur_log_prob

    def load_transitions_matrix_file(self, prob_file, k, marker_indexes):
        """
        This function loads the transitions probabilities matrix from the given prob_file
        :param k: the order of the Markov model
        :param marker_indexes: a list of indexes of the markers that related to the self.region (chromosome)
        :param prob_file: a binary file that contains the transitions probabilities (size: #Tissues, 2**(k), #Blocks)
        :return: the transitions probabilities matrix. The [i,:,j] element is the probabilities of each transition from
                some k sites to some specific site in the i'th tissue and the j'th block - P(Xi|Xi-k, ... , X_i-1)
                For example, the order for k=2 is: P(T|TT), P(C|TT), P(T|TC), P(C|TC), P(T|CT), P(C|CT), P(T|CC), P(C|CC).
                the order for k=3 is: P(T|TTT), P(C|TTT), P(T|TTC), P(C|TTC), P(T|TCT), P(C|TCT), P(T|TCC), P(C|TCC),
                                    P(T|CTT), P(C|CTT), P(T|CTC), P(C|CTC), P(T|CCT), P(C|CCT), P(T|CCCT), P(C|CCC)
        """
        is_file_exists(prob_file)
        data = np.fromfile(prob_file, dtype=np.uint16).reshape((self.num_of_tissues, 2**k, -1))[:,:,marker_indexes]
        if self.tissue_indexes_list is None:
            new_data = ((data / self.LOAD_FACTOR).astype(np.float32))
        else:
            new_data = ((data/self.LOAD_FACTOR).astype(np.float32))[self.tissue_indexes_list, :, :]
        data = None
        gc.collect()
        assert (np.all(new_data < 1))
        j1, j2, j3 = np.where(new_data == 0)
        new_data[j1, j2, j3] = self.ZERO_CORRECT
        return new_data

    def generate_transitions_matrix_array(self, trans_file):
        """
        :param trans_file: a binary file that contains the transitions probabilities (size: #Tissues, 2**(self.k), #Blocks)
        :return: array of size (self.k). the m element is a matrix of size
                (#Tissues, 2**(k), #Blocks). The [i,:,j] element of the m element is the probabilities of each
                 transition from some m sites to some specific site in the i'th tissue and the j'th block.
        """
        if self.markers_of_region_df.empty:
            marker_indexes = list()
        else:
            marker_indexes = list(range(self.first_marker_index, self.last_marker_index+1))  # +1 to include the last index
        all_transitions_matrices = [None]*(self.k+1)
        k_transition_matrix = self.load_transitions_matrix_file(trans_file, self.k, marker_indexes)
        loc = self.k
        all_transitions_matrices[loc] = k_transition_matrix
        loc -= 1
        cur_k = self.k
        while loc >= 0:
            trans_file_pref, trans_file_suff = trans_file.split(self.TRANSITION_FILE_SEPARATOR)
            trans_file = trans_file_pref + self.TRANSITION_FILE_SEPARATOR + str(cur_k-1) + trans_file_suff[1:]
            cur_trans_matrix = self.load_transitions_matrix_file(trans_file, cur_k-1, marker_indexes)
            all_transitions_matrices[loc] = cur_trans_matrix
            loc -= 1
            cur_k -= 1
        return all_transitions_matrices


