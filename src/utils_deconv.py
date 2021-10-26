import os.path as op
import pandas as pd
import numpy as np
import subprocess
from io import StringIO
import os

PAT_FILE_COLUMNS = ['chr', 'start', 'pat', 'count']
BED_FILE_COLUMNS = ['chr', 'start_bp', 'end_bp']
SAVE_FACTOR = 50000  # the maximum is (2**16)-1
PAT_GZ_SUFFIX = ".pat.gz"
TRANSITIONS_FILE_PATTERN = "transitions_matrix_binary_MC{}_g_{}"

class IllegalArgumentError(ValueError):
    pass


def create_directory_if_does_not_exist(out_dir):
    """
    check if the given output directory exists, if not- create this directory
    :param out_dir: a path to the output directory
    """
    if not (os.path.exists(out_dir)):
        os.mkdir(out_dir)


def is_file_exists(fpath):
    """
    Check if the given file exists, if not - raise error
    """
    if not op.isfile(fpath):
        raise IllegalArgumentError('No such file: {}'.format(fpath))


def load_np_array(path, shape=None, factor=SAVE_FACTOR, load_type=np.float32):
    """
    load the given array from the given binary file of type uint16
    :param path: input file
    :param shape: tuple of the shape of the loaded data
    :param factor: the save-load factor
    :param load_type - the type of the data to load
    :return: the loaded np-array of probabilities
    """
    is_file_exists(path)
    if shape is None:
        data = np.fromfile(path, dtype=np.uint16)
    else:
        data = np.fromfile(path, dtype=np.uint16).reshape(shape)
    distribution = (data / factor).astype(load_type)
    return distribution


def save_np_array(np_array, output_file, factor=SAVE_FACTOR):
    """
    save the given array to binary file in type of uint16
    :param np_array: numpy array to be saved
    :param output_file: the path to save in the given np-array
    :param factor: the save-load factor
    """
    file_obj = open(output_file, mode='wb')
    (np_array * factor).astype(dtype=np.uint16).tofile(file_obj)
    file_obj.close()
    return


def get_tissues_by_col_idx(tissues_file, group_idx):
    """
    :param tissues_file: a path to the file that contains the names of the tissues (sorted), according to
                            their order in self.probabilities_matrix and self.transitions_matrix.
    :param group_idx: the index of the column that represent the relevant tissues division.
    :return: numpy array that contains the names of the tissues (sorted) that appear in the given column index.
    """
    is_file_exists(tissues_file)
    df = pd.read_csv(tissues_file, header=0, sep=',').iloc[:, group_idx]
    tissues = sorted(np.unique(df))
    return tissues


def get_tissue_index(tissue_name, relevant_tissues, tissue_file, group_idx):
    """
    :param tissue_name: a name of tissue
    :param group_idx: the index of the column that represent the relevant tissues division.
    :param relevant_tissues: a list of tissues (part/all of the tissues that appears in the given tissue_file)
    :param tissue_file: path to a file that contains tissues name
    :return: if the relevant_tissues list and the list in the given file has the same size - return the index
                of the given tissue_name in the list of the tissues from the given file.
            else- sort the relevant_tissues according to there order in the given file and then return the index
                of the given tissue_name in the sorted list.
    """
    tissues_list = get_tissues_by_col_idx(tissue_file, group_idx)
    tissues_list = [x.lower() for x in tissues_list]
    tissue_idx = tissues_list.index(tissue_name.lower())
    if len(relevant_tissues) != len(tissues_list):
        tissue_idx_list = [tissues_list.index(t_name.lower()) for t_name in relevant_tissues]
        tissue_idx_list.sort()
        tissue_idx = tissue_idx_list.index(tissue_idx)
    return tissue_idx


def read_shell(command, shell=True, **kwargs):
    """
    Takes a shell command as a string and and reads the result into a Pandas DataFrame.
    Additional keyword arguments are passed through to pandas.read_csv.
    :param command: a shell command that returns tabular data
    :type command: str
    :param shell: passed to subprocess.Popen
    :type shell: bool
    :return: a pandas data frame. if the shell command returns an empty output - return None.
    :rtype: :class:`pandas.dataframe`
    """

    if command is None:
         return
    proc = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = proc.communicate()

    if len(output) == 0:
        err = error.decode()
        if err.strip():
            print(err)
        return None

    if proc.returncode == 0:
        with StringIO(output.decode()) as buffer:

            return pd.read_csv(buffer, **kwargs)
    else:
        message = ("Shell command returned non-zero exit status: {0}\n\n"
                   "Command was:\n{1}\n\n"
                   "Standard error was:\n{2}")
        raise IOError(message.format(proc.returncode, command, error.decode()))
