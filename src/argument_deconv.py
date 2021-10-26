from argparse import ArgumentParser
from utils_deconv import get_tissues_by_col_idx, TRANSITIONS_FILE_PATTERN
import numpy as np
import os

DEFAULT_GROUP = 1
ALL = "all"

TRANSITIONS_FILE_DIR = r"/home/ubuntu/data/wgbs_tools/deconvolution/kidney.400.may/train_trans_params/"
INPUT_TISSUES_FILE = r"groups.csv"


def parse_args(add_pat_path=False, add_markers_file=False, add_min_len=True, add_delay=False, add_target=False, add_plabels=False):
    p = ArgumentParser()
    g = p.add_argument_group('I/O')

    g.add_argument('--genome',  help='genome [hg19]', default='hg19') # model_order , k, default=4
    g.add_argument('--model_order', '-model', type=int, help='the order of the markov chain (int)') # model_order , k, default=4

    g.add_argument('--directory_name', '-output', type=str, help="name of the output directory", default='.')  # out_dir , o , default:"."

    g.add_argument('--tissues', '-tissues_col', type=int,  default=1,
                   help='the index of the column that contains the names of the tissues that we want to '
                        'check the probabilities for them in groups_file')  # remove - always 1 group

    g.add_argument('--prior_tissues_distribution', '-pai', nargs='+', type=float, default=None,
                   help="The distribution over the tissues (prior), corresponding to lexicography order "
                        "of the tissues that appears in groups_file in tissues_col column. "
                        "\n default: uniform distribution")  # prior, -p , uniform

    g.add_argument('--column_index', '-g', type=int, default=DEFAULT_GROUP,
                   help='the index of the column that represent the relevant tissues division in '
                        '/cs/cbio/netanel/tools/markers/groups.csv file')  # # remove

    g.add_argument('--transitions_dir', '-trans', type=str, default=TRANSITIONS_FILE_DIR,
                   help="The directory of the transitions matrix files."
                        " \n default: " + TRANSITIONS_FILE_DIR)  # work_dir, -w, '.'

    g.add_argument('--groups_file', '-groups_f', type=str, default=INPUT_TISSUES_FILE,
                   help="a path to a file in which the samples files are mapping to the tissues"
                        " \n default: " + INPUT_TISSUES_FILE)  # groups_file, g , group.csv

    if add_min_len:
        add_read_min_len_arg(g)
    if add_markers_file:
        add_bed_file_arg(g)
    if add_pat_path:
        add_pat_path_arg(g)
    if add_delay:
        add_delay_arg(g)
    if add_target:
        add_target_arg(g)
    if add_plabels:
        g.add_argument('--labels_path', help='path for predicted labels column dump')

    args = p.parse_args()
    return args


def add_read_min_len_arg(obj):
    obj.add_argument('--read_length', '-len', type=int, help='The minimum length of read for the analysis')  # min_sites, -l,  default =3


def add_pat_path_arg(obj):
    obj.add_argument('--pat_path', '-input', type=str, default=None,
                   help="A path to a directory that contains pat.gz files, or a path to a single pat.gz file")  # input, i


def add_bed_file_arg(obj):
    obj.add_argument('--markers_bed_path', '-L', type=str,
                     help="A path to a bed.gz file that contains markers (DMRs)")  # markers, L, no default


def add_delay_arg(obj):
    obj.add_argument('--delay', '-d', type=int, default=0,
                     help='if true- adding delay in order to prevent a situation of loading the '
                          'transition file in all processes simultaneously in markov_model_analyzer.py')  # remove - no delay


def add_target_arg(obj):
    obj.add_argument('--tissue_name', '-target', type=str, help="name of the target tissue", default='liver')  # -tissue_name', '-t


def get_transitions_path(transitions_files_dir, model, group_idx):
    return os.path.join(transitions_files_dir, TRANSITIONS_FILE_PATTERN.format(model, group_idx))


def parse_tissues_distribution(tissues_distribution, relevant_tissues):
    if not tissues_distribution:
        return np.full(len(relevant_tissues), 1 / len(relevant_tissues), dtype=float)
    else:
        return tissues_distribution


def validate_prior_input(prior_pai, selected_tissues_lst):
    assert (np.round(np.sum(prior_pai), 4) == float(1))
    print('priors: ', prior_pai)
    print('selected tissues lst: ', selected_tissues_lst)
    assert (len(prior_pai) == len(selected_tissues_lst))
