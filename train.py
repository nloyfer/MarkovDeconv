#!/usr/bin/python3 -u

import os
import sys
import os.path as op
import argparse
import subprocess
import pandas as pd
from multiprocessing import Pool


def eprint(*args,  **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_executable(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(os.path.join(p, cmd), os.X_OK):
            return True
    eprint(f'executable {cmd} not found in PATH')
    return False

def get_wgbs_tools_path(wtpath):
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


def pat_csi(pat, wgbs_tools):
    if not op.isfile(pat):
        eprint(f'File {pat} does not exists')
        exit()
    if not op.isfile(pat + '.csi'):
        eprint(f'Warning: file {pat} is not indexed. attempting to index it...')
        cmd = f'{wgbs_tools} index {pat}'
        subprocess.check_call('{wt} index {p}'.format(wt=wgbs_tools, p=pat), shell=True)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('markers_path', help='path to markers file. columns 4-5 must be the CpG-Index of the markers')
    parser.add_argument('--out_dir', '-o', help='path to output directory. Default is the name of the markers file' )
    parser.add_argument('--counter_tool', help='path to counter tools. Default: ./counter/counter', default='counter/counter')
    parser.add_argument('--deconv_src', help='path to source code directory with the deconv python files. Default: ./src', default='src')
    parser.add_argument('--groups_file', '-g', help='path to groups.csv file. default: ./groups.csv', default='groups.csv')
    parser.add_argument('--genome', help='genome reference name. default: hg19', default='hg19')
    parser.add_argument('--reference_data', '-r', help='directory of pat files of the reference data', default='./reference_data/')
    parser.add_argument('--wgbstools', '-wt', help='path to wgbstools, if it\'s not in PATH')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads for the counter [4]')
    parser.add_argument('--force', '-f', action='store_true', help='delete existing output directory if exists')
    return parser.parse_args()


def prefix2filepath(input_dir, prefix):
    res = None
    for f in sorted(os.listdir(input_dir)):
        if f.startswith(prefix) and f.endswith('.pat.gz'):
            if res is not None:
                eprint(f'Error ambigous file name {prefix} in groups file')
                exit()
            res = op.join(input_dir, f)
    if res is None:
        eprint('Couldn\'t find file {prefix} in dir {input_dir}')
        exit()
    return res


class MarkovTrain:
    def __init__(self, args):
        self.args = args
        self.verbose = args.verbose
        self.wgbs_tools = get_wgbs_tools_path(args.wgbstools)
        self.markers_path = args.markers_path
        self.groups_file = args.groups_file
        self.pats = []
        self.cmds = []
        self.ss_dir = None
        self.tt_dir = None
        self.deconv_src = None
        self.outdir = None

    def run(self):
        # validate arguments
        self.validate_args()
        # parse pat files and groups file:
        self.parse_pats()
        # Build directory and subdirectories:
        self.outdir = self.create_directory_tree()
        # Create markov counts (m*.gz files)
        self.create_counts()
        # Create transition matrices (Deconvolution code binary file)
        self.wrap_transition()

    def parse_pats(self):

        # validate reference data directory:
        refdir = self.args.reference_data
        if not op.isdir(refdir):
            eprint('Invalid reference data directory (of pat files):', refdir)
            exit()

        # test groups file is a legal csv, with the same files as in the reference_dir
        df = pd.read_csv(self.groups_file, usecols=[0, 1], comment='#').dropna().reset_index(drop=True)
        df.columns = ['name', 'group']
        self.pats = [prefix2filepath(refdir, name) for name in df['name']]

        # make sure pat.gz.csi files are present:
        for pat in self.pats:
            pat_csi(pat, self.wgbs_tools)
        nr_all_pats = len([f for f in os.listdir(refdir) if f.endswith('.pat.gz')])
        eprint('Processing {} pat files (out of {} files in {})...'.format(len(self.pats), nr_all_pats, refdir))

    def validate_args(self):
        # validate markers path:
        if not op.isfile(self.markers_path):
            eprint('Invalid markers path:', self.markers_path)
            exit()

        # validate groups file:
        if not op.isfile(self.groups_file):
            eprint('Invalid groups file (--groups_file):', self.groups_file)
            exit()

        # markers file must have at leaset 5 columns: 
        peek_df = pd.read_csv(self.markers_path, sep='\t', nrows=1, header=None)
        if len(peek_df.columns) < 5:
            eprint(f'Invalid markers path: {self.markers_path}. Should be tab separated with >= 5 columns')
            exit()
        if not (str(peek_df.iloc[0, 4]).isdigit() and str(peek_df.iloc[0, 3]).isdigit()):
            eprint('Invalid markers file. There should be no header, and columns 4-5 should be CpG-Index.')
            eprint('First line of markers file is:\n{}'.format(peek_df.iloc[0, :].values))
            exit()

        # validate counter tool path
        counter_tool = op.abspath(self.args.counter_tool)
        if not op.isfile(counter_tool):
            eprint('counter tool path is wrong:', counter_tool)
            exit()

        # validate deconvolution code src dir:
        dec_src = self.args.deconv_src
        self.gen_tran_tool = op.join(dec_src, 'generate_transition_matrix_MC_k_order.py')
        if not (op.isdir(dec_src) and op.isfile(self.gen_tran_tool)):
            eprint('Invalid deconvolution code directory (--deconv_src):', dec_src)
            exit()

    def gen_outdir(self):
        outdir = self.args.out_dir
        if outdir is not None:
            return outdir

        # If no output directory is given, generate a default one:
        name = op.basename(self.markers_path)
        if name.endswith('.gz'):
            name = name[:-3]
        name = op.splitext(name)[0]
        return name

    def create_directory_tree(self):
        outdir = self.gen_outdir()
        # If directory already exists, either delete it or abort
        if op.isdir(outdir):
            if self.args.force:
                subprocess.check_call('rm -rf ' + outdir, shell=True)
            else:
                eprint('directory already exists. Use -f to overwriteit. Aborting.')
                exit()
        # Else, just create it
        os.mkdir(outdir)

        # create subdirectories:
        self.ss_dir = op.join(outdir, 'train_suff_stats')
        os.mkdir(self.ss_dir)
        self.tt_dir = op.join(outdir, 'train_trans_params')
        os.mkdir(self.tt_dir)
        os.mkdir(op.join(outdir, 'results'))

        # copy markers file and groups file:
        subprocess.check_call('cp {} {}'.format(op.abspath(self.markers_path), op.join(outdir, 'markers.tsv')), shell=True)
        subprocess.check_call('cp {} {}'.format(op.abspath(self.groups_file), op.join(outdir, 'groups.csv')), shell=True)
        eprint('Output directory:', outdir)
        return outdir

    def create_counts(self):

        # run counter tool on all reference pat files:
        for pat in self.pats:
            prefix = op.join(self.ss_dir, op.basename(pat)[:-7])
            cmd = f'{self.wgbs_tools} cview -L {self.markers_path} {pat} --strict --strip --genome {self.args.genome}'
            cmd += f' | {self.args.counter_tool} - {prefix} -k 5 -b {self.markers_path} --pigz'
            self.cmds.append(cmd)
        eprint(f'Processing {len(self.cmds)} pat files...')
        self.perform_counter_commands()

        # make sure countes were created:
        nr_m_files = len([f for f in os.listdir(self.ss_dir) if '.m0' in f])
        if nr_m_files != len(self.pats):
            eprint('Failed creating counts files ({}/{} created). Aborting.'.format(nr_m_files, len(self.pats)))
            exit()
        # gzip m files
        for f in os.listdir(self.ss_dir):
            cmd = 'gzip ' + op.join(self.ss_dir, f)
            subprocess.check_call(cmd, shell=True)
        #subprocess.check_call('gzip {}/*'.format(self.ss_dir), shell=True)

    def perform_counter_commands(self):
        # Debug mode on: don't run commands, only print them
        if self.args.debug:
            for c in self.cmds:
                print(c)
            return

        # If a single thread is allowed, perform in a queue
        if self.args.threads ==  1:
            for cmd, pat in zip(self.cmds, self.pats):
                single_thread(cmd, pat, self.verbose)
        else:
            # else, multiprocess:
            with Pool(self.args.threads) as p:
                for cmd, pat in zip(self.cmds, self.pats):
                    p.apply_async(single_thread, (cmd, pat, self.verbose))
                p.close()
                p.join()

    def wrap_transition(self):
        for k in range(6):
            cmd = 'python3 {} -g 1 -k {} -otrans {} -in_tissue '.format(self.gen_tran_tool, k, self.tt_dir)
            cmd += '{} -i {}'.format(self.groups_file, self.ss_dir)
            if self.verbose:
                eprint('k = ', k)
                eprint(cmd)
            subprocess.check_call(cmd, shell=True)

        #subprocess.check_call('rm -rf ' + self.ss_dir, shell=True)
        return

def single_thread(cmd, pat, verbose):
    if verbose:
        eprint(op.basename(pat))
    subprocess.check_call(cmd, shell=True, stderr=None if verbose else subprocess.PIPE)

def main():
    args = parse_args()
    MarkovTrain(args).run()


if __name__ == '__main__':
    main()

