#!/usr/bin/python3 -u

import os
import sys
from io import StringIO
import numpy as np
import os.path as op
import argparse
import subprocess
import pandas as pd
from multiprocessing import Pool
from train import get_wgbs_tools_path, eprint, pat_csi


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('working_dir', help='path to working directory, created by train.py script.')
    parser.add_argument('--pats', '-p', nargs='+', help='path to pat file[s] to deconvolve', required=True)
    parser.add_argument('--prior', type=float, help='Prior for target cell type, in range [0, 1]. default: 0.05', default=0.05)
    parser.add_argument('--target', help='Target name, from the groups file. ' \
                        'Default is the groups of the first line in this file')
    parser.add_argument('--genome', help='genome reference name. default: hg19', default='hg19')
    parser.add_argument('--deconv_src', help='path to source code directory with the deconv python files. Default: ./src', default='src')
    parser.add_argument('--wgbstools', '-wt', help='path to wgbstools, if it\'s not in PATH')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    return parser.parse_args()


def validate_dir(d):
    if not op.isdir(d):
        eprint('Invalid directory:', d)
        exit()


class MarkovDecon:
    def __init__(self, args, wgbs_tools='wgbstools'):
        self.args = args
        self.verbose = args.verbose
        self.wgbs_tools = get_wgbs_tools_path(self.args.wgbstools)
        self.pats = args.pats
        self.res_dir = None
        self.tt_dir = None
        self.dec_script = None
        self.results_script = None
        self.wd = self.args.working_dir
        self.gpath = op.join(self.wd, 'groups.csv')
        self.markers_path = op.join(self.wd, 'markers.tsv')
        self.target, self.priors = self.set_target_and_priors()

    def set_target_and_priors(self):
        groups = list(pd.read_csv(self.gpath, usecols=[1]).iloc[:, 0].unique())
        target = self.args.target
        if target is None:
            target = groups[0]
            if target in ('background', 'bg', 'plasma', 'blood', 'other'):
                target = groups[1]
        elif target not in groups:
            eprint(f'Invalid target {target}: Does not appear in group file')
            exit()
        eprint('Target:', target)

        p = self.args.prior
        if not 0 <= p <= 1:
            eprint('Invalid prior: Not in range [0, 1]')
            exit()
        if target == groups[0]:
            return target, [p,  1 -  p]
        elif target == groups[1]:
            return target, [1 - p, p]
        else:
            eprint('Invalid target or groups file: this script supports at most 2 groups')
            exit()

    def run(self):
        # validate arguments
        self.validate_args()
        self.deconvolve()

    def validate_args(self):
        # validate working directory:
        cwd = self.wd
        validate_dir(cwd)
        self.res_dir = op.join(cwd, 'results')
        validate_dir(self.res_dir)
        self.tt_dir = op.join(cwd, 'train_trans_params')
        validate_dir(self.tt_dir)

        # validate deconvolution code src dir:
        dec_src = self.args.deconv_src
        self.dec_script = op.join(dec_src, 'single_save_of_probabilistic_data.py')
        self.results_script = op.join(dec_src, 'create_deconv_summary.py')
        if not (op.isdir(dec_src) and op.isfile(self.dec_script) and op.isfile(self.results_script)):
            eprint('Invalid deconvolution code directory (--deconv_src):', dec_src)
            exit()

        # validate pat file[s]:
        for pat in self.pats:
            pat_csi(pat, self.wgbs_tools)

    def single_pat(self, inresdir, pat, sbout):

        tmp_labels_path = op.join(inresdir, 'tmp.{}'.format(np.random.randint(10000)))
        tmp_labels_path = op.basename(pat)[:-7]
        eprint(pat)
        eprint(tmp_labels_path)
        cmd = f'python3 {self.dec_script} -input {pat} -model 4 -len 3 -tissues 1'
        cmd += f' -output {inresdir} '
        cmd += ' -pai {} {} -g 1 '.format(*self.priors)
        cmd += f' -L {self.markers_path} -groups_f {self.gpath} '
        cmd += f' -trans {self.tt_dir} -d 0'
        cmd += f' --labels_path {tmp_labels_path}'
        cmd += f' --genome {self.args.genome}'
        if self.args.debug:
            print(cmd)
        else:
            subprocess.check_call(cmd, shell=True, stderr=sbout, stdout=sbout)

        # generate dpat file
        dpat_path = op.join(inresdir, op.basename(pat)[:-7] + '.pat')
        if self.args.debug:
            self.per_marker(dpat_path)
            return
        cmd = '{wt} cview --strip -L {markers} --strict --min_len 3 --genome {g} {pat}'.format(wt=self.wgbs_tools, markers=self.markers_path, g=self.args.genome, pat=pat)
        ptxt = subprocess.check_output(cmd, shell=True).decode()
        pdf = pd.read_csv(StringIO(ptxt), sep='\t', header=None)
        pdf['label'] = pd.read_csv(tmp_labels_path, sep='\t', header=None).values
        pdf.dropna(axis=1, inplace=True)
        pdf.to_csv(dpat_path, sep='\t', header=None, index=None)
        subprocess.check_call(f'{self.wgbs_tools} index {dpat_path}', shell=True)
        print(dpat_path + '.gz')

        if op.isfile(tmp_labels_path):
            os.remove(tmp_labels_path)
        self.per_marker(dpat_path)

    def per_marker(self, dpat_path):
        dpat_path += '.gz'
        mdf = pd.read_csv(self.markers_path, sep='\t', usecols=range(5), names=['chr', 'start', 'end', 'startCpG', 'endCpG'])
        mdf[self.target] = -1
        mdf['total'] = -1
        nr_empty_markers = 0
        for i, row in mdf.iterrows():
            if i > 0 and i % 20 == 0:
                print('finished {} markers'.format(i))
            cmd = '{wt} cview --strict --strip --min_len 3 -s {s}-{e} --genome {g} {dp}'.format(wt=self.wgbs_tools, s=row['startCpG'], e=row['endCpG'], g=self.args.genome, dp=dpat_path)
            txt = subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE).decode()
            if txt.strip() in ('', 'empty'):
                nr_empty_markers += 1
                target_count = 0
                total_count = 0
            else:
                df = pd.read_csv(StringIO(txt), sep='\t', header=None)
                nr_cols = len(list(df.columns))
                df = df.iloc[:, [3, nr_cols - 1]]
                total_count = df.iloc[:, 0].sum()
                target_count = df[df.iloc[:, 1] == self.target].iloc[:, 0].sum()
            #print('{} / {}'.format(target_count, total_count))
            mdf.loc[i, [self.target, 'total']] = target_count, total_count
        print(mdf)
        print(f'{nr_empty_markers} empty markers out of {mdf.shape[0]}')
        mdf_path = dpat_path.replace('.pat.gz', '.per_marker.tsv')
        mdf.to_csv(mdf_path, sep='\t', index=None)
        print(mdf_path)

    def deconvolve(self):
        sbout = None if self.verbose else subprocess.PIPE
        inresdir = op.join(self.res_dir, f'len_3_MC4_priorL_{self.priors[0]}')
        for pat in self.pats:
            self.single_pat(inresdir, pat, sbout)
        cmd = f'python3 {self.results_script} -model 4 -len 3 -tissues 1 -target {self.target}'
        cmd += f' -output {inresdir}/ -g 1 -groups_f {self.gpath}'
        if self.args.debug:
            print(cmd)
            return
        subprocess.check_call(cmd, shell=True, stderr=sbout, stdout=sbout)

        # remove binary files:
        for f in os.listdir(inresdir):
            if f.endswith('binary'):
                os.remove(op.join(inresdir, f))
        print(op.join(inresdir, 'deconv_summary.csv'))


def main():
    args = parse_args()
    MarkovDecon(args).run()


if __name__ == '__main__':
    main()

