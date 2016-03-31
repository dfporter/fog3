import cluster_combine
import config

import pandas
import HTSeq
import argparse
import os
import sys
import numpy as np
import pandas
import collections
import glob
import re


def sum_reads_in_gene(full, lib=None):
    if lib is not None:
        beds = set([lib[x] for x in lib if re.match('exp.*', x)])
        exp_counts = [x.rstrip('.bed') + '_counts.txt' for x in beds]
        beds = set([lib[x] for x in lib if (
            re.match('control.*', x) or re.match('.*n2.*', x))])
        control_counts = [x.rstrip('.bed') + '_counts.txt' for x in beds]
    else:
        counts = [x for x in full.columns if re.match('.*_counts\.txt', x)]
        control_counts = []
        exp_counts = []
        for x in counts:
            print x
            if re.match('.*control.*', x) or re.match('.*n2.*', x):
                control_counts.append(x)
                print 'contr'
            else:
                exp_counts.append(x)
                print 'ex'
    for i, row in full.iterrows():
        full.loc[i, 'exp'] = sum([row[x] for x in exp_counts])
        full.loc[i, 'control'] = sum([row[x] for x in control_counts])
        full.loc[i, 'ratio'] = float(full.loc[i, 'exp'])/max([1., full.loc[i, 'control']])
    return full


def apply_cutoff(full, ratio_cutoff=10, min_reads_cutoff=200,
                 only_mrna=False):
    past = full[full['ratio']>=ratio_cutoff]
    past = past[past['exp']>min_reads_cutoff]
    if only_mrna:
        past = past[past['biotype']=='protein_coding']
    return past


def read_args():
    parser = argparse.ArgumentParser(description='''
        Convert cluster tables to peaks files by comparing replicates/p values.''')
    parser.add_argument('-c', '--config_ini',
                        default='config.ini',
                        help='''config.ini file.''')
    parser.add_argument('-i', '--peaks',
                        help='Peaks file.')
    parser.add_argument('-t', '--counts',
                        help='Counts per gene.')
    parser.add_argument('-r', '--ratio_cutoff',
                        help='Minimum RPKM enrichment in gene over \
negative control (default 10).',
                        default=10)
    parser.add_argument('-n', '--min_reads_cutoff',
                        help='Minimum RPKM in gene (selects for \
frequent targts, default 200 RPKM summed accross replicates).',
                        default=200)
    args = parser.parse_args()
    args.ratio_cutoff = float(args.ratio_cutoff)
    args.min_reads_cutoff = float(args.min_reads_cutoff)
    return args


def run(args, lib):
    if not os.path.exists('./tables/'): os.system('mkdir tables/')
    table = pandas.read_csv(args.peaks, sep='\t')
    table = cluster_combine.get_rpkm(
        table, lib, counts_fname=args.counts)
    table = sum_reads_in_gene(table, lib=None)
    table.to_csv('tmp', sep='\t')
    print 'ratio: min {z} max {zz} mean {aa}'.format(
        z=min(table['ratio'].tolist()),
        zz=max(table['ratio'].tolist()),
        aa=np.mean(table['ratio'].tolist())
    )
    print 'exp: min {z} max {zz} mean {xx}'.format(
        z=min(table['exp'].tolist()),
        zz=max(table['exp'].tolist()),
        xx=np.mean(table['exp'].tolist()),
    )
    table = apply_cutoff(table, ratio_cutoff=args.ratio_cutoff,
                         min_reads_cutoff=args.min_reads_cutoff)
    table.to_csv('tables/past_cutoff.txt', sep='\t')


if __name__ == '__main__':
    args = read_args()
    lib = config.config(args.config_ini)
    run(args, lib)

