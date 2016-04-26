"""
Load txpt objects and read locations within txpts.

For each txpt: call peaks by read shuffling (in negative and positive) and by CWT.

Rank by p value for the shuffle.

Apply a BH FDR.

Compare replicates and output peaks.

INPUT:
Config file.

OUTPUT:
Outputs to lib['clusters_dir'] from the .ini config file.
Outputs called peaks to lib['permutation_peaks'].

Run as;
python find_peaks_by_permutations.py -c config.ini

From config.ini, this uses the following variables:
exp_bed1, ...
control_bed1, ...
gtf_raw
read_beds
gtf
clusters_dir
permutation_peaks_dir

"""

import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import numpy as np
import pandas
import collections
import re

import cluster_reads
import p_values_of_clusters
import cluster_combine
import filter_by_reads_per_gene
from get_bed_files import get_bed_files, read_a_bed_file
import config


def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument('-c', '--config_ini',
                        default='config/',
                        help='''config.ini file.''')
    parser.add_argument('-v', '--no_ui',
                        help='No user input.',
                        default=False, action='store_true')
    parser.add_argument('-t', '--filter_by_counts',
                        help='''
(Optional) Filter the output using read count values supplied by
 the file passed to this argument. Default: False''',
                        default='')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = read_args()
    lib = config.config(args.config_ini)
    if not os.path.exists('logs/'): os.system('mkdir logs')
    if not os.path.exists('figs/'): os.system('mkdir figs')
    logging.basicConfig(
        filename='logs/%s_find_peaks_by_permutations.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    # Load raw reads.
    print lib['read_beds']
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    ga, reads_by_gene, counts_by_gene = get_bed_files(
        bed_folder=lib['read_beds'], gtf=gtf,
        args=args, lib=lib)
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf'])
    if args.no_ui:
        cluster_reads.run(
            reads_by_gene, counts_by_gene, exons_as_rows,
            lib=lib, args=args)
        cluster_combine.run(args, lib, gtf, exons_as_rows)
        if args.filter_by_counts != '':
            args.counts = args.filter_by_counts
            args.peaks = 'permutation_peaks/combined_exp.txt'
            os.system('python annotate_peaks.py -c {a} -p permutation_peaks/combined_exp.txt'.format(
                a=args.config_ini))
            os.system('python annotate_peaks.py -c {a} -p permutation_peaks/combined_control.txt'.format(
                a=args.config_ini))
            filter_by_reads_per_gene.run(args, lib)
        print "Successfully finished."
        sys.exit()
    while True:
        try:
            cluster_reads.run(
                reads_by_gene, counts_by_gene, exons_as_rows,
                lib=lib, args=args)
            cluster_combine.run(args, lib, gtf, exons_as_rows)
            os.system('python annotate_peaks.py -c {a} -p permutation_peaks/combined_exp.txt'.format(
                a=args.config_ini))
            os.system('python annotate_peaks.py -c {a} -p permutation_peaks/combined_control.txt'.format(
                a=args.config_ini))
            if args.filter_by_counts != '':
                import filter_by_reads_per_gene
                args.counts = args.filter_by_counts
                args.peaks = 'permutation_peaks/combined_exp.txt'
                filter_by_reads_per_gene.run(args, lib)
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(config)
                reload(cluster_reads)
                reload(p_values_of_clusters)
                reload(cluster_combine)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

