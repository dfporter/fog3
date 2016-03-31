"""
Load the clusters output by find_peaks_by_permutations.

"""
import HTSeq
import argparse
import os
import traceback
import sys
import pandas
import glob

import cluster_reads
import p_values_of_clusters
import cluster_combine
import config


def read_args():
    parser = argparse.ArgumentParser(description='''
        Convert cluster tables to peaks files by comparing replicates/p values.''')
    parser.add_argument('-c', '--config_ini',
                        default='config.ini',
                        help='''config.ini file.''')
    parser.add_argument('-v', '--no_ui',
                        help='No user input.',
                        default=False, action='store_true')
    parser.add_argument('-s', '--max_padj',
                    help='Maximum adjusted p value to accept (default 0.01).',
                    default=0.01,)
    parser.add_argument('-r', '--min_num_reps',
            help='Miniumum number of replicates with a cluster (default 2).',
            default=2,)
    args = parser.parse_args()
    args.max_padj = float(args.max_padj)
    args.min_num_reps = int(args.min_num_reps)
    return args


if __name__ == '__main__':
    args = read_args()
    lib = config.config(args.config_ini)
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    # get_txpts() uses the peaks argument to select what txpts to load by a call to
    #  set(peaks['gene_name'].tolist()). That is, peaks is a dataframe.
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf'])
    if args.no_ui:
        cluster_combine.run(
            args, lib, gtf, exons_as_rows, max_padj=args.max_padj,
            min_num_reps=args.min_num_reps)
        sys.exit()
    while True:
        try:
            cluster_combine.run(
                args, lib, gtf, exons_as_rows, max_padj=args.max_padj,
                min_num_reps=args.min_num_reps)
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
                reload(cluster_combine)
                reload(cluster_reads)
                reload(p_values_of_clusters)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

