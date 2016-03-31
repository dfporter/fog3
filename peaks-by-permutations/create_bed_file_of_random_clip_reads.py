__author__ = 'dp'
import HTSeq
import numpy as np
import bisect
import pandas
# from collections import Counter
import collections
import copy
import os
import sys
import argparse
import randomize_reads
import find_peaks_by_permutations
import p_values_of_clusters
import cluster_combine
import traceback
import assign_to_genes
import config

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_ini',
                        help="""Directory holding config.py""")
    args = parser.parse_args()
    lib = config.config(args.config_ini)
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    ga, reads_by_gene, counts_by_gene = find_peaks_by_permutations.get_bed_files(
        bed_folder=lib['read_beds'], gtf=gtf,
        args=args, lib=lib)
    exons_as_rows = p_values_of_clusters.get_exonic_ranges(lib['gtf'])
    while True:
        try:
            randomize_reads.run(
                reads_by_gene, counts_by_gene, exons_as_rows,
                lib=lib)
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
                reload(randomize_reads)
                reload(p_values_of_clusters)
                reload(cluster_combine)
                reload(assign_to_genes)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

