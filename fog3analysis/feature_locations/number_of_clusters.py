import collections
import numpy as np
import HTSeq
import logging
import datetime
import argparse
import os
import time
import traceback
import sys
import pandas
import matplotlib.pyplot as plt

import config
import subpeaks
import determine_feature_locations



def analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, coverage, txpts, chr_lens,
        output_dirname='figs/'):
    # Analysis by applying CWT to the total coverage in a gene.
    subpeak_list = []
    for n, gene_name in enumerate(gene_targets):
        #if n > 100: break
        if gene_name not in txpts:
            print "Missing " + gene_name
            continue
        gene = txpts[gene_name]
#        if gene.strand == '-': continue
        sp = sbGene(HTSeq.GenomicInterval(
            gene.chrom, gene.txpt_left, gene.txpt_right, gene.strand))
        sp.gene_name = gene_name
#        sp.determine_coverage_accross_exons(coverage, txpts[gene_name].exons_dict)
        rc_corrected_exons_dict = {}
        if gene.strand == '-':
            num_exons = len(txpts[gene_name].exons_dict)
            for exon_num in txpts[gene_name].exons_dict:
                (chrm, left, right, strand) = txpts[gene_name].exons_dict[exon_num]
                _left = chr_lens[chrm] - right
                _right = chr_lens[chrm] - left
                rc_corrected_exons_dict[num_exons - exon_num + 1] = (chrm, _left, _right, '-')
                sp.find_subpeaks(coverage, rc_corrected_exons_dict)
        else:
            sp.find_subpeaks(coverage, txpts[gene_name].exons_dict)
        subpeak_list.append(sp)


    num_of_peaks = collections.defaultdict(int)
    for sp in subpeak_list:
        num_of_peaks[sp.gene_name] = sp.n_peaks
    ranked = sorted(num_of_peaks.keys(), key=lambda x: num_of_peaks[x])
    peaks_filename = 'number_of_clusters'
    with open('tables/peaks_per_gene_{a}'.format(
            a=os.path.basename(peaks_filename)), 'w') as f:
        li = 'gene_name\tnumber_of_peaks\n'
        for gene in ranked[::-1]:
            li += '{a}\t{b}\n'.format(a=gene, b=num_of_peaks[gene])
        f.write(li)
    with open('tables/peaks_per_gene.txt', 'w') as f:
        for g in ranked: f.write('%s\t%i\n' % (g, num_of_peaks[g]))
    plt.clf()
    fig, ax1 = plt.subplots(1, 1)
    bins_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    num, bins, p =ax1.hist(num_of_peaks.values(), bins=bins_list, color='black', )
    print(num, bins)
    print("total: {a}. 1 peak {b}. >1 peak {v}. % with more than one peak {zz}".format(
        a=sum(num), b=num[0], v=sum(num[1:]), zz=float(sum(num[1:]))/float(sum(num))
    ))
    # ax1.set_xlim(1, 12)
    ax1.set_ylabel('Number of genes')
    ax1.set_xlabel('Number of peaks')
    ax1.set_xlim([min(bins_list), max(bins_list)])
    ax1.set_xticks([x + 0.5 for x in bins_list])
    ax1.set_xticklabels([str(x) for x in bins_list])
    plt.savefig('figs/peaks_per_gene_hist.pdf', format='pdf')
    plt.clf()
    plt.close()

def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument(
        '-n', '--no_user_input', action='store_true', default=False,
        help='Do not request user input or run as a loop (Default: False).'
    )
    parser.add_argument('-i', '--input',
                        default=None,
                        help='''Input peaks file (Required).''')
    parser.add_argument('-c', '--config_ini',
        help='config.ini file',
        default='/groups/Kimble/Aman Prasad/clip/analysis/src/config.ini')
    parser.add_argument('-l', '--load_txpts',
                        help='Load marked transcript information.',
                        default=False, action='store_true')
    parser.add_argument('-g', '--load_gtf',
                        help='Load objects created by parsing a gtf file. Implies load_txpts is False.',
                        default=False, action='store_true')

    args = parser.parse_args()
    return args

def process_file(filename, coverage, txpts, sequences, chr_lens, args=None,
                 lib=None):
    print 'process_file():'
    output_dirname = 'figs/{d}/'.format(d=os.path.basename(os.path.dirname(filename)))
    print "Writing to {s}.".format(s=output_dirname)
    gene_targets = set(
        pandas.read_csv(filename, sep='\t')['gene_name'].tolist())
    gene_targets = filter(lambda x: isinstance(x, str), gene_targets)
    analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, coverage, txpts,  chr_lens, 
        output_dirname=output_dirname)

if __name__ == '__main__':
    args = read_args()
    lib = config.config(filepath=args.config_ini)
    # Logging.
    if not os.path.exists('logs'): os.system('mkdir logs')
    logging.basicConfig(
        filename='logs/%s_subpeaks.log' % datetime.datetime.now().strftime('%Hh%Mm'),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    logger = logging.getLogger(__name__)
    if not os.path.isfile(args.input):
        raise IOError('Input peaks file (-i) was not found.')
#    input_dir = os.path.dirname(args.input) + '/' #'../clip/pypeaks_fdr1_negip_local/'
    print "Loading bedgraphs..."
    plus_file = lib['bedgraph_exp_plus']
    minus_file = lib['bedgraph_exp_minus']
    ga = subpeaks.get_a_bedgraph(plus_file, minus_file)
    (peaks, txpts, chr_lens) = determine_feature_locations.get_data(
        args, lib=lib)
    (sequences, rc_sequences, chr_lens) = subpeaks.get_sequences(lib)
    process_file(args.input, ga, txpts, sequences,
                 chr_lens, args=args, lib=lib)
