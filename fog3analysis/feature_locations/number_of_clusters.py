from __future__ import division
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

sys.path.insert(
    0, os.path.dirname(os.path.realpath(__file__)) + '../cliputil/')
import cliputil
from cliputil.sbGene import sbGene

def add_to_ga(ga, ga_to_add):
    if type(ga_to_add) == type([]):
        for _iv in ga_to_add:
            ga[_iv] += 1
        return
    if type(ga_to_add) == type(HTSeq.GenomicPosition('I', 1, '+')):
        ga[ga_to_add] += 1
        return
    for _iv, val in ga_to_add.steps():
        ga[_iv] += val

def write_bedgraph(bedgraph, filename):
    bedgraph.write_bedgraph_file(
        filename.partition('.wig')[0] + '_+.wig', strand='+')
    bedgraph.write_bedgraph_file(
        filename.partition('.wig')[0] + '_-.wig', strand='-')





def analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, coverage, _gtf,
        output_dirname='figs/'):
    # Analysis by applying CWT to the total coverage in a gene.
    print 'Analysis by applying CWT to the total coverage in a gene...'
    subpeak_list = []
    clust_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    exons_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    peak_pos_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    peak_maxima_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    #coverage_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    for n, gene_name in enumerate(gene_targets):
        if not( n % 100):
            cliputil.printProgress(
                n, len(gene_targets), prefix='Identifying subpeaks... ')
            print '',
        if gene_name not in _gtf.sbgenes:
            print "Missing " + gene_name
            continue
        sp = _gtf.sbgenes[gene_name]
#        if gene.strand == '-': continue
        #sp = sbGene(HTSeq.GenomicInterval(
        #    gene.chrom, gene.txpt_left, gene.txpt_right, gene.strand))
        sp.gene_name = gene_name
        sp.chr_lens = _gtf.chr_lens
        #sp.determine_coverage_accross_exons(coverage, txpts[sp.gene_name].exons_dict)
        #exons_dict = sp.rc_exons_dict_if_minus(
        #    sp.strand, txpts[sp.gene_name].exons_dict)
        if gene_name == 'mpk-1':
            print 'nofc'
            print sp.exons_dict
        sp.find_subpeaks(
            coverage, sp.exons_dict,
            rc_exons_dict_if_minus=sp.minus_strand_flipped)
        subpeak_list.append(sp)
        add_to_ga(clust_bedgraph, sp.clusters_as_ga)
        add_to_ga(exons_bedgraph, sp.exon_ivs)
        add_to_ga(peak_pos_bedgraph, sp.peak_pos_iv)
        add_to_ga(peak_maxima_bedgraph, sp.peak_maxima_iv)
    write_bedgraph(clust_bedgraph, 'pattern_bedgraphs/clusters.wig')
    write_bedgraph(exons_bedgraph, 'pattern_bedgraphs/exons.wig')
    write_bedgraph(peak_pos_bedgraph, 'pattern_bedgraphs/peak_pos.wig')
    write_bedgraph(peak_maxima_bedgraph, 'pattern_bedgraphs/peak_maxima.wig')
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
        for g in ranked[::-1]: f.write('%s\t%i\n' % (g, num_of_peaks[g]))
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
    #parser.add_argument(
    #    '-n', '--no_user_input', action='store_true', default=False,
    #    help='Do not request user input or run as a loop (Default: False).'
    #)
    parser.add_argument('-i', '--input',
                        default=None,
                        help='''Input peaks file (Required).''')
    parser.add_argument('-c', '--config_ini',
        help='config.ini file',
        default='auto.ini')
    #parser.add_argument('-l', '--load_txpts',
    #                    help='Load marked transcript information.',
    #                    default=False, action='store_true')
    args = parser.parse_args()
    args.no_user_input = False
    args.load_txpts = False
    return args


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
    ga = cliputil.get_a_bedgraph(plus_file, minus_file)
    _gtf = cliputil.load_gtf(lib=lib)
    peaks = _gtf.flip_minus_strand_peak_features(
        pandas.read_csv(args.input, sep='\t'))
    _gtf.add_peak_locations_to_transcripts(peaks)
    #(sequences, rc_sequences, chr_lens) = cliputil.get_sequences(lib)
    output_dirname = 'figs/{d}/'.format(
        d=os.path.basename(os.path.dirname(args.input)))
    print "Writing to {s}.".format(s=output_dirname)
    gene_targets = set(peaks['gene_name'].tolist())
    analyze_by_applying_cwt_to_all_reads_in_gene(
        gene_targets, ga, _gtf,
        output_dirname=output_dirname)
