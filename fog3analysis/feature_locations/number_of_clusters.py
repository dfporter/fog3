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
import locale
locale.setlocale(locale.LC_ALL, 'en_US')

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
    #peak_pos_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    peak_maxima_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    #coverage_bedgraph = HTSeq.GenomicArray('auto', stranded=True)
    for n, gene_name in enumerate(gene_targets):
        if not(n % 100):
            cliputil.printProgress(
                n, len(gene_targets), prefix='Identifying subpeaks... ')
            print ''
        if gene_name not in _gtf.sbgenes:
            print "Missing " + gene_name
            continue
        sp = _gtf.sbgenes[gene_name]
        sp.gene_name = gene_name
        sp.chr_lens = _gtf.chr_lens
        sp.find_subpeaks(
            coverage, sp.exons_dict,
            rc_exons_dict_if_minus=sp.minus_strand_flipped)
        subpeak_list.append(sp)
        add_to_ga(clust_bedgraph, sp.clusters_as_ga)
        add_to_ga(exons_bedgraph, sp.exon_ivs)
        #add_to_ga(peak_pos_bedgraph, sp.peak_pos_iv)
        add_to_ga(peak_maxima_bedgraph, sp.peak_maxima_ga)
    #write_bedgraph(clust_bedgraph, 'pattern_bedgraphs/clusters.wig')
    write_bedgraph(exons_bedgraph, 'pattern_bedgraphs/exons.wig')
    #write_bedgraph(peak_pos_bedgraph, 'pattern_bedgraphs/peak_pos.wig')
    write_bedgraph(peak_maxima_bedgraph, 'pattern_bedgraphs/peak_maxima.wig')
    num_of_peaks = {'clust': collections.defaultdict(int),
                    'maxima': collections.defaultdict(int)
                    }
    for sp in subpeak_list:
        num_of_peaks['clust'][sp.gene_name] = sp.n_peak_clusters
        num_of_peaks['maxima'][sp.gene_name] = sp.n_peak_maxima
    def write_peaks_per_gene_file(_dict, outfile):        
        ranked = sorted(_dict.keys(), key=lambda x: _dict[x])
        with open(outfile, 'w') as f:
            li = 'gene_name\tnumber_of_peaks\n'
            for gene in ranked[::-1]:
                li += '{a}\t{b}\n'.format(a=gene, b=_dict[gene])
            f.write(li)
    # Peaks per gene: csv and xls.
    write_peaks_per_gene_file(
        num_of_peaks['clust'], 'tables/peaks_per_gene_clusters.txt')
    write_peaks_per_gene_file(
        num_of_peaks['maxima'], 'tables/peaks_per_gene_maxima.txt')
    writer = pandas.ExcelWriter('tables/peaks_per_gene_histogram.xls',
                                engine='xlsxwriter')
    hist_df = dict_values_as_hist_df(num_of_peaks['maxima'])
    hist_df.to_excel(
        writer, sheet_name='Peak maxima')
    hist_df = dict_values_as_hist_df(num_of_peaks['clust'])
    hist_df.to_excel(
        writer, sheet_name='Peak clusters')
    writer.save()
    # Histograms.
    plot_hist(num_of_peaks['clust'], 'figs/peaks_per_gene_clusters.pdf')
    plot_hist(num_of_peaks['maxima'], 'figs/peaks_per_gene_maxima.pdf')


def dict_values_as_hist_df(_dict):
    hist = collections.defaultdict(int)
    for v in _dict.values(): hist[v] += 1
    hist = pandas.DataFrame(
        hist.items(), columns=['Number of FOG-3 peaks in a gene',
                               'Number of FOG-3 targets']) 
    return hist


def plot_hist(num_of_peaks, outfile):
    plt.clf()
    fig, ax1 = plt.subplots(1, 1)
    bins_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    num, bins, p =ax1.hist(
        num_of_peaks.values(), bins=bins_list, color='black', )
    print(num, bins)
    print("total: {a}. 1 peak {b}. >1 peak {v}. % with more than one peak {zz}".format(
        a=sum(num), b=num[0], v=sum(num[1:]), zz=float(sum(num[1:]))/float(sum(num))
    ))
    # ax1.set_xlim(1, 12)
    ax1.set_ylabel('Number of genes')
    ax1.set_xlabel('Number of peaks')
    ax1.set_xlim([min(bins_list), max(bins_list) + 1])
    ax1.set_xticks([x + 0.5 for x in bins_list])
    ax1.set_xticklabels([str(x) for x in bins_list])
    plt.savefig(outfile, format='pdf')
    plt.clf()
    plt.close()
    x = [x for x in num_of_peaks.values() if type(x) in [type(0), type(0.)]]
    print outfile
    print np.mean(x)
    print np.median(x)


def read_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument('-i', '--input',
                        default=None,
                        help='''Input peaks file (Required).''')
    parser.add_argument('-c', '--config_ini',
        help='config.ini file', default='auto.ini')
    parser.add_argument('-n', '--negative_control',
                        default='',
        help='(Optional) Negative control peaks file.')
    args = parser.parse_args()
    args.no_user_input = False
    args.load_txpts = False
    return args


if __name__ == '__main__':
    args = read_args()
    lib = config.config(filepath=args.config_ini)
    # Logging.
    #if not os.path.exists('logs'): os.system('mkdir logs')
    #logging.basicConfig(
    #    filename='logs/%s_subpeaks.log' % datetime.datetime.now().strftime('%Hh%Mm'),
    #    level=logging.DEBUG)
    #logging.info('Module called %s' % str(time.localtime()))
    #logger = logging.getLogger(__name__)
    if not os.path.isfile(args.input):
        raise IOError(
            'Input peaks file (-i) {0} was not found.'.format(args.input))
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
