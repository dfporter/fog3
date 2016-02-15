import pandas
import csv
import re
import numpy as np
import glob
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import sys
import matplotlib
import os
import argparse


def locate_in_gene(gtf_sep_cols, combined, cf_utr_elements=False):
    print combined.head()
    for index, peak_row in combined.iterrows():
        if not index % 100:
            print "Locating peak %i..." % index
        gene = str(combined.loc[index, 'gene_name'])
        rows = gtf_sep_cols[gtf_sep_cols['gene_name']==gene]
        combined.loc[index, 'biotype'] = rows['biotype'].tolist()[0] 
        if len(rows) == 0:
            combined.loc[index, 'location'] = "Unknown"
            continue
        rows = rows[rows['2']=='CDS']
        if len(rows) == 0:
            combined.loc[index, 'location'] = "ncRNA"
            continue
        gene_left = min(rows['3'])
        gene_right = max(rows['4'])
        gene_strand = rows['6']
        dist = get_distance((gene_left, gene_right),
                            (combined.loc[index, 'left'],
                             combined.loc[index, 'right']))
        if combined.loc[index, 'strand'] == '-':
            dist = -1 * dist
        combined.loc[index, 'dist_to_CDS'] = dist
        if not dist:
            combined.loc[index, 'location'] = "CDS"
        if dist < 0:
            combined.loc[index, 'location'] = '''5'UTR'''
        if dist > 0:
            combined.loc[index, 'location'] = '''3'UTR'''


def dist_between_elements(a, b):
    if a[0] > b[1]: return a[0] - b[1]
    elif a[1] < b[0]: return a[1] - b[0]
    else: return 0

      
def get_distance(geneiv, peakiv):
    if geneiv[1] < peakiv[0]:
        return peakiv[0] - geneiv[1]
    if peakiv[1] < geneiv[0]:
        return peakiv[1] - geneiv[0]
    if geneiv[0] < peakiv[0] < geneiv[1]:
        if peakiv[1] < geneiv[1]:
            return 0  # Entirely in CDS.
        else:
            in_cds = geneiv[1] - peakiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[1]
    if geneiv[0] < peakiv[1] < geneiv[1]:
        if peakiv[0] > geneiv[0]:
            return 0  # Entirely in CDS.
        else:
            in_cds = peakiv[1] - geneiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[0]
    return 0


def pie_chart_of_peak_locations(peaks, gtf_sep_cols,
                                label='no_label', subdir='',
                                discard_ncRNA=True,
                                already_located=False):
    if not already_located: locate_in_gene(gtf_sep_cols, peaks)
    locations = dict(peaks['location'].value_counts())
    if 'Unknown' in locations:
        del locations['Unknown']
    if discard_ncRNA and 'ncRNA' in locations: del locations['ncRNA']
    make_a_pie_chart(locations, label=label, subdir=subdir)


def make_a_pie_chart(locations, label='no_label',
                     subdir=''):
    total_locs = sum(locations.values())
    for pos in locations:
        locations[pos] = 100.0 * float(locations[pos])/float(total_locs)
    plt.clf()
    matplotlib.rcParams['lines.linewidth'] = 0
    _fig = plt.figure(1, figsize=(6,6))
    #ax = plt.axes([0.1,0.9])
    _labels = tuple(locations.keys())
    fracs = locations.values()
    # I have no idea how this color stuff works.
    # Revisit this to learn how matplotlib handles colors.
    num_values = len(_labels)
    color_vals = []
    for label_num in range(0, num_values):
        fraction = float(label_num)/float(num_values)
        color_vals.append(fraction)
    my_norm = matplotlib.colors.Normalize(-1, 1)
    my_cmap = matplotlib.cm.get_cmap('Set1')
    #plt.pie(fracs, labels=_labels, colors=my_cmap(my_norm(color_vals)),
    #        autopct='%1.1f%%')
    # Following trick to get rid of the lines is from:
    # http://stackoverflow.com/questions/1915871/matplotlib-controlling-pie-chart-font-color-line-width
    pieWedgesCollection = plt.pie(
        fracs,
        labels=_labels,
        colors=my_cmap(my_norm(color_vals)),
        autopct='%1.1f%%')[0]
    for wedge in pieWedgesCollection:
        wedge.set_lw(0)
    plt.axes().set_aspect('equal')
    #plt.show()
    if not os.path.exists('figs/'): os.system('mkdir figs/')
    if not os.path.exists('figs/%s/' % subdir): os.system('mkdir figs/%s' % subdir)
    print "\tSaving to %s" % str('figs/%s/pie_chart_%s.pdf' % (subdir, label))
    plt.savefig('figs/%s/pie_chart_%s.pdf' % (subdir, label))
    plt.clf()


def pie_chart_of_biotype(peaks, gtf_sep_cols,
                                label='biotype', subdir='',
                         already_located=False):
    if not already_located: locate_in_gene(gtf_sep_cols, peaks)
    locations = dict(peaks['biotype'].value_counts())
    if 'Unknown' in locations:
        del locations['Unknown']
    if len(locations.keys()) == 0:
        return
    make_a_pie_chart(locations, label=label+'_all_types', subdir=subdir)
    total_ncRNA = 0
    for _label in locations:
        if _label == 'protein_coding':
            continue
        else:
            total_ncRNA += locations[_label]
    simplified_locations = {'protein_coding': locations['protein_coding'],
                            'ncRNA': total_ncRNA}
    make_a_pie_chart(simplified_locations, label=label+'_coding_vs_nc', subdir=subdir)

def run(gtf_sep_cols, lib, peaks_dir, already_located=False):
    gtf_sep_cols = pandas.read_csv(lib['gtf'], sep='\t')
    top_level_dir = peaks_dir + '/'
    peaks = {}
    for filename in glob.glob(top_level_dir + '/*.txt'):
        print "\tReading peaks from {p}".format(p=filename)
        df = pandas.read_csv(filename, sep='\t')
        if 'location' in df.columns:
            peaks[filename] = df
    for exp in peaks:
        print '~~~'
        print exp
        print peaks[exp].head(1)
        pie_chart_of_peak_locations(peaks[exp],
                                    gtf_sep_cols,
                                    label='all_peaks_in_gene_%s' % os.path.basename(exp),
                                    subdir='figs/',
                                    already_located=already_located)
        pie_chart_of_biotype(
            peaks[exp], gtf_sep_cols,
            label='all_peaks_in_gene_bioytpe_%s' % os.path.basename(exp),
            subdir='figs/',
            already_located=True)
        peaks[exp].sort('height', inplace=True, ascending=False)
        peaks[exp].drop_duplicates('gene_name', inplace=True)
        pie_chart_of_peak_locations(
            peaks[exp], gtf_sep_cols,
            label='highest_peak_per_gene_%s' % os.path.basename(exp),
            subdir='figs/',
            already_located=True)
        pie_chart_of_biotype(
            peaks[exp], gtf_sep_cols,
            label='highest_peak_per_gene_bioytpe_%s' % os.path.basename(exp),
            subdir='figs/',
            already_located=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Compares CLIP-seq data to RIP.''')
    parser.add_argument(
        '-d', '--peaks_dir', default='all_peaks.txt',
        help='''A directory of peaks with counts of reads for each dataset.'''
    )
    parser.add_argument(
        '-o', '--output', default='annotated_peaks.txt',
        help='''Output filename.'''
    )
    parser.add_argument(
        '-c', '--config',
        help='''Directory holding config.py file.'''
    )
    args = parser.parse_args()
    sys.path.insert(0, args.config)
    import config
    lib = config.config()
    run(gtf_sep_cols, args.peaks_dir, lib)
