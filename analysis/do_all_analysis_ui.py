import re
import glob
import sys
import argparse
sys.path.insert(0,
    '/groups/Kimble/Common/fog_iCLIP/calls/analysis/src/')
import traceback
# import rip_comparison
import pandas
import os
import HTSeq

import clean_up_peaks_file
import get_args
import cf_fog1_rip
import cf_fog3_rip
import cf_sperm_oocyte
import peak_locations
import replicate_venns
import write_fastas
sys.path.insert(0,
    '/groups/Kimble/Common/fog_iCLIP/calls/peaks-by-permutations/')
import find_peaks_by_permutations
import annotate_peaks
del sys.path[0]

def read_collapsed_bed(bed_file):
    with open(bed_file, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            num_reads = int(re.match('.*n=(\d+).*', li).groups()[0])
            yield (
                HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5]),
                num_reads)


def read_uncollapsed_bed(bed_file):
    with open(bed_file, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            yield (
                HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5]), 1)


def create_bedgraphs_from_bed_files(
        overwrite=False,
        bed_folder='beds/',
        bedgraph_folder='bedgraphs/', lib=None):
    have_all = True
    for fname in ['controls_+.wig', 'controls_-.wig',
                  'exp_+.wig', 'exp_-.wig']:
        fname = bedgraph_folder + '/' + fname
        if not os.path.exists(fname): have_all = False
    if have_all and not overwrite:
        return
    print "Creating bedgraphs from bed files..."
    control_beds, exp_beds = get_bed_filenames(bed_folder, lib)
    ga = {}
    ga['control'] = HTSeq.GenomicArray('auto', stranded=True)
    for bed in control_beds:
        ga[bed] = HTSeq.GenomicArray('auto', stranded=True)
        for almnt, num_reads in read_uncollapsed_bed(bed):
            ga['control'][almnt] += num_reads
            ga[bed][almnt] += num_reads
        write_bedgraphs(ga[bed], bedgraph_folder, bed)
    ga['exp'] = HTSeq.GenomicArray('auto', stranded=True)
    for bed in exp_beds:
        ga[bed] = HTSeq.GenomicArray('auto', stranded=True)
        for almnt, num_reads in read_uncollapsed_bed(bed):
            ga['exp'][almnt] += num_reads
            ga[bed][almnt] += num_reads
        write_bedgraphs(ga[bed], bedgraph_folder, bed)
    write_bedgraphs(ga['control'], bedgraph_folder, 'controls')
    write_bedgraphs(ga['exp'], bedgraph_folder, 'exp')
    #
    # ga['control'].write_bedgraph_file(
    #     bedgraph_folder + '/controls_+.wig', '+')
    # ga['control'].write_bedgraph_file(
    #     bedgraph_folder + '/controls_-.wig', '-')
    # ga['exp'].write_bedgraph_file(
    #     bedgraph_folder + '/exp_+.wig', '+')
    # ga['exp'].write_bedgraph_file(
    #     bedgraph_folder + '/exp_-.wig', '-')


def write_bedgraphs(ga, bedgraph_folder, bed_name):
    ga.write_bedgraph_file(
        "{d}/{n}_+.wig".format(
            d=bedgraph_folder, n=os.path.basename(bed_name).partition('.bed')[0]),
        '+')
    ga.write_bedgraph_file(
        "{d}/{n}_-.wig".format(
            d=bedgraph_folder, n=os.path.basename(bed_name).partition('.bed')[0]),
        '-')


def get_bed_filenames(bed_folder, lib):
    control_beds = []
    exp_beds = []
    name = ('control_bed', 1)
    print lib
    while ''.join([str(x) for x in name]) in lib:
        control_beds.append(bed_folder + '/' + lib[''.join([str(x) for x in name])])
        name = (name[0], name[1] + 1)
    name = ('exp_bed', 1)
    while ''.join([str(x) for x in name]) in lib:
        exp_beds.append(bed_folder + '/' + lib[''.join([str(x) for x in name])])
        name = (name[0], name[1] + 1)
    print control_beds
    print exp_beds
    return control_beds, exp_beds


def add_reads_in_peak(args, lib):
    # sys.path.insert(0, '/groups/Kimble/Common/fog_iCLIP/calls/src/')
    # import peak_calling_tools
    # import callpeaks
    # sys.path.insert(0, '/groups/Kimble/Aman Prasad/clip/analysis/src/feature_locations/')
    # sys.path = sys.path[1:]
    peak_table = pandas.read_csv(args.peaks, sep='\t')
    # peak_table['gene_name'] = peak_table['gene_id']
    labels = [x.partition('.bed')[0] for x in [
        lib['control_bed1'], lib['control_bed2'],
        lib['control_bed3'], lib['control_bed4'],
        lib['exp_bed1'], lib['exp_bed2'],
        lib['exp_bed3'], lib['exp_bed4']
        ]
    ]
    gtf = HTSeq.GFF_Reader(lib['gtf_raw'], end_included=True)
    ga, reads_by_gene, counts_by_gene = find_peaks_by_permutations.get_bed_files(
        bed_folder=lib['read_beds'], gtf=gtf,
        args=args, lib=lib)
    ga_rename = {}
    for k in ga:
        ga_rename[os.path.basename(k).partition('.bed')[0]] = ga[k]
    for index, row in peak_table.iterrows():
        peak_iv = HTSeq.GenomicInterval(
            row['chrm'], row['left'], row['right'], row['strand'])
        for label in labels:
            read_ids = set()
            max_cov = max([x[1] for x in ga_rename[label][peak_iv].steps()])
            peak_table.loc[index, label] = max_cov
    # peak_table['num_reps_with_peak'] = peak_table['max_coverage']
    # del peak_table['max_coverage']
    # output_filename=lib['permutation_peaks_dir'] + os.path.basename(args.peaks)
    peak_table.to_csv(args.peaks, sep='\t', index=False)
    # ga_raw = peak_calling_tools.load_bed_folder(reads_dir)
    # callpeaks.fill_in_peaks(ga_raw, peak_table, output_filename=output_filename)


def hist_of_peaks_per_gene(peaks_filename):
    peaks = pandas.read_csv(peaks_filename, sep='\t')
    import collections
    num_of_peaks = collections.defaultdict(int)
    for gene in set(peaks['gene_name'].tolist()):
        num_of_peaks[gene] = len(peaks[peaks['gene_name']==gene])
    ranked = sorted(num_of_peaks.keys(), key=lambda x: num_of_peaks[x])
    if not os.path.exists('tables'): os.system('mkdir tables/')
    tablename = 'tables/{b}/'.format(b=os.path.dirname(peaks_filename))
    if not os.path.exists(tablename):
        os.system('mkdir ' + tablename)
    with open(tablename + 'peaks_per_gene_{a}'.format(
            a=os.path.basename(peaks_filename)), 'w') as f:
        li = 'gene_name\tnumber_of_peaks\n'
        for gene in ranked[::-1]:
            li += '{a}\t{b}\n'.format(a=gene, b=num_of_peaks[gene])
        f.write(li)
    with open('tables/peaks_per_gene.txt', 'w') as f:
        for g in ranked: f.write('%s\t%i\n' % (g, num_of_peaks[g]))
    import matplotlib.pyplot as plt
    import numpy as np
    plt.clf()
    fig, ax1 = plt.subplots(1, 1)
    bins_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    num, bins, p =ax1.hist(num_of_peaks.values(), bins=bins_list, color='black', )
    print num, bins
    print "total: {a}. 1 peak {b}. >1 peak {v}. % with more than one peak {zz}".format(
        a=sum(num), b=num[0], v=sum(num[1:]), zz=float(sum(num[1:]))/float(sum(num))
    )
    # ax1.set_xlim(1, 12)
    ax1.set_ylabel('Number of genes')
    ax1.set_xlabel('Number of peaks')
    ax1.set_xlim([min(bins_list), max(bins_list)])
    ax1.set_xticks([x + 0.5 for x in bins_list])
    ax1.set_xticklabels([str(x) for x in bins_list])
    # ax2.set_xlim(5, 12)
    # ax2.set_ylim(0, 50)
    # ax2.set_ylabel('Number of genes')
    # ax2.set_xlabel('Number of peaks')
    # ax2.hist(num_of_peaks.values(), color='black')
    plt.savefig('figs/peaks_per_gene_hist.pdf', format='pdf')
    plt.clf()
    plt.close()

def report_location(peaks_fname):
    peaks = pandas.read_csv(peaks_fname, sep='\t')
    cds = peaks[peaks['location']=='CDS']
    cds = len(cds['gene_id'].tolist())  # Allow multiple peaks per gene.
    fiveprime = peaks[peaks['location']=="5'UTR"]
    fiveprime = len(fiveprime['gene_id'].tolist())
    threeprime = peaks[peaks['location']=="3'UTR"]
    threeprime = len(threeprime['gene_id'].tolist())
    ncrna = peaks[peaks['location']=="ncRNA"]
    ncrna = len(ncrna['gene_id'].tolist())
    sum_of_mrna_locs = cds + fiveprime + threeprime
    gene_to_location = zip(
        [x for x in peaks['gene_id'].tolist()],
        [x for x in peaks['location'].tolist()])
    mrna_location = [
        x for x in gene_to_location if x[1] in ["5'UTR", "3'UTR", "CDS"]]
    non_mrna_location = [
        x for x in gene_to_location if x[1] not in ["5'UTR", "3'UTR", "CDS"]]
    print """
Location of peaks:
Peaks: {n} by len(.index), {nn} by len(genes)
Locations:
\tCDS: {a} ({aa}%) 5'UTR: {b} ({bb}%) 3'UTR: {c} ({cc}%)
\tSum of the above: {s}
\tncRNA {d} ({dd}%) mRNA {q} ({qq}%)
\tSum of the above: {allsum}
""".format(
    n=len(peaks.index), nn=len(peaks['gene_id'].tolist()),
    a=cds, b=fiveprime, c=threeprime, d=ncrna, q=sum_of_mrna_locs,
    aa="%.2f" % float(100. * float(cds)/float(sum_of_mrna_locs)),
    bb="%.2f" % float(100. * float(fiveprime)/float(sum_of_mrna_locs)),
    cc="%.2f" % float(100. * float(threeprime)/float(sum_of_mrna_locs)),
    s=cds + fiveprime + threeprime,
    dd="%.2f" % float(100. * float(ncrna)/float(ncrna + sum_of_mrna_locs)),
    qq="%.2f" % float(
        100. * float(sum_of_mrna_locs)/float(ncrna + sum_of_mrna_locs)),
    allsum=ncrna + sum_of_mrna_locs,
    )

if __name__ == '__main__':
    args = get_args.get_args()
    sys.path.insert(0, args.config)
    import config
    print 'config.__file__='
    print config.__file__
    lib = config.config()
    del sys.path[0]
    if not os.path.exists('figs/'): os.system('mkdir figs/')
    print lib['spermatogenic_program']
    reads_dir = lib['read_beds']
    peaks = pandas.read_csv(args.peaks, sep='\t')
    if 'location' not in peaks.columns:
        annotate_peaks.run(lib, args.peaks)
    if not os.path.exists('counts/'):
        assign_to_gene.run(use_this_lib=lib)
        os.system('cp combined_counts.txt counts/')
    if 'gene_id' in peaks.columns:
        target_wb_ids = set(peaks['gene_id'].tolist())
        if 'gene' not in peaks.columns:
            peaks['gene'] = peaks['gene_id']
    elif 'gene' in peaks.columns:  target_wb_ids = set(peaks['gene'].tolist())
    else: raise Exception
    hist_of_peaks_per_gene(args.peaks)
    print report_location(args.peaks)
    # peaks_l = peaks.to_dict('records')
    gtf_sep_cols = pandas.read_csv(lib['gtf'], sep='\t')
    # add_reads_in_peak(args, lib)
    while True:
        try:
            replicate_venns.run(
                lib['clusters_dir'], lib['permutation_peaks_dir'], lib=lib,
                only_use_these_targets=target_wb_ids
            )
            cf_fog1_rip.run(args.peaks, gtf_sep_cols, lib=lib)
            cf_fog3_rip.run(args.peaks, gtf_sep_cols, lib=lib)
            cf_sperm_oocyte.run(lib=lib, peaks_fname=args.peaks,
                                deseq_fname='/groups/Kimble/Common/fog_iCLIP/pre_calls/fog_clip_deseq.txt')
            peak_locations.run(gtf_sep_cols,
                               lib, os.path.dirname(args.peaks),
                               already_located=True)
            write_fastas.write_fastas(args.peaks)
            # create_bedgraphs_from_bed_files(
            #    overwrite=True,
            #    bed_folder=args.bed,
            #    bedgraph_folder=lib['bedgraphs_folder'],
            #    lib=lib)
            #
            print "Successfully finished."
        except:
            print traceback.format_exc()
            print "Failed."
        print "Hit enter to reload, CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(clean_up_peaks_file)
                reload(write_fastas)
                reload(get_args)
                reload(peak_locations)
                reload(replicate_venns)
                reload(cf_fog1_rip)
                reload(cf_fog3_rip)
                reload(cf_sperm_oocyte)
                reload(annotate_peaks)
                print "Successfully recompiled."
                reloaded = True
            except:
                print "Crash when trying to compile."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()
