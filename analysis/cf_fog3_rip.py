import re
import glob
import sys
import os
import pandas
import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scipy as sp
import scipy.stats


def load_file_into_dict(fname, use_header_val=None):
    _table = {}
    if use_header_val is not None:
        key = use_header_val
        df = pandas.read_csv(fname, sep='\t')
        print "header value: {h}. Read {n} peaks into {g} genes.".format(
            h=use_header_val, n=len(df.index), g=len(df[key].value_counts()))
        dfd = df.to_dict('records')
        for row in dfd:
            _table[row[key]] = row
        return _table
    if len(_table) == 0:
        print "%s was zero-length" % fname
    return _table


def build_array_of_ranks_in_each_dataset(rip, clip):
    all_genes = set(rip.keys()) | set(clip.keys())
    common_genes = set(rip.keys()) & set(clip.keys())
    array = []
    try_keys = ['no_antibody_ratio', 'n2_ratio', 'height', 'log2FoldChange']
    for gene in all_genes:
        if gene in common_genes:
            array.append([rip[gene]['Rank']] + [clip[gene][k] for k in try_keys if k in clip[gene]])
    return array, ['SAM rank', 'CLIP no-antibody ratio',
                   'CLIP N2 ratio', 'CLIP height',
                   'CLIP log2FoldChange vs no-antibody']


def correlations(array, header):
    if len(array) == 0:
        print '0 len array'
        return
    done = set()
    for i in range(0, 1):
        for j in range(0, len(array[0])):
            if i == j: continue
            if set([i, j]) in done: continue
            done |= set([i, j])
            print header[i], header[j]
            r = scipy.stats.spearmanr(
                [float(x[i]) for x in array],
                [float(x[j]) for x in array])
            print r


def compare_low_sam_ranks(rip, clip, header):
    print 'compare_low...' + '@' * 7
    quarts = separate_quartiles(rip, clip)
    print quarts.keys()
    for q_num in set([x[0] for x in quarts.keys()]):
        a = str(quarts[(q_num, 'rip')])
        b = str(quarts[(q_num, 'clip')])
        print "*" * 10
        print 'Quartile: %i' % q_num
        print a[0:100]
        print b[0:100]
        ar, header = build_array_of_ranks_in_each_dataset(quarts[(q_num, 'rip')], quarts[(q_num, 'clip')])
        for n, row in enumerate(ar):
            ar[n] = [float(x) for x in ar[n]]
        print q_num
        correlations(ar, header)


def separate_quartiles(rip, clip):
    all_genes = set(rip.keys()) | set(clip.keys())
    common_genes = set(rip.keys()) & set(clip.keys())
    s_rip = list(common_genes)

    q1_cutoff = (len(s_rip) * 0.33, 0)
    q2_cutoff = (len(s_rip) * 0.67, len(s_rip) * 0.33)
    q3_cutoff = (len(s_rip), len(s_rip) * 0.67)
    #q4_cutoff = (len(rip), len(rip) * 0.75)
    quarts = {}
    for gene in common_genes:
        for n, q in enumerate([q1_cutoff, q2_cutoff, q3_cutoff]):
            #print n
            #print 'cf %f < %f < %f' % (q[0],float(rip[gene]['Rank']), q[1])
            if (q[0] >= float(rip[gene]['Rank'])) and (
                float(rip[gene]['Rank']) >= q[1]
            ):
                label = (n, 'clip')
                quarts.setdefault(label, {})
                quarts[label][gene] = clip[gene]
                label = (n, 'rip')
                quarts.setdefault(label, {})
                quarts[label][gene] = rip[gene]
    # Check.
    for n, q in enumerate([q1_cutoff, q2_cutoff, q3_cutoff]):#, q4_cutoff]):
        label = (n, 'rip')
        print "=" * 10
        print label
        print min([quarts[label][x]['Rank'] for x in quarts[label]])
        print max([quarts[label][x]['Rank'] for x in quarts[label]])
        print 'len'
        print len(quarts[label])
        print "*"
        #print gene
    return quarts


def compare_overlap_as_function_of_rank(rip, clip):
    bins, in_two = build_bins(rip, clip)
    n_groups = len(bins)
    overlaps = [len(list(bins[x]['shared'])) for x in bins]
    non_overlaps = [len(list(bins[x]['unique'])) for x in bins]
    overlaps = [float(overlaps[n])/float(non_overlaps[n] + overlaps[n]) for n in range(len(overlaps))]
    #fig, (ax1, ax2) = plt.subplots(2, 1)
    fig = plt.figure(figsize=(8,6))
    #gs = gridspec.GridSpec(1, 2, width_ratios=[3,1])
    #ax1 = plt.subplot(gs[0])
    index = np.arange(n_groups)
    bar_width = 0.7
    opacity = 0.9
    plt.ylim([0, 0.7])
    rects1 = plt.bar(index, overlaps,
                     bar_width, alpha=opacity,
                     linewidth=0,
                     color='k', label='Overlapping with iCLIP targets')
#    rects2 = plt.bar(index + bar_width, non_overlaps,
#                     bar_width, alpha=opacity,
#                     linewidth=0,
#                     color='k', label='Not overlapping with iCLIP targets')
    plt.xlabel('Bins by SAM rank (100 genes per bin, highest ranked targets on the left)')
    plt.ylabel('Overlap with iCLIP')
    plt.xticks(index + bar_width/2., [str(x) for x in range(0, n_groups)])
    plt.tight_layout()
    plt.savefig('figs/overlapping_with_clip_by_sam_rank.pdf', format='pdf')
    plt.clf()
    plt.close()
    # Second figure.
    try:
        in_two_overlaps = [
            float(len(list(in_two[x]['shared'])))/float(
                len(list(in_two[x]['shared'])) + len(list(in_two[x]['unique']))) for x in in_two
        ]
    except ZeroDivisionError:
        return
    plt.ylim([0, 0.7])
    plt.bar([0,1], in_two_overlaps, 0.7, alpha=0.9,
            color='k')
    plt.xticks([0.35, 1.35], ['Above 532', 'Below 532'])
    plt.ylabel('Overlap with iCLIP')

    plt.tight_layout()
    plt.savefig('figs/overlap_with_clip_above_and_below_532_sam_rank.pdf', format='pdf')
    plt.clf()
    plt.close()


def build_bins(rip, clip):
    bins = {}
    in_two = {'above': {'shared': set(), 'unique': set()},
              'below': {'shared': set(), 'unique': set()}}
    bin_num = -1
    for n, gene in enumerate(rip):
        if not n % 100:
            bin_num += 1
            bins[bin_num] = {'shared': set(), 'unique': set()}
        add_to = 'above'
        if n >= 532:
            add_to = 'below'
        if gene in clip:
            bins[bin_num]['shared'] |= set([gene])
            in_two[add_to]['shared'] |= set([gene])
        else:
            bins[bin_num]['unique'] |= set([gene])
            in_two[add_to]['unique'] |= set([gene])
    return bins, in_two


def significance_of_overlap(rip, clip, gtf):
    """
    :param rip: set
    :param clip: set
    :param gtf: set
    :return:
    """
    # clips = set(clip.keys())
    # rips = set(rip.keys())
    both = clip & rip
    clip_onlys = clip - rip
    rip_onlys = rip - clip
    eithers = clip | rip
    neither = gtf - eithers
    td = {
        'neither': len(neither),
        'clip_only': len(clip_onlys),
        'rip_only': len(rip_onlys),
        'both': len(both)
    }
    table = [
        [td['neither'], td['clip_only']],
        [td['rip_only'], td['both']]
    ]
    import scipy.stats as sps
    odds, pval = sps.fisher_exact(table)
    print "\nTable for fisher test: {t}.\np value of independance {p}.\n".format(
        t=td, p=pval
    )


def load_mrnas_from_file(filename, all_mrna, use_header_val=None):
    """
    :param filename: filename
    :param all_mrna: set
    :return: dict
    """
    strip_end = re.compile(r'[^\d\.]+$')
    rip = load_file_into_dict(filename, use_header_val=use_header_val)
    print "*" * 50
    print "{v}:\nGene IDs: {a}".format(v=filename, a=len(set(rip.keys())))
    # print str(rip)[:1000]
    for gene in set([x for x in rip.keys() if remove_end(x, strip_end) not in all_mrna]):
        del rip[gene]
    print '\t...done loading file as {z} mrna rows.'.format(z=len(rip))
    return rip


def remove_end(name, strip_end):
    # print name
    name = name.split('.')
    if len(name) <= 1:
        # print ">" + strip_end.sub('', "".join(name))
        return strip_end.sub('', "".join(name))
    else:
        # print ">" + strip_end.sub('', ".".join(name[:2]))
        return strip_end.sub('', ".".join(name[:2]))


def run(clipfname, gtf_df, lib=None):
    """
    :param clipfname: Peaks filename.
    :param gtf_df: pandas dataframe of gtf.
    :return:
    """
    # strip_end = re.compile(r'[^\d\.]+$')
    # wb_to_biotype = zip(gtf_df['gene_id'].tolist(), gtf_df['biotype'].tolist())
    # locus_to_biotype = zip(gtf_df['gene_id'].tolist(), gtf_df['biotype'].tolist())
    # locus_to_biotype = [(remove_end(x[0], strip_end), x[1]) for x in locus_to_biotype]
    # all_mrna_wb = set([x[0] for x in wb_to_biotype if x[1] == 'protein_coding'])
    # all_mrna_locus =set([x[0] for x in locus_to_biotype if x[1] == 'protein_coding'])
    # # rip = load_mrnas_from_file(lib['fog3_rip_ranks'], all_mrna_locus, use_header_val='WB ID')
    # # clip = load_mrnas_from_file(clipfname, all_mrna_locus, use_header_val='gene_id')
    # rip = load_mrnas_from_file(lib['fog3_rip_ranks'], all_mrna_locus, use_header_val='Seq ID')
    # clip = load_mrnas_from_file(clipfname, all_mrna_wb, use_header_val='gene_id')

    strip_end = re.compile(r'[^\d\.]+$')
    # Create maps of names to biotype.
    wb_to_biotype = zip(gtf_df['gene_id'].tolist(), gtf_df['biotype'].tolist())
    locus_to_biotype = zip(gtf_df['transcript_name'].tolist(), gtf_df['biotype'].tolist())
    locus_to_biotype = [(remove_end(x[0], strip_end), x[1]) for x in locus_to_biotype]
    # Create sets of protein coding gene names.
    all_mrna_wb = set([x[0] for x in wb_to_biotype if x[1] == 'protein_coding'])
    all_mrna_locus = set([x[0] for x in locus_to_biotype if x[1] == 'protein_coding'])
    # Map of WB ID to locus ID.
    wb_to_locus = dict(zip(gtf_df['gene_id'], gtf_df['transcript_name']))
    # Load RIP and CLIP data.
    rip = load_mrnas_from_file(lib['fog3_rip_ranks'], all_mrna_locus, use_header_val='Seq ID')
    clip = load_mrnas_from_file(clipfname, all_mrna_wb, use_header_val='gene_id')
    clip = dict([(remove_end(wb_to_locus[x], strip_end), clip[x]) for x in clip if x in wb_to_locus])



    print "Protein coding genes only:\n\tgenome: {a} CLIP: {c} RIP: {r}".format(
        a=len(all_mrna_locus), c=len(clip), r=len(rip)
    )
    significance_of_overlap(set(rip.keys()), set(clip.keys()), all_mrna_locus)
#    compare_overlap_as_function_of_rank(rip, clip)
    array, header = build_array_of_ranks_in_each_dataset(rip, clip)
    correlations(array, header)
#    compare_low_sam_ranks(rip, clip, header)
    set_labels = ('FOG-3 CLIP', 'FOG-3 RIP')
    output_name = 'figs/targets_vs_fog3rip.pdf'
    make_venn(set(clip.keys()), set(rip.keys()),
              set_labels=set_labels, output_name=output_name)


def make_venn(
        targsk, targsj, output_name='figs/venn.pdf',
        set_labels=('k', 'j')):
    """
    :param targsk: set
    :param targsj: set
    :param output_name: filename
    :param set_labels: labels for diagram
    :return:
    """
    from matplotlib_venn import venn2, venn3
    import matplotlib.pyplot as plt
    both = targsk & targsj
    only_k = targsk - targsj
    only_j = targsj - targsk
    plt.clf()
    venn2([targsk, targsj], set_labels=set_labels)
    plt.savefig(output_name, format='pdf')
    plt.clf()
    plt.close()


if __name__ == '__main__':
    print "~" * 21
    # onlytargsfname = 'for_comp/both_ratios_and_deseq.txt'
    # onlytargsfname = 'for_comp/only_no_antibody_and_deseq.txt'
    # clipfname = 'deseq_counts_in_region/both_controls.txt'
    # clipfname = onlytargsfname
    import get_args
    args = get_args.get_args()
    clipfname = args.peaks
    run(clipfname)

