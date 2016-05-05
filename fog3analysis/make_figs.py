'''
Do all analysis on FOG-3 iCLIP.
Input:
permutation_peaks/ folder of peaks files.
bedgraph_unnorm/all_exp wig files.

'''
import pandas
import numpy as np
import re
import sys
import collections
import HTSeq
import argparse
import os
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
import scipy.stats as scs

from peaksList import *
import load_bedgraph

def load_peaks(infile):
    pk = peaksList(name='name')
    pk.read_csv(infile)
    pk.read_sp_vs_oo()
    pk.annotate_sp_vs_oo()
    pk.read_sp_vs_oo_as_programs()
    pk.add_rip()
    return pk


def add_locations(pk_list):
    need_loc = [pk for pk in pk_list if ('location' not in pk.df.columns)]
    if len(need_loc) == 0: return pk_list
    ga = load_bedgraph.load_bedgraph('bedgraph_unnorm/all_exp.wig')
    gtf = pandas.read_csv('/opt/lib/gtf_with_names_column.txt',
                          index_col=False, sep='\t').to_dict('records')
    gtf_d = collections.defaultdict(list)
    for row in gtf: gtf_d[row['gene_name']].append(row)
    for pk in pk_list:
        pk.add_location_from_integral(gtf_d, ga)
    return pk_list

def make_pie_biotype(pk, fname='figs/biotype_pie.pdf', simplify=False):
    biotypes = dict(pk.df['biotype'].value_counts())
    if simplify:
        n_bio = collections.defaultdict(int)
        for k, v in biotypes.items():
            if k != 'protein_coding': k = 'ncRNA'
            n_bio[k] += v
        make_pie_from_dict(n_bio, fname)
    else:
        make_pie_from_dict(biotypes, fname)


def make_pie_gender(pk, fname='figs/gender_pie.pdf'):
    programs = dict(pk.df['Program'].value_counts())
    make_pie_from_dict(programs, fname)


def make_pie_location(pk, fname='figs/location_pie.pdf'):
    locs = dict(pk.df['location'].value_counts())
    if 'ncRNA' in locs: del locs['ncRNA']
    make_pie_from_dict(locs, fname)


def make_pie_from_dict(adict, fname):
    labels = adict.keys()
    values = adict.values()
    labels = set(labels)
    print labels
    print tuple(values)
    colors = ['gold', 'brown', 'lightcoral', 'yellowgreen',
              'red', 'blue', 'green']
    colors = colors[:len(labels)]
    plt.clf()
    plt.rcParams['patch.linewidth'] = 0  
    plt.pie(adict.values(), labels=adict.keys(),
            colors=colors, autopct='%1.1f%%')
    plt.axis('equal')
    plt.savefig(fname, format='pdf')
    plt.close()
    plt.clf()


def overlap_rip(pk, rippk, fname='figs/overlap_with_rip_venn.pdf',
                only_mrna=True):
    #'Gene public ID', '
    if only_mrna:
        rip = set(rippk.df['WB ID'].tolist())
        clip = set(pk.df[pk.df['biotype']=='protein_coding']['gene_id'].tolist())
    else:
        rip = set(rippk.df['WB ID'].tolist())
        clip = set(pk.df['gene_id'].tolist())
    v = venn2([rip, clip], ('RIP-chip', 'iCLIP'))
    plt.savefig(fname, format='pdf')
    plt.clf()
    plt.close()
    oo = peaksList('Oo')
    sp = peaksList('Sp')
    rip_oo = set(
        rippk.df[rippk.df['Program']!='Spermatogenic only']['WB ID'].tolist())
    rip_sp = set(
        rippk.df[rippk.df['Program']=='Spermatogenic only']['WB ID'].tolist())
    v = venn3([rip_oo, rip_sp, clip], ('RIP-chip, Not spermatogenic',
                                       'RIP-chip, Spermatogenic',
                                       'iCLIP'))
    plt.savefig(fname.partition('.pdf')[0] + '_sp_split.pdf')
    plt.clf()
    plt.close()
    oo.read_csv('/opt/lib/ortiz/TableS1_oogenic.txt')
    oo.df[oo.df['Ortiz Gonad']!='Not Present']
    sp.read_csv('/opt/lib/ortiz/TableS2_spermatogenic.txt')
    sp.df[sp.df['Ortiz Gonad']!='Not Present']
    all_gl_rnas = set(oo.df['Wormbase ID'].tolist())
    all_gl_rnas |= set(sp.df['Wormbase ID'].tolist())
    table_pat = """
     [[(-R, -C), (+R, -C)],
      [(-R, +C), (+R, +C)]]"""
    fisher_table = [
        [len(all_gl_rnas - rip - clip), len(rip - clip)],
        [len(clip - rip), len(rip & clip)]]
    print "Statistics for %s overlap with RIP:" % fname
    print table_pat
    print "Genome size %i (all germline RNAs)." % len(all_gl_rnas)
    res = scs.fisher_exact(fisher_table)
    print fisher_table
    print res
    fisher_table = [
        [2e4 - len(rip | clip), len(rip - clip)],
        [len(clip - rip), len(rip & clip)]]
    print "Genome size %i (20k gene genome)." % len(all_gl_rnas)
    res = scs.fisher_exact(fisher_table)
    print fisher_table
    print res

"""
Needed config.ini parameters for determine_feature_locations:
lib['bedgraphs_folder']  # In auto.
lib['txpt_obj_path']  # In auto.
lib['coverage_wigs']  # In auto. = bedgraphs_folder.
lib['bedgraph_exp_plus']  # In auto.
lib['bedgraph_exp_minus']  # In auto.
lib['gtf_one_txpt_per_gene']  # In auto.
lib['gtf_pickle']  # In auto.
lib['feat_loc_dir']  # In auto.
lib['motif']  # In auto.
lib['gtf']  # In auto.
lib['gtf_one_txpt_per_gene']  # In auto.
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        default='permutation_peaks/')
    parser.add_argument('-r', '--ratio_is_raw_reads',
                        default=False, action='store_true')
    parser.add_argument('-c', '--config_ini',
                        default='auto.ini')
    args = parser.parse_args()
    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) \
                    + '/feature_locations/')
    import config
    lib = config.config(args.config_ini)
    ##################
    # Make line graph (Fig 1C).
    # Used by determine_feature_locations:
    args.load_txpts = False
    args.heat_map = False
    args.load_gtf = False
    hold_original = args.input
    args.input = args.input + '/combined_exp.txt'
    import determine_feature_locations
    import peak_in_gene_line_plot
    import build_image_objects_for_heatmaps
    (peaks, txpts, chr_lens) = determine_feature_locations.get_data(
        args, lib=lib)
#    build_image_objects_for_heatmaps.heatmap_of_raw_signal(
#        peaks, txpts, output_dirname='figs/Fig1E_raw_signal_heatmap.pdf')
    peak_in_gene_line_plot.normalize_distances(txpts)
    (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_peak_v_polyA,
     ave_fbe_v_polyA, ave_highest_peak, ave_secondary_peaks
     ) = peak_in_gene_line_plot.get_ave_pos(txpts)
    peak_in_gene_line_plot.plot_features(
       (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_highest_peak,
        ave_secondary_peaks),
       output_filename='figs/Fig1Cmaybe_features_in_normalized_gene.pdf')
        #determine_feature_locations.make_figs(
        #    peaks, txpts, chr_lens, args, hold_original)
    sys.exit()
    args.input = hold_original
    ##################
    cims_file = 'tables/cims_annotated.txt'
    pk_cims = load_peaks(cims_file)
    #if os.path.isdirectory(args.input):
    pk_pos = load_peaks(os.path.join(args.input, 'combined_exp.txt'))
#        pk_neg = load_peaks(os.path.join(args.input, 'combined_control.txt'))
    #pk.add_permutation_peaks(clip_fname)
    pk_cims = add_locations([pk_cims])[0]
    rip = peaksList('RIP') 
    rip.read_csv('/opt/lib/fog3rip/TableS4_fog3_rip_targets.txt')
    rip.read_sp_vs_oo()
    #rip.annotate_sp_vs_oo()
    rip.read_sp_vs_oo_as_programs()
    pk_cims.df['gene_id'] = pk_cims.df['Wormbase ID']
    overlap_rip(pk_cims, rip, fname='figs/overlap_with_rip_venn_cims.pdf')
    overlap_rip(pk_pos, rip, fname='figs/overlap_with_rip_venn_clip.pdf')
    make_pie_location(pk_pos, fname='figs/Fig1D_location_pie_exp.pdf')
    make_pie_biotype(pk_cims, fname='figs/biotype_cims.pdf')
    make_pie_gender(pk_cims, fname='figs/gender_cims_pie.pdf')
    make_pie_biotype(pk_pos, fname='figs/Fig1A_alt_biotype_exp.pdf')
    make_pie_biotype(pk_pos, fname='figs/Fig1A_biotype_exp.pdf',
                     simplify=True)
    make_pie_gender(pk_pos, fname='figs/Fig1B_gender_clip_pie.pdf')
    pk_cims.to_csv('tables/ann.txt')
