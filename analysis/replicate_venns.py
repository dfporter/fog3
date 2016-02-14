from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
import collections
import os
import pandas
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage
import matplotlib.colors as mcolors
import sys
import glob
import cluster_combine
sys.path.insert(0, '/groups/Kimble/Aman Prasad/clip/analysis/src/feature_locations/')
import filter_by_reads_per_gene
del sys.path[0]

def make_venn3(
        targs1, targs2, targs3,
        output_name='figs/venn.pdf', set_labels=('1', '2', '3')):
    plt.clf()
    venn2([targs1, targs2, targs3], set_labels=set_labels)
    plt.savefig(output_name, format='pdf')
    plt.clf()
    plt.close()


def make_venn2(
        targsk, targsj, output_name='figs/venn.pdf',
        set_labels=('k', 'j')):
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


def run(reps_folder, combined_folder,
        only_use_these_targets=None,
        min_rpkm_cutoff=None,
            min_rpkm_enrichment_cutoff=None, lib=None):
    # Clusters.
    list1 = glob.glob(reps_folder + '/exp*.txt')
    if len(list1) < 2:
        list1 = glob.glob(reps_folder + '/fog*.txt')
        if len(list1) < 2:
            print "Error - could not find experimental files in clusters directory."
            return False
    overlap_matrix(list1, only_use_these_targets=only_use_these_targets,
        pvalue_cutoff=None,
        min_rpkm_cutoff=min_rpkm_cutoff,
        min_rpkm_enrichment_cutoff=min_rpkm_enrichment_cutoff, lib=lib)
    # expsf = glob.glob(reps_folder + '/fog*.txt')
    # exps, _ = list_of_targs_in_reps(expsf)
    # controlsf = glob.glob(reps_folder + '/control*.txt')
    # controls, _ = list_of_targs_in_reps(controlsf)
    # Combined peaks.
    exp, _ = list_of_targs_in_reps(
        [combined_folder + '/combined_exp.txt'], lib=lib)
    exp = exp[0]
    control, _ = list_of_targs_in_reps(
        [combined_folder + '/combined_control.txt'], lib=lib)
    control = control[0]
    make_venn2(exp, control, output_name='exp_vs_control_venn.pdf',
               set_labels=('FOG-3', 'No anti-flag antibody IP'))


def overlap_matrix(file_list, only_use_these_targets=None,
            pvalue_cutoff=None,
            min_rpkm_cutoff=None,
            min_rpkm_enrichment_cutoff=None, lib=None):
    _, gene_to_reps = list_of_targs_in_reps(
        file_list,
        pvalue_cutoff=pvalue_cutoff, min_rpkm_cutoff=min_rpkm_cutoff,
        min_rpkm_enrichment_cutoff=min_rpkm_enrichment_cutoff, lib=lib)
    known_rep = set()
    if only_use_these_targets is not None:
        for gene in [x for x in gene_to_reps.keys() if x not in only_use_these_targets]:
            del gene_to_reps[gene]
    for gene in gene_to_reps:
        known_rep |= gene_to_reps[gene]
    rep_order = list(known_rep)
    as_number = collections.defaultdict(list)
    for gene in gene_to_reps:
        for rep in rep_order:
            if rep in gene_to_reps[gene]:
                as_number[gene].append(1)
            else:
                as_number[gene].append(0)
    list2d = as_number.values()
    arr2d = np.array([np.array(x) for x in list2d])
    # for r in list2d:
    #     if max(r) == 0:
    #         print r
    dists_image = scipy.spatial.distance.pdist(arr2d, metric='euclidean')
    Z = sch.linkage(dists_image)
    heatmap_order = sch.leaves_list(Z)
    arr2d = arr2d[heatmap_order,:]
    # print arr2d
    # for row in arr2d:
    #     if max(row) == 0:
    #         print row
    # import matplotlib.colors as mcolors
    with open('replicate_raw_data.txt', 'w') as f: f.write(
        "\n".join(["\t".join(str(x) for x in row) for row in arr2d])
    )
    make_fig_from_file('replicate_raw_data.txt')


def make_fig_from_file(fname):
    arr2d = []
    with open(fname, 'r') as f:
        for li in f:
            arr2d.append(np.array([int(x) for x in li.rstrip('\n').split('\t')]))
    arr2d = np.array(arr2d)
    plt.clf()
    plt.xlim([0, len(arr2d[0])])
    plt.ylim([0, len(arr2d)])
    fig, ax1 = plt.subplots()
    # palette = plt.cm.get_cmap('gray')
    # palette.set_bad('white')
    cmap, norm = mcolors.from_levels_and_colors([-1, 0.5, 2], ['white', 'black'])
    cax = ax1.imshow(arr2d, interpolation='None', #aspect='auto',
                cmap=cmap)#, norm=norm)#, origin='lower')
    n_thousands = len(arr2d)/1e3
    ax1.set_yticks([x for x in range(len(arr2d) - 500, 1, -500)])
    ax1.set_yticklabels([str(len(arr2d) - x) for x in range(len(arr2d) - 500, 1, -500)])
    ax1.xaxis.set_major_formatter(plt.NullFormatter())
    ax1.xaxis.set_minor_formatter(plt.NullFormatter())
    ax1.set_aspect(2./ax1.get_data_ratio())
    plt.tick_params(
        axis='x',
        which='both',
        top='off',
        bottom='off'
    )
    plt.tick_params(
        axis='y',
        which='both',
        right='off',
    )
    ax1.set_ylabel('Gene')
    ax1.set_xlabel('FOG-3 iCLIP replicate')
    if not os.path.exists('figs/'): os.system('mkdir figs')
    plt.savefig(
        'figs/replicate_overlap_heatmap_bw_total_genes_%s_max_pvalue_%s.pdf' % (
            str(len(arr2d)), 'from_file'
    ), format='pdf')


def list_of_targs_in_reps(file_list, pvalue_cutoff=None,
                          min_rpkm_cutoff=None,
                          min_rpkm_enrichment_cutoff=None,
                          lib=None):
    db = {}
    # targs = collections.defaultdict(list)
    targs = []
    gene_to_reps = collections.defaultdict(set)
    for fname in file_list:
        bname = os.path.basename(fname).partition('.txt')[0]
        counts_name = bname + '_counts.txt'
        db[bname] = pandas.read_csv(fname, sep='\t', index_col=False)
        if (min_rpkm_cutoff is not None) or (min_rpkm_enrichment_cutoff is not None):
            db[bname] = cluster_combine.get_rpkm(
                db[bname], lib)
            if counts_name not in db[bname].columns:
                print "Could not find in {u} column in {a}".format(
                    u=counts_name, a=db[bname].columns
                )
                sys.exit()
        if min_rpkm_cutoff is not None:
            counts_name = bname + '_counts.txt'
            print "applying min_rpkm at {o}: before {g} rows".format(
                o=min_rpkm_cutoff, g=len(db[bname].index)
            )
            db[bname] = db[bname][db[bname][counts_name]>min_rpkm_cutoff]
            print "after: {o}".format(o=len(db[bname].index))
        if min_rpkm_enrichment_cutoff is not None:
            print "applying min_rpkm_enrichment_cutoff at {o}: before {g} rows.".format(
                o=min_rpkm_enrichment_cutoff, g=len(db[bname].index)
            )
            db[bname] = filter_by_reads_per_gene.sum_reads_in_gene(db[bname])
            db[bname]['keep'] = 0
            db[bname]['ave_control'] = [float(x)/4.0 for x in db[bname]['control']]
            for i, row in db[bname].iterrows():
                if row[counts_name] > row['ave_control']:
                    db[bname].loc[i, 'keep'] = 1
            db[bname] = db[bname][db[bname]['keep']>0]
            print "after: {o}".format(o=len(db[bname].index))
        if ('gene_id' in db[bname].columns) and ('gene_name' not in db[bname].columns):
            db[bname]['gene_name'] = db[bname]['gene_id']
        if ('gene' in db[bname].columns) and ('gene_name' not in db[bname].columns):
            # print db[bname]['gene_name']
            # print db[bname]['gene']
            aseries = db[bname]['gene']
            print'len of gene col {l}'.format(l=len(aseries))
            db[bname]['gene_name'] = 0
            print 'len of empty col {l}, db.index {p}'.format(
                l=len(db[bname]['gene_name']),
                p=len(db[bname].index))
            db[bname]['gene_name'] = aseries
        # if 'padj' in db[bname].columns:
        #     db[bname] = db[bname][db[bname]['padj']<pvalue_cutoff]
        for gene in set(db[bname]['gene_name'].tolist()):
            gene_to_reps[gene] |= set([bname])
        targs.append(set(db[bname]['gene_name'].tolist()))
    return targs, gene_to_reps


if __name__ == '__main__':
    run('fog3/clusters/', 'fog3/permutation_peaks/', lib={'gtf': 'lib/gtf_with_names_column.txt'})
