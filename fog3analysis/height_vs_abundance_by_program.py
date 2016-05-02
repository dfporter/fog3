import pandas
import re
import os
import sys
import glob
import numpy as np
from peaksList import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scs
import matplotlib.lines as mlines
import matplotlib

def height_vs_sp_exp_from_df(df, xcol='Sp. expression', ycol='height'):
    y = df[ycol].tolist()
    y = [np.log10(_y) for _y in y]
    x_sp = df[xcol].tolist()
    x_sp = [np.log10(_x) for _x in x_sp] 
    return x_sp, y

def define_plot(ax, x, y, subplotx, subploty):
    ax[subplotx][subploty].hexbin(x, y, gridsize=30)
    ax[subplotx][subploty].set_ylabel('CLIP peak height log10')
    ax[subplotx][subploty].set_xlabel('Sp. RPKM log10')
#    ax[0][0].axis([ 0.25, 1.5, -1, 3])

def get_id_to_height(mapping):
    def id_to_height(gene_id):
        if gene_id in mapping:
            return mapping[gene_id]
        else:
            return 0
    return id_to_height

def get_id_to_val(mapping, null_return=0):
    def id_to_val(gene_id):
        if gene_id in mapping:
            return mapping[gene_id]
        else:
            return null_return
    return id_to_val


def split_by_program(ortiz_df):
    a = ortiz_df[ortiz_df['Program']=='Oogenic only'].copy()
    b = ortiz_df[ortiz_df['Program']=='Spermatogenic only'].copy()
    c = ortiz_df[ortiz_df['Program']=='Oogenic and Spermatogenic'].copy()
    return a, b, c

def both_ok(a, b, i):
    if (np.isfinite(a[i]) and np.isfinite(b[i])) and not (
        np.isnan(a[i]) or np.isnan(b[i])):
        return True
    else:
        return False

def get_oo_sp_exp(sub, col_a='Sp. RPKM',
                  col_b='Expression in fem-3 gonads (RPKM)'):
    cols = sub[[col_a, col_b]].values
    cols = np.apply_along_axis(np.log10, 0, cols)
    a = cols[:,0]
    b = cols[:,1]
    x, y = [], []
    for i in range(len(cols[:,0])):
        if both_ok(a, b, i):
            x.append(a[i])
            y.append(b[i]) 
    v = np.array([x, y]).T
    return np.array([x, y]).T

def close_fig(outname):
    plt.savefig(outname, format='pdf') 
    plt.clf()
    plt.close()


def aspect_equal(ax1):
    x0,x1 = ax1.get_xlim()
    y0,y1 = ax1.get_ylim()
    ax1.set_aspect((x1-x0)/(y1-y0))
    
def cor_with_abundance_by_program(pk):
    blackline = mlines.Line2D([], [], color='black', linewidth=0.01) 
    redline = mlines.Line2D([], [], color='red', linewidth=0.01)
    pinkline = mlines.Line2D([], [], color='pink', linewidth=0.01)
    blueline = mlines.Line2D([], [], color='blue', linewidth=0.01)
    orangeline = mlines.Line2D([], [], color='orange', linewidth=0.01)
    greenline = mlines.Line2D([], [], color='green', linewidth=0.01)
    df = pk.df
    ortiz_df = pandas.read_csv(
        '/opt/lib/ortiz/DESeq_genes_in_gonad.txt', sep='\t')
    id_to_fog3_height = dict(zip(df['gene_name'].tolist(), df['height'].tolist()))
    id_to_height = get_id_to_val(id_to_fog3_height, null_return=0)
    id_to_program = get_id_to_val(pk.program, null_return='')
    ortiz_df['Sp. RPKM'] = ortiz_df['Expression in  fog-2 gonads (RPKM)']
    ortiz_df['fog3 height'] = [id_to_height(x) for x in ortiz_df['Gene name'].tolist()]
    ortiz_df['Program'] = [id_to_program(x) for x in ortiz_df['WormBase ID (WS240)'].tolist()]
    # Subset to FOG-3 targets.
    fog3_targets = get_oo_sp_exp(ortiz_df[ortiz_df['fog3 height'] != 0].copy())
    oo, sp, ne = split_by_program(ortiz_df)
    oo_exp = rotate(get_oo_sp_exp(oo))
    sp_exp = rotate(get_oo_sp_exp(sp))
    ne_exp = rotate(get_oo_sp_exp(ne))
    ort = get_oo_sp_exp(ortiz_df)
    gl_exp = rotate(ort)
    tup_fog3_targets = rotate(fog3_targets)
    fig, ax = plt.subplots(ncols=2, nrows=2)
    abundance_kde = scs.gaussian_kde(gl_exp[:,0])
    abundance_kde_targ = scs.gaussian_kde(tup_fog3_targets[:,0])
    ahist = np.histogram(tup_fog3_targets[:,0])
    print '*********'
    ax[0][0].scatter(gl_exp[:,0], gl_exp[:,1],
                     linewidths=0, alpha=0.05, c='k')
    ax[0][0].scatter(tup_fog3_targets[:,0], tup_fog3_targets[:,1],
                     linewidths=0, alpha=0.25, c='r')
    font = {'size':7}
    plt.rc('font', size=7)
    def add_xylabel(_ax):
        _ax.set_xlabel('Overall abundance log10 RPKM')
        _ax.set_ylabel('Enrichment in Sp GL (log10 SP/OO)')
    add_xylabel(ax[0][0])
    aspect_equal(ax[0][0])
    ax[0][0].legend([blackline, redline], ['All GL RNAs', 'FOG-3 targets'],
                 fontsize='xx-small', loc='best')
    #ax[0].set(adjustable='box-forced', aspect='equal')
    ax[0][1].scatter(oo_exp[:,0], oo_exp[:,1],
                     linewidths=0, alpha=0.25, c='pink')
    ax[0][1].scatter(sp_exp[:,0], sp_exp[:,1],
                     linewidths=0, alpha=0.25, c='b')
    ax[0][1].scatter(ne_exp[:,0], ne_exp[:,1],
                     linewidths=0, alpha=0.05, c='k')
    add_xylabel(ax[0][1])
    ax[0][1].legend([pinkline, blueline, blackline],
                    ['Oo. program', 'Sp. program', 'Ne. program',],
                 fontsize='xx-small', loc='best')
    aspect_equal(ax[0][1])
    import annotate_cims
    cims_pk = annotate_cims.cims_peaks(
        'tables/cims_annotated.txt', 'permutation_peaks/combined_exp.txt')
    cims = get_oo_sp_exp(
        cims_pk.df.copy(), col_a='Sp. expression', col_b='Oo. expression')
    cims_exp = rotate(cims)
    ax[1][0].scatter(gl_exp[:,0], gl_exp[:,1],
                     linewidths=0, alpha=0.05, c='k')
    ax[1][0].scatter(cims_exp[:,0], cims_exp[:,1],
                  linewidths=0, alpha=0.5, c='r')
    ax[1][0].legend([blackline, redline],
                    ['All GL RNAs', 'FOG-3 CIMS'],
                 fontsize='xx-small', loc='best')
    add_xylabel(ax[1][0])
    aspect_equal(ax[1][0])
    rip = peaksList(name="RIP", gene_name_col='WB ID')
    rip.read_csv('/opt/lib/fog3rip/TableS4_fog3_rip_targets.txt')
    rip.read_sp_vs_oo()
    rip.annotate_sp_vs_oo()
    rip.add_permutation_peaks('permutation_peaks/combined_exp.txt')
    print "**"
    print rip.df.columns
    print rip.df.head(5)
    rip_exp_shared = rotate(get_oo_sp_exp(
        rip.df[rip.df['Is CLIP target?']].copy(),
        col_b='Sp. expression', col_a='Oo. expression'))
    rip_exp_unshared = rotate(get_oo_sp_exp(
        rip.df[~rip.df['Is CLIP target?']].copy(),
        col_b='Sp. expression', col_a='Oo. expression'))
    # Ortiz peaksList.
    ortiz_pk = peaksList(name='RNA-seq', gene_name_col='Gene name')
    ortiz_pk.read_csv('/opt/lib/ortiz/DESeq_genes_in_gonad.txt')
    ortiz_pk.read_sp_vs_oo()
    ortiz_pk.annotate_sp_vs_oo()
    ortiz_pk.read_sp_vs_oo_as_programs()
    ortiz_pk.add_rip()
    ortiz_pk.add_permutation_peaks('permutation_peaks/combined_exp.txt')
    clip_targ = ortiz_pk.df[ortiz_pk.df['Is CLIP target?']].copy()
    clip_exp_unshared = rotate(get_oo_sp_exp(
        clip_targ[~clip_targ['Is RIP target?']].copy(),
        col_b='Sp. expression', col_a='Oo. expression'))

    gl_exp = get_oo_sp_exp(
        ortiz_pk.df.copy(), col_b='Sp. expression', col_a='Oo. expression')
    gl_exp = rotate(gl_exp)
    #ax[1][1].scatter(gl_exp[:,0], gl_exp[:,1],
    #                 linewidths=0, alpha=0.01, c='k')
    ax[1][1].scatter(clip_exp_unshared[:,0], clip_exp_unshared[:,1],
                     linewidths=0, alpha=0.1, c='green')
    ax[1][1].scatter(rip_exp_unshared[:,0], rip_exp_unshared[:,1],
                     linewidths=0, alpha=0.1, c='b')
    ax[1][1].scatter(rip_exp_shared[:,0], rip_exp_shared[:,1],
                     linewidths=0, alpha=0.1, c='r')

    ax[1][1].legend([redline, blueline, greenline],
                    ['FOG-3 RIP-chip and CLIP targets',
                     'FOG-3 RIP-chip only', 'FOG-3 CLIP target only'],
                 fontsize='xx-small', loc='best')
    add_xylabel(ax[1][1])
    aspect_equal(ax[1][1])
    #ax[1].set(adjustable='box-forced', aspect='equal')
    #plt.tight_layout()
    close_fig('height_v_exp_by_prog_gauss.pdf')
    fig, ax = plt.subplots(2,2)
    x_grid = np.linspace(-2, 6, num=1e3)
#    ax[0][1].plot(x_grid, abundance_kde.pdf(x_grid), 'k--')
#    ax[0][1].plot(x_grid, abundance_kde_targ.pdf(x_grid), 'r--')
    abundance_kde = scs.gaussian_kde(gl_exp[:,1].ravel())
    abundance_kde_targ = scs.gaussian_kde(tup_fog3_targets[:,1].ravel())
    x_grid = np.linspace(-4, 5, num=1e3)
    ax[0][0].plot(x_grid, abundance_kde.pdf(x_grid), 'k--')
    ax[0][0].plot(x_grid, abundance_kde_targ.pdf(x_grid), 'r--')
    ax[0][0].set_xlabel('Enrichment in Sp. GLs (log10 RPKM)')
    ax[0][0].set_ylabel('Frequency')
    ax[0][0].legend([blackline, redline],
                    ['All GL RNAs', 'FOG-3 targets'], fontsize='xx-small')
    ax[0][1].plot(x_grid,
            scs.gaussian_kde(oo_exp[:,1].ravel()).pdf(x_grid), 'pink')
    ax[0][1].plot(x_grid,
            scs.gaussian_kde(sp_exp[:,1].ravel()).pdf(x_grid), 'b-')
    ax[0][1].plot(x_grid,
            scs.gaussian_kde(ne_exp[:,1].ravel()).pdf(x_grid), 'k-')
    ax[0][1].plot(x_grid,
            scs.gaussian_kde(tup_fog3_targets[:,1].ravel())(x_grid), 'r--')
    ax[0][1].set_xlabel('Enrichment in Sp. GLs (log10 RPKM)')
    ax[0][1].set_ylabel('Frequency')
    ax[0][1].legend([pinkline, blueline, blackline, redline],
                    ['Oo', 'Sp', 'Ne', 'FOG-3 target'], fontsize='xx-small')
    ax[1][0].plot(x_grid,
                  scs.gaussian_kde(gl_exp[:,0].ravel())(x_grid), 'k--')
    ax[1][0].plot(x_grid,
                  scs.gaussian_kde(tup_fog3_targets[:,0].ravel())(x_grid),
                  'r--')
    ax[1][0].legend([blackline, redline],
                    ['All GL RNAs', 'FOG-3 targets'], fontsize='xx-small',
                    loc='best')
    ax[1][0].set_xlabel('Overall RNA abundance (log10 RPKM)')
    ax[1][0].set_ylabel('Frequency')
    ax[1][1].plot(x_grid,
            scs.gaussian_kde(oo_exp[:,0].ravel()).pdf(x_grid), 'pink')
    ax[1][1].plot(x_grid,
            scs.gaussian_kde(sp_exp[:,0].ravel()).pdf(x_grid), 'b-')
    ax[1][1].plot(x_grid,
            scs.gaussian_kde(ne_exp[:,0].ravel()).pdf(x_grid), 'k-')
    ax[1][1].plot(x_grid,
            scs.gaussian_kde(tup_fog3_targets[:,0].ravel())(x_grid), 'r--')
    ax[1][1].legend([pinkline, blueline, blackline, redline],
                    ['Oo', 'Sp', 'Ne', 'FOG-3 target'], fontsize='xx-small',
                    loc='upper left')
    ax[1][1].set_ylabel('Frequency')
    ax[1][1].set_xlabel('Overall RNA abundance (log10 RPKM)')
    plt.tight_layout()
    close_fig('gaussian_cfs.pdf')
    # TO-DO
    # Scatter (or heatmap using contourf)
    # X axis RPKM in Sp. germline
    # Y axis negative IP gene RPKM (second panel: FOG-3 gene RPKM)
    # Color = spermatogenic or not (or target or not?)
    # Idea is to see how abundance, program and targetting are related.

    fig, ax = plt.subplots(2, 3)
    x_sp, y = height_vs_sp_exp_from_df(df[df['Program']=='Oogenic only'])
    define_plot(ax, x_sp, y, 0, 0)
    print scs.spearmanr(x_sp, y)
    x_sp, y = height_vs_sp_exp_from_df(df[df['Program']=='Spermatogenic only'])
    define_plot(ax, x_sp, y, 0, 1)
    print scs.spearmanr(x_sp, y)
    x_sp, y = height_vs_sp_exp_from_df(df[df['Program']=='Oogenic and Spermatogenic'])
    define_plot(ax, x_sp, y, 0, 2)
    print scs.spearmanr(x_sp, y)
    sub = ortiz_df[ortiz_df['Program']=='Oogenic only']
    y_oo = [np.log10(y) for y in sub['fog3 height'].tolist() if y != 0]
    sub = ortiz_df[ortiz_df['Program']=='Spermatogenic only']
    y_sp = [np.log10(y) for y in sub['fog3 height'].tolist() if y != 0]
    sub = ortiz_df[ortiz_df['Program']=='Oogenic and Spermatogenic']
    y_ne = [np.log10(y) for y in sub['fog3 height'].tolist() if y != 0]

    ax[1][0].boxplot([y_oo, y_sp, y_ne], [1, 3, 5])
    ax[1][0].set_xticks([1, 2, 3])
    ax[1][0].set_xticklabels(['Oo.', 'Sp.', 'Ne.'])
    ax[1][0].set_ylabel('FOG-3 height log10')
    sub = ortiz_df[ortiz_df['Program']=='Oogenic only']
    y_oo = [np.log10(y) for y in sub['Sp. RPKM'].tolist() if y != 0]
    sub = ortiz_df[ortiz_df['Program']=='Spermatogenic only']
    y_sp = [np.log10(y) for y in sub['Sp. RPKM'].tolist() if y != 0]
    sub = ortiz_df[ortiz_df['Program']=='Oogenic and Spermatogenic']
    y_ne = [np.log10(y) for y in sub['Sp. RPKM'].tolist() if y != 0]
    ax[1][1].boxplot([y_oo, y_sp, y_ne], [1, 3, 5])
    ax[1][1].set_ylabel('Sp. RPKM log10')
    ax[1][1].set_xticklabels(['Oo.', 'Sp.', 'Ne.'])
    n_df = df[df['Classification OO/SP GL']=='Gender Neutral']
    oo_df = df[df['Classification OO/SP GL']=='Oogenic']
    sp_df = df[df['Classification OO/SP GL']=='Spermatogenic']
    sub = ortiz_df[ortiz_df['fog3 height'] != 0].copy()
    y_sp_fog3_targets = [np.log10(y) for y in sub['Sp. RPKM'].tolist() if y != 0]
    y_oo_fog3_targets = [np.log10(y) for y in sub['Expression in fem-3 gonads (RPKM)'].tolist() if y != 0]
    y_sp = [np.log10(y) for y in ortiz_df['Sp. RPKM'].tolist() if y != 0]
    y_oo = [np.log10(y) for y in ortiz_df['Expression in fem-3 gonads (RPKM)'].tolist() if y != 0]
    tup = zip(y_sp, y_oo)
    tup_fog3_targets = zip(y_sp_fog3_targets, y_oo_fog3_targets)
    ax[1][2].scatter([x[0] for x in tup], [x[1] for x in tup],
                     linewidths=0, alpha=0.1, c='k')
    ax[1][2].scatter([x[0] for x in tup_fog3_targets],
                     [x[1] for x in tup_fog3_targets],
                     linewidths=0, alpha=0.4, c='r')
    ax[1][2].plot([-4, 5], [-4, 5], 'b--')
    ax[1][2].set_xlabel('Sp. RPKM log10')
    ax[1][2].set_ylabel('Oo. RPKM log10')
    ax[1][2].set_xlim([-4, 5])
    ax[1][2].set_ylim([-4, 5])
    print scs.ttest_ind(
        [x[0]/x[1] for x in tup], [x[0]/x[1] for x in tup_fog3_targets])
    plt.tight_layout()
    plt.savefig('height_v_exp_by_prog.pdf', format='pdf')
    plt.clf()
    plt.close()


def rotate(arr2d, deg=-45):
    theta = np.radians(deg)
    c, s = np.cos(theta), np.sin(theta)
    R = np.matrix([[c, -s], [s, c]])
    rotated_x = []
    rotated_y = []
    for t in arr2d:
        new = R * np.matrix(t).T
        new = np.array(new.T)[0]
        rotated_x.append(new[0])
        rotated_y.append(new[1])
    print "ROTATED:"
    print np.array([rotated_x, rotated_y]).T
    return np.array([rotated_x, rotated_y]).T

def clean(y_fc, y_oo, y_sp):
    tup = zip(y_fc, y_oo, y_sp)
    tup = [x for x in tup if np.all(np.isfinite(x))]
    return (
        [x[0] for x in tup], [x[1] for x in tup], [x[2] for x in tup])


def get_limits(x, y):
    x = np.array(x)
    y = np.array(y)
    a = np.nanmin(x)
    b = np.nanmax(x)
    c = np.nanmin(y)
    d = np.nanmax(y)
    return [a, b, c, d]

if __name__ == '__main__':
    indir = sys.argv[1]
    for fname in glob.glob(indir + '/*exp*'):
        pk = peaksList(name='FOG-3')
        pk.read_csv(fname)
        pk.read_sp_vs_oo()
        pk.annotate_sp_vs_oo()
        # This next line creates a 'Program' column using the Noble method.
        pk.read_sp_vs_oo_as_programs()
        cor_with_abundance_by_program(pk)
        pk.read_csv(fname)
