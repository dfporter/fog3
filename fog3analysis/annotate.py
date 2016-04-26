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


def cor_with_abundance(df):
    print df[df['gene_name']=='fog-3']
    print df.head(1)
    x = df['height'].tolist()
    x = [np.log10(_x) for _x in x]
    fig, ax = plt.subplots(2, 3)
    y_oo = df['Oo. expression'].tolist()
    y_oo = [np.log10(_y) for _y in y_oo] 
    ax[0][0].hexbin(x, y_oo, gridsize=30)
    ax[0][0].set_xlabel('CLIP peak height log10')
    ax[0][0].set_ylabel('Oo. RPKM log10')
    print scs.spearmanr(x, y_oo)
    ax[0][0].axis([ 0.25, 1.5, -1, 3])
    y_sp = df['Sp. expression'].tolist()
    y_sp = [np.log10(_y) for _y in y_sp]
    print scs.spearmanr(x, y_sp)
    ax[0][1].hexbin(x, y_sp, gridsize=30)
    ax[0][1].axis([0.25, 1.5, -1, 3])
    ax[0][1].set_xlabel('CLIP peak height log10')
    ax[0][1].set_ylabel('Sp. RPKM log10')
    y_fc = df['Fold change SP/OO GL'].tolist()
    y_fc = [np.float(_y) for _y in y_fc]
#    y = [np.log10(_y) for _y in y]
    ax[0][2].hexbin(x, y_fc, gridsize=20)
    ax[0][2].set_xlabel('CLIP peak height log10')
    ax[0][2].set_ylabel('Fold change SP/OO expression')
    print scs.spearmanr(x, y_fc)
    ax[1][0].hexbin(y_sp, y_oo, gridsize=30)
    #ax[1][1].hexbin(y_fc, y_oo, gridsize=30)
    #ax[1][1].axis([-3, 3, -1, 3])
    #ax[1][2].hexbin(y_fc, y_sp, gridsize=30)
    #ax[1][2].axis([-3, 3, -1, 3])
    (y_fc, y_oo, y_sp) = clean(y_fc, y_oo, y_sp)
    print "SP/OO log2 fold change, Oo. RPKM. Spearman"
    print scs.spearmanr([x for x in y_fc], y_oo)
    print "PEARSON"
    print "SP/OO log2 fold change, Oo. RPKM. Pearson"
    print scs.pearsonr([x for x in y_fc], y_oo)
    print "SP/OO log2 fold change, Sp. RPKM. Spearman"
    print scs.spearmanr(y_fc, y_sp)
    # Are more abundant RNAs more SP or OO?
    n_df = df[df['Classification OO/SP GL']=='Gender Neutral']
    oo_df = df[df['Classification OO/SP GL']=='Oogenic']
    sp_df = df[df['Classification OO/SP GL']=='Spermatogenic']
    n_heights = [np.log10(x) for x in n_df.height]
    oo_heights = [np.log10(x) for x in oo_df.height]
    sp_heights = [np.log10(x) for x in sp_df.height]
    ex = {}
    ex_lens = {}
    ex_medians = {}
    ex['n_rna_oo_exp'] = [np.log10(x) for x in n_df['Oo. expression'].tolist()]
    ex['oo_rna_oo_exp'] = [np.log10(x) for x in oo_df['Oo. expression'].tolist()]
    ex['sp_rna_oo_exp'] = [np.log10(x) for x in sp_df['Oo. expression'].tolist()]
    ex['n_rna_sp_exp'] = [np.log10(x) for x in n_df['Sp. expression'].tolist()]
    ex['oo_rna_sp_exp'] = [np.log10(x) for x in oo_df['Sp. expression'].tolist()]
    ex['sp_rna_sp_exp'] = [np.log10(x) for x in sp_df['Sp. expression'].tolist()]
    print "Expressions"
    for y in ex:
        ex_lens[y] = len(ex[y])
    print ex_lens
    table = [ex['n_rna_oo_exp'], ex['n_rna_sp_exp'],
         ex['oo_rna_oo_exp'], ex['oo_rna_sp_exp'],
         ex['sp_rna_oo_exp'], ex['sp_rna_sp_exp']]
    ax[1][1].boxplot(
        table, [1, 2, 4, 5, 7, 8])
    print "N vs OO peak height"
    print scs.ttest_ind(n_heights, oo_heights)
    print "N vs SP peak height"
    print scs.ttest_ind(n_heights, sp_heights)
    print "OO vs SP peak height"
    print scs.ttest_ind(oo_heights, sp_heights)
    ax[1][2].violinplot([n_heights, oo_heights, sp_heights],
                        [1, 3, 5], showmedians=True)
#    ax[1][2].set_ylim([0,10])
    plt.tight_layout()
    plt.savefig('corwithabundance.pdf', format='pdf')

"""

"""
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
    for fname in glob.glob(indir + '/*exp*txt'):
        pk = peaksList(name='FOG-3')
        pk.read_csv(fname)
        print pk
        pk.read_sp_vs_oo()
        pk.annotate_sp_vs_oo()
        pk.read_sp_vs_oo_as_programs()
        print pk
        cor_with_abundance(pk.df)
