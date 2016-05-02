import pandas
from peaksList import *
import os
import sys


from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles

def make_venn(pk, fname='figs/cims_rip_clip_venn.pdf'):
    rip = pandas.read_csv(
        '/opt/lib/fog3rip/TableS4_fog3_rip_targets.txt',
        sep='\t', index_col=False)
    rip = set(rip['WB ID'].tolist())
    clip = set(pandas.read_csv(clip_fname, sep='\t', index_col=False
        )['gene_id'].tolist())
    cims = set(pk.df['Wormbase ID'].tolist())
    v = venn3([rip, clip, cims], ('RIP-chip', 'iCLIP peaks', 'CIMS'))
    plt.savefig(fname, format='pdf')


def cims_peaks(fname, clip_fname):
    pk = peaksList(name='CIMS')
    pk.read_csv(fname)
    pk.read_sp_vs_oo()
    pk.annotate_sp_vs_oo()
    pk.read_sp_vs_oo_as_programs()
    pk.add_rip()
    pk.add_permutation_peaks(clip_fname)
    return pk

if __name__ == '__main__':
    clip_fname = 'permutation_peaks/combined_exp.txt'
    if not os.path.isfile(sys.argv[1]):
        print "python annotate_cims.py INPUT_CIMS_PEAKS_FILE"
        sys.exit()
    if not os.path.exists('tables/'): os.system('mkdir tables/')
    pk = cims_peaks(sys.argv[1], clip_fname)
    pk.to_csv('tables/cims_annotated.txt')
    make_venn(pk)
    # Make tables.
    pk.df[pk.df['Is RIP target?']].to_csv('tables/cims_rip.txt')
    pk.df[pk.df['Is CLIP target?']].to_csv('tables/cims_clip.txt')
    filt = [(x['Is RIP target?'] and x['Is CLIP target?']) \
            for x in pk.df.to_dict('records')]
    pk.df[filt].to_csv('tables/cims_rip_and_clip.txt')
    
