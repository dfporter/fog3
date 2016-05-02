import sys
import numpy as np
import pandas
import collections
import glob
import HTSeq
import logging
import argparse
import os
import re

from geneList import *


def get_rpkm():
    import get_rpkm


def build_geneList(args):
    gtf_longest_path = '/opt/lib/gtf_one_txtp_per_gene.txt'
    base_gtf_path = '/opt/lib/gtf_with_names_column.txt'
    gl = geneList(name='CLIP', language='wbid')
    if not os.path.exists(gtf_longest_path) or (
        args.make_gtf_of_longest_txpt_per_gene):
        make_gtf_of_longest_txpt_per_gene(
            base_gtf_path, output_filename=gtf_longest_path)
    print gl
    gl.read_peaks_file(
        args.peaks, col_to_use_as_name='gene_id',
        col_to_use_as_value='height', language='wbid')
    print gl.as_dataframe()
    gl.add_attribute(
        args.counts, name='gene',  value_name=[])
    gl.add_attribute('/opt/lib/ortiz/TableS1_oogenic.txt',
                    name='Wormbase ID',
                    value_name=['Ortiz Gonad', 'Program'],
                    col_names=['Oo. file Ortiz Gonad', 'Oo. file Program'])
    gl.add_attribute('/opt/lib/ortiz/TableS2_spermatogenic.txt',
                    name='Wormbase ID',
                    value_name=['Ortiz Gonad', 'Program'],
                    col_names=['Sp. file Ortiz Gonad', 'Sp. file Program'])
    gl.add_attribute('/opt/lib/ortiz/DESeq_genes_in_gonad.txt',
                     name='WormBase ID (WS240)',
                     value_name=['Expression in  fog-2 gonads (RPKM)',
                                 'Expression in fem-3 gonads (RPKM)',
                                 'log2 fold change (of normalized reads)'],
                     col_names=['Oo. gl RPKM', 'Sp. gl RPKM',
                                'log2 SP/OO gl RPKM']
    )
    print ")()()("
    gl.convert_to_rpkm(
        gtf=gtf_longest_path,
        files=glob.glob(args.bed + '/*.bed'))
    #print gl
    gl.as_dataframe().to_csv(args.peaks, sep='\t', index=False)
    print "Finished."


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks',
        default='../permutation_peaks/combined_exp.txt')
    parser.add_argument('-u', '--counts',
        default='counts/combined_counts.txt')
    parser.add_argument('-b', '--bed',
        default='../bed_collapsed/',
        help='Folder of bed files of all reads.')
    parser.add_argument('-g', '--make_gtf_of_longest_txpt_per_gene',
                        default=False,
                        action='store_true')
    args = parser.parse_args()
    build_geneList(args)
