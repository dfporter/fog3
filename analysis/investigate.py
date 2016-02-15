import pandas
import sys
import os
import glob
import re
import collections


def investigate(fname):
    df = pandas.read_csv(fname, sep='\t')
    biotypes = df['biotype'].value_counts()
    if 'gene_name' in df.columns:
        n_targets = set(df['gene_name'])
    elif 'gene_id' in df.columns:
        n_targets = set(df['gene_id'])
    else:
        print 'No gene_name or gene_id column...'
        n_targets = set([])
    print "\n{f}\nPeaks {n}. genes {g}. biotypes:\n{b}.".format(
        f=fname,
        n=len(df.index), g=len(n_targets), b=biotypes
    )
def make_table_of_replicate_stats(indir, lib):
    data = collections.defaultdict(dict)
    for fname in glob.glob(lib['read_beds'] + '/*.bed'):
        bfname = os.path.basename(fname).partition('.bed')[0]
        # data[bfname]['num_fastq_reads'] = 0
        data[bfname]['num_collapsed_bed_reads'] = sum([1 for li in open(fname, 'r')])
    for fname in glob.glob('../pre_calls/trimmed_cliped/*fastq'):
        bfname = os.path.basename(fname).partition('.fastq')[0]
        data[bfname]['num_fastq_reads'] = sum([1 for li in open(fname, 'r')])
        data[bfname]['num_fastq_reads'] /= 4
    for fname in glob.glob(lib['clusters_dir'] + '/*'):
        bfname = os.path.basename(fname).partition('.txt')[0]
        data[bfname]['clusters'] = sum([1 for li in open(fname, 'r')]) - 1
        df = pandas.read_csv(fname, sep='\t')
        df = df[df['padj']<0.01]
        data[bfname]['significant clusters'] = len(df.index)
    wd(data)
def wd(d):
    print '*' * 14
    for k in d:
        print k
        for k2 in d[k]:
            print "   {a:30}: {b:>10}".format(
                a=k2, b=d[k][k2]
            )
    cols = set()
    for k in d:
        cols |= set(d[k].keys())
    exp_order = [x for x in d.keys() if not (cols - set(d[x].keys()))]
    print 'exp order'
    print exp_order
    li = 'Category\t' + '\t'.join(exp_order) + '\n'
    for col in cols:
        li += col + '\t' + '\t'.join([str(d[x][col]) for x in exp_order]) + '\n'
    print li


sys.path.insert(0, 'analysis/')
import config
lib = config.config()
make_table_of_replicate_stats('', lib)

def investigate_file_list():
    fname1 = sys.argv[1]
    print "**\n" * 2
    i = 1
    while True:
        if len(sys.argv) > i:
            fname = sys.argv[i]
            investigate(fname)
        else: break
        i += 1
