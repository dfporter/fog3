import pandas
import glob
import os
import sys

in_dir = sys.argv[1]
for f in glob.glob(in_dir + '/*'):
    df = pandas.read_csv(f, sep='\t')
    total_clusters = len(df.index)
    sig_clusters = df[df['padj']<0.05]
    li = '{name:>20}\t{total_clusters:>10}\t{sig_clusters:>20}\t{n_targs}'.format(
        name=os.path.basename(f), total_clusters=total_clusters,
        sig_clusters=len(sig_clusters.index),
        n_targs=len(set(df['gene_id'].tolist()))
    )
    print li