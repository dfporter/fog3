import pandas
import os
import config
import re
import glob
import argparse

class stats(object):

    def __init__(self, name='Unnamed'):
        self.name = name
        self.stats = {}

    def fastq_files(self, dirname):
        fsizes = [(os.path.basename(x).partition('.fastq')[0],
               len(open(x).readlines())) for x in \
              glob.glob(dirname + '/*fastq')]
        self.stats['Number of reads'] = dict([(x[0], x[1]/4) for x in fsizes])

    def uniquely_mapping(self, dirname):
        self.stats['Number of reads mapping uniquely at 20 AS'] = \
                dict(self.sam_sizes_in_dir(glob.glob(dirname + '/*sam')))
        self.as_perc('Number of reads mapping uniquely at 20 AS',
                'Number of reads')

    def as_perc(self, numer, denom):
        self.stats[numer + ' (%)'] = {}
        for j in self.stats[numer]:
            if self.stats[denom][j] == 0:
                self.stats[numer + ' (%)'][j] = 0
                continue
            self.stats[numer + ' (%)'][j] = \
                100*float(self.stats[numer][j])/float(self.stats[denom][j])

    def map_to_two_places(self, dirname):
        self.stats['Number of alignments mapping to exactly two places'] = \
                dict(self.sam_sizes_in_dir(glob.glob(dirname + '/*sam')))
        self.as_perc('Number of alignments mapping to exactly two places',
                'Number of reads')

    def total_reads_in_sam_files(self, dirname):
        self.stats['Number of alignments mapped or unmapped by STAR'] = \
                dict(self.sam_sizes_in_dir(glob.glob(dirname + '/*sam')))

    def unmapped(self, dirname):
        self.stats['Number of reads unmapped by STAR'] = \
                dict(self.sam_sizes_in_dir(glob.glob(dirname + '/*sam')))
        self.as_perc('Number of reads unmapped by STAR',
                'Number of reads')

    def uncollapsed_reads(self, dirname):
        self.stats['Number of reads mapping uniquely before collapsing'] = \
                dict(self.bed_sizes_in_dir(glob.glob(dirname + '/*bed')))
        self.as_perc('Number of reads mapping uniquely before collapsing',
                'Number of reads')

    def collapsed_reads(self, dirname):
        self.stats['Number of reads mapping uniquely after collapsing duplicates'] = \
                dict(self.bed_sizes_in_dir(glob.glob(dirname + '/*bed')))
        self.as_perc('Number of reads mapping uniquely after collapsing duplicates',
                'Number of reads')

    def collapsed_reads_rrna(self, dirname):
        self.stats['Number of rRNA reads mapping uniquely after collapsing duplicates'] = \
                dict(self.bed_sizes_in_dir(glob.glob(dirname + '/*bed')))
        self.stats['rrna'] = self.stats['Number of rRNA reads mapping uniquely after collapsing duplicates']
        if 'non_rrna' in self.stats:
            self.set_total_bed()

    def set_total_bed(self):
        self.stats['total bed'] = {}
        for k in self.stats['rrna']:
            self.stats['total bed'][k] = \
                self.stats['rrna'][k] + self.stats['non_rrna'][k]
        self.stats['Number of reads mapping uniquely after collapsing duplicates'] = self.stats['total bed']
        self.as_perc('Number of reads mapping uniquely after collapsing duplicates',
                'Number of reads')

    def collapsed_reads_no_rrna(self, dirname):
        self.stats['Number of non-rRNA reads mapping uniquely after collapsing duplicates'] = \
                dict(self.bed_sizes_in_dir(glob.glob(dirname + '/*bed')))
        self.stats['non_rrna'] = self.stats['Number of non-rRNA reads mapping uniquely after collapsing duplicates']
        if 'rrna' in self.stats:
            self.set_total_bed()
        self.as_perc('Number of non-rRNA reads mapping uniquely after collapsing duplicates',
                'Number of reads mapping uniquely after collapsing duplicates')

    def __str__(self):
        li = '\n==={a}===\n'.format(a=self.name)
        for k in self.stats:
            li += "{k}:\n".format(k=k)
            for j in self.stats[k]:
                li += "\t{a}\t{b:>14}\n".format(a=j, b=self.stats[k][j])
        return li
    
    def bed_sizes_in_dir(self, file_list):
        return [(os.path.basename(x).partition('.bed')[0],
                 len(open(x).readlines())) for x in file_list]

    def sam_sizes_in_dir(self, file_list):
        return [(os.path.basename(x).partition('.sam')[0],
               self.sam_size(x)) for x in file_list]
    
    def sam_size(self, fname):
        n = 0
        with open(fname) as f:
            for li in f:
                if li[0] != '@':
                    n += 1
        return n

    def as_dataframe(self):
        self.df = pandas.DataFrame(self.stats)
        print self.df
        return self.df

class statsWithPeaksByPermutations(stats):

    def read_clusters_dir(self, dirname):
        self.stats['Peaks, before FDR cutoff'] = {}
        self.stats['Peaks, after FDR cutoff'] = {}
        self.stats['Genes, before FDR cutoff'] = {}
        self.stats['Genes, after FDR cutoff'] = {}
        for fname in glob.glob(dirname + '/*txt'):
            f = os.path.basename(fname).partition('.txt')[0]
            df = pandas.read_csv(fname, sep='\t')
            self.stats['Peaks, before FDR cutoff'][f] = len(df.index)
            self.stats['Genes, before FDR cutoff'][f] = len(
                set(df['gene_id'].tolist()))
            df = df[df['past_fdr']]
            self.stats['Peaks, after FDR cutoff'][f] = len(df.index)
            self.stats['Genes, after FDR cutoff'][f] = len(
                set(df['gene_id'].tolist()))

    def read_peaks_dir(self, dirname):
        print 'asdfasdf'
        self.stats['Peaks, final'] = {}
        self.stats['Genes, final'] = {}
        for fname in glob.glob(dirname + '/*txt'):
            f = os.path.basename(fname).partition('.txt')[0]
            df = pandas.read_csv(fname, sep='\t')
            self.stats['Peaks, final'][f] = len(df.index)
            self.stats['Genes, final'][f] = len(
                set(df['gene_id'].tolist()))
            print self.stats
    
def run(args):
    top = 'v2/'
    lib = config.config(args.config_ini)
    exp_bed = [lib[x] for x in lib if re.search('\Aexp_bed', x)]
    control_bed = [lib[x] for x in lib if re.search('\Acontrol_bed', x)]
    s = statsWithPeaksByPermutations(name='FOG-3')
    s.read_clusters_dir('./clusters/')
    s.read_peaks_dir('./permutation_peaks/')
    s.fastq_files(os.path.join(top, 'temp_clipped'))
    print s
    s.uniquely_mapping(os.path.join(top, 'uniquely_mapping'))
    print s
    s.map_to_two_places(os.path.join(top, 'map_to_two'))
    s.unmapped(os.path.join(top, 'unmapped'))
    s.as_dataframe()
    s.uncollapsed_reads(os.path.join(top, 'bed_uncollapsed3'))
 #   s.collapsed_reads(os.path.join(top, 'bed_collapsed/'))
    s.collapsed_reads_rrna(os.path.join(top, 'bed_collapsed/rrna'))
    s.collapsed_reads_no_rrna(os.path.join(top, 'bed_collapsed/no_rrna/'))
    s.total_reads_in_sam_files(os.path.join(top, 'unfiltered_star_sams_output'))
    print s
    s.as_dataframe().to_csv('stats.txt', sep='\t', index=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''File stats.''')
    parser.add_argument('-c', '--config_ini',
                        default='config.ini',
                        help='''config.ini file.''')
    args = parser.parse_args()
    run(args)
