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
import random


class gene(object):

    def __init__(self, name, language, value=None, origin=None):
        self.name = name
        self.language = language
        self.names = collections.defaultdict(str)
        self.names[language] = name
        self.vals = {}
        if (value is not None) and (origin is not None):
            self.vals[origin] = value

    def add_name(self, name, language):
        self.names[language] = name

    def __str__(self):
        return "{k}:\n{a}, values:{v}\n".format(
            k=self.name, a=self.names, v=self.vals)

    def add_value(self, origin, val):
        self.vals[origin] = val

    def as_dict(self):
        d = {'Name': self.name, 'Language': self.language}
        d.update(self.vals)
        return d

    def apply_rpkm(self, dataset_sizes, gene_lengths):
        for col in [x for x in self.vals if x in dataset_sizes]:
            self.vals[col] = float(1e6 * self.vals[col])/float(
                dataset_sizes[col])


class geneList(object):

    def __init__(self, name='Unnamed', language='Unknown'):
        self.name = name
        self.language = language
        self.genes = {}

    def __str__(self):
        li = '{a}, {b}\n'.format(a=self.name, b=self.language)
        li += ''.join(str(x) for x in self.genes.values())
        return li

    def add_gene(self, gene):
        if gene.language != self.language:
            print "Expected language %s, found %s. Input: %s" % (
                gene.language, self.language, gene)
        if gene.name in self.genes:
            return True
        self.genes[gene.name] = gene

    def add_vals_from_dict(self, adict):
        for gene in adict:
            if gene in self.genes:
                self.genes[gene] = adict[gene]

    def as_list(self):
        return [gene.as_dict() for gene in self.genes.values()]

    def as_dataframe(self):
        return pandas.DataFrame(self.as_list())
            
    def read_peaks_file(
        self, fname, col_to_use_as_name='gene_id',
        col_to_use_as_value='height', language='wbid'):
        df = pandas.read_csv(fname, sep='\t')
        dfr = df.to_dict('records')
        for row in dfr:
            if row[col_to_use_as_name] in self.genes:
                continue
            self.add_gene(
                gene(row[col_to_use_as_name], language,
                     value=row[col_to_use_as_value], origin=col_to_use_as_value)
                )
        print "*\n" * 14

    def add_attribute(
        self, fname, name='gene_id', value_name='height',
        col_names=None):
        df = pandas.read_csv(fname, sep='\t', index_col=False,)
        # Deal with the case of default column names:
        if type(value_name) == type([]) and col_names is None:
            if len(value_name) > 0:
                col_names = [None] * len(value_name)
            else:
                col_names = [None] * len([x for x in df.columns if x != name])
        print df.head()
        # Add attributes.
        if type(value_name)==type([]):
            for (col_name, attr) in zip(col_names, value_name):
                self.add_an_attribute(df, name, attr, col_name=col_name)
            if len(value_name) == 0:
                for col_name, attr in \
                    zip(col_names, [x for x in df.columns if x != name]):
                    print attr
                    self.add_an_attribute(df, name, attr, col_name=col_name)
        elif type(value) == type(''):
            self.add_an_attribute(df, name, value, col_name=col_name)

    def add_an_attribute(self, df, name_col, value_col, col_name=None):
        name_to_val = dict(zip(df[name_col].tolist(), df[value_col].tolist()))
        if col_name is None:
            col_name = value_col.partition('.')[0]
            col_name = re.sub('_counts', '', col_name)
        print "&&&" * 20
        if random.randrange(2) == 1:
            print "add_an_attribute():"
            print name_col, value_col
            print str(name_to_val)[:100]
            print str(self.genes)[:100]
        for gene_name in self.genes:
            if gene_name in name_to_val:
                self.genes[gene_name].add_value(col_name, name_to_val[gene_name])
            else:
                self.genes[gene_name].add_value(col_name, 0)
        
    def convert_to_rpkm(self, gtf='', files=[]):
        lens = get_gene_len(
            self.genes.keys(), gtf, use='gene_id')
        print files
        print "***"
        dataset_sizes = dict([(os.path.basename(x).partition('.')[0],
                               len(open(x).readlines())) for x in files])
        print dataset_sizes
        print str(lens)[:99]
        for name in self.genes:
            self.genes[name].apply_rpkm(dataset_sizes, lens)


def read_as_table_of_lists(fname, use_header_val=None):
    print "looking for header value %s" % use_header_val
    if use_header_val is not None:
        _table = collections.defaultdict(list)
        df = pandas.read_csv(fname, sep='\t')
        dfd = df.to_dict('records')
        for row in dfd:
            _table[row[use_header_val]].append(row)
        return _table


def make_gtf_of_longest_txpt_per_gene(
    filename, output_filename='/opt/lib/gtf_one_txtp_per_gene.txt'):
    gtf = read_as_table_of_lists(filename, use_header_val='gene_name')
    header = open(filename).readline()
    outrows = header
    header_list = header.rstrip('\n').split('\t')
    for gene_name in gtf:
        bytxpt = collections.defaultdict(list)
        for dict_row in gtf[gene_name]:
            bytxpt[dict_row['transcript_id']].append(dict_row)
        txptlen = {}
        for _id in bytxpt:
            txptlen[_id] = sum([int(x['4']) - int(x['3']) for \
                                x in bytxpt[_id] if x['2']=='exon'])
        longest_txpt = sorted(txptlen.keys(), key=lambda x: txptlen[x])[-1]
#        print bytxpt
#        print longest_txpt
#        print bytxpt[longest_txpt]
        rowlist = ["\t".join([str(x[y]) for y in header_list]) for x in bytxpt[longest_txpt]]
#        print rowlist
        outrows += "\n".join(rowlist)
        outrows += '\n'
    with open(output_filename, 'w') as f:
        f.write(outrows)


def get_gene_len(genes, filename, use='gene_id'):
    gtf = read_as_table_of_lists(filename, use_header_val=use)
    lens = {}
    print str(gtf)[:500]
    for k in gtf:
        print "K:{k}".format(k=k)
        for row in gtf[k]:
            print row.items()
        break
    found, missing = set(), set()
    for gene in genes:
        if gene not in gtf:
            lens[gene] = 1e3
            print 'Missing ' + gene
            missing.add(gene)
            if random.randrange(100) == 1: sys.exit()
            continue
        found.add(gene)
        txpts = set([x['transcript_id'] for x in gtf[gene]]) - set([np.nan]) - set(['.'])
        if len(txpts) > 1:
            print "more than one txpt"
            print txpts
        try:
            this_txpt = list(txpts)[0]  # Expect only one txpt, but just in case.
        except:
            print "get_gene_len() error 119"
            print gene
            print txpts
            print gtf[gene]
        lens[gene] = sum([
            int(row['4']) - int(row['3']) for row in gtf[gene] if (
                (row['transcript_id'] == this_txpt) and (row['2'] == 'exon')
            )
        ])
        if lens[gene] == 0:
            print "No gene len?"
            print this_txpt
            print set([x['transcript_id'] for x in gtf[gene]])
            print [x['2'] for x in gtf[gene]]
    print "Found {a}/{b} gene ids in {name}. Failed to find {c}.".format(
        a=len(found), b=len(genes), name=filename, c=len(missing))
    return lens

