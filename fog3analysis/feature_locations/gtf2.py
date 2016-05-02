import os
import pandas

class gtf(object):

    def __init__(self, filename, name='Unnamed'):
        self.name = name
        if filename is not None:
            if os.path.exists(filename):
                self.read_csv(filename)

    def read_csv(self, in_filename):
        utrs = {}
        exon = {}
        cds = {}
        other = {}
        lt = pandas.read_csv(in_filename, sep='\t')
        ltr = lt.to_dict('records')
        txpt_to_gene = dict(zip(lt['transcript_id'].tolist9),
                            lt['gene_name'].tolist())
        to_line = collections.defaultdict(dict)
        for t in zip(lt['gene_name'].tolist(),
                     lt['transcript_id'].tolist(), ltr):
            to_line[t[0]] = collections.defaultdict(list)
            to_line[t[0]][t[1]].append(t[1])
        for g in to_line:
            utrs[g] = collections.defaultdict(list)
            exon[g] = collections.defaultdict(list)
            cds[g] = collections.defaultdict(list)
            for t in to_line[g]:
                exon[g][t] = [x for x in txpt_to_line[t] if x['2'] == 'exon']
                utrs[g][t] = [x for x in txpt_to_line[t] if x['2'] == 'UTR']
                cds[g][t] = [x for x in txpt_to_line[t] if x['2'] == 'CDS']
                other[g][t] = [x for x in txpt_to_line[t] if \
                              x['2'] not in ['exon', 'UTR', 'CDS']]
        self.utrs = utrs
        self.exon = exon
        self.cds = cds
        self.other = other
        print "load_gtf(): created dicts of genes size: utrs: %i exons: %i cds: %i other: %i" % (
            len(utrs), len(exon), len(cds), len(other))
        print "load_gtf(): total txpt numbers are %i %i %i %i" %(
            sum([len(utrs[gene_name] )for gene_name in utrs]),
            sum([len(exon[gene_name]) for gene_name in exon]),
            sum([len(cds[gene_name]) for gene_name in cds]),
            sum([len(other[gene_name]) for gene_name in other])
            )
        return (utrs, exon, cds, other)

    def as_list(self):
        return (self.utrs, self.exon, self.cds, self.other)
        no = """
        self.data_list = convert_to_list(self.data)
        self.utrsl = [x for x in self.data_list if x['2']=='UTR']
        self.exonl = [x for x in self.data_list if x['2']=='exon']
        self.cdsl = [x for x in self.data_list if x['2']=='cds']
        self.utrs = collections.defaultdict(dict)
        self.exon = collections.defaultdict(dict)
        self.cds = collections.defaultdict(dict)
        for _feat, d_feat in [
            (self.utrsl, self.utrs), (self.exonl, self.exon),
            (self.cdsl, self.cds)]:
            for row in _feat:
                if row['gene_name'] in d_feat:
                    if row['transcript_id'] in d_feat[row['gene_name']]:
                        d_feat[row['gene_name']][row['transcript_id']].append(row)
                    else:
                        d_feat[row['gene_name']][row['transcript_id']] = [row]
                else:
                    d_feat[row['gene_name']] = {row['transcript_id']:[row]}

    def convert_to_list(self, df):
        _list = df.to_dict('records')
        _list = [dict([(t, x[t]) for t in self.order]) for x in _list]
        return _list"""
