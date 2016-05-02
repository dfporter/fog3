import HTSeq
import pandas
import os
import sys
import collections
import time

from flocs import *

class gtf(object):

    def __init__(self, filename, name='Unnamed', debug=False):
        self.name = name
        self.order = [
            '0', '1', '2', '3', '4', '5', '6', '7', '8',
            'gene_name', 'transcript_id', 'transcript_name',
            'exon_number', 'gene_id', 'biotype']
        if filename is not None:
            if os.path.exists(filename):
                self.read_csv(filename)

    def read_csv(self, in_filename):
        self.utrs = {}
        self.exon = {}
        self.cds = {}
        self.other = {}
        lt = pandas.read_csv(in_filename, sep='\t')
        ltr = lt.to_dict('records')
        to_line = collections.defaultdict(dict)
        for gene, tx, row in zip(lt['gene_name'].tolist(),
                     lt['transcript_id'].tolist(), ltr):
            if gene not in to_line:
                to_line[gene] = collections.defaultdict(list)
            to_line[gene][tx].append(row)
        for g in to_line:
            self.utrs[g] = collections.defaultdict(list)
            self.exon[g] = collections.defaultdict(list)
            self.cds[g] = collections.defaultdict(list)
            self.other[g] = collections.defaultdict(list)
            for t in to_line[g]:
                self.exon[g][t] = [x for x in to_line[g][t] if x['2'] == 'exon']
                self.utrs[g][t] = [x for x in to_line[g][t] if x['2'] == 'UTR']
                self.cds[g][t] = [x for x in to_line[g][t] if x['2'] == 'CDS']
                self.other[g][t] = [x for x in to_line[g][t] if \
                              x['2'] not in ['exon', 'UTR', 'CDS']]
        print "load_gtf(): created dicts of genes size: utrs: %i exons: %i cds: %i other: %i" % (
            len(self.utrs), len(self.exon), len(self.cds), len(self.other))
        print "load_gtf(): total txpt numbers are %i %i %i %i" %(
            sum([len(self.utrs[gene_name] )for gene_name in self.utrs]),
            sum([len(self.exon[gene_name]) for gene_name in self.exon]),
            sum([len(self.cds[gene_name]) for gene_name in self.cds]),
            sum([len(self.other[gene_name]) for gene_name in self.other])
            )
        return (self.utrs, self.exon, self.cds, self.other)

    def as_list(self):
        return (self.utrs, self.exon, self.cds, self.other)

    def flip_minus_strand_features(self, chr_len):
        for _feat in [self.utrs, self.exon, self.cds]:
            for gene_name in _feat:
                for txpt_id in _feat[gene_name]:
                    try:
                        if _feat[gene_name][txpt_id][0]['6'] == '+':
                            continue
                        chrm = _feat[gene_name][txpt_id][0]['0']
                    except: continue
                    for _row in _feat[gene_name][txpt_id]:
                        old_right = int(_row['4'])
                        _row['4'] = chr_len[chrm] - (int(_row['3']) - 0) + 1
                        _row['3'] = chr_len[chrm] - (old_right - 0)

    def flip_minus_strand_peak_features(self, peaks, chr_len):
        _ivs = zip(peaks.chrm, peaks.left, peaks.right, peaks.strand)
        def left(_iv):
            if _iv[3] == '+': return _iv[1]
            return chr_len[_iv[0]] - _iv[2]
        def right(_iv):
            if _iv[3] == '+': return _iv[2]
            return chr_len[_iv[0]] - _iv[1] + 1
        peaks['left'] = [left(_iv) for _iv in _ivs]
        peaks['right'] = [right(_iv) for _iv in _ivs]
#        for index, row in peaks.iterrows():
#            if row['strand'] == '+': continue
#            _new_right = chr_len[row['chrm']] - row['left'] + 1
#            peaks.loc[index, 'left'] = chr_len[row['chrm']] - peaks.loc[
#                index, 'right']
#            peaks.loc[index, 'right'] = _new_right
        return peaks

    def longest_txpt(self):
        self.longest_tx = collections.defaultdict(str)
        for gene in self.exon:
            tx_len = {}
            for tx in self.exon[gene]:
                tx_len[tx] = sum([
                    row['4'] - row['3'] for row in self.exon[gene][tx]])
            self.longest_tx[gene] = sorted(tx_len, key=lambda x: tx_len[x])[-1]

    def make_genes(self, sequences, rc_sequences):
        if not hasattr(self, 'longest_tx'): self.longest_txpt()
        self.flocs_map = {}
        for gene in self.cds:
            if gene not in self.longest_tx: continue
            if self.longest_tx[gene] == '': continue
            txpt = self.longest_tx[gene]
            if (txpt not in self.cds[gene]) or (
                len(self.cds[gene][txpt])==0): continue
            cds = sorted(
                self.cds[gene][txpt], key=lambda x: int(x['3']))
            left_cds = sorted(cds, key=lambda x: x['3'])[0]['3']
            right_cds = sorted(cds, key=lambda x: x['4'])[-1]['4']
            left_utr_span = [cds[0]['3'], cds[0]['3']]
            right_utr_span = [cds[-1]['4'], cds[-1]['4']]
            given_left_utr = False
            given_right_utr = False
            if (gene in self.utrs) and (txpt in self.utrs[gene]) and \
               (len(self.utrs[gene][txpt]) > 0):
                utrs = self.utrs[gene][txpt]
                left_utr_iv = sorted(utrs, key=lambda x: x['3'])
                right_utr_iv = sorted(utrs, key=lambda x: x['4'])
                if int(left_utr_iv[0]['3']) < left_cds:
                    left_utr_span = [left_utr_iv[0]['3'], left_cds]
                    given_left_utr = True
                if int(right_utr_iv[-1]['4']) > right_cds:
                    right_utr_span = [right_cds, right_utr_iv[-1]['4']]
                    given_right_utr = True
            strand, chrom = (cds[0]['6'], cds[0]['0'])
            exon_dict, exon_borders_in_seq, exon_num = ({}, {}, 1)
            # Start with the left utr region, if it exists.
            left_utr_iv = [chrom, left_utr_span[0]-1, left_utr_span[1]-1, '+']
            _seq = sequences if strand == '+' else rc_sequences
            seq = seq_from_iv(left_utr_iv, _seq)
            exon_dict[exon_num] = left_utr_iv
            exon_borders_in_seq[exon_num] = [0, len(seq)]
            for ex in cds:
                slice_coord = [ex['0'], ex['3'] - 1, ex['4'], '+']
                init_seq_len = len(seq)
                exon_dict[exon_num] = slice_coord
                seq += seq_from_iv(slice_coord, _seq)
                exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
                exon_num += 1
            # Add right utr.
            init_seq_len = len(seq)
            right_utr_iv = [chrom, right_utr_span[0]-1, right_utr_span[1], '+']
            seq += seq_from_iv(right_utr_iv, _seq)
            exon_dict[exon_num] = right_utr_iv
            exon_borders_in_seq[exon_num] = [init_seq_len, len(seq)]
            txpt_span = [left_utr_span[0]-1, right_utr_span[1]]
            if given_left_utr:
                txpt_span[0] = left_utr_span[0]
            if given_right_utr:
                txpt_span[1] = right_utr_span[1]
            seq = extend_to_stop(
                [left_cds, right_cds], txpt_span, seq, cds, strand,
                sequences, rc_sequences, exon_borders_in_seq, gene)
            self.flocs_map[gene] = flocs(
                txpt, chrom, txpt_span[0], txpt_span[1],
                left_cds, right_cds, strand, seq, exon_dict, gene,
                exon_borders_in_seq, given_left_utr, given_right_utr)
            self.flocs_map[gene].given_left_utr = given_left_utr
            self.flocs_map[gene].given_right_utr = given_right_utr


def extend_to_stop(cds_span, txpt_span, seq, sorted_cds, strand, sequences,
                   rc_sequences,
                   exon_borders_in_seq, gene_name):
    verbose=True
    if (cds_span[1] == txpt_span[1]) and (seq[-3:] not in ['TAA', 'TAG', 'TGA']):
        #if verbose: print "Extending %s to get stop..." % gene_name
        _add = seq_from_iv(
            (sorted_cds[0]['0'], cds_span[1], cds_span[1] + 3, '+'),
            sequences if strand == '+' else rc_sequences)
        if _add in ['TAA', 'TAG', 'TGA']:
            txpt_span[1] += 3
            seq = seq + _add
            #if verbose: print "Did it work? Expected stop: %s" % seq[-3:]
            try:
                exon_borders_in_seq[len(exon_borders_in_seq)][1] = len(seq)
            except:
                loggern.warn(
                    "%s: Missing an exon? %s" % (gene_name, str(exon_borders_in_seq)))
        else:
            if verbose: print "%s: Did it work? Expected stop: %s" % (
                gene_name, _add)
    return seq


def seq_from_iv(iv, sequences):
    """Returns a slice from 0-based, slice coordinates.
    """
    if iv[2] > len(sequences[iv[0]]):
        end = len(sequences[iv[0]])
    else:
        end = iv[2]
    a_seq = sequences[iv[0]][iv[1]:end]
    if iv[3] == "-":
        a_seq = rc(a_seq)
    return a_seq


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join([basecomplement[base] for base in list(s)])


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def get_sequences():
    fasta_filename = '/scratch/indexes/WS235.fa'
    sequences = dict(
        (p.name.split(' ')[0], p.seq) for p in HTSeq.FastaReader(fasta_filename))
    rc_sequences = dict(
        (p.name.split(' ')[0], rc(p.seq)) for p in HTSeq.FastaReader(fasta_filename))
    chr_lens = dict(
        [(name, len(sequences[name])) for name in sequences])
    return (sequences, rc_sequences, chr_lens)

def _time(start):
    return time.time() - start

if __name__ == '__main__':
    print "Debug mode for gtf.py."
    print sys.version
    start = time.time()
    seqs, rc_seqs, chr_lens = get_sequences()
    print "get_sequneces(): {0} sec".format(_time(start))
    start = time.time()
    g = gtf('/opt/lib/gtf_with_names_column.txt', debug=True)
    print "gtf init(): {0} sec".format(_time(start))
    start = time.time()
    g.flip_minus_strand_features(chr_lens)
    g.longest_txpt()
    print "gtf flip strand and longest txpt(): {0} sec".format(_time(start))
    g.make_genes(seqs, rc_seqs)
