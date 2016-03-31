__author__ = 'dp'
import p_values_of_clusters
import random
import collections
import numpy as np
import os
import HTSeq
import sys


def is_good_gene(rows):
    # print rows
    for row in rows:
        if row['biotype'] == 'rRNA':
            return False
        if row['biotype'] == 'snoRNA':
            return False
        # if row['biotype'] == 'protein_coding':
        #     return False
        if row['biotype'] != 'protein_coding':
        # # if row['gene_id'] != 'WBGene00010957':
            return False
    return True


def run(reads_by_gene_by_fname, counts_by_gene, exons, lib=None):
    read_table = collections.defaultdict(list)
    exons = dict((k, v) for k, v in exons.iteritems() if is_good_gene(v))
    gene_set = set(exons.keys())
    exon_ga = HTSeq.GenomicArray('auto', stranded=True)
    for fname in reads_by_gene_by_fname:
        print "*" * 7 + '\nRandomizing reads in ' + fname
        reads_by_gene = reads_by_gene_by_fname[fname]
        obs_gene_list = set(reads_by_gene.keys()) & set(exons.keys())
        for obs_gene in reads_by_gene:
            rig = list(reads_by_gene[obs_gene])
            read_lengths = [x.end - x.start for x in rig]
            rand_gene = random.choice(list(gene_set))
            # gene_set -= set([rand_gene])
            # print rand_gene
            # print rig
            txpt_id, txpt_exons = p_values_of_clusters.find_longest_txpt(
                exons[rand_gene])
            # print txpt_exons[0]
            txpt_exons = sorted(txpt_exons, key=lambda x: x[0])
            for exnum, ex in enumerate(txpt_exons, start=1):
                # print exnum, ex
                exon_ga[HTSeq.GenomicInterval(
                    exons[rand_gene][0]['0'], ex[0], ex[1] + 1, exons[rand_gene][0]['6']
                )] += exnum
            reads = randomize_reads(txpt_exons, read_lengths, exons[rand_gene])
            read_table[fname].extend([
                [exons[rand_gene][0]['0'], str(x[0]), str(x[1]), '.',
                 '100.00', exons[rand_gene][0]['6']
                 ] for x in reads
            ])
    if not os.path.exists('random_bed/'): os.system('mkdir random_bed/')
    for fname in read_table:
        with open('random_bed/' + os.path.basename(fname), 'w') as f:
            li = ''
            for read in read_table[fname]:
                li += '\t'.join(read) + '\n'
            f.write(li)
    import bed_to_wig
    bed_to_wig.run('random_bed/', 'random_bedgraph/')
    exon_ga.write_bedgraph_file('exons_+.wig', strand='+')
    exon_ga.write_bedgraph_file('exons_-.wig', strand='-')
    import assign_to_genes
    assign_to_genes.run(use_this_lib=lib, bed='random_bed/')


def randomize_reads(exons, read_lengths, exons_of_gene):
    exons_of_gene = sorted(exons_of_gene, key=lambda x: x['3'])
    total_len_of_ranges = 0
    range_starts = {}
    longest_txpt, _ = p_values_of_clusters.find_longest_txpt(exons_of_gene)
    exons_of_longest_txpt = [
        x for x in exons_of_gene if x['transcript_id']==longest_txpt]
    # print exons_of_longest_txpt
    ranges = sorted(exons, key=lambda x: x[0])
    for _range in ranges:
        for pos in range(
                total_len_of_ranges,
                total_len_of_ranges + _range[1] - _range[0] + 1, 1):
            range_starts[pos] = _range[0] + pos - total_len_of_ranges
        total_len_of_ranges += _range[1] - _range[0]
        # print "pos {i} > {g}. total len = {t}.".format(
        #     i=total_len_of_ranges, g=range_starts[total_len_of_ranges],
        #     t=total_len_of_ranges)
    if min(range_starts.values()) + 50 < min([x['3'] for x in exons_of_longest_txpt]):
        print "Too far left?"
        sys.stdin.readline()
    if max(range_starts.values()) - 50 > max([x['4'] for x in exons_of_longest_txpt]):
        print "Too far right?"
        sys.stdin.readline()
    reads = []
    rand_starts = np.random.random_integers(
        0, total_len_of_ranges, size=len(read_lengths))
    for i, n in enumerate(rand_starts):
        reads.append([range_starts[n], range_starts[n] + read_lengths[i]])
        # reads.append([n, n + np.random.random_integers(15, high=30)])
        # print "rand start {o} -> range start {v}  = {s}".format(
        #     o=n, v=range_starts[n],  s=reads[-1])
    return reads

