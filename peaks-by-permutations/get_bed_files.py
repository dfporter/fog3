import HTSeq
import collections
import re
import cluster_reads
import time
import os

def get_bed_files(bed_folder='beds/', gtf=None, args=None, lib=None):
    ht_exons = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    gtfd = collections.defaultdict(set)
    for feature in gtf:
        if feature.type == 'exon':
            ht_exons[feature.iv] += feature.name
            gtfd[feature.name].add(feature)
    control_names = [lib[x] for x in lib.keys() if re.match('control_bed.*', x) is not None]
    exp_names = [lib[x] for x in lib.keys() if re.match('exp_bed.*', x) is not None]
    bed_file_list = [
        "/".join([bed_folder, x]) for x in control_names + exp_names
    ]
    print "get_bed_files:"
    print bed_file_list
    reads_by_gene = {}
    counts_by_gene = {}
    ga = {}
    for fname in bed_file_list:
        (_ga, reads_by_gene[fname], counts_by_gene[fname]) = read_a_bed_file(fname, ht_exons)
        ga[fname] = _ga
    return ga, reads_by_gene, counts_by_gene


def read_a_bed_file(fname, ht_exons):
    import time
    start_time = time.time()
    counts = collections.defaultdict(int)
    reads_by_gene = collections.defaultdict(set)
    _ga = HTSeq.GenomicArray('auto', stranded=True)
    if not os.path.exists(fname):
        print "bed file not found: {o}".format(o=fname)
        return _ga, counts, reads_by_gene
    for n, line in enumerate(open(fname)):
        # if n > 1e6: break
        if (n > 0) and (not n % 1e6):
            min_per_line = float(time.time() - start_time)/float(60 * max([1., n]))
            print "Line {n}: Elapsed: {te} s. Minutes for 10 million reads: {f}".format(
                n=n, te=float(time.time() - start_time),
                f=min_per_line * 1e7/60.
            )
        s = line.rstrip('\n').split('\t')
        read_iv = HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5])
        _ga[read_iv] += 1
        gene_ids = set()
        for iv, val in ht_exons[read_iv].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[gene_id] += 1
            reads_by_gene[gene_id].add(read_iv)
    print "Stats on {a}:".format(a=fname)
    cluster_reads.stats(reads_by_gene)
    return _ga, reads_by_gene, counts
