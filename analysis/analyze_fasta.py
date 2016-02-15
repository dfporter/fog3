import os
import sys
import pandas
import collections
import re
import random
import numpy as np

def analyze_fasta(fasta_filename, motif=None):
    motif_list = ['ctcac', 'ctctct', 'ctca', 'atcc', 'ccta']
    motif_per = collections.defaultdict(list)
    motif_binary = collections.defaultdict(int)
    n_seqs = 0
    counts = collections.defaultdict(dict)
    num_shuffles = 1000
    rand_counts = {}
    for motif in motif_list:
        counts[motif]['per'] = []
        for index in range(0, num_shuffles):
            rand_counts.setdefault(index, {})
            rand_counts[index][motif] = {'per': []}
    with open(fasta_filename) as f:
        while True:
            _id = f.readline()
            seq = f.readline().rstrip('\n')
            if not seq: break
            seq = re.sub('u', 't', seq.lower())
            n_seqs += 1
            for motif in motif_list:
                motif_regex = re.compile(motif)
                counts[motif]['per'].append(0)
                for m in motif_regex.finditer(seq):
                    counts[motif]['per'][-1] += 1
                for index in range(0, num_shuffles):
                    seqarr = list(seq)
                    random.shuffle(seqarr)
                    rand_seq = ''.join(seqarr)
                    rand_counts[index][motif]['per'].append(0)
                    for m2 in re.finditer(motif_regex, rand_seq):
                        rand_counts[index][motif]['per'][-1] += 1
    for pat in counts:
        counts[pat] = add_frac_with(counts[pat])
        print pat
        n_experimental = sum(counts[pat]['with'])
        for index in rand_counts:
            rand_counts[index][pat] = add_frac_with(rand_counts[index][pat])
        mean_random_with = np.mean([rand_counts[x][pat]['frac_with'] for x in rand_counts])
        mean_random_per = np.mean([
                np.mean(rand_counts[x][pat]['per']) for x in rand_counts
            ])
        n_above = [sum(rand_counts[x][pat]['with']) for x in rand_counts]
        n_above = [1 for x in n_above if x>=n_experimental]
        p_value = float(sum(n_above))/float(len(rand_counts))
        print """{ai}/{aii} sequences have the motif {pat} ({asdf}). The\
 average sequence has {qq} motif instances.
{i}/{ii} ({mrw}) of shuffled sequences have {a} or more sequences with the motif\
 (p value {ff}). The average shuffled sequence has {conqq} instances.
 The average enrichment is {enr}.""".format(
            ai=n_experimental, aii=len(counts[pat]['with']),
            pat=pat, asdf=counts[pat]['frac_with'],
            qq=np.mean(counts[pat]['per']),
            i=sum(n_above), ii=len(rand_counts),
            mrw=mean_random_with, a=n_experimental,
            ff=p_value, conqq=mean_random_per,
            enr=np.mean(counts[pat]['per'])/mean_random_per
        )

def add_frac_with(_counts):
    _counts['with'] = [max(0, min(x, 1)) for x in _counts['per']]
    _counts['frac_with'] = np.mean(_counts['with'])
    return _counts


if __name__ == '__main__':
    try:
        fa_filename = sys.argv[1]
    except:
        print "run python analyze_fasta.py fasta_file"# <motif> (motif optional)"
        sys.exit()
    analyze_fasta(fa_filename)