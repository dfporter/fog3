import collections
import numpy as np
import HTSeq
import os
import sys
import pandas
import matplotlib.pyplot as plt

class sbGene:

    def __init__(self, iv):
        self.iv = iv

    def find_subpeaks(self, coverage, exons_dict):
        self.determine_coverage_accross_exons(coverage, exons_dict)
        self.n_peaks = self.find_clusters()
        n = """
        self.peak_ivs_genomic_coord = [HTSeq.GenomicPosition(
            self.iv.chrom,
            self.position_in_exon_coverage_mapped_to_genomic_position[pos],
            self.iv.strand,
            ) for pos in self.peaks]
        self.peak_ivs_relative = [HTSeq.GenomicPosition(
            self.iv.chrom, pos, self.iv.strand,
            ) for pos in self.peaks]"""

    def filter_by_intersection(self, list_of_ivs):
        self.filtered_peak_ivs_genomic_coord = []
        for iv in self.peak_ivs_genomic_coord:
            for alt_iv in list_of_ivs:
                if alt_iv[3] != iv.strand: continue
                if alt_iv[0] != iv.chrom: continue
                if alt_iv[1] <= iv.pos <= alt_iv[2]:
                    self.filtered_peak_ivs_genomic_coord.append(iv)
        return self.filtered_peak_ivs_genomic_coord

    def add_to_wig(self, ga):
        for iv in self.peak_ivs_genomic_coord:
            ga[iv] += 1

    def find_clusters(self):
        try:
            n_clusters = len(
                filter(
                    lambda x: np.nanmax(x)>1.01,
                    self.consecutive(np.nonzero(self.exon_coverage)[0])
                    )
                )
            arr = self.consecutive(np.nonzero(self.exon_coverage)[0])
            arr = np.array([self.exon_coverage[i] for i in arr])
            after = filter(lambda x: np.nanmax(x)>1, arr)
            n_clusters = len(after)
        except:
            print "...Failure"
            n_clusters = 0
        return n_clusters

    def consecutive(self, data):
        return np.split(data, np.where(np.diff(data) != 1)[0]+1)

    def determine_coverage_accross_exons(
            self, coverage, exons_dict):
        '''
        exon    0      1       2
        genomic a      b       c
        rela    012345601234567012345
        '''
        self.exon_coverage = np.array([])
        self.position_in_exon_coverage_mapped_to_genomic_position = {}
        first_exon_left = exons_dict[sorted(exons_dict.keys(), key=lambda x: int(x))[0]][1]
        strand = '+'
        reverse = False
        if exons_dict.values()[0][-1] == '-':
            strand = '-'
            reverse = True
        for exon_num in sorted(
                exons_dict.keys(), key=lambda x: int(x), reverse=reverse):
            (chrm, left, right, strand) = exons_dict[exon_num]
            iv = HTSeq.GenomicInterval(chrm, left, right, strand)
            this_exon_coverage = np.fromiter(coverage[iv], dtype='f')
            for pos in range(0, len(this_exon_coverage)):
                self.position_in_exon_coverage_mapped_to_genomic_position[pos] \
                        = left + pos 
            self.exon_coverage = np.concatenate([
                self.exon_coverage, this_exon_coverage])
            #sorted_pos = sorted(
            #    self.position_in_exon_coverage_mapped_to_genomic_position.keys(),
            #key=lambda x: int(x))
