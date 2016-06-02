import pandas
import numpy as np
import re
import sys
import collections
import HTSeq


from locatedPeak import locatedPeak

class peaksList(object):

    def __init__(self, dataframe=pandas.DataFrame(), name="Unnamed",
                 gene_name_col='gene_name'):
        self.name = name
        self.df = dataframe
        self.gene_name_col = gene_name_col

    def __str__(self):
        li = """
{a}:
{b} peaks
{c} genes
""".format(
    a=self.name,
    b=len(self.df.index),
    c=len(set(self.df[self.gene_name_col].tolist()))
    )
        if 'Classification OO/SP GL' in self.df.columns:
            li += """\nClassifications  (from Ortiz Deseq file):\n{d}""".format(
                d=str(self.df['Classification OO/SP GL'].value_counts()))
        if 'Program' in self.df.columns:
            li += """\nPrograms:\n{d}""".format(
                d=str(self.df['Program'].value_counts()))
        return li
    
    def read_csv(self, fname):
        self.df = pandas.read_csv(fname, sep='\t', index_col=False)
        if ('Wormbase ID' not in self.df.columns) and (
            'WB ID' in self.df.columns):
            self.df['Wormbase ID'] = self.df['WB ID']
        if ('WB ID' not in self.df.columns) and (
            'Wormbase ID' in self.df.columns):
            self.df['WB ID'] = self.df['Wormbase ID']
        if ('Wormbase ID' not in self.df.columns) and (
            'WB Gene ID' in self.df.columns):
            self.df['Wormbase ID'] = self.df['WB Gene ID']

    def to_csv(self, fname):
        self.df.to_csv(fname, sep='\t', index=False)

    def dict_by_first_instance(self, list_of_tups):
        if len(list_of_tups) < 1: return collections.defaultdict(str)
        if len(list_of_tups[0]) < 2: return collections.defaultdict(str)
        to_y = collections.defaultdict(type(list_of_tups[1]))
        for x in list_of_tups:
            if x[0] in to_y: continue
            to_y[x[0]] = x[1]
        return to_y

    def add_rip(self, fname='/opt/lib/fog3rip/TableS4_fog3_rip_targets.txt'):
        rip = pandas.read_csv(fname, sep='\t', index_col=False)
        rip['SAM rank'] = [int(x) for x in rip['FOG-3 Rank']]
        #rip_top722 = rip[rip['SAM rank']<723]
        to_sam = self.dict_by_first_instance(
            zip(rip['WB ID'].tolist(), rip['SAM rank'].tolist()))
        rip_wb = set(rip['WB ID'].tolist())
        self.df['Is RIP target?'] = [
            (x in rip_wb) for x in self.df['Wormbase ID'].tolist()]
        self.df['SAM rank'] = [
            to_sam[x] for x in self.df['Wormbase ID'].tolist()]

    def add_permutation_peaks(self, fname):
        peaks = pandas.read_csv(fname, sep='\t', index_col=False)
        to_peaks = self.dict_by_first_instance(
            zip(peaks['gene_id'].tolist(), peaks['height'].tolist()))
        clip_wb = set(peaks['gene_id'].tolist())
        print "clip wb IDs: {a}".format(a=len(clip_wb))
        
        self.df['Is CLIP target?'] = [
            (x in clip_wb) for x in self.df['Wormbase ID']]
        self.df['CLIP height'] = [
            to_peaks[x] for x in self.df['Wormbase ID']]

    def read_sp_vs_oo_as_programs(
            self, fname_oo='/opt/lib/ortiz/TableS1_oogenic.txt',
            fname_sp='/opt/lib/ortiz/TableS2_spermatogenic.txt',
            gtf_fname='/opt/lib/gtf_with_names_column.txt'):
        print "self.gene_name_col="
        print self.gene_name_col
        if not hasattr(self, 'wbid_to_name'):
            self.transl(gtf_fname)
        oo_df = pandas.read_csv(fname_oo, sep='\t', index_col=False)
        sp_df = pandas.read_csv(fname_sp, sep='\t', index_col=False)
        sp_df['Program'] = [
            re.sub('Spermatogenic and Oogenic', 'Oogenic and Spermatogenic', x) \
            for x in sp_df['Program'].tolist()]
        wb_program_oo = dict(zip(oo_df['Wormbase ID'].tolist(),
                                oo_df['Program'].tolist()))
        wb_program_sp = dict(zip(sp_df['Wormbase ID'].tolist(),
                                sp_df['Program'].tolist()))
        for k in (set(wb_program_oo) & set(wb_program_sp)):
            if wb_program_oo[k] != wb_program_sp[k]:
                print "Contradiction between %s and %s: %s (%s vs %s)." % (
                    fname_oo, fname_sp, k, wb_program_oo[k], wb_program_sp[k])
                sys.exit()
        wb_program_oo.update(wb_program_sp)
        self.program = wb_program_oo
        if ('Wormbase ID' not in self.df.columns):
            self.df['Wormbase ID'] = [self.name_to_wbid[x] for x in \
                                      self.df[self.gene_name_col].tolist()]
        self.df['Program'] = [self.what_program(x) for x in \
                              self.df['Wormbase ID'].tolist()]

    def what_program(self, x):
        if x in self.program:
            return self.program[x]
        else:
            return ''

    def transl(self, gtfname):
        gtf = pandas.read_csv(gtfname, sep='\t')
        self.wbid_to_name = collections.defaultdict(str)
        self.wbid_to_name.update(dict(zip(
            gtf['gene_id'].tolist(), gtf['gene_name'].tolist())))
        self.name_to_wbid = collections.defaultdict(str)
        self.name_to_wbid.update(dict(zip(
            gtf['gene_name'].tolist(), gtf['gene_id'].tolist())))
        self.name_to_biotype = dict(zip(
            gtf['gene_name'].tolist(), gtf['biotype'].tolist()))

    def add_gene_name(self, gtfname):
        if not hasattr(self, 'wbid_to_name'):
            self.transl(gtfname)
        self.df['gene_name'] = [self.wbid_to_name[x] for x in self.df.gene_id]

    def add_biotype(self, gtfname):
        if not hasattr(self, 'name_to_biotype'):
            self.transl(gtfname)
        self.df['biotype'] = [self.name_to_biotype[x] for x in self.df.gene_name]

    def read_sp_vs_oo(self, fname='/opt/lib/ortiz/DESeq_genes_in_gonad.txt'):
        """/opt/lib/ortiz/DESeq_genes_in_gonad.txt
        """
        self.ortiz = pandas.read_csv(fname, sep='\t')
        fc = self.ortiz['Fold change (of normalized reads)']
        self.loci_to_fc = dict(zip(self.ortiz['Gene ID'].tolist(),
                  self.ortiz['log2 fold change (of normalized reads)'].tolist()))
        self.gene_name_to_fc = dict(zip(self.ortiz['Gene name'].tolist(),
                  self.ortiz['log2 fold change (of normalized reads)'].tolist()))
        self.gene_name_to_class = dict(zip(
            self.ortiz['Gene name'].tolist(),
            self.ortiz['Gene expression'].tolist()))
        self.loci_to_class = dict(zip(
            self.ortiz['Gene ID'].tolist(),
            self.ortiz['Gene expression'].tolist()))
        self.gene_name_to_oo_rpkm = dict(zip(
            self.ortiz['Gene name'].tolist(),
            self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.loci_to_oo_rpkm = dict(zip(
            self.ortiz['Gene ID'].tolist(),
            self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.gene_name_to_sp_rpkm = dict(zip(
            self.ortiz['Gene name'].tolist(),
            self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))
        self.loci_to_sp_rpkm = dict(zip(
            self.ortiz['Gene ID'].tolist(),
            self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))
        self.wb_id_to_oo_rpkm = dict(zip(
            self.ortiz['WormBase ID (WS240)'].tolist(),
            self.ortiz['Expression in  fog-2 gonads (RPKM)'].tolist()))
        self.wb_id_to_sp_rpkm = dict(zip(
            self.ortiz['WormBase ID (WS240)'].tolist(),
            self.ortiz['Expression in fem-3 gonads (RPKM)'].tolist()))

    def name_to_fc(self, name):
        if name in self.gene_name_to_fc:
            return self.gene_name_to_fc[name]
        elif name in self.loci_to_fc:
            return self.loci_to_fc[name]
        else:
            return np.nan

    def name_to_sp_oo_class(self, name):
        if name in self.gene_name_to_class:
            return self.gene_name_to_class[name]
        elif name in self.loci_to_class:
            return self.loci_to_class[name]
        else:
            return np.nan

    def name_to_oo_rpkm(self, name):
        if name in self.gene_name_to_oo_rpkm:
            return self.gene_name_to_oo_rpkm[name]
        elif name in self.loci_to_oo_rpkm:
            return self.loci_to_oo_rpkm[name]
        elif name in self.wb_id_to_oo_rpkm:
            return self.wb_id_to_oo_rpkm[name]
        else:
            return np.nan

    def name_to_sp_rpkm(self, name):
        if name in self.gene_name_to_sp_rpkm:
            return self.gene_name_to_sp_rpkm[name]
        elif name in self.loci_to_sp_rpkm:
            return self.loci_to_sp_rpkm[name]
        elif name in self.wb_id_to_sp_rpkm:
            return self.wb_id_to_sp_rpkm[name]
        else:
            return np.nan

    def annotate_sp_vs_oo(self):
        self.read_sp_vs_oo()
        print '%s annotate_sp_vs_oo, df.columns:' % self.name
        print self.df.columns
        self.df['Fold change SP/OO GL'] = [
            self.name_to_fc(x) for x in self.df[self.gene_name_col]]
        self.df['Classification OO/SP GL'] = [
            self.name_to_sp_oo_class(x) for x in self.df[self.gene_name_col]]
        self.df['Oo. expression'] = [
            self.name_to_oo_rpkm(x) for x in self.df[self.gene_name_col]]
        self.df['Sp. expression'] = [
            self.name_to_sp_rpkm(x) for x in self.df[self.gene_name_col]]
        if 'fog_TGGC' in self.df.columns:
            t = zip(
                    self.df['fog_TGGC'].tolist(),
                    self.df['fog_GGTT'].tolist(),
                    self.df['fog_GGCA'].tolist(),
                    self.df['fog_CGGA'].tolist()
              )
            self.df['height'] = [np.mean(x) for x in t]
            print self.df['height']
        
    def add_reads(
        self, ga=HTSeq.GenomicArray('auto'), name='unnamed array'):
        if self.df is None:
            print "Asked to set reads on an empty peaks object!"
            return False
        peaks = self.df  # Shorter.
        try:
            ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(),
                peaks['right'].tolist(), peaks['strand'].tolist())
        except:
            print "Could not create ivs from %s." % self.name
            return False
        self.df[name] = [
            self.get_val(ga, HTSeq.GenomicInterval(*iv)) for iv in ivs]
#        self.df = self.df
    
    def get_val(self, ga, iv):
        return np.max(np.fromiter(ga[iv], dtype=np.float))
    
    def add_location_from_integral(self, gtf, ga, use_name=None):
        if use_name is None: use_name = self.gene_name_col
        #assert(type(gtf) == type(pandas.DataFrame()))
        if 'chrm' not in self.df.columns:
            self.df['chrm'] = self.df['chrom']
        ivs = zip(self.df.chrm, self.df.left, self.df.right,
                  self.df.strand, self.df[use_name].tolist())
        located_peaks = [locatedPeak(*iv) for iv in ivs]
        locations = [x.add_location_from_integral(gtf, ga) for \
                     x in located_peaks]
        self.df['location'] = locations
