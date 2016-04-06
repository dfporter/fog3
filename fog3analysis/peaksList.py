import pandas
import numpy as np


class peaksList(object):

    def __init__(self, dataframe=pandas.DataFrame(), name="Unnamed"):
        self.name = name
        self.df = dataframe

    def __str__(self):
        li = """
{a}:
{b} peaks
{c} genes
""".format(
    a=self.name,
    b=len(self.df.index),
    c=len(set(self.df['gene_name'].tolist()))
    )
        if 'Classification OO/SP GL' in self.df.columns:
            li += """
Classifications:
{d}""".format(d=str(self.df['Classification OO/SP GL'].value_counts()))
        return li
    

    def read_csv(self, fname):
        self.df = pandas.read_csv(fname, sep='\t')

    def read_sp_vs_oo(self, fname='lib/ortiz/DESeq_genes_in_gonad.txt'):
        """lib/ortiz/DESeq_genes_in_gonad.txt
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
        else:
            return np.nan

    def name_to_sp_rpkm(self, name):
        if name in self.gene_name_to_sp_rpkm:
            return self.gene_name_to_sp_rpkm[name]
        elif name in self.loci_to_sp_rpkm:
            return self.loci_to_sp_rpkm[name]
        else:
            return np.nan

    def annotate_sp_vs_oo(self):
        self.read_sp_vs_oo()
        self.df['Fold change SP/OO GL'] = [
            self.name_to_fc(x) for x in self.df['gene_name']]
        self.df['Classification OO/SP GL'] = [
            self.name_to_sp_oo_class(x) for x in self.df['gene_name']]
        self.df['Oo. expression'] = [
            self.name_to_oo_rpkm(x) for x in self.df['gene_name']]
        self.df['Sp. expression'] = [
            self.name_to_sp_rpkm(x) for x in self.df['gene_name']]
        t = zip(
                self.df['fog_TGGC'].tolist(),
                self.df['fog_GGTT'].tolist(),
                self.df['fog_GGCA'].tolist(),
                self.df['fog_CGGA'].tolist()
          )
        self.df['height'] = [np.mean(x) for x in t]
        print self.df['height']
        
    
