import sys
from peaksList import *
import pandas

def load_peaks(infile):
    pk = peaksList(name='name')
    pk.read_csv(infile)
    pk.read_sp_vs_oo()
    pk.annotate_sp_vs_oo()
    pk.read_sp_vs_oo_as_programs()
    pk.add_rip()
    return pk

xls = pandas.read_excel(sys.argv[1])
xls.to_csv('temp.txt', sep='\t', index=False)
print xls
pk = load_peaks('temp.txt')
print pk
pkne = pk.df[pk.df['Program']=='Oogenic and Spermatogenic']
pkne.to_csv('Shared_both_oo_and_sp_program.txt', sep='\t', index=False)

pkne2 = pk.df[pk.df['Classification OO/SP GL']=='Gender Neutral']
pkne2.to_csv('Shared_gender_neutral_class_from_ortiz_deseq_file.txt',
            sep='\t', index=False)
