import os
import sys
import pandas


def write_fastas(peaks_filename, output_filename=None):
    peaks = pandas.read_csv(peaks_filename, sep='\t')
    if not os.path.exists('fastas/'):
        os.system('mkdir fastas/')
    if 'biotype' in peaks.columns:
        no_ncrna = peaks[peaks['biotype']=='protein_coding']
        print peaks['biotype'].value_counts()
        no_ncrna_filename = 'fastas/{b}{x}_{v}_seqs_{y}_genes.fa'.format(
            b=os.path.basename(peaks_filename).rstrip('.txt'),
            x='_no_ncrna',
            v=len(no_ncrna.index),
            y=len(set(no_ncrna['gene_id']))
        )
        with open(no_ncrna_filename, 'w') as f:
            f.write(create_fasta_lines(no_ncrna))
    if output_filename is None:
        output_filename = 'fastas/{b}_{v}_seqs_{x}_genes.fa'.format(
            b=os.path.basename(peaks_filename).rstrip('.txt'),
            v=len(peaks.index),
            x=len(set(peaks['gene_id']))
        )
    with open(output_filename, 'w') as f:
        f.write(create_fasta_lines(peaks))


def create_fasta_lines(peaks):
    seqs = peaks['seq'].tolist()
    names = peaks['left'].tolist()
    genes = peaks['gene_id'].tolist()
    li = ''
    for i in range(len(seqs)):
        li += '>{a}{v}\n{b}\n'.format(a=names[i], v=genes[i], b=seqs[i])
    return li


if __name__ == '__main__':
    try:
        peaks_filename = sys.argv[1]
    except:
        print "run python write_fastas.py peaks_file"
        sys.exit()
    if not os.path.exists(peaks_filename):
        print "run python write_fastas.py peaks_file"
    write_fastas(peaks_filename)