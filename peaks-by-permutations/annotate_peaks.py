__author__ = 'dp'
import pandas
import argparse
import re
import HTSeq
import sys
import collections
import numpy as np
import glob
import os
import config



def read_as_table(fname, use_col='WB ID', use_header_val=None):
    if use_header_val is not None:
        _table = {}
        key = use_header_val
        df = pandas.read_csv(fname, sep='\t')
        dfd = df.to_dict('records')
        for row in dfd:
            _table[row[key]] = row
        return _table
    with open(fname, 'r') as f:
        table = {}
        header = next(f).rstrip('\n').split('\t')
        for li in f:
            s = li.rstrip('\n').split('\t')
            row = {}
            for i, val in enumerate(header):
                row[val] = s[i]
            table[row[use_col]] = row  # By gene name.
    return table


def wb_to_public_name(gtf_fname):
    public_names = {}
    with open(gtf_fname, 'r') as f:
        for li in f:
            m = re.match('.*gene_id "([^"]+)";', li)
            if m is not None:
                name_m = re.match('.*gene_name "([^"]+)";', li)
                if name_m is not None:
                    public_names[m.groups()[0]] = name_m.groups()[0]
    return public_names


def read_as_table_of_lists(fname, use_col='WB ID', use_header_val=None):
    if use_header_val is not None:
        _table = collections.defaultdict(list)
        key = use_header_val
        df = pandas.read_csv(fname, sep='\t')
        dfd = df.to_dict('records')
        for row in dfd:
            _table[row[key]].append(row)
        return _table
    with open(fname, 'r') as f:
        table = {}
        header = next(f).rstrip('\n').split('\t')
        for li in f:
            s = li.rstrip('\n').split('\t')
            row = {}
            for i, val in enumerate(header):
                row[val] = s[i]
            table.setdefault(row[use_col], [])  # By gene name.
            table[row[use_col]].append(row)
    return table


def get_sequences(peaks, lib):
    fasta_filename = lib['fasta']
    # fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for row in peaks:
        start = int(row['left'])
        end = int(row['right'])
        chrm = row['chrm']
        seq = sequences[chrm][start:end]
        if row['strand'] == '-':
            seq = rc(seq)
        row['seq'] = seq

def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def load_gtf_as_dict_and_locate_in_gene(peaks, filename, use_header_val='gene_name'):
    gtf = read_as_table_of_lists(filename, use_header_val=use_header_val)
#    with open(filename, 'r') as f:
#        for li in f:
#            s = li.rstrip('\n').split('\t')
    for gene in peaks:
        for p in peaks[gene]:
            gene = p['gene_name']
            if gene not in gtf:
                p['location'] = 'Unknown'
                p['biotype'] = 'Unknown'
                p['gene_len'] = 'Unknown'
                print 'no'
                continue
            this_txpt = gtf[gene][0]['transcript_id']  # Expect only one txpt, but just in case.
            p['gene_len'] = sum([
                int(row['4']) - int(row['3']) for row in gtf[gene] if (
                    (row['transcript_id'] == this_txpt) and (row['2'] == 'exon')
                )
            ])
            is_nc = True
            p['biotype'] = gtf[gene][0]['biotype']
            for row in gtf[gene]:
                if row['2'] == 'CDS':
                    is_nc = False
            if is_nc:
                p['location'] = 'ncRNA'
                continue
            gene_left = min([int(x['3']) for x in gtf[gene] if x['2']=='CDS'])
            gene_right = max([int(x['4']) for x in gtf[gene] if x['2']=='CDS'])
            #gene_strand = rows['6']
            dist = get_distance((gene_left, gene_right),
                                (int(p['left']),
                                 int(p['right'])))
            if p['strand'] == '-':
                dist = -1 * dist
            p['dist_to_CDS'] = dist
            if not dist:
                p['location'] = "CDS"
            if dist < 0:
                p['location'] = '''5'UTR'''
            if dist > 0:
                p['location'] = '''3'UTR'''
    return peaks

# Duplicated from peak_locations.py.
def locate_in_gene(
        peak_list, gtf_filename):  #'lib/gtf_with_id_column.txt'):
    print "\tLocating peaks..."
    gtf_sep_cols = pandas.read_csv(gtf_filename, sep='\t')
    for n, peak in enumerate(peak_list):
        if not n % 200: print n
        #try:
        #    if not index % 1000:
        #        print "\tLocating peak %i..." % index
        #except:
        #    pass
        gene = str(peak['gene_name'])
        rows = gtf_sep_cols[gtf_sep_cols['gene_name']==gene]
        if len(rows) == 0:
            peak['location'] = "Unknown"
            continue
        peak['biotype'] = rows['biotype'].tolist()[0]
        rows = rows[rows['2']=='CDS']
        if len(rows) == 0:
            peak['location'] = "ncRNA"
            continue
        gene_left = min([int(x) for x in rows['3']])
        gene_right = max([int(x) for x in rows['4']])
        #gene_strand = rows['6']
        dist = get_distance((gene_left, gene_right),
                            (int(peak['left']),
                             int(peak['right'])))
        if peak['strand'] == '-':
            dist = -1 * dist
        peak['dist_to_CDS'] = dist
        if not dist:
            peak['location'] = "CDS"
        if dist < 0:
            peak['location'] = '''5'UTR'''
        if dist > 0:
            peak['location'] = '''3'UTR'''

# Duplicated from peak_locations.py
def get_distance(geneiv, peakiv):
    if geneiv[1] < peakiv[0]:
        return peakiv[0] - geneiv[1]
    if peakiv[1] < geneiv[0]:
        return peakiv[1] - geneiv[0]
    if geneiv[0] < peakiv[0] < geneiv[1]:
        if peakiv[1] < geneiv[1]:
            return 0  # Entirely in CDS.
        else:
            in_cds = geneiv[1] - peakiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[1]
    if geneiv[0] < peakiv[1] < geneiv[1]:
        if peakiv[0] > geneiv[0]:
            return 0  # Entirely in CDS.
        else:
            in_cds = peakiv[1] - geneiv[0]
            frac_in_cds = float(in_cds)/float(peakiv[1] - peakiv[0])
            if frac_in_cds > 0.5:
                return 0
            else:
                midpoint = float(peakiv[1] + peakiv[0])/2.0
                return  midpoint - geneiv[0]
    return 0


def get_comparisons(comparisons={}, peaks={}):
    comparisons['n2'] = ['fbf1n2_acaa', 'fbf1n2_ccgg', 'fbf1n2_tgcc', 'fbf2n2_acaa', 'fbf2n2_ccgg', 'fbf2n2_tgcc']
    comparisons['no_antibody'] = ['control_AATA',  'control_CCGG',
        'control_TTAA', 'control_TTGT']
    columns = peaks[peaks.keys()[0]].keys()
    del_keys = set()
    for key in comparisons:
        for name in comparisons[key]:
            if name not in columns:
                del_keys.add(key)
    for key in del_keys: del comparisons[key]
    return comparisons


def add_ratio_column(
        input_filename='fog3_peaks.txt', label='n2_ratio',
        positives_include='exp_.*\.txt', negatives_include='control_.*\.txt'):
    df = pandas.read_csv(input_filename, sep='\t')
    dfd = df.to_dict('records')
    positive_labels = [x for x in df.columns if (re.search(positives_include, x))]
    negative_labels = [x for x in df.columns if (re.search(negatives_include, x))]
    dfd = [add_ratio_to_row(row, positive_labels, negative_labels, label) for row in dfd]
    df = pandas.DataFrame(dfd)
    out_cols = get_cols(dfd)
    df.to_csv(input_filename, index=False, sep='\t')
    return df


def add_ratio_to_row(row, positive_labels, negative_labels, label):
    numerator = sum([row[x] for x in positive_labels])
    denominator = sum([row[x] for x in negative_labels])
    row[label] = int(float(numerator)/float(max([1., denominator])))
    return row


def _as_list(df):
    df_l = []
    for p in df:
        df_l.append(df[p])
    return df_l


def add_gene_name(peak_df, lib):  #, fname='./lib/gtf_with_names_column.txt'):
    fname = lib['gtf']
    gtf = read_as_table(fname, use_col='gene_id')
    for i, row in peak_df.iterrows():
        if row['gene_id'] in gtf:
            peak_df.loc[i, 'gene_name'] = gtf[row['gene_id']]['gene_name']
    return peak_df


def add_wb_name(
        peaks_fname, fname='./lib/gtf_with_names_column.txt'):
    print "Loading gtf..."
    gtf = read_as_table(fname, use_header_val='gene_name')
    print "...finished loading gtf."
    peaks_df = pandas.read_csv(peaks_fname, sep='\t')
    def name_to_id(name, gtf):
        if name in gtf: return gtf[name]['gene_id']
        else: return '***'
    peaks_df['gene_id'] = [name_to_id(x, gtf) for x in peaks_df['gene_name'].tolist()]
    peaks_df.to_csv(peaks_fname, sep='\t', index=False)


def fasta_from_peaks(peaks, filename='peaks.fa'):
    li = ''
    with open(filename, 'w') as f:
        for index, row in peaks.iterrows():
            li += '>%i\n%s\n' % (index, peaks.loc[index, 'seq'])
        f.write(li)


def get_cols(clip_list):
    try:
        obs_cols = clip_list[0].keys()
    except: 
        try: obs_cols = clip_list.columns  # Allow dataframes to be passed.
        except:
            print "get_cols(): Did not recognize type: %s" % str(clip_list)
    out_cols = [
        'chrm', 'left', 'right', 'gene_name', 'strand']
    ratio_cols = []
    for name in obs_cols:
        if re.match('.*ratio.*', name) is not None:
            print 'found %s' % name
            ratio_cols.append(name)
    out_cols.extend(ratio_cols)
    out_cols.extend(['unnorm_' + x.partition('.bed')[0] for x in get_file_list(lib)])
    extras = ['gene_id', 'biotype', 'location', 'RIP_target', 'RIP_rank', 'height', 'padj', 'pvalue', 'log2FoldChange', 'seq',
         #'fog_CGGA',  'fog_GGCA',  'fog_GGTT',  'fog_TGGC',
         #'control_AATA',  'control_CCGG',  'control_TTAA', 'control_TTGT',
              'chrom', 'gene_len', 'Seq ID', 'pval',# 'dist_to_CDS',
         #'fbf1n2_acaa', 'fbf1n2_ccgg', 'fbf1n2_tgcc', 'fbf2n2_acaa',
         #'fbf2n2_ccgg', 'fbf2n2_tgcc',
              'max_coverage']
    for name in extras:
        if name in obs_cols:
            out_cols.append(name)
    extra_unexpected = set(obs_cols) - set(out_cols)
    #for col in extra_unexpected: out_cols.append(col)
    print "Extra unexpected "
    print extra_unexpected
    return out_cols


def add_height(clip_list, positives_include='exp'):
    file_list = ['unnorm_' + x.partition('.bed')[0] for x in get_file_list(lib, just='positives')]
    file_list = set(file_list)
    print file_list
    print "***"
    for gene in clip_list:
        for p in clip_list[gene]:
            print p
            #if 'height' in p and p['height'] != 0: continue
            p['height'] = 0
            for key in [x for x in p.keys()]:# if x in file_list]:
                print ">" + key
                print "in %s?" % str(file_list)
                if str(key) in file_list:
                #if re.match('.*%s.*' % positives_include, key) is not None:
                    p['height'] += p[key]
                    print 'y'
                else: print 'n'


def read_bedgraph(fname, strand, bed=None):
    if bed is None:
        bed = HTSeq.GenomicArray('auto', stranded=True)
    else:
        assert(isinstance(bed, type(HTSeq.GenomicArray('auto', stranded=True))))
    with open(fname) as f:
        next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            bed[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)
                ] += int(float(s[3]))
    return bed


def get_file_list(lib, just=''):
    file_list = []
    if (just == '') or (just == 'positives'):
        found_all = False
        try_num = 1
        while not found_all:
            try_name = 'exp_bed{i}'.format(i=try_num)
            if try_name in lib:
                file_list.append(lib[try_name])
            else:
                found_all = True
            try_num += 1
    if (just == '') or (just == 'control') or (just == 'negatives'):
        found_all = False
        try_num = 1
        while not found_all:
            try_name = 'control_bed{i}'.format(i=try_num)
            if try_name in lib:
                file_list.append(lib[try_name])
            else:
                found_all = True
            try_num += 1
    print "get_file_list(lib): %s" % str(file_list)
    return file_list


def add_reads_in_peak_from_bedgraph(peak_df, bedgraph_folder, lib):
    expected = [x.rstrip('.bed') for x in get_file_list(lib)]
    # expected = [
    #              'fog_CGGA',  'fog_GGCA',  'fog_GGTT',  'fog_TGGC',
    #      'control_AATA',  'control_CCGG',  'control_TTAA', 'control_TTGT',
    # ]
    if 'chrom' in peak_df.columns and 'chrm' not in peak_df.columns:
        peak_df['chrm'] = peak_df['chrom']
    gas = read_bedgraph_folder(peak_df, bedgraph_folder, lib)
    peaks = peak_df.to_dict('records')
    for row in peaks:
        for label in expected:
            row['unnorm_' + label] = max_depth(
                gas[label], HTSeq.GenomicInterval(row['chrm'], row['left'],
                                                  row['right'], row['strand']))
    peak_df = pandas.DataFrame(peaks)
    return peak_df


def max_depth(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype='i'))


def read_bedgraph_folder(peak_df, bedgraph_folder, lib):
    gas = {}
    print bedgraph_folder
    print os.system('ls %s' % bedgraph_folder)
    # expected = [
    #              'fog_CGGA',  'fog_GGCA',  'fog_GGTT',  'fog_TGGC',
    #      'control_AATA',  'control_CCGG',  'control_TTAA', 'control_TTGT',
    # ]
    expected = [x.rstrip('.bed') for x in get_file_list(lib)]
    print expected
    for fname in glob.glob(bedgraph_folder + '/*_+.wig'):
        print "Found %s" % fname
        bname = os.path.basename(fname).partition('_+.wig')[0]
        if bname not in expected:
            print "Unexpected file: |{a}| Expected {aa}".format(a=bname,
                                                    aa=expected)
            continue
        gas[bname] = read_bedgraph(fname, '+')
    for fname in glob.glob(bedgraph_folder + '/*_-.wig'):
        bname = os.path.basename(fname).partition('_-.wig')[0]
        if bname not in expected:
            print "Unexpected file: {a}".format(a=bname)
            continue
        gas[bname] = read_bedgraph(fname, '-', bed=gas[bname])
    return gas


def gene_name_from_wb_id(clip_df, lib):
    wb = wb_to_public_name(lib['gtf'])
    if 'gene_id' in clip_df.columns:
        if 'gene_name' not in clip_df.columns:
            clip_df['gene_name'] = clip_df['gene_id']
        for i, row in clip_df.iterrows():
            if row['gene_id'] in wb.keys():
                clip_df.loc[i, 'gene_name'] = wb[row['gene_id']]
            else:
                print 'not found'
    elif clip_df.loc[0, 'gene_name'][:6] == 'WBGene':
        clip_df['gene_id'] = clip_df['gene_name']
        for i, row in clip_df.iterrows():
            if row['gene_id'] in wb.keys():
                clip_df.loc[i, 'gene_name'] = wb[row['gene_id']]
            else:
                print 'not found'
    elif ('gene' in clip_df.columns) and (
                clip_df.loc[0, 'gene'][:6] == 'WBGene'):
        clip_df['gene_id'] = clip_df['gene']
        for i, row in clip_df.iterrows():
            if row['gene_id'] in wb.keys():
                clip_df.loc[i, 'gene_name'] = wb[row['gene_id']]
            else:
                print 'not found'

    return clip_df


def run(lib, peaks_filename):
    peaks_df = pandas.read_csv(peaks_filename, sep='\t')
    peaks_df = add_reads_in_peak_from_bedgraph(peaks_df, lib['bedgraphs_folder'], lib)
    print peaks_df.head(1)
    if 'gene_name' not in open(peaks_filename, 'r').readline().rstrip('\n').split('\t'):
        print "Adding a gene_name column..."
        peak_df = add_gene_name(peaks_df, lib)
        print peaks_df.columns
    peaks_df.to_csv(peaks_filename, sep='\t', index=False)
    if 'gene_id' not in open(peaks_filename).readline().rstrip('\n').split('\t'):
        print "Adding a gene ID column..."
        add_wb_name(peaks_filename)
    print "Calculating ratios..."
    add_ratio_column(input_filename=peaks_filename,
                     label='n2_ratio', positives_include='fog')
    clip = read_as_table_of_lists(peaks_filename, use_header_val='gene_id')
    if 'fog3_rip_ranks' in lib:
        rip = read_as_table(lib['fog3_rip_ranks'], use_header_val='WB ID')
        print "Peaks %i, ex: %s" % (len(clip.keys()), str(clip.keys()[0]))
        for gene in clip:
            if gene in rip:
                extra = {'RIP_target': 1, 'RIP_rank': rip[gene]['Rank'],
                              'Seq ID': rip[gene]['Seq ID']}
            else:
                extra = {'RIP_target': 0, 'RIP_rank': -1, 'Seq ID': ''}
            for row in clip[gene]:
                row.update(extra)
    clip = load_gtf_as_dict_and_locate_in_gene(clip, lib['gtf'],
                                               use_header_val='gene_name')
    add_height(clip, positives_include='exp')
    flattened = []
    for gene in clip: flattened.extend(clip[gene])
    get_sequences(flattened, lib)
    clip_df = pandas.DataFrame(flattened)
    out_cols = get_cols(clip_df)
    fasta_from_peaks(clip_df, filename='peaks.fa')
    clip_df = clip_df[get_cols(clip_df)]
    clip_df = gene_name_from_wb_id(clip_df, lib)
    print clip_df.loc[0]
    clip_df.sort(columns=['height'], inplace=True, ascending=False )
    clip_df.to_csv(peaks_filename,
                   # cols=out_cols,
                   sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Annotates a peaks file (IN PLACE!).''')
    parser.add_argument(
        '-p', '--peaks', default='all_peaks.txt',
        help='''A file or folder of peaks.'''
    )
    parser.add_argument(
        '-c', '--config_ini', default='analysis/',
        help='config.ini file.'
    )
    args = parser.parse_args()
    lib = config.config(args.config_ini)
    if os.path.isfile(args.peaks):
        run(lib, args.peaks)
    elif os.path.isdir(args.peaks):
        for peaks_fname in glob.glob(args.peaks +'/*'):
            run(lib, peaks_fname)
