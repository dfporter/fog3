"""
Determines the average locations of an FBF peak, FBE, poly(A),
txpt bounds and CDS bounds.

1. Loads the entire gtf.
2. Load the peaks.
3. Subset the gtf to the targets.

Feature locations will be held in a class, flocs.
4. Marks the CDS and txpt bounds with mark_txpts.
"""
import HTSeq
from build_image_objects_for_heatmaps import *
import peak_in_gene_line_plot
import scatterplot_correlation_by_wig
import output_heatmap_figures
import logging
import datetime
import cliputil
from cliputil.gtf import *
from cliputil import rc
from cliputil.flocs import *

norm_left_len = 200.
norm_right_len = 200.
norm_cds_len = 1000.
utr_len = 800
verbose = False

logger = logging.getLogger(__name__)



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




def has_intron(txpt_id, gtf_df):
    if len(gtf_df[(gtf_df['transcript_id']==txpt_id) & (gtf_df['2']=='CDS')]) > 1:
        return True
    return False


def make_cdna_gtf(in_name, out_name):
    outf = open(out_name, 'w')
    with open(in_name, 'r') as f:
        header = next(f)
        outf.write(header)
        for li in f:
            if re.search('protein_coding', li) is not None:
                outf.write(li)


def load_signal_from_bedgraph(
        peaks, txpts, chr_lens,
        load_existing_data=False, do_combine_bedgraphs=False,
        bedgraphs_folder='data/wigs_coverage/',
        lib=None):
    if lib is not None:
        bedgraphs_folder = lib['bedgraphs_folder']
        txpts_obj_path = lib['txpt_obj_path']
    else:
        txpts_obj_path = 'data/for_heatmap.p'
    if load_existing_data:
        logger.info('%s: Loading data/for_heatmap.p holding the txpts object.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        with open(txpts_obj_path, 'r') as f:
            txpts = pickle.load(f)
    else:
        logger.info('%s: Building a genomic array from %s.' % (
            datetime.datetime.now().strftime('%Hh%Mm'), bedgraphs_folder))
        ga = get_bedgraph(do_combine_bedgraphs=do_combine_bedgraphs,
                          bedgraphs_folder=bedgraphs_folder, lib=lib)
        logger.info('%s: Adding raw reads to txpts.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        for gene_name in set(peaks['gene_name']):
            if gene_name not in txpts: continue
            txpts[gene_name].add_raw_reads_to_utr(ga, chr_lens)
            txpts[gene_name].add_raw_reads_to_all_peak_regions(ga)
        logger.info('%s: Saving txpts object to data/for_heatmap.' % (
            datetime.datetime.now().strftime('%Hh%Mm')))
        with open(txpts_obj_path, 'w') as f:
            pickle.dump(txpts, f)
    logger.info('Number of peaks: %i. Number of transcripts: %i.' % (
        len(peaks.index), len(txpts)
    ))
    locs = []
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("after loading, average number of peak locs {mn}".format(mn=np.mean(locs)))
    return peaks, txpts


def get_bedgraph(
        do_combine_bedgraphs=False, bedgraphs_folder='data/wigs/',
        lib=None):
    if lib is not None:
        bedgraphs_folder = lib['coverage_wigs']
        bedgraph_exp_plus = lib['bedgraph_exp_plus']
        bedgraph_exp_minus = lib['bedgraph_exp_minus']
    else:
        bedgraph_exp_plus = bedgraphs_folder + 'both_fbfs_plus.bed'
        bedgraph_exp_minus = bedgraphs_folder + 'both_fbfs_minus.bed'
    bedgraphs_folder = bedgraphs_folder.rstrip('/') + '/'
    if do_combine_bedgraphs: combine_bedgraphs(bedgraphs_folder=bedgraphs_folder)
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    with open(bedgraph_exp_plus, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '+')] = float(s[3])
    with open(bedgraph_exp_minus, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), '-')] = float(s[3])
    return ga


def combine_bedgraphs(bedgraphs_folder='data/wigs_five_prime/'):
    ga = {}
    bedgraphs_folder = bedgraphs_folder.rstrip('/') + '/'
    for filename_list in [
            (bedgraphs_folder + 'fbf1_reads_plus.bed',
             bedgraphs_folder + 'fbf1_reads_minus.bed',
             'combined_fbf1.txt'),
            (bedgraphs_folder + 'fbf2_reads_plus.bed',
             bedgraphs_folder + 'fbf2_reads_minus.bed',
             'combined_fbf2.txt')]:
        peaks_filename = filename_list[2]
        scatterplot_correlation_by_wig.load_bedgraph(filename_list, ga, use_key=peaks_filename)
    ga['combined'] = HTSeq.GenomicArray(chroms='auto', stranded=True)
    for iv, score in ga['combined_fbf1.txt'].steps():
        ga['combined'][iv] += score
    for iv, score in ga['combined_fbf2.txt'].steps():
        ga['combined'][iv] += score
#    with open('temp_ga.p', 'w') as f:
#        pickle.dump(ga, f)
    ga['combined'].write_bedgraph_file(bedgraphs_folder + 'both_fbfs_plus.bed', '+')
    ga['combined'].write_bedgraph_file(bedgraphs_folder + 'both_fbfs_minus.bed', '-')
    return ga


def get_args():
    parser = argparse.ArgumentParser(description='''
        Determine the locations of features in genes across the genome.''')
    parser.add_argument('-l', '--load_txpts',
                        help='Load marked transcript information.',
                        default=False, action='store_true')
    parser.add_argument('-g', '--load_gtf',
                        help='Load objects created by parsing a gtf file. Implies load_txpts is False.',
                        default=False, action='store_true')
    parser.add_argument('-c', '--config',
                        help='Directory holding config.py',
                        default='None')
#    parser.add_argument('-w', '--continuous_signal',
#                        default=False, action='store_true',
#                        help='Load true mapping data to display peak signal (slow).')
    args = parser.parse_args()
    return args


def create_one_txpt_per_gene(gtf_in, gtf_out):
    gtf = pandas.read_csv(gtf_in, sep='\t')
    gene_to_txpt = zip(gtf.gene_name, gtf.transcript_id)
    #txpt_to_left = zip(gtf.transcript_id, gtf.3)
    #txpt_to_right =zip(gtf.transcript_id, gtf.4)
    min_left = collections.defaultdict(list)
#    for t in gtf.transcript_id:
        

def get_txpts(args, sequences, rc_sequences, peaks, utr, exon, cds, other,
              lib=None):
    print 'get_txpts'
    if lib is None:
        feat_loc = 'data/feat_loc'
    else:
        feat_loc = lib['feat_loc_dir']
    if not os.path.exists(feat_loc):
        os.system('mkdir ' + feat_loc)
    if not args.load_txpts:
        print 'marking txpts...'
        txpts = mark_txpts(sequences, rc_sequences, peaks, utr, exon, cds, other)
        print len(txpts)
        print 'get_txpts after mark_txpts'
        print lib
        for t in txpts:
            txpts[t].check_coherency()
            txpts[t].motif = lib['motif']
            # flocs_map[gene_name].find_features()
            print "finding features"
            txpts[t].find_features()
        #add_peak_locations_to_transcripts(peaks, txpts, chr_lens)
        with open(feat_loc + '/%s.p' % os.path.basename(args.input).partition('.txt')[0], 'w') as f:
            pickle.dump(txpts, f)
    else:
        with open(feat_loc + '/%s.p' % os.path.basename(args.input).partition('.txt')[0], 'r') as f:
            txpts = pickle.load(f)
    logger.info("\tLoaded %i transcripts with features." % (len(txpts)))
    return txpts


def add_peak_locations_to_transcripts(peaks, txpts, chr_lens):
    peaks_d = to_dict(peaks)
    all_added_peak_locs = []
    for gene_name in peaks_d:
        if gene_name in txpts:
            txpts[gene_name].find_peaks(peaks_d, chr_lens)
    return txpts


def get_data(args, lib=None):
    if lib is None:
        lib = {'gtf': 'lib/gtf_with_names_column.txt',
               'gtf_one_txpt_per_gene': 'lib/gtf_one_txtp_per_gene.txt'}
    if not os.path.exists(lib['gtf_one_txpt_per_gene']):
        print "Creating a GTF-like of the longest txpt per gene..."
        cliputil.create_gtf_of_longest_txpt_per_gene(
            in_filename=lib['gtf'],
            out_filename=lib['gtf_one_txpt_per_gene'])
    print "Getting sequences..."
    (sequences, rc_sequences, chr_lens) = cliputil.get_sequences(lib)
    print "Getting GTF..."
    g = gtf(lib['gtf_one_txpt_per_gene'])
    g.flip_minus_strand_features(chr_lens)
    g.make_genes(sequences, rc_sequences, chr_lens)
    txpts = g.flocs_map
    (utr, exon, cds, other) = g.as_list()
    peaks = g.flip_minus_strand_peak_features(
        pandas.read_csv(args.input, sep='\t'), chr_lens)
    print '---- finished getting txpts. Adding peaks to transcripts...'
    txpts = add_peak_locations_to_transcripts(peaks, txpts, chr_lens)
    print txpts['mpk-1'].exons_dict
    print "Loading bedgraph signal..."
    peaks, txpts = load_signal_from_bedgraph(
        peaks, txpts, chr_lens, load_existing_data=args.load_txpts,
        lib=lib)
    print txpts['mpk-1'].exons_dict
    return peaks, txpts, chr_lens


def to_dict(df):
    _li = {}
    for index, row in df.iterrows():
        _li.setdefault(row['gene_name'], [])
        _li[row['gene_name']].append(row.to_dict())
    return _li


def make_figs(peaks, txpts, chr_lens, args, input_dir, lib=None):
    locs = []
    #pandas.DataFrame.sort(columns='height', ascending=False, inplace=True)
    peaks.sort(columns='height', ascending=False, inplace=True)
    output_dirname = 'figs/%s/' % os.path.basename(
        os.path.dirname(os.path.realpath(args.input)))
    for g in txpts:
        locs.append(len(txpts[g].peak_locs))
    logger.info("Assigned %i peaks and %i transcripts." % (len(peaks.index), len(txpts)))
#    logger.info("before loading, average number of peak locs {mn}".format(mn=np.mean(locs)))
#    args.load_txpts = True
    print "Assigned %i peaks and %i transcripts." % (len(peaks.index), len(txpts))
    logger.info('%s: Peak in gene line plot.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    peak_in_gene_line_plot.normalize_distances(txpts)
    logger.info('%s: Getting average positions.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_peak_v_polyA, ave_fbe_v_polyA,
     ave_highest_peak, ave_secondary_peaks) = peak_in_gene_line_plot.get_ave_pos(txpts)
    #logger.info('%s: Peak vs FBE line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    #peak_in_gene_line_plot.peak_vs_fbe(txpts)
    #logger.info('%s: Features in normalized gene line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_of_raw_signal(peaks, txpts, output_dirname=output_dirname,
                          include_motif=True)
    peak_in_gene_line_plot.plot_features(
       (ave_peaks, ave_fbes, ave_negs, ave_polyA, ave_highest_peak, ave_secondary_peaks),
       output_filename='figs/%s/features_in_normalized_gene.pdf' % os.path.dirname(input_dir))
    #logger.info('%s: Plotting UTR line plot.' % (
    #    datetime.datetime.now().strftime('%Hh%Mm')))
    peak_in_gene_line_plot.plot_utr(ave_peak_v_polyA, ave_fbe_v_polyA,
              output_filename='figs/%s/features_in_utr.pdf' % os.path.dirname(input_dir))
    logger.info('%s: Heatmap by gene length, only peak ranges.' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_by_gene_length_just_peak_ranges(
        peaks, txpts, output_dirname='figs/%s/' % os.path.dirname(input_dir))
    logger.info('\n***\n%s: Creating heatmaps.\n***' % (
        datetime.datetime.now().strftime('%Hh%Mm')))
    heatmap_of_peak_region(peaks, txpts, output_dirname=output_dirname)
    output_heatmap_figures.plot_intrapeak_distances(peaks, txpts, output_dirname=output_dirname)

