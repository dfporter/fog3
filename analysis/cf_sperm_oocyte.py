import pandas
import sys
import matplotlib.pyplot as plt
import numpy as np
import re
import collections


def read_programs(lib):
    oogenic = {}
    spermatogenic = {}
    with open(lib['spermatogenic_program'], 'r') as f:
        for li in f:
            s = li.rstrip('\t').split('\t')
            spermatogenic[s[0]] = s[5]
    with open(lib['oogenic_program'], 'r') as f:
        for li in f:
            s = li.rstrip('\t').split('\t')
            oogenic[s[0]] = s[5]
    return oogenic, spermatogenic


def get_data(clip_fname, peaks_fname, lib):
    print locals()
    # Get maps.
    wb_id_to_gene_name = wb_to_public_name(lib['gtf_raw'])
    gene_name_to_wb_id = collections.defaultdict(list)
    for wbid in wb_id_to_gene_name:
        gene = wb_id_to_gene_name[wbid]
        gene_name_to_wb_id[gene].append(wbid)
    # Get RNA-seq
    gl_deseq = pandas.read_csv(lib['gl_expression'], sep='\t')
    print gl_deseq.columns
    gl_deseq = gl_deseq.to_dict('records')
    gl_by_gene = {}
    for row in gl_deseq:
        gl_by_gene[row['WormBase ID (WS240)']] = row
    # Get Clip.
    peaks = pandas.read_csv(peaks_fname, sep='\t')
    peaks = peaks.to_dict('records')
    peaks_by_gene = collections.defaultdict(dict)
    for row in peaks:
        if row['gene_id'] in peaks_by_gene:
            peaks_by_gene[row['gene_id']]['height'] += float(row['height'])
            continue
        peaks_by_gene[row['gene_id']] = row
        peaks_by_gene[row['gene_id']]['height'] = float(row['height'])
    return gl_by_gene, peaks_by_gene


def cf_deseq_w_sperm_oocyte_expression(clip_deseq_fname, peaks_fname, lib):
    # Build combined table.
    gl_by_gene, peaks = get_data(
        clip_deseq_fname, peaks_fname, lib)
    understood_genes = set(gl_by_gene.keys()) & set(peaks.keys())
    print """
    Genes in Noble et al. germline expression dataset (DESeq2) {i}
    Genes in FOG-3 iCLIP (called peaks): {pgs}
    Genes in both FOG-3 clip peaks and Noble datasets: {bpgs}
    Genes in FOG-3 clip peaks and not in germline: {npgs}
    """.format(i=len(gl_by_gene.keys()),
               pgs=len(peaks),
               bpgs=len(list(understood_genes)),
               npgs=len(set(peaks.keys()) - understood_genes)
    )
    split_by_noble_cat(gl_by_gene, peaks)
#    table = {}
#    for gene in understood_genes:
#        rnaseq = gl_by_gene[gene]['Fold change (of normalized reads)']
#        if rnaseq == 'Infinite':
#            rnaseq = 1e3
#        clip_val = np.log2(float(clip_by_gene[gene]['baseMean']))
#        rnaseq = np.log2(float(rnaseq))
        #if float(clip_by_gene[gene]['pvalue']) < 0.2:
#        table[gene] = [clip_val, rnaseq]
#    as_arr = table.values()
    #print table.values()
#    make_fig(table, output_filename=lib['figs'] + '/cf_w_sperm_oocyte.pdf')


def split_by_noble_cat(gl, peaks):
    by_cat = collections.defaultdict(dict)
    for gene in gl:
        cat = gl[gene]['Gene expression']
        by_cat[cat][gene] = gl[gene]
    print """
    From the dataset of germline genes from Noble, Ortiz et al.:
    Genes in the germline: {n}
    Categories: {c}
    """.format(n=len(gl), c=len(by_cat))
    for cat in by_cat:
        print "\tThere are {n} genes in the category {c}.".format(
            n=len(by_cat[cat]), c=cat
        )
    peaks_in_gl = [x for x in peaks if x in gl]
    peaks_by_category = {}
    gl_gene_by_category = {}
    for cat in by_cat:
        peaks_by_category[cat] = len([
            x for x in by_cat[cat] if x in peaks])
        gl_gene_by_category[cat] = len(by_cat[cat].keys())
    pie_of_categories(gl_gene_by_category.values(), gl_gene_by_category.keys(),
                       out_filename='figs/pie_of_gamete_categories_in_germline.pdf')
    pie_of_categories(peaks_by_category.values(), peaks_by_category.keys(),
                       out_filename='figs/pie_of_gamete_categories_in_peaks.pdf')
    for cat in ['Spermatogenic', 'Oogenic', 'Gender Neutral']:
        exp_in_fem3 = [
            float(by_cat[cat][gene]['Expression in fem-3 gonads bash(normalized reads)']) for gene in by_cat[cat]]
        exp_in_fog2 = [
            float(by_cat[cat][gene]['Expression in fog-2 gonads (normalized reads)']) for gene in by_cat[cat]]
        fold_change = [
            by_cat[cat][gene]['Fold change (of normalized reads)'] for gene in by_cat[cat]
        ]
        for i in range(len(fold_change)):
            if fold_change[i] == 'Infinite':
                fold_change[i] = 1.e3
        fold_change = [float(x) for x in fold_change]
        fold_change = np.array(fold_change)
        peak_heights = [
            float(peaks[gene]['height']) for gene in by_cat[cat] if gene in peaks
        ]
        lclip_peaks = [
            x for x in by_cat[cat] if x in peaks
        ]
        frac_clip = float(len(lclip_peaks))/float(len(peaks))
        frac_clip_in_gl = float(len(lclip_peaks))/float(len(peaks_in_gl))
        frac_gl = float(len(by_cat[cat].keys()))/float(len(gl.keys()))
        outside_cat = set()
        for other_cat in by_cat:
            if other_cat == cat: continue
            outside_cat = outside_cat | set(by_cat[other_cat].keys())
        import scipy.stats as sps
        ov_t = {'in_cat_and_clip': len(
            [x for x in by_cat[cat] if x in peaks_in_gl]
            ),
            'in_cat_and_not_clip': len(
                [x for x in by_cat[cat] if x not in peaks_in_gl]
            ),
            'not_in_cat_and_in_clip': len(
                [x for x in outside_cat if x in peaks_in_gl]
            ),
            'not_in_cat_and_not_in_clip': len(
                [x for x in outside_cat if x not in peaks_in_gl]
            )
                }
        table_for_stats = [
            [ov_t['in_cat_and_clip'], ov_t['in_cat_and_not_clip']],
            [ov_t['not_in_cat_and_in_clip'], ov_t['not_in_cat_and_not_in_clip']]
        ]
        odds, pval = sps.fisher_exact(table_for_stats)
        print """
        {c}: Median expression in fem-3 worms {n}
        {c}: Median expression in fog-2 worms {fog2}
        {c}: Median fold change spermatogenic/oogenic {q}
\t\t({percgl}% of germline genes)
        {c}: Median FOG-3 iCLIP height: {peaksh}
        {c}: Number of FOG-3 iCLIP targets: {numpeaks}\n\t\t({perc}% of iCLIP targets)\n\t\t({perc_in_gon}% of iCLIP targets in germline database)
        {c}: Significance of overlap with FOG-3 iCLIP peaks within gonad-expressed genes (Fisher): {pval}
        """.format(
            c=cat, n=np.median(exp_in_fem3),
            fog2=np.median(exp_in_fog2),
            q=np.median(fold_change),
            percgl="%.3f" % float(100 * frac_gl),
            peaksh=np.median(peak_heights),
            numpeaks=len(lclip_peaks),
            perc="%.3f" % float(100 * frac_clip),
            perc_in_gon="%.3f" % float(100 * frac_clip_in_gl),
            pval=pval
        )


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


def cf_w_sperm_oocyte_expression(peaks_fname):
    lib = config.config()
    deseq = pandas.read_csv(lib['gl_expression'], sep='\t')
    gl = deseq.to_dict('records')
    peaks = pandas.read_csv(peaks_fname, sep='\t')
    peaks = peaks.to_dict('records')
    gl_by_gene = {}
    for row in gl:
        gl_by_gene[row['Gene name']] = row
    for peak in peaks:
        if peak['gene_name'] in gl_by_gene:
            row = gl_by_gene[peak['gene_name']]
            peak.update(row)
    #make_fig(peaks)
    #print peaks

def make_fig(table, output_filename='figs/cf_w_sperm_oocyte.pdf'):
    points = np.array(table.values())
    print np.max([x[0] for x in points])
    print np.median([x[0] for x in points])
    plt.clf()
    plt.scatter([x[0] for x in points], [x[1] for x in points],
                c='black', marker='.', alpha=0.2)
    plt.xlabel('CLIP log2 coverage')
    plt.ylabel('sperm/oocyte log2')
    plt.savefig(output_filename, format='pdf')
    plt.clf()
    plt.close()


def pie_of_categories(vals, labels, out_filename='figs/pie.pdf'):
    print len(vals)
    print labels
    print vals
    if not isinstance(vals, type([])):
        vals = list(vals)
    labels = tuple(labels)
    import matplotlib.pyplot as plt
    plt.clf()
    plt.axis('equal')
    fig = plt.figure(1, figsize=(2,2))
    ax = fig.gca()
    plt.rcParams['patch.linewidth'] = 0
    wedges, texts, autotext = ax.pie(
        vals, labels=labels, colors=['gold', 'red', 'lightskyblue', 'green'],
        autopct='%.1f%%')
    ax.set_aspect('equal')
    plt.savefig(out_filename, format='pdf')
    for w in wedges:
        w.set_linewidth(0)
    plt.clf()
    plt.close()


def run(lib=None, deseq_fname='/groups/Kimble/Common/fog_iCLIP/pre_calls/fog_clip_deseq.txt',
        peaks_fname='/groups/Kimble/Common/fog_iCLIP/calls/for_comp/both_ratios_and_deseq.txt'):
    oogenic, spermatogenic = read_programs(lib)
    peaks_df = pandas.read_csv(peaks_fname, sep='\t')
    cat_oo = collections.defaultdict(set)
    cat_sp = collections.defaultdict(set)
    peaks_mrna = peaks_df[peaks_df['biotype']=='protein_coding']
    for gene in peaks_mrna['gene_id'].tolist():
        if gene in spermatogenic:
            cat_sp[spermatogenic[gene]].add(gene)
        else:
            cat_sp['Not spermatogenic'].add(gene)
        if gene in oogenic:
            cat_oo[oogenic[gene]].add(gene)
        elif gene in spermatogenic:
            cat_oo['Spermatogenic only\n'].add(gene)
        else:
            cat_oo['Not in germline'].add(gene)
    cats = {}
    # cats['Oogenic'] = cat_oo['Oogenic and Spermatogenic\n'] | cat_oo['Oogenic only\n']
    # cats['Spermatogenic'] = cat_sp['Spermatogenic only\n']
    oo_total = sum([float(len(cat_oo[x])) for x in cat_oo])
    sp_total = sum([float(len(cat_sp[x])) for x in cat_sp])
    cats_total = sum([float(len(cats[x])) for x in cats])
    for val in cat_oo:
        cat_oo[val] = 100 * len(list(cat_oo[val]))/oo_total
    for val in cat_sp:
        cat_sp[val] = 100 * len(list(cat_sp[val]))/sp_total
    pie_of_categories(cat_oo.values(), tuple(cat_oo.keys()), out_filename='figs/new_oogenic.pdf')
    pie_of_categories(cat_sp.values(), tuple(cat_sp.keys()), out_filename='figs/new_spermatogenic.pdf')

    cf_deseq_w_sperm_oocyte_expression(
        deseq_fname,
        peaks_fname, lib)
    print "end of run."
        #'/groups/Kimble/Common/fog_iCLIP/pre_calls/counts/')
    #cf_w_sperm_oocyte_expression(peaks_fname)


if __name__ == '__main__':
    peaks_fname = sys.argv[1]
    import config
    lib = config.config()
    run(peaks_fname=peaks_fname, lib=lib)
