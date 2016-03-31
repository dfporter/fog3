import sys
import os
import re
import glob


if len(sys.argv) > 1:
    top = sys.argv[1]
else:
    if raw_input(
        "Use current directory as top directory? [y/n]>").lower() \
        == 'y':
        top = '.'
    else:
        print "Exiting..."
        sys.exit()

# Defaults.
bed_dir = './bed_collapsed/'
bedgraphs_dir = './bedgraphs_unnorm/'
control_beds = []
exp_beds = []

if os.path.exists(top + '/bed_collapsed/'):
    bed_dir = top + '/bed_collapsed/'
    for rep in glob.glob(bed_dir + '/*.bed'):
        if re.search('control', rep):
            control_beds.append(rep)
        else:
            exp_beds.append(rep)
if os.path.exists(top + '/bedgraph_unnorm/'):
    bedgraphs_dir = top + '/bedgraph_unnorm'
elif os.path.exists(top + '/bedgraphs/'):
    if raw_input(
        "Looked for bedgraphs_unnorm but only found bedgraphs/ folder. \
Use that? [y/n]>").lower() == 'y':
        bedgraphs_dir = top + '/bedgraphs'

outli = """
[library]

top: {k}""".format(k=top)
outli += """
lib: %(top)s/lib/
fog3_lib: /groups/Kimble/Common/fog_iCLIP/calls/lib/
gtf_raw: %(lib)s/Caenorhabditis_elegans.WBcel235.78.noheader.gtf
fasta: %(lib)s/c_elegans.WS235.genomic.fa
fai: %(lib)s/c_elegans.WS235.genomic.fa.fai
chr_sizes: %(lib)s/chr_sizes.txt
gtf_one_txpt_per_gene: %(lib)s/gtf_one_txpt_per_gene.txt
gtf_pickle: %(top)s/data/feat_loc/gtf_data.p
gtf: %(lib)s/gtf_with_names_column.txt
fog3_rip_ranks: %(fog3_lib)s/fog3_rip.txt
fog1_rip_ranks: %(fog3_lib)s/fog1_2x_rip.txt

spermatogenic_program: %(fog3_lib)s/TableS2_spermatogenic.txt
oogenic_program: %(fog3_lib)s/TableS1_oogenic.txt
gl_expression: %(fog3_lib)s/ortiz/DESeq_genes_in_gonad.txt

#
txpt_obj_path: %(top)s/data/txpt_objects.p
feat_loc_dir:  %(top)s/data/feat_loc
#
"""
outli +="""
bedgraphs_folder: {v}
# Used by subpeaks.py
coverage_wigs:  {v}

# Used by subpeaks.py
bedgraph_exp_plus:  %(bedgraphs_folder)s/exp_+.wig
bedgraph_exp_minus: %(bedgraphs_folder)s/exp_-.wig
""".format(
    v=bedgraphs_dir)
outli += """
read_beds:  {a}
#
figs: %(top)s/figs/
#
""".format(a=bed_dir)
for i, x in enumerate(control_beds, start=1):
    outli += "control_bed{z}: {q}\n".format(
        z=i, q=os.path.basename(x))
for i, x in enumerate(exp_beds, start=1):
    outli += "exp_bed{z}: {q}\n".format(
        z=i, q=os.path.basename(x))

outli += """
clusters_dir: %(top)s/clusters/
permutation_peaks_dir: %(top)s/permutation_peaks/
motif: aatt
motif_control: ttaa
"""
print outli
"""
#control_bed1: control_AATA.bed
#control_bed2: control_CCGG.bed
#control_bed3: control_TTAA.bed
#control_bed4: control_TTGT.bed
#exp_bed1: fog_CGGA.bed
#exp_bed2: fog_GGCA.bed
#exp_bed3: fog_GGTT.bed
#exp_bed4: fog_TGGC.bed
control_bed1: n2_oo_lane1_rt15.bed
control_bed2: n2_oo_lane1_rt16.bed
control_bed3: n2_oo_lane1_rt3.bed

exp_bed1: fbf1_oo_lane2_rt1.bed
exp_bed2: fbf1_oo_lane2_rt6.bed
exp_bed3: fbf1_oo_lane2_rt9.bed
exp_bed4: fbf2_oo_lane2_rt11.bed
exp_bed5: fbf2_oo_lane2_rt13.bed
exp_bed6: fbf2_oo_lane2_rt2.bed

#exp_bed1: fbf1_sp_lane1_rt1.bed
#exp_bed2: fbf1_sp_lane1_rt6.bed
#exp_bed3: fbf1_sp_lane1_rt9.bed
#exp_bed4: fbf2_sp_lane1_rt13.bed
#exp_bed5: fbf2_sp_lane1_rt14.bed
#exp_bed6: fbf2_sp_lane1_rt2.bed


#
"""
