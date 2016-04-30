FOG-3 iCLIP analysis
====

Pre-processing
----
The initial fastq zip file and the result of splitting by barcode are not included in the repository, but are currently in pre_calls/fastq/.

Processing from fastq to uncollapsed bed is done with clip-preprocess, followed by collapsing reads (handled with the output from write_commands_for_collapsing.py in order to use multiple CPUs).
Bedgraphs are then made by running bed_to_wig.py, which also makes the normalized bedgraphs (coverage per million reads).

Scripts are in analysis/.
Most output is in peaks-by-permutations/.


```bash
python split_fastq_to_bed.py
# Split rRNA reads.
python write_commands_for_collapsing.py
# Ran these commands.
# Make bedgraphs.
python bed_to_wig.py
# Edit .ini files.
# Call peaks.
python peaks-by-permutations.py
```

The mapping and filtering can be summarized as:

```AsciiDoc
Reads were aligned to the WS235 genome using STAR with previously described
 parameters (CSEQ http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf),
 except for the parameters --alignEndsType Local and
 --outFilterMultimapScoreRange 0.
Uniquely mapping reads were identified by selecting MAPQ=255, high-confidence
 mappings were selected as those with alignment scores (STAR "AS" field) of
 at least 20.
 ```

```TeX
Reads were aligned to the WS235 genome using STAR with previously described
 parameters (CSEQ http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf),
 except for the parameters --alignEndsType Local and
 --outFilterMultimapScoreRange 0.
Uniquely mapping reads were identified by selecting MAPQ=255, high-confidence
 mappings were selected as those with alignment scores (STAR "AS" field) of
 at least 20.
 ```

```Tea
Reads were aligned to the WS235 genome using STAR with previously described
 parameters (CSEQ http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf),
 except for the parameters --alignEndsType Local and
 --outFilterMultimapScoreRange 0.
Uniquely mapping reads were identified by selecting MAPQ=255, high-confidence
 mappings were selected as those with alignment scores (STAR "AS" field) of
 at least 20.
 ```

```Twig
Reads were aligned to the WS235 genome using STAR with previously described
 parameters (CSEQ http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf),
 except for the parameters --alignEndsType Local and
 --outFilterMultimapScoreRange 0.
Uniquely mapping reads were identified by selecting MAPQ=255, high-confidence
 mappings were selected as those with alignment scores (STAR "AS" field) of
 at least 20.
 ```

```Vue
Reads were aligned to the WS235 genome using STAR with previously described
 parameters (CSEQ http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf),
 except for the parameters --alignEndsType Local and
 --outFilterMultimapScoreRange 0.
Uniquely mapping reads were identified by selecting MAPQ=255, high-confidence
 mappings were selected as those with alignment scores (STAR "AS" field) of
 at least 20.
 ```

Here's the mapping code:

```python

local_paths={
     'indexes': '/opt/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/',
     'gtf': '/opt/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf',
     'sjdb': '/opt/lib/sjdb.txt'}

def call_star(args, paths=None):
    in_dir = args.input_dir
    for fastq_filename in glob.glob(in_dir + '/*.fastq'):
        print fastq_filename
        bname = os.path.basename(fastq_filename).partition('.fastq')[0]
        # CSEQ parameters. http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf
        # Defaults changed:
        # --outFilterMultimapScoreRange is default 1
        # --alignEndsType EndToEnd for CSEQ.
        # We've added the --runThreadN for multithreading (faster speed, same output).
        #
        cmd = '''
STAR --alignIntronMax 1 --sjdbGTFfile {sjdb} \
--genomeDir {indexes} --readFilesIn {rin} --outSAMunmapped "Within" \
--outFilterMultimapNmax 3 --outFilterMismatchNmax 2 \
--seedSearchStartLmax 6 --winAnchorMultimapNmax 10000 --alignEndsType Local \
--sjdbGTFtagExonParentTranscript transcript_id \
--outFilterMultimapScoreRange 0  --runThreadN 8 \
--outFileNamePrefix {prefix} \
--alignTranscriptsPerReadNmax 100000
'''.format(sjdb=paths['sjdb'],
        indexes=paths['indexes'], rin=fastq_filename, prefix=bname + '_')
```

Here's where unique/multimapping reads are split:

```python
        with open('{d}/{a}'.format(d=unique_dir, a=os.path.basename(sam)), 'w') as uniquef:
            unique = get_lines_with_mapq_val(samn=sam, mapq='255')
            uniquef.write(unique)
            with open('{d}/{a}'.format(d=multi_dir, a=os.path.basename(sam)), 'w') as multif:
                multif.write(unique)
        with open('{d}/{a}'.format(d=multi_dir, a=os.path.basename(sam)), 'a') as multif:
            multi = get_lines_with_mapq_val(samn=sam, mapq='3', include_header=False)
            multif.write(multi)
```

And here's the filtering out of < 20 quality reads.
```python

    for sam in glob.glob('uniquely_mapping/*.sam'):
        header = ''
        outli = ''
        for li in open(sam).readlines():
            if li[0] == '@':
                header += li
                continue
            s = li.rstrip('\n').split('\t')
            as_flag = re.match('.*AS:i:(\d+).*', str(s[11:]))
            if as_flag is not None:
                as_value = int(as_flag.group(1))
            else:
                print "No AS value?"
            if as_value < 20: continue
            outli += li
        open('uniquely_mapping_20_AS/{a}'.format(a=os.path.basename(sam)), 'w').write(header + outli)
```

We will keep track of where the reads have gone in tables/read_stats.xls.

|    Folder     | Description |
|    ------     | ----------- |
|mapping        | Fastq files, beds, and sams|
|cims           | CIMS analysis and novoalign mapping|
|peaks-by-permutations | Peak calling |
|fog3analysis | Analyis scripts |
|analyis | Analysis scripts |
|fbf_analysis | Analysis scripts ported from FBF |
|calls | Old peak calling scripts |
|permutation_peaks | Peaks called from bedgraph_unnorm/ and bed_collapsed/ |
|bedgraph_unnorm | Produced from bed_collapsed. Not normalized to dataset size. |
|bedgraph_norm | Bedgraphs normalized to reads/million. |
|bed_collapsed | A copy of mapping/bed_collapsed |
|pre_calls | Old pre-processing. |
|old_analysis_files | Old analysis files. |


General analysis
----

Running find_peaks_by_permutations.py will make a call to peaks_by_pemutations/annotate_peaks.py.

This will load the -c .ini file, and add reads in peaks from lib['bedgraphs_folder'], which will be the unnormalized bedgraphs.
The added values will represent the max depth in the peak region in absolute read number.


Supp Dataset 1. FOG-3 iCLIP mapping statistics. 
Made by fog3/stats/file_stats.py

Supp Dataset 2. FOG-3 target RNAs.
Supp Dataset 3. FOG-3 CIMS.

```bash
python stats/peaks_file.py -i permutation_peaks

# To get Dataset 3 as well:
$ python stats/peaks_file.py -i permutation_peaks -c tables/cims_annotated.txt
```

From text
----

After normalization to a negative control (see Methods),
most of the enriched sequences were from mRNAs (Fig 4A).

```bash
make_figs.py
# Outputs to figs/location_*.pdf
```

Most of the 605  [larger list: 1388] FOG-3 mRNA targets are associated with
oogenesis (Fig 4B), as previously observed in microarray studies (cit XXX). 

```bash
make_figs.py
# Outputs to figs/gender_*.pdf
```

The iCLIP targets overlapped significantly with the microarray target list
(Supp Fig 5F-G, p value < 10-208 [larger list, 10^-227]), but was even more
depleted in genes associated with spermatogenesis
(Supp Fig 5H, p value < 10-24 [larger list < 10^49]).

```bash
make_figs.py
# Outputs to figs/*.pdf

# This is not discarding ncRNA. Doing so increases significance
# around 1e3 fold for both aligners.

# novoalign
Statistics for figs/overlap_with_rip_venn_clip.pdf overlap with RIP:
Genome size 12839 (all germline RNAs).
[[11220, 456], [952, 266]]
(6.875, 1.2939754179801647e-96)
Genome size 12839 (20k gene genome).
[[18326.0, 456], [952, 266]]
(11.229166666666666, 7.2166322496851259e-142)

# star
Statistics for figs/overlap_with_rip_venn_clip.pdf overlap with RIP:
     [[(-R, -C), (+R, -C)],
      [(-R, +C), (+R, +C)]]
Genome size 12839 (all germline RNAs).
[[11402, 492], [764, 230]]
(6.9767058272676969, 4.4869857195357939e-87)
Genome size 12839 (20k gene genome).
[[18514.0, 492], [764, 230]]
(11.328427616736901, 2.8738015736746662e-126)
```

A majority of the sequences mapped to the 3’ UTRs (Fig 4C-D),
with an apparent preference for the beginning versus the end of the region
(Fig 4E).

```bash
python fog3analysis/feature_locations/determine_feature_locations_ui.py \
       -c auto.ini -i permutation_peaks/combined_exp.txt
```

Furthermore,326/605 [larger list: 983/1388] mRNAs had two or more sequence
peaks spreading across the 3’UTR (Fig 4F, Supp Fig 5I-J). 

This binding pattern was reminiscent of other non-specific RNA binding
proteins, like ataxin and HuR (XXX cit). 

Few peaks were observed in the 5’ UTR and protein coding regions (Fig 4C-D),
implying FOG-3 RNA binding was restricted to the 3’ region. 

We did not observe a strong sequence preference in the identified peaks,
but in a subset of the peaks we did see a preference for CU repeats,
a UA rich motif and a CUCA motif (Supp Fig 5K-L). CU repeats may show be
evidence for sequence preference by FOG-3 or bias for UV protein-RNA
crosslinking (cit XXX).


Mapping
----

|    Folder     | Description |
|    ------     | ----------- | 
|mapping/fastq  | Raw fastq   |
|mapping/temp_clipped | Fastq files ready to be mapped|
|mapping/sams   | Filtered sams - uniquely mapping, 20 AS |
|mapping/bed_uncollapsed | Uncollapsed bed files made from mapping/sams|
|mapping/bed_collapsed | Beds collapsed from bed_uncollapsed |

Mapping statistics are obtained with:

```bash
# Looks for subfolders in the folders given by i, b and c arguments.
# The order does not matter.
python file_stats.py -i cims/ -b mapping/ -c ./
# This outputs Dataset 1 to tables/.

# With FBF:
# I put symlinks to the mapping folder in the cims folder
# so this script can find everything related to mapping, both for STAR
# and novoalign, in the cims/ folder.
dfporter at kimble-1 in /groups/Kimble/Common/fog_iCLIP/stats on master* b8caca3
$ python file_stats.py -i ../cims -b ../../fbf_celltype/cims/ -o fog3_file_stats
```

Gender analysis
----

Three ways to define gender:
    
    1. Use definitions from /opt/lib/ortiz/DESeq_genes_in_gonad.txt.
    2. Use fold changes from /opt/lib/ortiz/DESeq_genes_in_gonad.txt.
    3. Use the method from Noble et al.:
        * Read programs from /opt/lib/ortiz/TableS1_oogenic.txt and /opt/lib/ortiz/TableS2_spermatogenic.txt
        * Count oogenic as those only in oogenic dataset, ect.

In the case of 3., these files have the format:
    Wormbase ID	Ortiz Gonad	Sperm X	Maternal X	Oocyte X	Program
WBGene00235362	Not Present	Not Present	Oogenic	Not Present	Oogenic only

It seems quite clear that, on the whole, being a FOG-3 target does not
 correlate with SP/OO RPKM ratio - they are not enriched in either direction.
 In addition, SP/OO RPKM ratio does not correlate with FOG-3 peak height.
 Higher FOG-3 peaks are not more likely to be SP or OO enriched.
 Overall, the subset of SP or OO expressed RNAs that are FOG-3 targets
 is largely the set of RNAs that are abundant in both germlines.

The next question is - how does this abundance bias affect the SP/OO
 program definitions of Noble et al.? To answer this, we have to ask
 IF we have the observed abundance bias, THEN how likely is a given
 category size to be observed? So we take the observed range of abundances,
 and fit a histogram along each axis.

On a different note, we could take each gene's (x_sp, y_oo) coordinate,
 and rotate by 45 degrees.



CIMS analysis
----

|    Folder     | Description |
|    ------     | ----------- |
|cims/cims_tables | Results of cims analysis |
|cims/CIMS      | Zhang lab CIMS scripts with modifications|
|cims/fasta     | Fasta files from cims analysis |
|cims/novoaligned | Raw output from novoalign in .novo |
|cims/novo_tags | .novo files in novoaligned converted to .bed|
|cims/novo_tags_collapse | novo_tags after Zhang collapsing (still bed)|
|cims/collapsed_reformated | novo_tags_collapse reformated (still bed) for input into CIMS analysis |
|cims/mismatch_tags | bed-like, holds read and mutation |
|cims/mismatches_by_type | mismatch_tags split by ins/del/sub |
|cims_out | Raw output of CIMS scripts |

There is a cims github repo.

```bash
# This maps the fastq files and writes nohup commands for read collapsing.
python cims_master.py -i FOLDER_OF_FASTQ_FILES --lib fog_cims.ini --map
# Then run the read collapsing.
# Then do the cims analysis.
python cims_master.py -i FOLDER_OF_FASTQ_FILES --lib fog_cims.ini --cims
# Output is in cims_tables/

# To get dinucleotide frequencies in a fasta.
python dimers.py FASTA_FILE
# Outputs a heatmap.

# To get frequncy near a morif of interest:
pyhton pos_vs_motif.py -i cims_tables/ -m A_MOTIF
# Outputs a line plot and some seqlogos.
```

Novoalign peak calls
----

Could use the novoaligned reads for peak calls.
Work is in cims/.

```bash
# Generated bedgraphs.
# Split non-rRNA files.
# Generated auto.ini.
# Ran peak caller.

python peaks-by-permutations/annotate_peaks.py -p permutation_peaks/combined_control.txt -c auto.ini 
python peaks-by-permutations/annotate_peaks.py -p permutation_peaks/combined_exp.txt -c auto.ini
python ../fog3analysis/make_figs.py -i permutation_peaks/
mkdir data
mkdir data/feat_loc
mkdir figs/permutation_peaks
python ../fog3analysis/feature_locations/determine_feature_locations_ui.py -c auto.ini -i permutation_peaks/combined_exp.txt 

```
