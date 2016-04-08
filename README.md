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

Mapping
----
|    Folder     | Description |
|    ------     | ----------- |
|mapping/fastq  | Raw fastq   |
|mapping/temp_clipped | Fastq files ready to be mapped|
|mapping/sams   | Filtered sams - uniquely mapping, 20 AS |
|mapping/bed_uncollapsed | Uncollapsed bed files made from mapping/sams|
|mapping/bed_collapsed | Beds collapsed from bed_uncollapsed |

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

