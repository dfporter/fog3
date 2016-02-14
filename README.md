FOG-3 iCLIP analysis
====

The initial fastq zip file and the result of splitting by barcode are not included in the repository, but are currently in pre_calls/fastq/.

Processing from fastq to uncollapsed bed is done with clip-preprocess, followed by collapsing reads (handled with the output from write_commands_for_collapsing.py in order to use multiple CPUs).
Bedgraphs are then made by running bed_to_wig.py, which also makes the normalized bedgraphs (coverage per million reads).

Scripts are in analysis/.
