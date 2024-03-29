# Test demultiplexing of MiSeq data

This folder contains results from a test demultiplexing of BCL files obtained from the BC Genome Sciences
Center MiSeq Nano QC run of 44 TruSeq libraries prepared with primers designed by GIL with default parameters.
The purpose of this test was to ensure that the sample sheets generated by GIL work correctly with bcl2fastq,
and to assess the quality of the primers we ordered.

# Demultiplexing BCL files

The raw BCL files from the GSC were demultiplexed using bcl2fastq v2.19.0.316 using the `demultiplex.pbs` script
in the `demultiplexing` directory. The `RunInfo.xml` and sample sheet used for demultiplexing are also in
this same directory. The relevant demultiplexed fastq files can be found [here](data/output/).

# Assessing oligo quality

We ordered the indexing primers as standard 100 nmol desalted oligos from IDT with a single phosphorothioate
modification at the 3' end. We were warned that oligos of this quality may not perform as well as higher grades
of oligos, so we wanted to assess how many reads were not demultiplexed due to deletions in one or both indexes.
Deletions are the most common error in oligo synthesis, and a deletion in the index would "frame shift" the remainder
of the index, preventing bcl2fastq from successfully demultiplexing the read.

The analysis script can be found [here](analysis/demultiplexing_analysis.ipynb).