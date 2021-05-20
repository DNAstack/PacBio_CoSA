# Analysis of PacBio SARS-CoV-2 data using the CoSA pipeline

This repo provides a [WDL wrapper](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) for running [PacBio's CoSA pipeline](https://github.com/PacificBiosciences/CoSA). The CoSA pipeline produces an assembly, variants file, and lineage assignment using PacBio SARS-CoV-2 HiFi data. This pipeline uses DeepVariant for variant calling, but [other variant callers are available](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#4-variant-calling). Lineage assignment is run using [pangolin](https://github.com/cov-lineages/pangolin).

This workflow starts from primer-trimmed bam or fastq files; for instructions on how to convert a PacBio `subreads.bam` file to a primer-trimmed fastq or bam, see steps 1-3 [here](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#1-generate-ccs-data).


## Input

These values can be filled in in [the template inputs file](wdl/in.cosa.json).

- `String samplename`
- `File? primer_trimmed_bam` - a demultiplexed, single-sample bam file containing primer-trimmed reads; one of this or `primer_trimmed_fastq` must be provided
- `File? primer_trimmed_fastq` - a demultiplexed, single-sample fastq file containing primer-trimmed reads; one of this or `primer_trimmed_bam` must be provided
- `File reference` - [NC_045512.2.fasta](https://github.com/Magdoll/CoSA/blob/master/data/NC_045512.2.fasta)
- `File reference_index` - `${reference}.fai`; Generate using e.g. `samtools faidx NC_045512.2.fasta`


## Output

- `File aligned_bam`, `File aligned_bam_index` - reads aligned to the reference & index
- `File variants`, `File variants_index` - high quality variant calls & index
- `File consensus_sequence` - consensus sequence
- `File consensus_sequence_fragments` - consensus sequence broken up by 'N' regions
- `File lineage_metadata` - file containing lineage assignment (from `pangolin`) and tool version information


## Test data

A number of PacBio SARS-CoV-2 runs are available from [the SRA](https://www.ncbi.nlm.nih.gov/sra?term=(((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases))%20AND%20((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases)%20AND%20%22platform%20pacbio%20smrt%22%5BProperties%5D))%20AND%20%22pacbio%20smrt%22%5BPlatform%5D). Most of this data is available as primer-trimmed fastqs, which can be used as input directly.


## Running the workflow

WDL files may be run using a number of execution engines; common engines include [Cromwell](https://github.com/broadinstitute/cromwell) and [miniwdl](https://github.com/chanzuckerberg/miniwdl).

Before running any of the following commands, ensure that input files have been localized and their local paths filled in to the [in.cosa.json](wdl/in.cosa.json) file. Note that only one of `primer_trimmed_bam` and `primer_trimmed_fastq` should be provided; the other input line should be removed from the input file.


### Running with Cromwell

1. Retrieve the latest Cromwell jar file; `java` is required.
```
wget https://github.com/broadinstitute/cromwell/releases/download/63/cromwell-63.jar
```

2. Run the workflow.
```
java -jar cromwell-63.jar run wdl/cosa.wdl -i wdl/in.cosa.json
```


### Running with miniwdl

1. Install miniwdl via pip; `python3` is required.
```
python3 -m pip install miniwdl
```

2. Run the workflow.
```
miniwdl run wdl/cosa.wdl -i wdl/in.cosa.json
```
