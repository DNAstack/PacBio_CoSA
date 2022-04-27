# Analysis of PacBio SARS-CoV-2 data using the CoSA pipeline

This reposotiry provides a [WDL wrapper](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) for running [PacBio's CoSA pipeline](https://github.com/PacificBiosciences/CoSA) to process Pacific Biosciences SARS-CoV-2 long read HiFi data.

The CoSA pipeline produces an assembly and a variants file. This pipeline uses DeepVariant for variant calling, but [other variant callers are available](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#4-variant-calling).

This workflow starts from primer-trimmed fastq files; for instructions on how to convert a PacBio `subreads.bam` file to a primer-trimmed fastq or bam, see steps 1-3 [here](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#1-generate-ccs-data).


## Workflow inputs

An input template file with some defaults pre-defined can be found [here](./workflows/inputs.json).

| Input | Description |
|:-|:-|
| `accession` | Sample ID |
| `primer_trimmed_fastq` | A demultiplexed, single-sample fastq file containing primer-trimmed reads |
| `reference`, `reference index` | [The SARS-CoV-2 reference genome](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3) and its index |
| `container_registry` | Registry that hosts workflow containers. All containers are hosted in [DNAstack's Dockerhub](https://hub.docker.com/u/dnastack) [`dnastack`] |


## Workflow outputs

| Output | Description |
|:-|:-|
| `aligned_bam`, `aligned_bam_index` | Reads aligned to the SARS-CoV-2 reference genome |
| `high_quality_variants`, `high_quality_variants_index` | Variant calls and index |
| `consensus_fa` | Genome assembly |


## Test data

PacBio SARS-CoV-2 runs are available from [the SRA](https://www.ncbi.nlm.nih.gov/sra?term=(((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases))%20AND%20((txid2697049%5BOrganism%3Anoexp%5D%20NOT%200%5BMbases)%20AND%20%22platform%20pacbio%20smrt%22%5BProperties%5D))%20AND%20%22pacbio%20smrt%22%5BPlatform%5D). Most of this data is available as primer-trimmed fastqs, which can be used as input directly.
