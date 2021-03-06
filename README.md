# ampatchlab/nf-typeHLA

[![Build Status](https://codebuild.ap-southeast-2.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiYWI0clRvd3JIOUNRWFFrajFwTTYrOHZha1FHQWVNa01ET2p1bUZ4cmRXRVdnUzRqVmFLV0dQWk1jT0VzWUZXMUNFRVRka28yaFFrQnNDWUg4M0llcDlvPSIsIml2UGFyYW1ldGVyU3BlYyI6InRpQkhSVVFkRThiQTBEQnQiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)](https://ap-southeast-2.console.aws.amazon.com/codesuite/codebuild/projects/nf-typeHLA/history)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

HLA typing Nextflow pipeline

## Usage

```
Usage:
    nextflow run -profile <profile> -revision <revision> ampatchlab/nf-typeHLA [options]


Nextflow execution options:

    -profile STR
        Nextflow configuration profile to use. Available profiles include:
        'awsbatch', 'conda', 'docker' and 'singularity'

    -revision STR
        Git branch/tag (version) of this workflow to use

    -work-dir DIR
        Directory where intermediate result files are stored

    -help
        Show additional execution options and exit


Required params:

    --readgroup_csv FILE
        Comma-separated list of sample and readgroup inputs


Reference genome params:

    --genome STR
        Reference genome name [Default: hs38DH]

    --idxbase FILE
        Override the BWA indexed FASTA file with FILE [Default: null]

    --hla_resource DIR
        Override the HLA resource directory with DIR [Default: null]


Adapter trimming params:

    --adapters STR
        The adapters to trim [Either: TruSeq, NexteraTransposase, BGISeq; Default: null]

    --r1_adapter_file FILE
        Override the R1 adapter file with FILE [Default: null]

    --r2_adapter_file FILE
        Override the R2 adapter file with FILE [Default: null]


Qualimap params:

    --qualimap_feature_file FILE
        Feature file with regions of interest in GFF/GTF or BED format [Default: null]


Output params:

    --publish_dir DIR
        Path where the results will be published [Default: ./results]

    --publish_mode STR
        The mode used to publish files to the target directory [Default: copy]


Standard params:

    --help
        Show this message and exit

    --version
        Show the pipeline version and exit
```

## Inputs

For paired-end data, the input CSV must have the following required columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq1: Absolute path of the 'R1' FASTQ file (required)
 * fastq2: Absolute path of the 'R2' FASTQ file (required)

For single-end data, the CSV must have the following columns:

 * sample: Unique sample name or ID (required)
 * readgroup: Unique readgroup name or ID (optional)
 * fastq: Absolute path of the FASTQ file (required)

If a particular sample has multiple FASTQ files (or pairs of FASTQ files), then these may
be specified on additional lines with a unique readgroup identifier. All readgroups belonging
to a particular sample will be aligned and merged.
