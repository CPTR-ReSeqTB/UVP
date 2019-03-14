# UVP

Set of scripts for analyzing NGS data

## Software dependencies:

| Name          | Version   |
|---------------|-----------|
| `bedtools`    | `2.17.0`  |
| `bcftools`    | `1.2`     |
| `bwa`         | `0.7.12`  |
| `fastqc`      | `0.11.5`  |
| `fqtools`     | `2.0`     |
| `gatk`        | `3.6`     |
| `kraken`      | `0.10.5`  |
| `picard`      | `1.141`   |
| `prinseq`     | `0.20.4`  |
| `pigz`        | `2.3.4`   |
| `qualimap`    | `2.1.3`   |
| `samtools`    | `1.2`     |
| `snpeff`      | `4.1`     |
| `vcftools`    | `0.1.16`  |

## Installation

The UVP requires at least 100GB RAM and up to 100GB storage space to run locally. Insatlling the UVP on your local machine is straight forward. Clone the entire repository, then create a conda environment using the `environment.yml` file provided:

```
conda env create -f environment.yml
```

Then activate the environment:

```
conda activate reseqtb-uvp
```

You will need to edit the config.yml file to point to your kraken database.

You will run the UVP using command line prompts as follows:

```
uvp -q 'input fastq' -r 'path to H37Rv reference genome fasta file' -n 'sample name' -q2 'paired fastq file' -a -v 
```
