# UVP

Set of scripts for analyzing NGS data

## Software dependencies:

| Name                | Version   |
|---------------------|-----------|
| `bedtools`          | `2.17.0`  |
| `bcftools`          | `1.2`     |
| `bwa`               | `0.7.12`  |
| `fastqc`            | `0.11.5`  |
| `fqtools`           | `2.0`     |
| `gatk`              | `3.6`     |
| `kraken`            | `0.10.5`  |
| `picard`            | `1.141`   |
| `prinseq`           | `0.20.4`  |
| `pigz`              | `2.3.4`   |
| `qualimap`          | `2.1.3`   |
| `samtools`          | `1.2`     |
| `snpeff`            | `4.1`     |
| `perl-vcftools-vcf` | `0.1.16`  | 
| `vcftools`          | `0.1.16`  |

## Installation

The UVP requires at least 100GB RAM and up to 100GB storage space to run locally. Insatlling the UVP on your local machine is straight forward. Clone the entire repository, then create a conda environment using the `environment.yml` file provided:

```
conda env create -f environment.yml
```

Then activate the environment:

```
conda activate reseqtb-uvp
```

With the `reseqtb-uvp` conda environment active, run the following command from withing the top-level directory of the 'UVP' repository to install the UVP code into the conda environment:

```
pip install -e .
```

Note: the `-e` flag installs UVP in 'editable' mode so that any changes to the code-base will be automatically installed into the conda environment.

Due to the terms of the GATK v3 license, it will be necessary to download and install the GATK `.jar` file separately. Download the GATK `.tar.bz2` file from this link:

```
https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.6-0-g89b7209
```

...then decompress it:

```
tar -xjvf GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 
```

...and (with the reseqtb-uvp conda environment active) register the `GenomeAnalysisTK.jar` file:

```
gatk-register GenomeAnalysisTK.jar
```

The `GenomeAnalysisTK.jar` file will be copied into the appropriate directory inside your `reseqtb-uvp` conda environment.

## Download the M. tuberculosis H37Rv snpEff database

A M. tuberculosis snpEff database is necessary for annotating variants. With the `reseqtb-uvp` conda environment active, run the following command to download it:

```
snpEff download m_tuberculosis_H37Rv
```

## Configuration

You will need to edit the config.yml file to point to your kraken database.

## Running the UVP

Run the UVP using command line prompts as follows:

```
uvp -q 'input fastq' -r 'path to H37Rv reference genome fasta file' -n 'sample name' -q2 'paired fastq file' -a -v 
```
