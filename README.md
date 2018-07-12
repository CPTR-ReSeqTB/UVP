# UVP

Set of scripts for analyzing NGS data

##Software dependencies:

BEDtools Version 2.17.0

Bcftools Version 1.2 

BWA Version 0.7.12

FastQC Version 0.11.5

Fastqvalidator Version 1.0.5

GATK Version 3.4.0

Kraken Version 0.10.5

Picard Version 1.134

Prinseq-lite.pl Version 0.204

Pigz Version 2.3.3

Qualimap Version 2.1.1

Samtools Version 1.2

SnpEff Version 4.1

Vcftools Version 0.1.126

## Installation

Insatlling the UVP on your local machine is straight forward. Clone the entire repository, and download the specific version of each of the third party tool listed above into the <Directory Path>/uvp/bin folder . You will need to edit the config.yml file in the <Directory Path>/uvp/bin folder to point to the correct directory and file paths of all the scripts and tools listed there in.

You will run the UVP using command line prompts, by invoking the UVP module in the <Directory Path>/uvp/scripts directory:
  <Directory Path>/uvp/scripts/UVP -q <input fastq> -r <path to H37Rv reference genome fasta file> -n <sample name> -q2 [paired fastq file] -a -v [verbose] 
  




