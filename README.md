# HPV_EM_Pipeline

The HPV EM Pipeline is an alignment tool that utilizes an expectation maximization algorithm to identify the potential presence of different HPV types in a particular genome. 

## Prerequisites
  - python
  - bwa
  - STAR
    - to compile STAR, you may need to install CMAKE or gcc
  - samtools
  
## Installation
  TBD
  
## Setup
  To begin, assure that you know the path to your reference fasta file, as well as the human reference genome (as needed). If you plan to use the pipeline on multiple files of the same name, multiple folders will be created (as shown in the Examples section).
  
## Input
  Upon running the pipeline with the -h option (or --help), the following options are presented:
  
  Positional arguments
  - sampleName : The name of the sample to be aligned. This file can be in either .bam or .fq (fastq) format. If it's in the fastq format, a second fastq file can be aligned using the -2 option explained below.
  
  - refFasta : The name of the reference fasta file to be used for alignment. For most cases, the reference fasta (.fa) file should be the human genome (which needs to be installed separately from the pipeline).
  
  - path : This argument indicates the path to the human genome directory as needed by STAR, which includes files such as genomeParameters.txt and chrName.txt.  
  
  Optional Arguments 
  - -h, -help :
  
  - -@ (CPUS) :
  
  - -2 (otherSample) :

## Output
  The output is formatted as 

## Examples
