#!/bin/bash

#SBATCH --job-name=PMthreetophat.sh
#SBATCH --time=12:00:00
#SBATCH --mail-user=nevinsc@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=PMthreetophat.log
#SBATCH --error=PMthreetophat.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8Gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcb4325
#SBATCH --qos=mcb4325


module load ufrc
module load bowtie2
module load tophat

#bowtie2-build hg38.fa hg38

tophat --no-coverage-search hg38 SRR445721.fastq

