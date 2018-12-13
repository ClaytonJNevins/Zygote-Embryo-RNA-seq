#!/bin/bash

#SBATCH --job-name=PMthreetophat3.sh
#SBATCH --time=120:00:00
#SBATCH --mail-user=nevinsc@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=PMthreetophat3.log
#SBATCH --error=PMthreetophat3.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60Gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcb4325
#SBATCH --qos=mcb4325

module load ufrc
module load bowtie2
module load tophat

#bowtie2-build hg38.fa hg38

tophat2 -p 4 --no-coverage-search hg38 SRR445721.fastq

