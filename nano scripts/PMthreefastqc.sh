#!/bin/bash

#SBATCH --job-name=PMthreefastqc.sh
#SBATCH --time=12:00:00
#SBATCH --mail-user=nevinsc@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=PMthreefastqc.log
#SBATCH --error=PMthreefastqc.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcb4325
#SBATCH --qos=mcb4325
#SBATCH --partition=hpg2-dev

module load ufrc
module load fastqc

#gunzip /ufrc/mcb4325/share/ModuleIII/SRR445721.fastq.gz
#gunzip /ufrc/mcb4325/share/ModuleIII/SRR445722.fastq.gz
#gunzip /ufrc/mcb4325/share/ModuleIII/SRR445723.fastq.gz
#gunzip /ufrc/mcb4325/share/ModuleIII/SRR440985.fastq.gz
#gunzip /ufrc/mcb4325/share/ModuleIII/SRR440986.fastq.gz
#gunzip /ufrc/mcb4325/share/ModuleIII/SRR440987.fastq.gz

#fastqc SRR445721.fastq
fastqc SRR445722.fastq
fastqc SRR445723.fastq
fastqc SRR490985.fastq
fastqc SRR490986.fastq
fastqc SRR490987.fastq

#fastqc /ufrc/mcb4325/share/ModuleIII/SRR445721.fastq
#fastqc /ufrc/mcb4325/share/ModuleIII/SRR445722.fastq .
#fastqc /ufrc/mcb4325/share/ModuleIII/SRR445723.fastq
#fastqc /ufrc/mcb4325/share/ModuleIII/SRR490985.fastq
#fastqc /ufrc/mcb4325/share/ModuleIII/SRR490986.fastq
#fastqc /ufrc/mcb4325/share/ModuleIII/SRR490987.fastq



