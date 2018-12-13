#!/bin/bash

#SBATCH --job-name=PMthreecl.sh
#SBATCH --time=12:00:00
#SBATCH --mail-user=nevinsc@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=PMthreecl.log
#SBATCH --error=PMthreecl.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcb4325
#SBATCH --qos=mcb4325
#SBATCH --partition=hpg2-dev

module load cufflinks
module load samtools

#samtools sort -n /ufrc/strauss/nevinsc/PMthree/sample3/accepted_hits.sorted.bam -o /ufrc/strauss/nevinsc/PMthree/sample3/accepted_hits.sorted.n
cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample2_g_output/ /ufrc/strauss/nevinsc/PMthree/sample2/accepted_hits.sorted.n

#cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample3_g_output/ /ufrc/strauss/nevinsc/PMthree/sample3/accepted_hits.sorted.bam

#cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample4_g_output/ /ufrc/strauss/nevinsc/PMthree/sample4/accepted_hits.sorted.bam

#cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample5_g_output/ /ufrc/strauss/nevinsc/PMthree/sample5/accepted_hits.sorted.bam

#cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample6_g_output/ /ufrc/strauss/nevinsc/PMthree/sample6/accepted_hits.sorted.bam

#cufflinks -g Homo_sapiens.GRCh38.86.OK.gtf -p 8 -o sample2_g_output/ /ufrc/strauss/nevinsc/PMthree/sample2/accepted_hits.sorted.bam



