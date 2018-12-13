#!/bin/bash

#SBATCH --job-name=fullhtseq.sh
#SBATCH --time=12:00:00
#SBATCH --mail-user=nevinsc@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --output=fullhtseq.log
#SBATCH --error=fullhtseq.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcb4325
#SBATCH --qos=mcb4325
#SBATCH --partition=hpg2-dev

module load ufrc
module load htseq
module load samtools

#samtools sort -n /ufrc/strauss/nevinsc/MCBWeek14/accepted_hits2.bam -o accepted_hits2SBName.bam

#htseq-count -f bam -s no -i gene_id -m union accepted_hits.sorted.bam Homo_sapiens.GRCh38.86.OK.gtf > ref_counts_union2.txt
htseq-count -f bam -s no -i gene_id -m union /ufrc/strauss/nevinsc/PMthree/sample3/accepted_hits.sorted.bam Homo_sapiens.GRCh38.86.OK.gtf > ref_counts_union3.txt
htseq-count -f bam -s no -i gene_id -m union /ufrc/strauss/nevinsc/PMthree/sample4/accepted_hits.sorted.bam Homo_sapiens.GRCh38.86.OK.gtf > ref_counts_union4.txt
htseq-count -f bam -s no -i gene_id -m union /ufrc/strauss/nevinsc/PMthree/sample5/accepted_hits.sorted.bam Homo_sapiens.GRCh38.86.OK.gtf > ref_counts_union5.txt
htseq-count -f bam -s no -i gene_id -m union /ufrc/strauss/nevinsc/PMthree/sample6/accepted_hits.sorted.bam Homo_sapiens.GRCh38.86.OK.gtf > ref_counts_union6.txt

