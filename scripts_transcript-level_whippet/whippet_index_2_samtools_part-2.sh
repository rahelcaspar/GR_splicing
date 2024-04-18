#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=whippet_index_samtools_part-2.txt
#SBATCH --job-name=samtools_part2

echo sort
srun samtools sort -o /nfs/proj/gr_splicing/whippet-files/whippet_index_samples.sort.bam /nfs/proj/gr_splicing/whippet-files/whippet_index_samples.bam -m 2G

echo rmdup
srun samtools rmdup -S /nfs/proj/gr_splicing/whippet-files/whippet_index_samples.sort.bam /nfs/proj/gr_splicing/whippet-files/whippet_index_samples.sort.rmdup.bam

echo index
samtools index /nfs/proj/gr_splicing/whippet-files/whippet_index_samples.sort.rmdup.bam



