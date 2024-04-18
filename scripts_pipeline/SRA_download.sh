#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=SRA_download.txt
#SBATCH --job-name=SRA_download


srun /nfs/home/students/rahelcaspar/sratoolkit.3.0.7-ubuntu64/bin/prefetch --option-file /nfs/home/students/rahelcaspar/MA/SRR_Acc_List.txt -O /nfs/data3/gr_splicing


for i in $(cat /nfs/home/students/rahelcaspar/MA/SRR_Acc_List.txt) ; do
	echo ${i};
	srun /nfs/home/students/rahelcaspar/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump --split-files /nfs/data3/gr_splicing/${i} -O /nfs/data3/gr_splicing/GEO_FASTQ --temp /nfs/data3/gr_splicing/GEO_FASTQ

done;


