#!/bin/bash

n=41;

for i in $(cat /nfs/home/students/rahelcaspar/MA/Samples_Tabellen/whippet_samples_new_data.csv) 
do 
	echo "${i}_S${n}";
	julia Whippet.jl/bin/whippet-quant.jl /nfs/data3/gr_splicing/new_data/Sample_${i}/${i}_S${n}_L004_R1_001.fastq.gz /nfs/data3/gr_splicing/new_data/Sample_${i}/${i}_S${n}_L004_R2_001.fastq.gz -o /nfs/proj/gr_splicing/whippet-files/output_quant/${i} --biascorrect 
	((n++));

done;




