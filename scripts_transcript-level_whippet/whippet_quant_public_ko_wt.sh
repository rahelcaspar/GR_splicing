#!/bin/bash


for i in $(cat /nfs/home/students/rahelcaspar/MA/Samples_Tabellen/whippet_samples_public_ko_wt.csv) 
do 
	echo "${i}";

	R1=/nfs/data3/gr_splicing/GEO_FASTQ/${i}_1.fastq.gz;
	R2=/nfs/data3/gr_splicing/GEO_FASTQ/${i}_2.fastq.gz;

	julia Whippet.jl/bin/whippet-quant.jl \
		<(zcat $R1 | awk -v count="0" '++count==3{{$0="+";count=-1}} 1') \
		<(zcat $R2 | awk -v count="0" '++count==3{{$0="+";count=-1}} 1') \
		-o /nfs/proj/gr_splicing/whippet-files/output_quant/public_ko_wt/${i}

done;




