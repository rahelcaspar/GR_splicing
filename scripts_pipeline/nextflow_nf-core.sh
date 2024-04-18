#!/bin/bash


nextflow run nf-core/rnaseq --input /nfs/home/students/rahelcaspar/MA/samplesheet.csv --outdir /nfs/proj/gr_splicing --email rahel@the-caspar-team.de -profile daisybio,apptainer --fasta /nfs/data/references/ensembl110_GRCm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa --gtf /nfs/data/references/ensembl110_GRCm39/Mus_musculus.GRCm39.110.gtf --aligner star_rsem --save_reference --skip_pseudo_alignment -resume


