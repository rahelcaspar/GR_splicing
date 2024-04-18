#!/bin/bash


d="/nfs/proj/gr_splicing/whippet-files/output_quant/public_ko_wt/";

#a = ko, b = wt

julia Whippet.jl/bin/whippet-delta.jl -a ${d}SRR9969249.psi.gz,${d}SRR9969250.psi.gz,${d}SRR9969251.psi.gz,${d}SRR9969261.psi.gz,${d}SRR9969262.psi.gz,${d}SRR9969263.psi.gz,${d}SRR9969271.psi.gz,${d}SRR9969272.psi.gz,${d}SRR9969281.psi.gz,${d}SRR9969282.psi.gz,${d}SRR9969283.psi.gz,${d}SRR9969293.psi.gz,${d}SRR9969294.psi.gz,${d}SRR9969304.psi.gz,${d}SRR9969305.psi.gz,${d}SRR9969306.psi.gz -b ${d}SRR9969252.psi.gz,${d}SRR9969253.psi.gz,${d}SRR9969254.psi.gz,${d}SRR9969264.psi.gz,${d}SRR9969265.psi.gz,${d}SRR9969266.psi.gz,${d}SRR9969273.psi.gz,${d}SRR9969274.psi.gz,${d}SRR9969284.psi.gz,${d}SRR9969285.psi.gz,${d}SRR9969286.psi.gz,${d}SRR9969295.psi.gz,${d}SRR9969296.psi.gz,${d}SRR9969297.psi.gz,${d}SRR9969307.psi.gz,${d}SRR9969308.psi.gz,${d}SRR9969309.psi.gz -o /nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta_ko_wt -d /nfs/proj/gr_splicing/whippet-files/output_quant/public_ko_wt


