#!/bin/bash


d="/nfs/proj/gr_splicing/whippet-files/output_quant/new_data/";

#a = dex; b = control

julia Whippet.jl/bin/whippet-delta.jl -a ${d}21L006487.psi.gz,${d}21L006489.psi.gz,${d}21L006491.psi.gz -b ${d}21L006486.psi.gz,${d}21L006488.psi.gz,${d}21L006490.psi.gz -o /nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta_new_data -d /nfs/proj/gr_splicing/whippet-files/output_quant/new_data


