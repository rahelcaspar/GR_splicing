#!/bin/bash


d="/nfs/proj/gr_splicing/whippet-files/output_quant/public_dex/";

#a = dex, b = control

julia Whippet.jl/bin/whippet-delta.jl -a ${d}SRR9969324.psi.gz,${d}SRR9969325.psi.gz,${d}SRR9969327.psi.gz,${d}SRR9969328.psi.gz,${d}SRR9969329.psi.gz -b ${d}SRR9969314.psi.gz,${d}SRR9969315.psi.gz,${d}SRR9969316.psi.gz,${d}SRR9969317.psi.gz -o /nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta_public_dex -d /nfs/proj/gr_splicing/whippet-files/output_quant/public_dex


