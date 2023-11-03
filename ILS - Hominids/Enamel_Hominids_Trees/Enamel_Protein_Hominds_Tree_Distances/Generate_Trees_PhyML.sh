#!/bin/bash
for GENE in AHSG ALB AMBN AMELX AMELY AMTN COL17A1 COL1A1 COL1A2 ENAM MMP20 ODAM
do

        phyml -i $GENE/$GENE_aln_e.phy -b 100 -m Blosum62 -a e -s BEST -v e -o tlr -f m --rand_start --n_rand_starts 4  --print_site_lnl --print_trace --no_memory_check;

done;