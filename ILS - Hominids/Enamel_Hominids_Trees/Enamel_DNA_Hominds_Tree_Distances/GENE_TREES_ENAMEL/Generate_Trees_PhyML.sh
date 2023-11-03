#!/bin/bash
for GENE in AHSG ALB AMBN AMELX AMELY AMTN COL17A1 COL1A1 COL1A2 ENAM MMP20 ODAM
do

        seqmagick convert --output-format phylip Hominid_$GENE.fa Hominid_$GENE.phy
        phyml -i Hominid_$GENE.phy -b 100 -m GTR -a e -s BEST -v e -o tlr -f m --rand_start --n_rand_starts 4  --print_site_lnl --print_trace --no_memory_check;
        sed  -i 's/HG00463.fi/HUM/' Hominid_$GENE.phy_phyml_tree.txt
        sed  -i 's/SAMEA10436/PON/' Hominid_$GENE.phy_phyml_tree.txt
        sed  -i 's/SAMN019204/GOR/' Hominid_$GENE.phy_phyml_tree.txt
        sed  -i 's/SAMN019205/CHI/' Hominid_$GENE.phy_phyml_tree.txt

done;