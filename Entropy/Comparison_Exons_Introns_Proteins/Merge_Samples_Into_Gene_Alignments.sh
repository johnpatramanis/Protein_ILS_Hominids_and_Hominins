#!/bin/bash
for GENE in AHSG ALB AMBN AMELX AMELY AMTN COL17A1 COL1A1 COL1A2 ENAM MMP20 ODAM;
do
        # >"$GENE"_INTRONS.fa;
        # >"$GENE"_EXONS.fa;
        # >"$GENE"_PROTEIN.fa;
        for SAMPLE in SAMN01920479_Gorilla_gorilla_gorilla      HG00463.final   SAMEA104361527_Pongo_abelii     SAMN01920517_Pan_troglodytes_ellioti;
        do
                cat "$SAMPLE"_"$GENE".fa >> "$GENE"_INTRONS.fa;
                cat "$SAMPLE"_"$GENE"_spliced.fa >> "$GENE"_EXONS.fa;
                cat "$SAMPLE"_"$GENE"_translated.fa >> "$GENE"_PROTEIN.fa;

        done;

done;
