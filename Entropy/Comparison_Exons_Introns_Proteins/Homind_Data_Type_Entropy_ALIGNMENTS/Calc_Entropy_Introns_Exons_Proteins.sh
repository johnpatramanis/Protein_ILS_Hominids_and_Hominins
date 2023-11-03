#!/bin/bash


for GENE in AHSG ALB AMBN AMELX AMELY AMTN COL17A1 COL1A1 COL1A2 ENAM MMP20 ODAM
do

	for TYPE in INTRONS EXONS PROTEIN
	do
		>"$GENE"_"$TYPE".entr 
		Rscript Calc_Entropy_of_Fasta_IEP.r "$GENE"_"$TYPE".fa 
		
	done;
done;

cp ./*.entr ../Entropy_Measurements/