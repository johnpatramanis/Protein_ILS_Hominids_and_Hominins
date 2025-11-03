# Assessing the potential of ancient protein sequences in the study of hominid evolution

What is in the repository:
This repository contains scripts to reproduce some of the analysis and results of the manuscript "Assessing the potential of ancient protein sequences in the study of hominid evolution".
Specifically two main analyses, can be reprodecuded following the commands described in the Supplementary materials.

1) The ``Iterative Phylogenetic Analyses"
We investigated how the concatenation of different numbers and different combinations of the 12 proteins might affect the topology of the inferred “consensus”  tree, which is often taken as an estimate of the population or species tree. For this analysis we utilised a ``hominid dataset'', consisting of \textit{Homo sapiens}, \textit{Pan troglodytes}, \textit{Gorilla gorilla} and \textit{Pongo abelii} (as an outgroup) and a second ``hominin dataset'' consisting of \textit{Homo sapiens}, Neanderthals, Denisovans and \textit{Pan troglodytes} (outgroup). To assess how the recovery of additional proteins, from different tissues affects phylogenetic analyses, we expanded the ``hominin dataset'', creating a third ``bone-dentin dataset''. This dataset consisted in the protein sequences that are most often recovered from dentin or bone tissue. In choosing which proteins to include in this anaysis, we utilised the list provided by Ruther et al. 2022 \cite{ruther2022spin}, which includes 20 proteins utilised in species identification: COL1A1, COL1A2, COL2A1, COL3A1, COL4A4, AHSG, COL5A2, ALB, BGN, COL5A3, COL5A1, CHAD, COL22A1, COL11A2, SERPINF1, F2, COL11A1, LUM, COL12A1, POSTN. Four of these 20 proteins (COL1A1, COL1A2, AHSG and ALB) were already included in the original 12 proteins, leading to a final combined dataset of 28 proteins. 

In each iteration of this analysis, we carry out a concatenation using a subset of proteins sampled from the full set of proteins, reflecting the fact that not all proteins in the full set might be available in practice. The subset ranges in size from 1 (a single protein recovered) to all proteins recovered (either 12 or 28, depending on the tested dataset). One representative individual per population or species is randomly chosen and included in the alignment, as ancient protein studies are often limited to single individuals that are made to represent an entire species. For each concatenation, we build a phylogenetic tree and record the resulting topology. We then compare it to the underlying population tree, as inferred from past DNA studies. In total, we do this over 1,000 iterations per each number of proteins, sampling different sets of proteins and different representative individuals, in each turn.

We performed the same iterative analysis on each of the following 6 datasets (folders):
a) ILS_Hominids
b) ILS_Hominids_Only_Enamel
c) ILS_Hominins
d) ILS_Hominins_Dentine-Bone
e) ILS_Hominins_Dentine-Bone_Only_Africans
f) ILS_Hominins_Only Africans



3) The ``Introgression Investigation"






**For the full description and instruction on the reproduction of the analysis, please see the Supplementary materials PDF file of the same manuscript (also in this repository here).**


