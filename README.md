# Assessing the potential of ancient protein sequences in the study of hominid evolution

What is in the repository:
This repository contains scripts to reproduce some of the analysis and results of the manuscript "Assessing the potential of ancient protein sequences in the study of hominid evolution".
Specifically two main analyses, can be reprodecuded following the commands described in the Supplementary materials.


**1) The "Iterative Phylogenetic Analyses"**

In this analysis we investigated how the concatenation of different numbers and different combinations of 12 enamel and collagen type I proteins might affect the topology of the inferred consensus  tree. In each iteration of this analysis, we carry out a concatenation using a subset of proteins sampled from the full set of proteins, reflecting the fact that not all proteins in the full set might be available in practice. The subset ranges in size from 1 (a single protein recovered) to all proteins recovered (either 12 or 28, depending on the tested dataset). One representative individual per population or species is randomly chosen and included in the alignment, as ancient protein studies are often limited to single individuals that are made to represent an entire species. For each concatenation, we build a phylogenetic tree and record the resulting topology. We then compare it to the underlying population tree, as inferred from past DNA studies. In total, we do this over 1,000 iterations per each number of proteins, sampling different sets of proteins and different representative individuals, in each turn.

We performed the same iterative analysis on each of the following 6 datasets (folders within this repository):

**a)** ILS_Hominids: For this version of the analysis we utilised a dataset of 12 enamel and collagen type I proteins from 4 species: _Homo sapiens_, _Pan troglodytes_, _Gorilla gorilla_ and _Pongo abelii_ (outgroup). The list of proteins is the following:  AHSG, ALB, AMBN, AMELX, AMELY, AMTN, COL17A1, ENAM, MMP20, ODAM (enamel) COL1A1, COL1A2 (collagen type I). 
**b)** ILS_Hominids_Only_Enamel: Same analysis as above, but this time the dataset was limited to only the 10 enamel proteins.
**c)** ILS_Hominins: For this version of the analysis we utilised a dataset, consisting of _Homo sapiens_, Neanderthals, Denisovans and _Pan troglodytes_ (outgroup).
**d)** ILS_Hominins_Dentine-Bone: This dataset is the same as above (**c**), but expanded in the protein used, to incude 20 proteins that are most often recovered from dentin or bone tissue: COL1A1, COL1A2, COL2A1, COL3A1, COL4A4, AHSG, COL5A2, ALB, BGN, COL5A3, COL5A1, CHAD, COL22A1, COL11A2, SERPINF1, F2, COL11A1, LUM, COL12A1, POSTN. Four of these 20 proteins (COL1A1, COL1A2, AHSG and ALB) were already included in the original 12 proteins, leading to a final combined dataset of 28 proteins. 
**e)** ILS_Hominins_Dentine-Bone_Only_Africans:  This dataset is the same as  **d**, but by sampling only from four present day human populations: Yoruba, Mende, Luhya and Mandinka. These African populations consist of mostly un-admixed representatives of _Homo sapiens_. 
**f)** ILS_Hominins_Only Africans: This dataset is the same as  **c**, but by sampling only from four present day human populations: Yoruba, Mende, Luhya and Mandinka. These African populations consist of mostly un-admixed representatives of _Homo sapiens_. 


**2) The "Introgression Investigation"**






**For the full description and instruction on the reproduction of the analysis, please see the Supplementary materials PDF file of the same manuscript (also in this repository here).**


