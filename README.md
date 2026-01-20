# Assessing the potential of ancient protein sequences in the study of hominid evolution
<br>

What is in the repository:
This repository contains scripts to reproduce some of the analysis and results of the manuscript ["Assessing the potential of ancient protein sequences in the study of hominid evolution"](https://www.biorxiv.org/content/10.1101/2025.04.08.647730v2.abstract).
Specifically two main analyses, can be reproduced following the step-by-step commands described in the Supplementary materials.

<br>

**1) Files corresponding to "Iterative Phylogenetic Analyses" section of methods**

In this analysis, we investigated how the concatenation of different numbers and different combinations of 12 enamel and collagen type I proteins might affect the topology of the inferred consensus  tree. In each iteration of this analysis, we carry out a concatenation using a subset of proteins sampled from the full set of proteins, reflecting the fact that not all proteins in the full set might be available in practice. The subset ranges in size from 1 (a single protein recovered) to all proteins recovered (either 12 or 28, depending on the tested dataset). One representative individual per population or species is randomly chosen and included in the alignment, as ancient protein studies are often limited to single individuals that are made to represent an entire species. For each concatenation, we build a phylogenetic tree and record the resulting topology. We then compare it to the underlying population tree, as inferred from past DNA studies. In total, we do this over 1,000 iterations per each number of proteins, sampling different sets of proteins and different representative individuals, in each turn.

We performed the same iterative analysis on each of the following 6 datasets (folders within this repository):

**a)** **ILS_Hominids**: For this version of the analysis we utilised a dataset of 12 enamel and collagen type I proteins from 4 species: _Homo sapiens_, _Pan troglodytes_, _Gorilla gorilla_ and _Pongo abelii_ (outgroup). The list of proteins is as follows:  AHSG, ALB, AMBN, AMELX, AMELY, AMTN, COL17A1, ENAM, MMP20, ODAM (enamel) COL1A1, COL1A2 (collagen type I). <br>
**b)** **ILS_Hominids_Only_Enamel**: Same analysis as above, but this time the dataset was limited to only the 10 enamel proteins. <br>
**c)** **ILS_Hominins**: For this version of the analysis we utilised a dataset, consisting of _Homo sapiens_, Neanderthals, Denisovans and _Pan troglodytes_ (outgroup). <br>
**d)** **ILS_Hominins_Dentine-Bone**: This dataset is the same as above **c)**, but expanded in the proteins used, to include 20 proteins that are most often recovered from dentin or bone tissue: COL1A1, COL1A2, COL2A1, COL3A1, COL4A4, AHSG, COL5A2, ALB, BGN, COL5A3, COL5A1, CHAD, COL22A1, COL11A2, SERPINF1, F2, COL11A1, LUM, COL12A1, POSTN. Four of these 20 proteins (COL1A1, COL1A2, AHSG and ALB) were already included in the original 12 proteins, leading to a final combined dataset of 28 proteins. <br>
**e)** **ILS_Hominins_Dentine-Bone_Only_Africans**:  This dataset is the same as  **d)**, but by sampling only from four present day human populations: Yoruba, Mende, Luhya and Mandinka. These African populations consist of mostly un-admixed representatives of _Homo sapiens_. <br>
**f)** **ILS_Hominins_Only Africans**: This dataset is the same as  **c)**, but by sampling only from four present day human populations: Yoruba, Mende, Luhya and Mandinka. These African populations consist of mostly un-admixed representatives of _Homo sapiens_. <br>

Within each folder is a bash file named **Master_Script.sh** , which when executed will attempt to activate the necessary conda environmnets (see bottom of this page) and then run 3 snakemake scripts. Each [snakemake](https://snakemake.readthedocs.io/en/stable/) script is tasked with performing a part of the Iterative analysis (asssembling the phylogenetic dataset, running the phylogenetic analysis on it and recording the resulting topology).

<br>
<br>

**2) Files corresponding to "Introgression Investigation" section of methods**

We assessed the impact of admixture, as a contributor to apparent tree discordance in the protein sequences of _Homo sapiens_, Neanderthals and Denisovans. We first identified how often the proteins under investigation here can be found within archaic-introgressed regions of present-day human genomes. We used previously reported archaic haplotypes found within two present-day human datasets (Skov et al. 2018, Chen et al. 2020) to assess this. The details of our methodology can be found in the Supplementary Material-S4, which is also attached here as pdf. The analysis can be repeated by following the steps described in the supplementary, in the folder below.

**Introgression_Investigation**: This folder contains two main subfolders - **Skov_et_al_2018** and **Chen_et_al_2020**. Each of these two folders corresponds to a dataset/publication in which archaic segments were inferred on the genomes of modern humans. The data from the publications are not provided here, as they are large in size (a couple of gigabytes), but code and links to download them are available in the Supplementary PDF. Inside each folder, are python scripts (**Check_Overlap.py** & **Report_on_found_introgressed_segments.py**) and some parameter files (**Protein_Location_Sorted.txt**) that can be used to inspect whether a collection of protein-coding genes are overlapping with introgressed regions of modern human genomes. After downloading the published data, running **Check_Overlap.py** will produce the **Introgressed_Proteins.txt** (output example can be found in each dataset folder), which can then be used by **Report_on_found_introgressed_segments.py** to produce a small report in the form of a .txt file (**Report_on_Introgressed_Proteins.txt**). The user can also modify the **Protein_Location_Sorted.txt** file to search for alternative genes than the ones investigated here. Lastly, for the dataset of **Skov_et_al_2018**, the results can also be plotted inside the **Plotting** folder, by first downloading some additional data (see Supplementary) and then running the pythnon script **Get_data.py**, followed by the R script **plot.R**. This should generate a plot similar to the Main text figure 7.

<br>
<br>

**Prerequisites**
<br>

A Linux machine (or server access to one) with all the necessary packages installed. The exact version of all pre-requisites can be seen within the two YML files in the repository (Entropy.yml and Analyser.yml), but all packages can be easily installed through [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). If conda is installed and running on your machine, simply use the following commands to download and prepare all pre-requisites:


```conda create -n Analyser  -c bioconda -c conda-forge snakemake phyml mafft trimal bioconductor-shortread r-stringr r-data.table r-phyclust seqmagick```

```conda create -n Entropy -c conda-forge -c bioconda biopython r-bio3d snakemake biopython```

<br>


**For the full description and instructions on the reproduction of any analysis, please see the Supplementary materials PDF file (which can be found in this repository here).**

