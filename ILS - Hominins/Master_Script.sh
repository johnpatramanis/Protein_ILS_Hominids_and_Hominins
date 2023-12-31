#!/bin/bash
#### Pre install both conda environments using conda -f

### How many cores for snakemake to use
CORES=64



###### load first conda env, sample individuals, prepare datasets for trees

eval "$(conda shell.bash hook)"
conda activate Entropy

snakemake -F -j$CORES

conda deactivate




###### load second conda env, infere trees through Paleoprophyler


eval "$(conda shell.bash hook)"
conda activate Analyser

# rm -rf Dataset_Analysis/Workspace/2_DATASETS/*
# rm -rf Protein_Samples/*

cd Dataset_Analysis
snakemake -F -j$CORES
cd ..
conda deactivate






#### load third conda env Comapare hominid tree with output tree


eval "$(conda shell.bash hook)"
conda activate Entropy
cd Tree_Comparisons
snakemake -F -j$CORES
cd ..
conda deactivate