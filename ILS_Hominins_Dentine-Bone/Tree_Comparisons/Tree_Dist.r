args = commandArgs(trailingOnly=TRUE)
library('TreeDist')
library("TreeTools", quietly = TRUE)


### Run with Rguments: Rscript ../Newick_Files/EXAMP1/Newick_File.txt ../Dataset_Analysis/Workspace/2_DATASETS/EXAMP1_Dataset/CONCATINATED/CONCATINATED_aln_e.phy_phyml_tree.txt ./EXAMP1/Tree_Distances.txt

### Get paths for files

INITIAL_TREE_PATH=args[1]
INFERED_TREE_PATH=args[2]
OUTPUT_FILES=args[3]


###################### Import Trees ####################

Tree1=ape::read.tree(INITIAL_TREE_PATH) #### Input Tree
Tree2=ape::read.tree(INFERED_TREE_PATH) #### Infered Tree

#### Root Trees at Chimpansee
Tree1=RootTree(Tree1, 'CHI')
Tree2=RootTree(Tree2, 'CHI')


print('Trees loaded')	
##################### Calculate Distances #################


##### Normal Tree Distances	
distance_Generic <- TreeDistance(Tree1, Tree2)
distance_Nye <- NyeSimilarity(Tree1, Tree2, normalize = TRUE)
distance_JRF <- JaccardRobinsonFoulds(Tree1,Tree2)
distance_Bod_Giaro <- MatchingSplitDistance(Tree1,Tree2)
distance_KC <- KendallColijn(Tree1,Tree2)

### Node Labels, get bootstrap support for tree (1 number between 0-100)
BSSupport=Tree2$node.label[[2]] 

Case1 = ape::is.monophyletic(Tree2,c('NEA','DEN'))
Case2 = ape::is.monophyletic(Tree2,c('HUM','NEA'))
Case3 = ape::is.monophyletic(Tree2,c('HUM','DEN'))

if (Case1==TRUE){
Tree_Topology='T1'
}

if (Case2==TRUE){
Tree_Topology='T2'
}

if (Case3==TRUE){
Tree_Topology='T3'
}

if (Case1==FALSE & Case2==FALSE & Case3==FALSE){
Tree_Topology='T4'
}



##### For comparisons of tree with not - equal number of leaves
Clustering_Entropy1=ClusteringEntropy(Tree1)
Clustering_Entropy2=ClusteringEntropy(Tree2)
Mutual_Clustering_Info=MutualClusteringInfo(Tree1, Tree2)



#### Output Distances into file
METRICS=c(distance_Generic,distance_Nye,distance_JRF,distance_Bod_Giaro,distance_KC,Mutual_Clustering_Info,Tree_Topology,BSSupport)
write(paste(METRICS,collapse='\t'), file = OUTPUT_FILES)
