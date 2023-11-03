library('TreeDist')



### names of Genes
GENES=c("AHSG","ALB","AMBN","AMELX","AMELY","AMTN","COL1A1","COL1A2","COL17A1","ENAM","MMP20","ODAM")


for (GENE in GENES){
	INITIAL_TREE_PATH="Hominid_Tree.txt"
	INFERED_TREE_PATH=paste0('GENE_TREES_ENAMEL/Hominid_',GENE,".phy_phyml_tree.txt",sep="")
	OUTPUT_FILES=paste0('GENE_TREES_ENAMEL/Hominid_',GENE,"_Tree_Distances.txt",sep="")


	###################### Import Trees ####################

	Tree1=ape::read.tree(INITIAL_TREE_PATH) #### Input Tree	
	Tree2=ape::read.tree(INFERED_TREE_PATH) #### Infered Tree
	



	
	print('Trees loaded')	
	##################### Calculate Distances #################


	##### Normal Tree Distances	
	distance_Generic <- TreeDistance(Tree1, Tree2)
	distance_Nye <- NyeSimilarity(Tree1, Tree2, normalize = TRUE)
	distance_JRF <- JaccardRobinsonFoulds(Tree1,Tree2)
	distance_Bod_Giaro <- MatchingSplitDistance(Tree1,Tree2)
	distance_KC <- KendallColijn(Tree1,Tree2)




	##### For comparisons of tree with not - equal number of leaves
	Clustering_Entropy1=ClusteringEntropy(Tree1)
	Clustering_Entropy2=ClusteringEntropy(Tree2)
	Mutual_Clustering_Info=MutualClusteringInfo(Tree1, Tree2)



	#### Output Distances into file
	METRICS=c(distance_Generic,distance_Nye,distance_JRF,distance_Bod_Giaro,distance_KC,Mutual_Clustering_Info)
	file.create(OUTPUT_FILES)
	write(paste(METRICS,collapse='\t'), file = OUTPUT_FILES,append=FALSE)
}