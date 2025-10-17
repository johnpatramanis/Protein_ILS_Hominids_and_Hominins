### install.packages("bio3d", dependencies=TRUE)
### http://thegrantlab.org/bio3d/reference/entropy.html
### Run like this: R Calc_Entropy_of_Fasta.r A1BG_PROT_REFERENCE.fa PER_PROTEIN Entropies_Hominids

## Start, args ,packages

args = commandArgs(trailingOnly=TRUE)
library('bio3d')



###############################################################
###############################################################
######### Functions



##################
### Flat Entropy
Calculate_Entropy_Flat <- function(Fasta_alignment) {

	Protein_Flat_Entropy=0

	if (length(Fasta_alignment$ali)>0){
		
		Entropy_Metrics <- entropy(Fasta_alignment)
		Entropy_Pure <- Entropy_Metrics$H
		Protein_Flat_Entropy<- sum(Entropy_Pure)
		
	}


	print(Protein_Flat_Entropy)
}


##################
### Mean Entropy (divided by length)
Calculate_Entropy_Mean <- function(Fasta_alignment) {

	Protein_Mean_Entropy=0

	if (length(Fasta_alignment$ali)>0){
		
		Entropy_Metrics <- entropy(Fasta_alignment)
		Entropy_Pure <- Entropy_Metrics$H
		Protein_Mean_Entropy<-  (sum(Entropy_Pure)/length(Entropy_Pure))
		
	}


	print(Protein_Mean_Entropy)
}


##################
### Flat Entropy Normalised for Long Branches
Calculate_Entropy_For_Long_Branches_Flat <- function(Fasta_alignment) {

	Protein_Flat_Entropy_Long_Branch=0
	Temp_alignment=Fasta_alignment
	
	if (length(Fasta_alignment$ali)>0){
		
		### Create all possible pairs of aligned sequences
		Samples=1:length(Fasta_alignment$id)
		All_Pairs=t(combn(Samples,2))
		
		Average_Entropy=c() #empty to fill up
		
		for ( i in 1:dim(All_Pairs)[1] ){	
		
		
		
			Temp_alignment$ali=Fasta_alignment$ali[All_Pairs[i,],]
			Temp_Entropy_Metrics <- entropy(Temp_alignment)
			Temp_Entropy_Pure <- Temp_Entropy_Metrics$H    
			Temp_Flat_Protein_Flat_Entropy<- sum(Temp_Entropy_Pure)
			
			Average_Entropy=c(Average_Entropy,Temp_Flat_Protein_Flat_Entropy)
			
		}
		
		
		
		Protein_Flat_Entropy_Long_Branch= sum(Average_Entropy)/length(Average_Entropy)
	}


	print(Protein_Flat_Entropy_Long_Branch)
}



##################
### Mean Entropy Normalised for Long Branches
Calculate_Entropy_For_Long_Branches_Mean <- function(Fasta_alignment) {

	Protein_Flat_Entropy_Long_Branch=0
	Temp_alignment=Fasta_alignment
	
	if (length(Fasta_alignment$ali)>0){
		
		### Create all possible pairs of aligned sequences
		Samples=1:length(Fasta_alignment$id)
		All_Pairs=t(combn(Samples,2))
		
		Average_Entropy=c() #empty to fill up
		
		for ( i in 1:dim(All_Pairs)[1] ){	
		
		
		
			Temp_alignment$ali=Fasta_alignment$ali[All_Pairs[i,],]
			Temp_Entropy_Metrics <- entropy(Temp_alignment)
			Temp_Entropy_Pure <- Temp_Entropy_Metrics$H    
			Temp_Mean_Protein_Flat_Entropy<- ( sum(Temp_Entropy_Pure) / length(Temp_Entropy_Pure) )
			
			
			Average_Entropy=c(Average_Entropy,Temp_Mean_Protein_Flat_Entropy)
			
		}
		
		
		
		Protein_Flat_Entropy_Long_Branch= sum(Average_Entropy)/length(Average_Entropy)
	}


	print(Protein_Flat_Entropy_Long_Branch)
}

###############################################################
###############################################################












###############################################################
##### load file and names
Directory_Path <- args[2]
Fasta_alignment <- read.fasta(paste(args[2], args[1], sep="/"),rm.dup = FALSE)  ##### should not contain '?', replace with '-'
Name_of_Protein<- gsub('_PROT_REFERENCE.fa','',args[1])




### Output file
OUTPUT_FODLER_FULL_PATH <- args[3]
OUTPUT_FODLER_FULL_PATH <- paste(OUTPUT_FODLER_FULL_PATH,Name_of_Protein,sep='/')
OUTPUT_FODLER_FULL_PATH <- paste(OUTPUT_FODLER_FULL_PATH,'.entr',sep='')

OUTPUT_LENGTH <- args[3]
OUTPUT_LENGTH <- paste(OUTPUT_LENGTH,'Protein_Lengths.txt',sep='/')



#### 
print(Name_of_Protein)




#################################################################################################################################
##### Clean up columns with missing value
GAPS=which(Fasta_alignment$ali=='-',arr.ind=TRUE) ## Find columns with gaps
if (length(GAPS[,2])>0){
	Fasta_alignment$ali <- Fasta_alignment$ali[,-(GAPS[,2])] ## Remove them from alignment, use the rest for entropy
	}





####
#### Calculate Entropy etc

Protein_Flat_Entropy = Calculate_Entropy_Flat(Fasta_alignment)
Protein_Mean_Entropy = Calculate_Entropy_Mean(Fasta_alignment)
Protein_Flat_Entropy_Long_Branch = Calculate_Entropy_For_Long_Branches_Flat(Fasta_alignment)
Protein_Mean_Entropy_Long_Branch = Calculate_Entropy_For_Long_Branches_Mean(Fasta_alignment)




OUTPUT=paste0(Protein_Mean_Entropy,'\t',Protein_Flat_Entropy,'\t',Protein_Flat_Entropy_Long_Branch,'\t',Protein_Mean_Entropy_Long_Branch,'\n',sep='')
# OUTPUT=paste(Protein_Mean_Entropy)

##### Output
write(OUTPUT, file=OUTPUT_FODLER_FULL_PATH,append = TRUE)
