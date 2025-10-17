import os
########### File names
Prot_Loc_File=open('Protein_Location_Sorted_Lift_Over.txt','r')
Introgressed_Segments_Folders=["./Vindija/"]
####### Introgressed_Segments_Folders=["./Vindija/","./Chagyrskaya/","./Altai/"] ### For checking all three archaic genomes


Output_File=open("Introgressed_Proteins.txt",'w')



###############

Probability_Cutoff=0.00

########### Load Protein Locations
Prot_to_Loc={}
Proteins=[]
for line in Prot_Loc_File:
    line = line.strip().split()
    
    Protein = line[0]
    Chromosome = line[1]
    Start = int(line[2])
    End = int(line[3])
    if Start > End: ## make sure start < end, flip them if not
        TEMP=Start
        Start=End
        End=TEMP
    
    
    
    Proteins.append(Protein)
    
    
    Prot_to_Loc[Protein]=[Chromosome ,Start ,End]


#### Go through each Neanderthal Genome
for Introgressed_Segments_Folder in Introgressed_Segments_Folders:
    
    
    
    ########### Go through folder, population file, by population file

    Population_Files=[name for name in os.listdir(Introgressed_Segments_Folder)]

    for POP_FILE in Population_Files:
        Opened_POP_FILE=open(Introgressed_Segments_Folder+POP_FILE,'r')
        Opened_POP_FILE.readline() #### ignore first line
        
    ######### Go through line by line, checking each segment
        for line in Opened_POP_FILE:
            ######### Get relevant info of segment
            line = line.strip().split()
            ID = line[0]
            Chrom = line[1]
            Start = int(line[2])
            End = int(line[3])
            Prob = float(line[4])
            ArchaichHum = Introgressed_Segments_Folder.split("/")[1]
            Hum_Pop = POP_FILE.split(".")[0]
        
            
            if Start > End: ## make sure start < end, flip them if not
                TEMP=Start
                Start=End
                End=TEMP
        
        
        ##### Go through all Enamel Proteins and check if they overlap
            for Protein in Proteins:
                
                Prot_Chrom = Prot_to_Loc[Protein][0]
                Prot_Start = Prot_to_Loc[Protein][1]
                Prot_End = Prot_to_Loc[Protein][2]
                IntroGr_Trigger=0
                
                if (Chrom == Prot_Chrom) and (Prob >=Probability_Cutoff): #### if right chromosome, check
                    
                    if ( (Start <= Prot_Start) and (End >= Prot_End) ): #### if protein segment falls withing inrogressed segment
                        IntroGr_Trigger+=1
                        
                    if ((Prot_Start >= Start) and (Prot_Start <= End )): ##### if start of protein is within bounds of introgressed segment
                        IntroGr_Trigger+=1
                    
                    if ((Prot_End >= Start) and (Prot_End <= End )): ##### if end of protei  is within bound of introgressed segment
                        IntroGr_Trigger+=1
        
            
                if IntroGr_Trigger>0:
                    print(f"Woah! Looks like {Protein} is within Introgressed Segment at Chromosome {Chrom}, {Start}:{End}")
                    Output_File.write(F"{Protein}  {Prot_Chrom}  {Prot_Start}  {Prot_End}  {Chrom}  {Start}  {End}  {Hum_Pop}  {ArchaichHum}  {Prob}\n")       