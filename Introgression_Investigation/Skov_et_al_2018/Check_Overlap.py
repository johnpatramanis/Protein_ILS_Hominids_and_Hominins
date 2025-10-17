########### File names
Prot_Loc_File=open('Protein_Location_Sorted.txt','r')
OneKGenomes_File=open('hg38_1000g_segments.txt','r')
HGDP_File=open('hg38_HGDP_segments.txt','r')

Output_File=open("Introgressed_Proteins.txt",'w')



###############

Probability_Cutoff=0.85

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



########### 1000 Genomes Scan
Labels_1kg=OneKGenomes_File.readline() ### ignore first line

########## Go through line by line, checking each segment
for line in OneKGenomes_File:
    ### Get relevant info of segment
    line = line.strip().split()
    Chrom = line[4]
    Chrom = Chrom.split("chr")[1]
    Start = int(line[5])
    End = int(line[6])
    Prob = float(line[7])
    ArchaichHum = str(line[8])
    Hum_Pop = str(line[2])
    
    
    if Start > End: ## make sure start < end, flip them if not
        TEMP=Start
        Start=End
        End=TEMP
    
    
    ###### Go through all Enamel Proteins and check if they overlap
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
            
            
            

########### HGDP Genomes Scan, Do exactly the same!
Labels_1kg=HGDP_File.readline() ### ignore first line

########## Go through line by line, checking each segment
for line in HGDP_File:
    ### Get relevant info of segment
    line = line.strip().split()
    Chrom = line[4]
    Chrom = Chrom.split("chr")[1]
    Start = int(line[5])
    End = int(line[6])
    Prob = float(line[7])
    ArchaichHum = str(line[8])
    Hum_Pop = str(line[2])
    
    
    if Start > End: ## make sure start < end, flip them if not
        TEMP=Start
        Start=End
        End=TEMP
    
    
    ###### Go through all Enamel Proteins and check if they overlap
    for Protein in Proteins:
        
        Prot_Chrom = Prot_to_Loc[Protein][0]
        Prot_Start = Prot_to_Loc[Protein][1]
        Prot_End = Prot_to_Loc[Protein][2]
        IntroGr_Trigger=0
        if (Chrom == Prot_Chrom) and (Prob >=0.75): #### if right chromosome, check
            
            if ( (Start <= Prot_Start) and (End >= Prot_End) ): #### if protein segment falls withing inrogressed segment
                IntroGr_Trigger+=1
                
            if ((Prot_Start >= Start) and (Prot_Start <= End )): ##### if start of protein is within bounds of introgressed segment
                IntroGr_Trigger+=1
            
            if ((Prot_End >= Start) and (Prot_End <= End )): ##### if end of protei  is within bound of introgressed segment
                IntroGr_Trigger+=1
    
        
        if IntroGr_Trigger>0:
            print(f"Woah! Looks like {Protein} is within Introgressed Segment at Chromosome {Chrom}, {Start}:{End}")
            Output_File.write(F"{Protein}  {Prot_Chrom}  {Prot_Start}  {Prot_End}  {Chrom}  {Start}  {End}  {Hum_Pop}  {ArchaichHum}  {Prob}\n")            