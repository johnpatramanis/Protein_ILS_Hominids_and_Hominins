########### File names
Introgressed_File=open('Introgressed_Proteins.txt','r')

Output_File=open("Report_on_Introgressed_Proteins.txt",'w')


########### Load Introgressed Segs
Proteins={}

for line in Introgressed_File:
    line = line.strip().split()
    
    Protein = line[0]
    Human_Pop = line[7]
    Archaic_Pop = line[8]
    Prob = float(line[9])
    
    
    if Protein not in Proteins.keys():
        Proteins[Protein]=[1,[Human_Pop],[Archaic_Pop],[Prob]]
        
    else:
        Proteins[Protein][0]+=1
        Proteins[Protein][1].append(Human_Pop)
        Proteins[Protein][2].append(Archaic_Pop)
        Proteins[Protein][3].append(Prob)
    
    
for PRT in Proteins.keys():
    
    Proteins[PRT][1]=",".join(list(set(Proteins[PRT][1])))
    Proteins[PRT][2]=",".join(list(set(Proteins[PRT][2])))
    Proteins[PRT][3]=sum(Proteins[PRT][3])/len(Proteins[PRT][3])
    
    print(F"Protein {PRT}  was found in {Proteins[PRT][0]} Introgressed Segment, with an average probability of {Proteins[PRT][3]} in the human populations of {Proteins[PRT][1]}, matching the archaic populations of {Proteins[PRT][2]}.\n")
    Output_File.write(F"{PRT}\t{Proteins[PRT][0]} \t{Proteins[PRT][3]}\t{Proteins[PRT][1]}\t{Proteins[PRT][2]}\n")
