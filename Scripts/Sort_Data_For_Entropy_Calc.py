import os
import random
import sys

from Bio import SeqIO



FILES=os.listdir('./PER_PROTEIN/')
print(F'Output Fasta file will be placed in: {sys.argv[1]}')
path=sys.argv[1]

if (os.path.exists(path)==False):
    os.makedirs(path)


#### Check if only a subset of the fastas are desired
if (os.path.exists('Protein_Subset')==True):
    Protein_Subset=[]
    Protein_Subset_File=open('Protein_Subset','r')
    for line in Protein_Subset_File:
        line=line.strip()
        Protein_Subset.append(line)



    

#### Clean up if folder exists and is full of files 
PROTEIN_FILES=os.listdir(path)
for FILE in PROTEIN_FILES:
    if '.fa' in FILE:
        os.remove(path+'/'+FILE)







##### Go through original protein files
for FILE in FILES:
    if '.fa' in FILE:
        fasta_sequences = SeqIO.parse(open('./PER_PROTEIN/'+FILE),'fasta')
        PROTEIN_NAME=FILE.split('_PROT_REFERENCE.fa')[0]
        
        ##### How many of each genus/species in new fasta
        Human_Count=0
        Neand_Count=0
        Chimp_Count=0
        Gorilla_Count=0
        Pongo_Count=0
        
        
        
        ######## Shuffle sequence
        FASTA_NAMES=[]
        FASTA_SEQUENCES=[]
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            FASTA_NAMES.append(name)
            FASTA_SEQUENCES.append(sequence)
            
         
        arrangement=[x for x in range(0,len(FASTA_NAMES))]
        random.shuffle(arrangement)

        
        
        
        #### Go through fasta sequences
        for fasta in arrangement:
            check=False
            
            name = FASTA_NAMES[fasta]
            sequence = FASTA_SEQUENCES[fasta]
            
            name=name.split('/')[0]
            
            SAMPLE_NAME=name
          
            sequence=sequence.replace('X','-')
            sequence=sequence.replace('?','-')
            
            DASH_COUNT=sequence.count('-')

            if DASH_COUNT<len(sequence) and sequence!='' and sequence!=' ':
            
          
                if ( ('HG' in SAMPLE_NAME[0:2]) or ('NA' in SAMPLE_NAME[0:2]) or ('sapiens' in SAMPLE_NAME) or ('HUM' in SAMPLE_NAME) ) and (Human_Count==0):
                    check=True
                    Human_Count+=1

                # if ( ('Altai' in SAMPLE_NAME) or ('Vindija' in SAMPLE_NAME) or ('Chagyrskaya' in SAMPLE_NAME) or ('Denisova' in SAMPLE_NAME) ) and (Neand_Count==0):
                    # check=True
                    # Neand_Count+=1            
                
                if ( ('Pan' in SAMPLE_NAME) or ('pan' in SAMPLE_NAME) or ('CHI' in SAMPLE_NAME)  ) and (Chimp_Count==0):
                    check=True
                    Chimp_Count+=1    
                    
                if ( ('Gorilla' in SAMPLE_NAME) or ('gorilla' in SAMPLE_NAME) or ('GOR' in SAMPLE_NAME) ) and (Gorilla_Count==0):
                    check=True
                    Gorilla_Count+=1            

                if ( ('Pongo' in SAMPLE_NAME) or ('pongo' in SAMPLE_NAME) or ('PON' in SAMPLE_NAME) ) and (Pongo_Count==0):
                    check=True
                    Pongo_Count+=1

            if ( ('Protein_Subset' in locals()) or ('Protein_Subset' in globals())):
                if PROTEIN_NAME not in Protein_Subset:
                    check=False


            if check==True:
                
                
                PROTEIN_FASTA=open(F'./{path}/{PROTEIN_NAME}_PROT_REFERENCE.fa','a')
                PROTEIN_FASTA.write('>'+name+'\n')
                PROTEIN_FASTA.write(sequence+'\n')
                PROTEIN_FASTA.close()
