import os
import random
import sys

from Bio import SeqIO

#### LOAD FILE
FASTA_FILE=sys.argv[1]
TO_BE_REMOVED=sys.argv[2:]

fasta_sequences = SeqIO.parse(open(FASTA_FILE),'fasta')

print(fasta_sequences)
print(TO_BE_REMOVED)


### Organise data, exclude desired samples
FASTA_NAMES=[]
FASTA_SEQUENCES=[]
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    
    if name not in TO_BE_REMOVED:
    
        FASTA_NAMES.append(name)
        FASTA_SEQUENCES.append(sequence)
    







#### Look through alignment for variant sites

VARIANTS=0  
for POSITION_NUMBER in range(0,len(FASTA_SEQUENCES[0])):
    AA_HERE=[]
    for SAMPLE_NUMBER in range(0,len(FASTA_SEQUENCES)):
        AA_HERE.append(FASTA_SEQUENCES[SAMPLE_NUMBER][POSITION_NUMBER])

    
    AA_HERE=[x for x in AA_HERE if (x!='-' and x!='?' and x!='X')]
    POSSIBLE_VARIANT=set(AA_HERE)
    if len(POSSIBLE_VARIANT)>=2:
        VARIANTS+=1
        
print(VARIANTS)

OUTPUT_NAME='/'.join(FASTA_FILE.split('/')[0:(len(FASTA_FILE.split('/'))-1)])
OUTPUT_NAME=OUTPUT_NAME+'/Number_of_Variants.txt'
OUTPUT_FILE=open(OUTPUT_NAME,'w')
OUTPUT_FILE.write(str(VARIANTS))

