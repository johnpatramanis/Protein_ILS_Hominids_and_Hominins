import os
import random
import sys
import statistics

from Bio import SeqIO


Entropy_Folder=sys.argv[1]
OF=sys.argv[2]
Output_File=open(OF,'w')

FILES=os.listdir(Entropy_Folder)
FILES=[x for x in FILES if '.entr' in x]
PROTEIN_TO_ENTROPY={}
PROTEIN_TO_ENTROPY_FLAT={}
PROTEIN_TO_ENTROPY_FLAT_LONGB={}
PROTEIN_TO_ENTROPY_MEAN_LONGB={}

PROTEINS=[]

for FILE in FILES:
    PROTEIN_NAME=FILE.split('.entr')[0]
    PROTEINS.append(PROTEIN_NAME)
    PROTEIN_ENTROPY_FILE=open(F'{Entropy_Folder}/{FILE}','r')
    PROTEIN_NAME=FILE.split('.entr')[0]
    ENTROPY_HERE_MEAN=[]
    ENTROPY_HERE_FLAT=[]
    ENTROPY_HERE_FLAT_LONGB=[]
    ENTROPY_HERE_MEAN_LONGB=[]
    
    for line in PROTEIN_ENTROPY_FILE:
        line=line.strip()
        if (line!='') and (line!=' '):
            line=line.split()
            ENTROPY_HERE_MEAN.append(float(line[0]))
            ENTROPY_HERE_FLAT.append(float(line[1]))
            ENTROPY_HERE_FLAT_LONGB.append(float(line[2]))
            ENTROPY_HERE_MEAN_LONGB.append(float(line[3]))
    
    
    PROTEIN_TO_ENTROPY[PROTEIN_NAME]=statistics.median(ENTROPY_HERE_MEAN)
    PROTEIN_TO_ENTROPY_FLAT[PROTEIN_NAME]=statistics.median(ENTROPY_HERE_FLAT)
    PROTEIN_TO_ENTROPY_FLAT_LONGB[PROTEIN_NAME]=statistics.median(ENTROPY_HERE_FLAT_LONGB)
    PROTEIN_TO_ENTROPY_MEAN_LONGB[PROTEIN_NAME]=statistics.median(ENTROPY_HERE_MEAN_LONGB)


PROTEINS=list(set(PROTEINS))



for X in PROTEINS:
    PROTEIN=X
    ENTROPY=PROTEIN_TO_ENTROPY[X]
    ENTROPY_FLAT=PROTEIN_TO_ENTROPY_FLAT[X]
    ENTROPY_FLAT_LONGB=PROTEIN_TO_ENTROPY_FLAT_LONGB[X]
    ENTROPY_MEAN_LONGB=PROTEIN_TO_ENTROPY_MEAN_LONGB[X]
    
    
    print(F'{PROTEIN} : {ENTROPY} : {ENTROPY_FLAT} : {ENTROPY_MEAN_LONGB} : {ENTROPY_FLAT_LONGB}')
    Output_File.write(F'{PROTEIN} : {ENTROPY} : {ENTROPY_FLAT} : {ENTROPY_MEAN_LONGB} : {ENTROPY_FLAT_LONGB}\n')
    
    
   