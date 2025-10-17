import os
import random
import sys
import statistics




Rates_Folder=sys.argv[1]
OF=sys.argv[2]
Output_File=open(OF,'w')

FOLDERS=[x[0] for x in os.walk(Rates_Folder)]  ### Get all folders in folder
FOLDERS=FOLDERS[1:] ### Remove first item, because that is the directory itself

# print(Rates_Folder)
# print(FOLDERS)

PROTEIN_TO_RATES_FLAT={}
PROTEIN_TO_RATES_MEAN={}
PROTEIN_TO_LENGTH={}

PROTEINS=[]

for FOLDER in FOLDERS:
    
    PATH = FOLDER
    FILES = os.listdir(PATH)
    
    
    #########################################################
    ##### Per Protein
    
    PROTEIN = PATH.split('/')[len(PATH.split('/'))-1] ### Get file name as last from path //
    PROTEIN = PROTEIN.split('_PROT_REFERENCE')[0] ### Get protein name from file name / path
    
    PROTEIN_TO_RATES_FLAT[PROTEIN]=[]
    PROTEIN_TO_RATES_MEAN[PROTEIN]=[]
    PROTEIN_TO_LENGTH[PROTEIN]=[]
    PROTEINS.append(PROTEIN)
    
    for FILE in FILES:
        
        
        
        
        
        PROTEIN_RATE_FILE = open(F'{PATH}/{FILE}','r')
    
    
        #########################################################
        ##### Per Output File
        
        RATE_PER_FILE=[]
        length_counter=0
        
        for line in PROTEIN_RATE_FILE:
            
            ROW=line.split()
            
            
            if line[0]!='#' and ROW!=[]:
                
                line=line.strip()
                
                RATE=float(ROW[2])
                length_counter+=1
                if (RATE>-2.0) and (RATE<2.0):
                    RATE=0
                
                RATE_PER_FILE.append(RATE)

         
        SUM_OF_FILE = sum(RATE_PER_FILE)  
        MEAN_OF_FILE = SUM_OF_FILE/len(RATE_PER_FILE)  
        LENGTH_OF_FILE = length_counter
        
        PROTEIN_TO_RATES_FLAT[PROTEIN].append(SUM_OF_FILE)
        PROTEIN_TO_RATES_MEAN[PROTEIN].append(MEAN_OF_FILE)
        PROTEIN_TO_LENGTH[PROTEIN].append(LENGTH_OF_FILE)
        
        
        
        ##### Per Output File
        #########################################################
        
        
    ##### Per Protein    
    #########################################################
    
#### Calculate mean of sums and mean of means, drop it inside a folder for each protein
print(PROTEIN_TO_RATES_FLAT['AMELX'])
print(len(PROTEIN_TO_RATES_FLAT['AMELX']))









################# OUTPUT

Output_File.write('Protein_Name : Sum_of_Rates_all_Sites : Mean_Rate_Per_Site : Avrg_Protein_Length\n')

for PROTEIN in PROTEINS:
    
    FLAT = sum(PROTEIN_TO_RATES_FLAT[PROTEIN])
    MEAN = sum(PROTEIN_TO_RATES_MEAN[PROTEIN])
    AVRG_LENGTH = sum(PROTEIN_TO_LENGTH[PROTEIN]) / len(PROTEIN_TO_LENGTH[PROTEIN])
    
    print(F'{PROTEIN} : {FLAT} : {MEAN} : {AVRG_LENGTH}')

    Output_File.write(F'{PROTEIN} : {FLAT} : {MEAN} : {AVRG_LENGTH}\n')


