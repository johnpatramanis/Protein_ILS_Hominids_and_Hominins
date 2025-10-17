import os
import random
import sys
import statistics
import math




Rates_Folder=sys.argv[1]
OF=sys.argv[2]
Output_File=open(OF,'w')

FOLDERS=[x[0] for x in os.walk(Rates_Folder)]  ### Get all folders in folder
FOLDERS=FOLDERS[1:] ### Remove first item, because that is the directory itself


PROTEIN_TO_RATES={} ##### Dictionary. Each entry corresponds to a list of lists, with length equal to the length of the protein. Each sub-list has the evolutionary rate of each iteration (usually 1000 iters)
PROTEINS=[] ##### List of proteins (names)

Iteration_counter=0
for FOLDER in FOLDERS:
    ##### For each protein in folder
    
    PATH = FOLDER
    FILES = os.listdir(PATH)
    
    
    #########################################################
    ##### Per Protein
    
    PROTEIN = PATH.split('/')[len(PATH.split('/'))-1] ### Get file name as last from path //
    PROTEIN = PROTEIN.split('_PROT_REFERENCE')[0] ### Get protein name from file name / path
    PROTEINS.append(PROTEIN)
    PROTEIN_TO_RATES[PROTEIN]={}
    
    for FILE in FILES:
    #### For each iteration of this protein   
        
        
        
        
        PROTEIN_RATE_FILE = open(F'{PATH}/{FILE}','r')
    
    
        #########################################################
        ##### Per Output File
        RATE_PER_POSITION={} #### Dictionary for this iteration (keys are position number and values are evol. rate of position
        length_counter=1
        
        for line in PROTEIN_RATE_FILE:
            
            ROW=line.split()
            
            
            if line[0]!='#' and ROW!=[]:
                
                line=line.strip()
                
                RATE=float(ROW[2])
                
                if RATE<0.1:
                    RATE=0
                
                
                RATE_PER_POSITION[length_counter]=RATE ### add to dictionary
                
                length_counter+=1
        
        for X,Y in RATE_PER_POSITION.items():
            
            if X not in PROTEIN_TO_RATES[PROTEIN].keys():
                PROTEIN_TO_RATES[PROTEIN][X]=[Y]
            else:
                PROTEIN_TO_RATES[PROTEIN][X].append(Y)


        
        
        Iteration_counter+=1
        # print(PROTEIN,Iteration_counter,len(RATE_PER_POSITION))
        
        
        
        ##### Per Output File
        #########################################################
        
        
    ##### Per Protein    
    #########################################################
    
    print(F"\n\nProtein: {PROTEIN} has a total of {len(PROTEIN_TO_RATES[PROTEIN])} sites, each site has roughly {len(PROTEIN_TO_RATES[PROTEIN][1])} measurements of evolutionary rate \n\n")
    
    
    ##### Calculate Phylogenetic Informativeness!
    TIMES=[0.001,0.01,0.1,0.5,1,6,10]
    INFORMATIVENESS_PER_TIME=[]
    INFORMATIVENESS_PER_TIME_V2=[]
    
    
    for T in TIMES: ### For timepoint of interest
        
        Average_Lambda=[]
        
        INFORMATIVENESS_PER_THIS_TIME = [] 
        
        for POSITION in PROTEIN_TO_RATES[PROTEIN].keys():
            
            Lambda = sum(PROTEIN_TO_RATES[PROTEIN][POSITION]) / len(PROTEIN_TO_RATES[PROTEIN][POSITION]) ### Average evolutionary rate of position (Lambda)
            Average_Lambda.append(Lambda)
            # Lambda = Lambda/10
            
            LambdaSqr = Lambda * Lambda
            
            E_comp = 1/(math.exp(4*T*Lambda))
            
            INFORMATIVENESS_PER_POSITION= 16 * LambdaSqr * T * E_comp
            
            INFORMATIVENESS_PER_THIS_TIME.append(INFORMATIVENESS_PER_POSITION)
            
            # print(F"Lambda:{Lambda},Lambda Squared: {LambdaSqr}, E: {E_comp}\n")
            
        INFORMATIVENESS_PER_TIME.append(sum(INFORMATIVENESS_PER_THIS_TIME))
        
        
        
        
        ### Alternative Calculation
        Average_Lambda=sum(Average_Lambda)/len(Average_Lambda)
        Average_LambdaSqrt = Average_Lambda * Average_Lambda
        E_comp_V2 = 1/(math.exp(4*T*Average_Lambda))
        
        INFORMATIVENESS_PER_THIS_TIME_V2= 16 * Average_LambdaSqrt * T * E_comp_V2 * len(PROTEIN_TO_RATES[PROTEIN].keys())
        INFORMATIVENESS_PER_TIME_V2.append(INFORMATIVENESS_PER_THIS_TIME_V2)
        
    print(PROTEIN, INFORMATIVENESS_PER_TIME,'\n')
    print(PROTEIN, INFORMATIVENESS_PER_TIME_V2,'\n')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#### TEST Print
### print(PROTEIN_TO_RATES_FLAT['AMELX'])
### print(len(PROTEIN_TO_RATES_FLAT['AMELX']))









################# OUTPUT

# Output_File.write('Protein_Name : Sum_of_Rates_all_Sites : Mean_Rate_Per_Site : Avrg_Protein_Length\n')

# for PROTEIN in PROTEINS:
    
    # FLAT = sum(PROTEIN_TO_RATES_FLAT[PROTEIN])
    # MEAN = sum(PROTEIN_TO_RATES_MEAN[PROTEIN])
    # AVRG_LENGTH = sum(PROTEIN_TO_LENGTH[PROTEIN]) / len(PROTEIN_TO_LENGTH[PROTEIN])
    
    # print(F'{PROTEIN} : {FLAT} : {MEAN} : {AVRG_LENGTH}')

    # Output_File.write(F'{PROTEIN} : {FLAT} : {MEAN} : {AVRG_LENGTH}\n')


