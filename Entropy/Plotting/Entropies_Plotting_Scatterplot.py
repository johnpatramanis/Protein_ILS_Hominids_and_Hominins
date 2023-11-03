import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import pandas as pd
import random


########################################################################
######### LOAD DATA


Real_Entropies={}
Real_Entropies_Flat={}
Real_Entropies_Mean_LongB={}
Real_Entropies_Flat_LongB={}


RL_ENTR=open('../Entropies_Hominids/Combined_Entropies','r')

for line in RL_ENTR:
    line=line.strip().split(' : ')
    PROTEIN=line[0]
    PROT_ENTROPY=float(line[1])
    PROT_ENTROPY_FLAT=float(line[2])
    PROT_ENTROPY_MEAN_LONGB=float(line[3])
    PROT_ENTROPY_FLAT_LONGB=float(line[4])
    
    
    Real_Entropies[PROTEIN]=PROT_ENTROPY
    Real_Entropies_Flat[PROTEIN]=PROT_ENTROPY_FLAT
    Real_Entropies_Mean_LongB[PROTEIN]=PROT_ENTROPY_MEAN_LONGB
    Real_Entropies_Flat_LongB[PROTEIN]=PROT_ENTROPY_FLAT_LONGB


Sorted_Real_Entropies=sorted(Real_Entropies.items(), key=lambda x:x[1])
Sorted_Real_Entropies_Flat=sorted(Real_Entropies_Flat.items(), key=lambda x:x[1])
Sorted_Real_Entropies_Mean_LongB=sorted(Real_Entropies_Mean_LongB.items(), key=lambda x:x[1])
Sorted_Real_Entropies_Flat_LongB=sorted(Real_Entropies_Flat_LongB.items(), key=lambda x:x[1])


count=1
print('Mean Entropies Sorted')
for j in range(0,len(Sorted_Real_Entropies)):
    print(F'{count} - {Sorted_Real_Entropies[j]}')
    count+=1



count=1
print('\n\n\n')
print('Flat Entropies Sorted')
for j in range(0,len(Sorted_Real_Entropies_Flat)):
    print(F'{count} - {Sorted_Real_Entropies_Flat[j]}')
    count+=1



count=1
print('\n\n\n')

print('Mean-LongB Entropies Sorted')
for j in range(0,len(Sorted_Real_Entropies_Mean_LongB)):
    print(F'{count} - {Sorted_Real_Entropies_Mean_LongB[j]}')
    count+=1


count=1
print('\n\n\n')

print('Flat-LongB Entropies Sorted')
for j in range(0,len(Sorted_Real_Entropies_Flat_LongB)):
    print(F'{count} - {Sorted_Real_Entropies_Flat_LongB[j]}')
    count+=1




REAL_DATA_X=[]
REAL_DATA_Y=[]
REAL_DATA_LABELS=[]
REAL_DATA_TEXT=[]

ENAMEL_PROTEIN_NAMES=['ENAM','AMELX','AMELY','AMBN','COL17A1','MMP20','ALB','AHSG','AMTN','ODAM']
COLLAGEN_PROTEIN_NAMES=['COL1A1','COL1A2']
FIBROPROTEINS_PROTEIN_NAMES=['FGB','FGG']
UBIQUITINES_PROTEIN_NAME=['USP46']
HISTONES_PROTEIN_NAME=['H2BC9','H2BC3']


ENAMEL=[]
COLLAGEN=[]
FIBROPROTEINS=[]
UBIQUITINES=[]
HISTONES=[]
OTHER=[]

for x,j in Real_Entropies_Flat_LongB.items():

    REAL_DATA_Y.append(j)
    ### REAL_DATA_X.append(random.uniform(-0.1, 0.1))
    REAL_DATA_X.append(0)
    REAL_DATA_TEXT.append(x)
    
    SPECIAL_CATEGORY=False
    
    if x in ENAMEL_PROTEIN_NAMES:
        ENAMEL.append([x,j])
        SPECIAL_CATEGORY=True
    
    if x in COLLAGEN_PROTEIN_NAMES:
        COLLAGEN.append([x,j])
        SPECIAL_CATEGORY=True
 
    if x in FIBROPROTEINS_PROTEIN_NAMES:
        FIBROPROTEINS.append([x,j])
        SPECIAL_CATEGORY=True
 
    if x in UBIQUITINES_PROTEIN_NAME:
        UBIQUITINES.append([x,j])
        SPECIAL_CATEGORY=True

    if x in HISTONES_PROTEIN_NAME:
        HISTONES.append([x,j])
        SPECIAL_CATEGORY=True
    
    if SPECIAL_CATEGORY==False:
        OTHER.append([x,j])
        
#######################################################################

## PLOT DATA


fig = plt.figure()
ax = fig.add_subplot()

plt.figure(figsize=(10,14))



#### Enamel Proteins
if OTHER!=[]:
    sns.stripplot(x=[0 for x in OTHER], y=[x[1] for x in OTHER],color='black') #### Scatter plot with proteins that don't match a specific category 

sns.stripplot(x=[0 for x in ENAMEL], y=[x[1] for x in ENAMEL],color='dodgerblue',size=11.5) #### Scatter plot with proteins that are found exclusively in Enamel    

sns.stripplot(x=[0 for x in COLLAGEN], y=[x[1] for x in COLLAGEN],color='firebrick',size=11.5) #### Scatter plot with collagen proteins   


sns.stripplot(x=[0 for x in FIBROPROTEINS], y=[x[1] for x in FIBROPROTEINS],color='limegreen',size=11.5) #### Scatter plot with fibro proteins
sns.stripplot(x=[0 for x in UBIQUITINES], y=[x[1] for x in UBIQUITINES],color='darkorchid',size=11.5) #### Scatter plot with ubiquitines 
sns.stripplot(x=[0 for x in HISTONES], y=[x[1] for x in HISTONES],color='gold',size=11.5) #### Scatter plot with histones 


for j in ax.get_children():
    print(j)

fig.suptitle('Real Proteins - Entropies', fontsize=14, fontweight='bold')


ax.set_xlabel('Proteins', fontsize=11, fontweight='bold')
ax.set_ylabel('Calculated Entropy', fontsize=11, fontweight='bold')


plt.show()

