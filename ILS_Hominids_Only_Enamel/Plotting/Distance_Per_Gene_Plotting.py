import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import glob


########################################################################
######### LOAD SIMULATED DATA


Simulated_Subfolders = [ f.path for f in os.scandir('./Simulated/') if f.is_dir() ]

Simulated_Subfolders = [x for x in Simulated_Subfolders if 'EXAMP' in x]


GENE_NAMES=[]
GENE_TREE_DISTANCES=[]


for FOLDER in Simulated_Subfolders:
    
    GENE_NAME=FOLDER.split('/')[2]

    
    Distances_FILE=open(FOLDER+F'/Tree_Distances.txt')
    
    print(GENE_NAME,Distances_FILE)
    
    for line in Distances_FILE:
        line=line.strip().split()
        line=[float(x) for x in line]
        GENE_NAMES.append(GENE_NAME)
        GENE_TREE_DISTANCES.append(line)







########################################################################

### PLOT DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Number of Genes and Distance to original tree', fontsize=14, fontweight='bold')


ax.set_xlabel('Name of Gene', fontsize=11, fontweight='bold')
ax.set_ylabel('Distance to Original Tree', fontsize=11, fontweight='bold')

ax.scatter([x for x in GENE_NAMES], [x[0] for x in GENE_TREE_DISTANCES])

### for Z in range(0,len(REAL_DATA_Y)):
    ### ax.annotate(REAL_DATA_TEXT[Z], (REAL_DATA_X[Z], REAL_DATA_Y[Z]))



plt.show()
