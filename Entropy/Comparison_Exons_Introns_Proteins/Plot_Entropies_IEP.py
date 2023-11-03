import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import random


########################################################################
######### LOAD DATA


GENE_TO_ENTROPY_INTRONS={}
GENE_TO_ENTROPY_EXONS={}
GENE_TO_ENTROPY_PROTEINS={}

GENE_NAMES=['AHSG','ALB','AMBN', 'AMELX', 'AMELY', 'AMTN', 'COL17A1', 'COL1A1', 'COL1A2', 'ENAM', 'MMP20', 'ODAM']
DATA_TYPE=['INTRONS','EXONS','PROTEIN']

for J in GENE_NAMES:
    for K in DATA_TYPE:
        FILE=open(F'Entropy_Measurements/{J}_{K}.entr','r')
        ENTROPY=FILE.readline().split()[1] ### Flat [1] or Divided by length [0]

        if K=='INTRONS':
            GENE_TO_ENTROPY_INTRONS[J]=float(ENTROPY)
            
        if K=='EXONS':
            GENE_TO_ENTROPY_EXONS[J]=float(ENTROPY)
            
        if K=='PROTEIN':
            GENE_TO_ENTROPY_PROTEINS[J]=float(ENTROPY)            


print(GENE_TO_ENTROPY_INTRONS,GENE_TO_ENTROPY_EXONS,GENE_TO_ENTROPY_PROTEINS)        
## PLOT DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Entropy Comparison of 3 Data types: Introns & Exons, Exons, Proteins', fontsize=14, fontweight='bold')



ax.set_xlabel('Gene', fontsize=11, fontweight='bold',labelpad=20)
ax.set_ylabel('Entropy', fontsize=11, fontweight='bold',labelpad=10)



ax.scatter([x for x in GENE_NAMES], [GENE_TO_ENTROPY_INTRONS[x] for x in GENE_NAMES],c='darkviolet',s=180,zorder=2)
ax.scatter([x for x in GENE_NAMES], [GENE_TO_ENTROPY_EXONS[x] for x in GENE_NAMES],c='orange',s=150,zorder=3)
ax.scatter([x for x in GENE_NAMES], [GENE_TO_ENTROPY_PROTEINS[x] for x in GENE_NAMES],c='gold',s=150,zorder=4)

# ax.set_yticks([],)
# ax.set_yticks([], minor=True)

plt.show()


