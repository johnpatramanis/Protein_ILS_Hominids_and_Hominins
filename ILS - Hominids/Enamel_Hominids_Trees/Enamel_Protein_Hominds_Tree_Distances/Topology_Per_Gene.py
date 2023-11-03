import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


########################################################################
######### LOAD DATA


Protein_Subfolders = [ f.path for f in os.scandir('./') if f.is_dir() ]
DNA_Folder = '../Enamel_DNA_Hominds_Tree_Distances/GENE_TREES_ENAMEL/'




GENE_NAMES=[]
GENE_TREE_DISTANCES_PROTEIN=[]


for FOLDER in Protein_Subfolders:
    
    GENE_NAME=FOLDER.split('/')[1]
    Distances_FILE=open(FOLDER+F'/{GENE_NAME}_Tree_Distances.txt')
    
    
    for line in Distances_FILE:
        line=line.strip().split()
        line=[float(x) for x in line]
        GENE_NAMES.append(GENE_NAME)
        GENE_TREE_DISTANCES_PROTEIN.append(line)



GENE_TREE_DISTANCES_DNA=[]

for GENE in GENE_NAMES:
    
    GENE_NAME=GENE
    Distances_FILE=open(DNA_Folder+F'Hominid_{GENE}_Tree_Distances.txt')
    
    
    for line in Distances_FILE:
        line=line.strip().split()
        line=[float(x) for x in line]
        GENE_TREE_DISTANCES_DNA.append(line)





########################################################################
### PLOT DNA DATA

# fig = plt.figure()
# ax = fig.add_subplot()

# fig.suptitle('Tree Topology of enamel Genes', fontsize=14, fontweight='bold')



# ax.set_xlabel('Gene', fontsize=11, fontweight='bold',labelpad=20)
# ax.set_ylabel('Alternative Tree Topologies', fontsize=11, fontweight='bold',labelpad=10)



# plt.axhline(y = 1, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
# plt.axhline(y = 2, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
# plt.axhline(y = 3, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
# plt.axhline(y = 4, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)

# ax.scatter([x for x in GENE_NAMES], [x[1] for x in GENE_TREE_DISTANCES_DNA],c='red',s=60,zorder=2)

# ax.set_yticks([],)
# ax.set_yticks([], minor=True)

# plt.show()







########################################################################
### PLOT PROTEIN DATA

# fig = plt.figure()
# ax = fig.add_subplot()

# fig.suptitle('Tree Topology of enamel Proteins', fontsize=14, fontweight='bold')



# ax.set_xlabel('Protein', fontsize=11, fontweight='bold',labelpad=20)
# ax.set_ylabel('Alternative Tree Topologies', fontsize=11, fontweight='bold',labelpad=10)



plt.axhline(y = 1, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 2, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 3, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 4, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)

# ax.scatter([x for x in GENE_NAMES], [x[1] for x in GENE_TREE_DISTANCES_PROTEIN],c='blue',s=60,zorder=3)


# ax.set_yticks([],)
# ax.set_yticks([], minor=True)

# plt.show()







########################################################################
### PLOT COMBINED DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Tree Topology of enamel Genes and Proteins', fontsize=14, fontweight='bold')



ax.set_xlabel('Gene', fontsize=11, fontweight='bold',labelpad=20)
ax.set_ylabel('Resulting Tree Topology', fontsize=11, fontweight='bold',labelpad=10)



plt.axhline(y = 1, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 2, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 3, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)
plt.axhline(y = 4, color = 'grey', linestyle = '--', xmax = 16,zorder=1,linewidth=0.8)

ax.scatter([x for x in GENE_NAMES], [x[0] for x in GENE_TREE_DISTANCES_DNA],c='darkviolet',s=180,zorder=2)
ax.scatter([x for x in GENE_NAMES], [x[0] for x in GENE_TREE_DISTANCES_PROTEIN],c='orange',s=150,zorder=3)

ax.set_yticks([],)
ax.set_yticks([], minor=True)

plt.show()
