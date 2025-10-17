import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import glob
import statistics
import numbers
########################################################################
#### Conda env for running this: conda install conda-forge::pandas conda-forge::matplotlib conda-forge::numpy conda-forge::seaborn

######### LOAD DATA


Simulated_Subfolders = [ f.path for f in os.scandir('../Tree_Comparisons/') if f.is_dir() ]


NUMBER_TO_METRICS=[]

for FOLDER in Simulated_Subfolders:
    if 'REP-' in FOLDER:
        Distances_FILE=open(FOLDER+'/Tree_Distances.txt')
        Number_Of_Genes = FOLDER.split('-')[1]
        Number_Of_Genes = int(Number_Of_Genes.split('_')[1])
        
        for line in Distances_FILE:
            line=line.strip().split()
            
            if len(line)<8: #### If bootstrap is missing, tree distances were smaller than 0.001, thus turn bootstrap into 0
                line.append(0)
            
            NUMBER_TO_METRICS.append([Number_Of_Genes,line])
        






#########################################################################
### Load and Organise Variation Data and 
### Hominids


Folders_with_number_of_variants = [ f.path for f in os.scandir('../Dataset_Analysis/Workspace/2_DATASETS/') if f.is_dir() ]

Variants_Per_Protein_Number={} #### Should have 12 keys, each key is a 'number of proteins' (1-12) and corresponds to a list of 1000 items, each of which is the variant sites of that run

for FOLDER in Folders_with_number_of_variants:

    Number_of_proteins=FOLDER.split('/')[len(FOLDER.split('/'))-1]
    Number_of_proteins=Number_of_proteins.split('-')[1]
    Number_of_proteins=int(Number_of_proteins.split('_')[1])

    
    
    Variant_file=open(FOLDER+'/CONCATINATED/Number_of_Variants.txt','r')
    Variants=int(Variant_file.readline())
    
    if Number_of_proteins in Variants_Per_Protein_Number.keys():
        Variants_Per_Protein_Number[Number_of_proteins].append(Variants)
        
    if Number_of_proteins not in Variants_Per_Protein_Number.keys():
        Variants_Per_Protein_Number[Number_of_proteins]=[Variants]




Max_Proteins=len(Variants_Per_Protein_Number.keys())
Variants_Output_Folder=open('Variant_Table_Hominins.tsv','w')
Variants_Output_Folder.write('Number_of_Proteins\tMean_Number_of_Variants\tMedian_Number_of_Variants\tMinimum_Number_of_Variants\tMaximum_Number_of_Variants\n')

for nop in range(1,Max_Proteins+1):

    print(F'$$$$ Hominids - {nop} Proteins $$$$')
    print(F' Mean: {np.median(Variants_Per_Protein_Number[nop])}')
    print(F' Median: {np.median(Variants_Per_Protein_Number[nop])}')
    print(F' Minimum: {min(Variants_Per_Protein_Number[nop])} and Maximum: {max(Variants_Per_Protein_Number[nop])}')
    print('\n')
    
    Variants_Output_Folder.write(F'{np.median(Variants_Per_Protein_Number[nop])}\t{np.median(Variants_Per_Protein_Number[nop])}\t{min(Variants_Per_Protein_Number[nop])}\t{max(Variants_Per_Protein_Number[nop])}\n')
    
    
    
#######################################################################################################
####### PLOTTING
#### Variation Boxplot over number of proteins

data = [Variants_Per_Protein_Number[nop] for nop in range(1,Max_Proteins+1)]
fig, ax = plt.subplots(figsize=(8.6, 2.5))

plt.xticks(np.arange(1, len(data)+2, step=1))
plt.yticks(np.arange(0, 31, step=5))

bp = ax.boxplot(data,showfliers=False)

ax.set_title("Number of variants per number of proteins")
ax.set_xlabel("Number of Proteins")
ax.set_ylabel("Number of Variants")



plt.savefig('Variation.pdf')
plt.savefig('Variation.svg')













########################################################################
########################################################################
####### Calcualte Bootstrap support and relation with accrucacy and number of genes, ignore polytomies


NUMBER_TO_BOOTSTRAPS_AND_TREE={} ### Dictionary. Keys are number of proteins and values are lists of lists. Each sublist contains two elements: the topology supported by a tree and the bootstrap support of that topology
for k in NUMBER_TO_METRICS:
    NoP=k[0]
    Supported_Tree=k[1][6]
    Bootstrap_Support=k[1][7]
    
    if NoP not in NUMBER_TO_BOOTSTRAPS_AND_TREE.keys():
        NUMBER_TO_BOOTSTRAPS_AND_TREE[NoP]=[]
    
    if NoP in NUMBER_TO_BOOTSTRAPS_AND_TREE.keys():
        NUMBER_TO_BOOTSTRAPS_AND_TREE[NoP].append([Supported_Tree,Bootstrap_Support]) 



print('BOOTSTRAP SUPPORT FOR T1 TREE TOPOLOGY')

SORTED_N_PROTEINS=sorted(NUMBER_TO_BOOTSTRAPS_AND_TREE.keys())
BOOTSTRAPS_SUPPORTING_T1_PER_NUMBER=[]

for PROTEIN_N in SORTED_N_PROTEINS: #### cycle through all protein numbers
    
    TOPOLOGIES_AND_SUPPORTS=NUMBER_TO_BOOTSTRAPS_AND_TREE[PROTEIN_N]
    BOOTSTRAPS_SUPPORTING_T1=[]
    
    for TAS in TOPOLOGIES_AND_SUPPORTS:
        if TAS[0]=='T1':
            BOOTSTRAPS_SUPPORTING_T1.append(float(TAS[1]))
    
    AVG_SUP_T1=np.median(BOOTSTRAPS_SUPPORTING_T1)
    BOOTSTRAPS_SUPPORTING_T1_PER_NUMBER.append(BOOTSTRAPS_SUPPORTING_T1)
    print(F"Avergage Bootstrap support of T1: {AVG_SUP_T1} for {PROTEIN_N} number of proteins")
    # print(F"List of bootstrap support of TT: {BOOTSTRAPS_SUPPORTING_T1}")

print("\n\n")


########################################################################


print('BOOTSTRAP SUPPORT FOR T2 TREE TOPOLOGY')
BOOTSTRAPS_SUPPORTING_T2_PER_NUMBER=[]

for PROTEIN_N in SORTED_N_PROTEINS: #### cycle through all protein numbers
    
    TOPOLOGIES_AND_SUPPORTS=NUMBER_TO_BOOTSTRAPS_AND_TREE[PROTEIN_N]
    BOOTSTRAPS_SUPPORTING_T2=[]
    
    for TAS in TOPOLOGIES_AND_SUPPORTS:
        if TAS[0]=='T2':
            BOOTSTRAPS_SUPPORTING_T2.append(float(TAS[1]))
    
    if BOOTSTRAPS_SUPPORTING_T2==[]: #### If no trees are supporting these topologies, e.g. in higher number of proteins
        BOOTSTRAPS_SUPPORTING_T2=[0]
        
    BOOTSTRAPS_SUPPORTING_T2 = [x for x in BOOTSTRAPS_SUPPORTING_T2 if isinstance(x, numbers.Number)] ### if somehow Nans are here
    AVG_SUP_T2=np.median(BOOTSTRAPS_SUPPORTING_T2)
    BOOTSTRAPS_SUPPORTING_T2_PER_NUMBER.append(BOOTSTRAPS_SUPPORTING_T2)
    
    print(F"Avergage Bootstrap support of T2: {AVG_SUP_T2} for {PROTEIN_N} number of proteins")
    # print(F"List of bootstrap support of T2: {BOOTSTRAPS_SUPPORTING_T2}")

print("\n\n")




########################################################################


print('BOOTSTRAP SUPPORT FOR T3 TREE TOPOLOGY')
BOOTSTRAPS_SUPPORTING_T3_PER_NUMBER=[]

for PROTEIN_N in SORTED_N_PROTEINS: #### cycle through all protein numbers
    
    TOPOLOGIES_AND_SUPPORTS=NUMBER_TO_BOOTSTRAPS_AND_TREE[PROTEIN_N]
    BOOTSTRAPS_SUPPORTING_T3=[]
    
    for TAS in TOPOLOGIES_AND_SUPPORTS:
        if TAS[0]=='T3':
            BOOTSTRAPS_SUPPORTING_T3.append(float(TAS[1]))
    
    if BOOTSTRAPS_SUPPORTING_T3==[]: #### If no trees are supporting these topologies, e.g. in higher number of proteins
        BOOTSTRAPS_SUPPORTING_T3=[0]
        
    BOOTSTRAPS_SUPPORTING_T3 = [x for x in BOOTSTRAPS_SUPPORTING_T3 if isinstance(x, numbers.Number)] ### cleanup if somehow Nans are here
    AVG_SUP_T3=np.median(BOOTSTRAPS_SUPPORTING_T3)
    
    BOOTSTRAPS_SUPPORTING_T3_PER_NUMBER.append(BOOTSTRAPS_SUPPORTING_T3)
    
    print(F"Avergage Bootstrap support of T3: {AVG_SUP_T3} for {PROTEIN_N} number of proteins")
    # print(F"List of bootstrap support of T3: {BOOTSTRAPS_SUPPORTING_T3}")

print("\n\n")



########################################################################
#### Plot bootstraps, 3 for each number of proteins


BOOTSTRAPS=[BOOTSTRAPS_SUPPORTING_T1_PER_NUMBER ,BOOTSTRAPS_SUPPORTING_T2_PER_NUMBER ,BOOTSTRAPS_SUPPORTING_T3_PER_NUMBER]

fig, ax = plt.subplots(figsize=(8.6, 2.5))

# plt.xticks(np.arange(1, Max_Proteins+2, step=1))
plt.yticks(np.arange(0, 101, step=10))

## Plot boxplots in positions
bp1 = ax.boxplot(BOOTSTRAPS[0], positions = np.arange(1, (Max_Proteins)*2, step=2)-0.5, sym='', widths=0.4,showfliers=True,showmeans=False, patch_artist=True)
bp2 = ax.boxplot(BOOTSTRAPS[1], positions = np.arange(1, (Max_Proteins)*2, step=2), sym='', widths=0.4,showfliers=True,showmeans=False, patch_artist=True)
bp3 = ax.boxplot(BOOTSTRAPS[2], positions = np.arange(1, (Max_Proteins)*2, step=2)+0.5, sym='', widths=0.4,showfliers=True,showmeans=False, patch_artist=True)

c1="blueviolet"
c2="dodgerblue"
c3="deepskyblue"



## Colours
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp1[item], color=c1)
plt.setp(bp1['boxes'], color="black",linewidth=0.8)
plt.setp(bp1['boxes'], facecolor=c1)
plt.setp(bp1["medians"], color="black")



## Colours
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp2[item], color=c2)
plt.setp(bp2['boxes'], color="black",linewidth=0.8)
plt.setp(bp2['boxes'], facecolor=c2)
plt.setp(bp2["medians"], color="black")



## Colours
for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp3[item], color=c3)
plt.setp(bp3['boxes'], color="black",linewidth=0.8)
plt.setp(bp3['boxes'], facecolor=c3)
plt.setp(bp3["medians"], color="black")



ax.set_title("Bootstraps supporting each topology")
ax.set_xlabel("Number of Proteins")
ax.set_ylabel("Bootstrap support")


ax.set_xticks([])
ax.set_xticks([], minor=True)

plt.savefig('Bootstraps.pdf')
plt.savefig('Bootstraps.svg')























########################################################################
########################################################################
####### Calcualte Average accuracy for each Gene number

NUMBER_TO_DIST={}
NUMBER_TO_TOPOLOGIES={}

for k in NUMBER_TO_METRICS:
    
    if k[0] in NUMBER_TO_DIST.keys():
        NUMBER_TO_DIST[int(k[0])].append(float(k[1][0]))
        NUMBER_TO_TOPOLOGIES[int(k[0])].append(k[1][6])
        
        
    if k[0] not in NUMBER_TO_DIST.keys():
        NUMBER_TO_DIST[int(k[0])]=[float(k[1][0])]
        NUMBER_TO_TOPOLOGIES[int(k[0])]=[k[1][6]]


Max_Number_of_Genes=max(list(NUMBER_TO_TOPOLOGIES.keys()))


# print(NUMBER_TO_TOPOLOGIES[1])

AVG=[]
SUM=[]
PRECENTAGE_OF_ACCURATE_TREES=[]

for j,k in NUMBER_TO_DIST.items():
    AVG.append([int(j),statistics.median(k)])
    SUM.append([int(j),sum(k)])
    PRECENTAGE_OF_ACCURATE_TREES.append([int(j),(k.count(0)/(len(k)))])
    
    
    
    
# print(NUMBER_TO_DIST)    
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')    

AVG=sorted(AVG, key=lambda x: x[0])   
SUM=sorted(SUM, key=lambda x: x[0]) 
PRECENTAGE_OF_ACCURATE_TREES=sorted(PRECENTAGE_OF_ACCURATE_TREES, key=lambda x: x[0]) 
# print(AVG)
# print(SUM)  
# print(PRECENTAGE_OF_ACCURATE_TREES)  


Number_of_Genes=np.asarray([int(x[0]) for x in AVG])
Number_of_Genes_Label=[str(x[0]) for x in AVG]
Average_Distance=np.asarray([float(x[1]) for x in AVG])
Summed_Distance=np.asarray([float(x[1]) for x in SUM])
Percentage_of_Trees=np.asarray([float(x[1]) for x in PRECENTAGE_OF_ACCURATE_TREES])

# print(Number_of_Genes,Summed_Distance)








################################################################################################################################################
##### Generate Dictionary for Stacked Barplot


Topologies_Supported_By_Trees_Per_Gene={
"T1":[],
"T2":[],
"T3":[],
"T4":[]
}

Discordant_Topologies={
"NOT1":[]
}

for J in range(1,Max_Number_of_Genes+1):
    # print(NUMBER_TO_TOPOLOGIES[J].count("T4"))
    Topologies_Supported_By_Trees_Per_Gene["T1"].append(NUMBER_TO_TOPOLOGIES[J].count("T1"))
    Topologies_Supported_By_Trees_Per_Gene["T2"].append(NUMBER_TO_TOPOLOGIES[J].count("T2"))
    Topologies_Supported_By_Trees_Per_Gene["T3"].append(NUMBER_TO_TOPOLOGIES[J].count("T3"))
    Topologies_Supported_By_Trees_Per_Gene["T4"].append(NUMBER_TO_TOPOLOGIES[J].count("T4"))
    Discordant_Topologies["NOT1"].append(1000-NUMBER_TO_TOPOLOGIES[J].count("T1"))
    
# print(Topologies_Supported_By_Trees_Per_Gene)




PROTEINS=[x for x in range(1,Max_Number_of_Genes+1)]
plotdata = pd.DataFrame(Topologies_Supported_By_Trees_Per_Gene, 
    index=PROTEINS
)

plotdata2 = pd.DataFrame(Discordant_Topologies, 
    index=PROTEINS
)

# my_colors0=('blueviolet','royalblue', 'dodgerblue','deepskyblue')
my_colors0=('blueviolet','dodgerblue', 'deepskyblue','black')
# my_colors0=('blueviolet','lightskyblue', 'deepskyblue','black')

#### Fix Size
fig, axs = plt.subplots(2,figsize=(8.70, 2.35),sharex=True)

### Fix Ticks subplot 1
axs[0].set_ylim([0,1000])
axs[0].set_yticks([0,200,400,600,800,1000])
axs[0].set_yticklabels([0,200,400,600,800,1000],fontsize=6)
axs[0].tick_params(axis='x', which='both',bottom=False)

### Fix Ticks subplot 2
axs[1].set_ylim([0,1000])
axs[1].set_yticks([0,200,400,600,800,1000])
axs[1].set_yticklabels([0,200,400,600,800,1000],fontsize=6)

plotdata2.plot(kind='bar',legend=False,width=0.8,color='maroon',ax=axs[0])
plotdata.plot(kind='bar', stacked=True,legend=False,width=0.8,color=my_colors0,ax=axs[1])
### , 



plt.savefig("Percentage_of_Trees.pdf", format="pdf", bbox_inches="tight")
plt.savefig("Percentage_of_Trees.svg", format="svg", bbox_inches="tight")
plt.show()




