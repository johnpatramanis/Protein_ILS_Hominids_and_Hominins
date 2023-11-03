import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import glob
import statistics
from scipy.optimize import curve_fit
########################################################################
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
            NUMBER_TO_METRICS.append([Number_Of_Genes,line])




########################################################################
####### Calcualte Average accuracy for each Gene number

NUMBER_TO_DIST={}
NUMBER_TO_TOPOLOGIES={}

for k in NUMBER_TO_METRICS:
    
    if k[0] in NUMBER_TO_DIST.keys():
        NUMBER_TO_DIST[int(k[0])].append(float(k[1][0]))
        NUMBER_TO_TOPOLOGIES[int(k[0])].append(k[1][-1])
        
        
    if k[0] not in NUMBER_TO_DIST.keys():
        NUMBER_TO_DIST[int(k[0])]=[float(k[1][0])]
        NUMBER_TO_TOPOLOGIES[int(k[0])]=[k[1][-1]]


Max_Number_of_Genes=max(list(NUMBER_TO_TOPOLOGIES.keys()))


# print(NUMBER_TO_TOPOLOGIES[1])

AVG=[]
SUM=[]
PRECENTAGE_OF_ACCURATE_TREES=[]

for j,k in NUMBER_TO_DIST.items():
    AVG.append([int(j),statistics.mean(k)])
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
#### Generate Curve

# Define the Gaussian function
def Gauss(x, A, B):
    
    y = A*np.exp(-1*B*x**2)
    return y

parameters, covariance = curve_fit(Gauss, Number_of_Genes, Average_Distance)

fit_A = parameters[0]
fit_B = parameters[1]
fit_y = Gauss(Number_of_Genes, fit_A, fit_B)


################################################################################################################################################
#### Generate Better Curve
def cos_func(x, D, E):
    y = D*np.cos(E*x)
    return y
    
    
guess = [0.25, 0.1]
parameters, covariance = curve_fit(cos_func, Number_of_Genes, Average_Distance, p0=guess)

fit_D = parameters[0]
fit_E = parameters[1]

fit_cosine = cos_func(Number_of_Genes, fit_D, fit_E)


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


































###########  PLOT DATA ### Avrg distance dots with fitted line

# fig = plt.figure()
# ax = fig.add_subplot()
# fig.suptitle('Number of Genes and Distance to original tree', fontsize=14, fontweight='bold')


# ax.set_xlabel('Number of Genes/Proteins', fontsize=11, fontweight='bold')
# ax.set_ylabel('Distance to Original Tree', fontsize=11, fontweight='bold')

# ax.scatter(Number_of_Genes, Average_Distance,cmap="jet")
# plt.plot(Number_of_Genes, fit_y, '-', label='fit',color='orange')
##### plt.plot(Number_of_Genes, fit_cosine, '-', label='fit')

# plt.xticks(np.arange(0, len(Average_Distance)+1, step=1))
#####  plt.yticks(np.arange(0, 1, step=1))


# plt.show()









################################################################################################################################################
##### Third Plot #### Barplot with trees that agree with known species tree

# fig = plt.figure()
# ax = fig.add_subplot()

# fig.suptitle('Number of Genes and Percentage of Accurate Trees', fontsize=14, fontweight='bold')


# ax.set_xlabel('Number of Genes/Proteins', fontsize=11, fontweight='bold')
# ax.set_ylabel('Percentage of trees that agree with species tree', fontsize=11, fontweight='bold')

# plt.bar(Number_of_Genes_Label, Percentage_of_Trees, color ='maroon', width = 0.8)

# plt.xticks(np.arange(0, len(Summed_Distance)+1, step=1))
# plt.yticks(np.arange(0, 1, step=1))


# plt.show()



################################################################################################################################################
##### Second Plot  #### Barplot with alternative trees


fig, axs = plt.subplots(2)


# fig.suptitle('Number Trees in GTD per nubmer of proteins used ', fontsize=14, fontweight='bold')


# ax.set_xlabel('Number of Genes/Proteins', fontsize=11, fontweight='bold')
# ax.set_ylabel('Percentage of GTD', fontsize=11, fontweight='bold')

# axs[0].bar([x for x in range(1,Max_Number_of_Genes+1)], Summed_Distance, color ='maroon', width = 0.8)




# plt.xticks(np.arange(0, len(Summed_Distance)+1, step=1))
# plt.yticks(np.arange(0, 1, step=1))


# plt.show()



################################################################################################################################################
##### Fourth Plot #### Stacked - Barplot with percentage of different trees
#### Using example from here https://www.shanelynn.ie/bar-plots-in-python-using-pandas-dataframes/
#### stylistic options: https://python-charts.com/part-whole/stacked-bar-chart-matplotlib/





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

plotdata2.plot(kind='bar',legend=False,width=0.8,color='maroon',ax=axs[0])
plotdata.plot(kind='bar', stacked=True,legend=False,width=0.8,color=my_colors0,ax=axs[1])
### , 

plt.figure(figsize=(18,5))
plt.show()