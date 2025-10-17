OUT_PUT_FILE=open('Input_File.txt','w')

Number_of_max_Proteins=28
Number_of_max_Reps=1000

for K in range(1,Number_of_max_Proteins+1):
    OUT_PUT_FILE.write(F'{K}\t{Number_of_max_Reps}\n')