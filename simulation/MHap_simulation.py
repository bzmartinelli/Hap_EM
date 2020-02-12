import os
from sys import *
from numpy import *

#cwd = os.getcwd()


## Run simulations and edit post files

for i in range(10, 30, 10):
    i = str(i)
    os.system("perl simulation.pl 3windows5cpgs  "+i)

for i in range(10, 30, 10):
    i = str(i)
    os.system("head -1 D"+i+"_sim_1_Estimated_Frequencies > Estimated_Frequencies_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Estimated_Frequencies >> Estimated_Frequencies_depth"+i)
    os.system("head -1 D"+i+"_sim_1_Estimated_Entropy > Estimated_Entropy_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Estimated_Entropy >> Estimated_Entropy_depth"+i)
    os.system("head -1 D"+i+"_sim_1_Estimated_Meth_Proportion > Estimated_Meth_Proportion_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Estimated_Meth_Proportion >> Estimated_Meth_Proportion_depth"+i)
    os.system("head -1 D"+i+"_sim_1_Simulated_Frequencies > Simulated_Frequencies_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Simulated_Frequencies >> Simulated_Frequencies_depth"+i)
    os.system("head -1 D"+i+"_sim_1_Simulated_Entropy > Simulated_Entropy_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Simulated_Entropy >> Simulated_Entropy_depth"+i)
    os.system("head -1 D"+i+"_sim_1_Simulated_Meth_Proportion > Simulated_Meth_Proportion_depth"+i)
    os.system("tail -n +2 -q D"+i+"_*_Simulated_Meth_Proportion >> Simulated_Meth_Proportion_depth"+i)

os.system("rm D*")


## Sort hap frequencies
for i in range(10, 30, 10):
    i = str(i)
    os.system("python sort_hap_freqs.py Simulated_Frequencies_depth"+i+" Estimated_Frequencies_depth"+i+" "+i)


## Normalize Frequencies
for i in range(10, 30, 10):
    i = str(i)
    os.system("python normalize_frequencies.py Frequencies_sorted_depth"+i+" "+i) 

## Normalize Entropy
for i in range(10, 30, 10):
    i = str(i)
    os.system("python normalize_entropy.py Frequencies_sorted_normalized_depth"+i+" "+i)
