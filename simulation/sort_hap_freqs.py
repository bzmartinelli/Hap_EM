# python sort_hap_freqs.py simulated_freq_file estimated_freq_file ID_sim

# retrieves only the estimated haplotypes corresponding to the simulated haplotypes
from sys import *
from numpy import *


Sf = open(argv[1])
Ef = open(argv[2])

header = Sf.readline()
simfreq = {}   # {id:{hap:freq, hap:freq}}

for line in Sf:
    fields = line.replace('\n','').split('\t')
    #    simfreq.setdefault(fields[2]+str('_rand'), {}).setdefault(fields[0],[]).append(float(fields[1]))
    simfreq.setdefault(fields[2], {}).setdefault(fields[0],[]).append(float(fields[1]))

header = Ef.readline()
estfreq = {}
for line in Ef:
    fields = line.replace('\n','').split('\t')
    estfreq.setdefault(fields[2], {}).setdefault(fields[0],[]).append(float(fields[1]))

outALL = open("Frequencies_"+argv[3]+"_sorted", "w") # argv[3] is the id
outALL.write('Haplotype\tSimulated_Frequency\tEstimated_Frequency\tSim_id\n')

counter = 0
for k, v in simfreq.iteritems():
    for k2 in v.keys():
        if k2 in estfreq[k]:
            outALL.write(str(k2)+'\t'+str(v[k2]).strip("[]")+'\t'+str(estfreq[k][k2]).strip("[]")+'\t'+ str(k)+'\n')
        else:
            counter += 1

            outALL.write(str(k2)+'\t'+str(v[k2]).strip("[]")+'\t'+str(0).strip("[]")+'\t'+ str(k)+'\n')

#print 'Number of haplotypes not found: ', counter

















