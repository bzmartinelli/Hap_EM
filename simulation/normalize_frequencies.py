# python normalize_frequencies.py freqs_sorted ID

# normalize haplotype frequencies
from sys import *
from numpy import *


freqfile = open(argv[1])

header = freqfile.readline()
simfreq = {}   # {id:{hap:freq, hap:freq}}
estfreq = {}

for line in freqfile:
    fields = line.replace('\n','').split('\t')
    simfreq.setdefault(fields[3], {}).setdefault(fields[0],[]).append(float(fields[1]))
    estfreq.setdefault(fields[3], {}).setdefault(fields[0],[]).append(float(fields[2]))

for k, v in estfreq.iteritems():
    hapfreqs_sum = sum(estfreq[k].values())
    for k2 in v.keys():
#        if hapfreqs_sum != 0:
        v[k2] = float(str(v[k2]).strip("[]"))/hapfreqs_sum  # float(str(v[k2]).strip("[]")) transform [0.2] in 0.2

outALL = open("Frequencies_"+argv[2]+"_sorted_normalized", "w")
outALL.write('Haplotype\tSimulated_Frequency\tEstimated_Frequency\tSim_id\n')


for k, v in estfreq.iteritems():
    for k2 in v.keys():
        outALL.write(str(k2)+'\t'+str(simfreq[k][k2]).strip("[]")+'\t'+str(v[k2]).strip("[]")+'\t'+ str(k)+'\n')









