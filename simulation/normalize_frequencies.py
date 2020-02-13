# python normalize_frequencies.py freqs_sorted_file IDdepth

# normalize haplotype frequencies
from sys import *
from numpy import *


freqfile = open(argv[1])
depth = str(argv[2]) 

header = freqfile.readline()
simfreq = {}   # {id:{hap:freq, hap:freq}}
estfreq = {}

for line in freqfile:
    fields = line.replace('\n','').split('\t')
    simfreq.setdefault(fields[3], {}).setdefault(fields[0],[]).append(float(fields[1]))
    estfreq.setdefault(fields[3], {}).setdefault(fields[0],[]).append(float(fields[2]))

for k, v in estfreq.items():
    for k2 in v.keys():
        v[k2] = float(str(v[k2]).strip("[]"))



# for k, v in estfreq.items():
#     print("estfreq[k].values()",estfreq[k].values(), "\n")
#     hapfreqs_sum = sum(estfreq[k].values())
#     print("hapsfreqsum",hapfreqs_sum, "\n")
#     #print(k,"\t", v, "\n")
#     for k2 in v.keys():
# #        if hapfreqs_sum != 0:
#         v[k2] = float(str(v[k2]).strip("[]"))/hapfreqs_sum  # float(str(v[k2]).strip("[]")) transform [0.2] in 0.2

for k, v in estfreq.items():
    for k2 in v.keys():
        v[k2] = float(str(v[k2]).strip("[]"))

for k, v in estfreq.items():
    hapfreqs_sum = sum(list(estfreq[k].values()))
#    print("\nhapsfreqsum",hapfreqs_sum, "\n")
    for k2 in v.keys():
 #       print("k2", k2, "\n")
#        if hapfreqs_sum != 0:
        v[k2] = float(v[k2])/hapfreqs_sum  # float(str(v[k2]).strip("[]")) transform [0.2] in 0.2
 #       print("vk2", v[k2],"\n")

outALL = open("Frequencies_sorted_normalized_depth"+depth, "w")
outALL.write('Haplotype\tSimulated_Frequency\tEstimated_Frequency\tSim_id\n')

for k, v in estfreq.items():
    for k2 in v.keys():
        outALL.write(str(k2)+'\t'+str(simfreq[k][k2]).strip("[]")+'\t'+str(v[k2]).strip("[]")+'\t'+ str(k)+'\n')










