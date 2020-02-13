# python normalize normentropy.py hap_frequencies_norm depth

from sys import *
from numpy import *


hap_freq_file = open(argv[1])
depth = str(argv[2])

haps = {}
haps_dec = {}
id_haps_f = {}
header = hap_freq_file.readline()
for line in hap_freq_file:
    fields = line.replace('\n','').split('\t')
    if float(fields[2]) != 0:
        id_haps_f.setdefault(fields[3], {}).setdefault(fields[0], []).append(float(fields[2]))


out = open("Estimated_Entropy_normalized_depth"+depth, "w")
out.write('Estimated_entropy\tSim_id\n')

for id, h in id_haps_f.items():
    entropy = 0
    for f in h.keys():
        hf = float(str(h[f])[1:-1])
        log2f = float(log2(hf))
        entropy = entropy - float(hf) * log2f
    N = len(list(f)) #number of cpgs
    entropy = float((1.0 / N)) * entropy
    out.write(str(entropy)+'\t'+str(id)+'\n')
