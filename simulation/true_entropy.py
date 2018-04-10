# python entropy.py hap_freq

from sys import *
from numpy import *
from decimal import Decimal


hap_freq_file = open(argv[1])
SimID = (argv[2])

haps = {}
haps_dec = {}
for line in hap_freq_file:
    fields = line.replace('\n','').split('\t')
    haps[fields[0]] = float(fields[1])
    haps_dec[fields[0]] = Decimal(fields[1])

entropy = 0
for f in haps.keys():
    log2f = float(log(haps[f])) / log(2)
    entropy = entropy - haps[f] * log2f

N = len(list(f)) #number of cpgs
#print "-----N: ", N
entropy = float((1.0 / N)) * entropy

out = open(SimID+"_Simulated_Entropy", "w")
out.write('Simulated_entropy\tSim_id\n')
out.write(str(entropy)+'\t'+str(SimID)+'\n')
#print '\nEntropy:', entropy, '\n'


#####

for hap in sorted(haps_dec.keys()):
    len_hap = len(hap)
    fd = {key: 0 for key in range(len_hap)}

for hap in sorted(haps_dec.keys()):
    len_hap = len(hap)
    get_hap = list(hap)
    for position in range(len_hap):
        if get_hap[position] == '1':
            fd[position] += haps_dec[hap]

for k, v in fd.iteritems():
    fd[k] = float(v)


out2 = open(SimID+"_Simulated_Meth_Proportion", 'w')
out2.write('Position\tProportion\tSim_id\n')
for k, v in fd.iteritems():
    out2.write(str(k)+'\t'+str(v)+'\t'+str(SimID)+'\n')
##

out3 = open(SimID+"_Simulated_Frequencies", "w")
out3.write('Haplotype\tSimulated_Frequency\tSim_id\n')
for h, f in sorted(haps.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    out3.write(str(h)+'\t'+str(f)+'\t'+str(SimID)+'\n')



