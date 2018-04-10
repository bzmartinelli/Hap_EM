#  python call_EM.py $outfile 0 $pos_cpgs_hap sim$i



import os
from sys import *


reads_file = argv[1]
pos_start = argv[2]
pos_end = argv[3]
simID = str(argv[4])


o = open(simID+"EM_temp","w")
o2 = open(simID+"_EntropyTemp","w")
o5 = open(simID+"_MethTemp","w")

print "Window:", pos_start, ' ', pos_end
cLine = "python MHap.py -f "+ reads_file +" -start_w "+ pos_start +" -len_w  "+ pos_end +" -id "+ str(simID)
os.system(cLine)

print '-'*45, '\n'
o.close()

o = open(simID+"EM_temp")
o3 = open(simID+"_Estimated_Frequencies","w")
o3.write('Haplotype\tFrequency\tSim_id\n')
outputFile = o.readlines()
for line in outputFile:
    elem = line.replace('\n','').split('\t')
    print >> o3, line.replace('\n','')
o3.close()

o2 = open(simID+"_EntropyTemp")
o4 = open(simID+"_Estimated_Entropy","w")
o4.write('Estimated_Entropy\tSim_id\n')
outputFile2 = o2.readlines()
for line in outputFile2:
    elem = line.replace('\n','').split('\t')
    print >> o4, line.replace('\n','')
o4.close()


o5 = open(simID+"_MethTemp")
o6 = open(simID+"_Estimated_Meth_Proportion","w")
o6.write('Position\tFrom_reads\tFrom_haplotypes\tSim_id\n')
outputFile5 = o5.readlines()
for line in outputFile5:
    elem = line.replace('\n','').split('\t')
    print >> o6, line.replace('\n','')
o5.close()


os.remove(simID+"EM_temp")
os.remove(simID+"_EntropyTemp")
os.remove(simID+"_MethTemp")

