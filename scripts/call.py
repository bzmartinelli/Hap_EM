from __future__ import print_function
import os
from sys import *
import argparse
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', required=True)
parser.add_argument('-wf', '--windows_file', required=True )
parser.add_argument('-gr', '--genomic_regions', required=True )
parser.add_argument('-data', '--cpg_reads_data', required=True)
parser.add_argument('-data_from', '--data_from', required=True)
parser.add_argument('-w', '--window_size', type=int, default=1000)
parser.add_argument('-ncpgs', '--number_of_cpgs', type=int, default=5)
parser.add_argument('-mmin', '--meth_min', type=float, default=0)
parser.add_argument('-mmax', '--meth_max', type=float, default=100)
parser.add_argument('-max_iter', '--max_iterations', type=int, default=1000)
parser.add_argument('-conv', '--convergence_criteria', type=float, default=0.000001)
parser.add_argument('-freq_cutoff', '--threshold_min_freq', type=float, default=0.0001)
parser.add_argument('-initial_freq', '--initial_freq_em', type=str, required=False)
parser.add_argument('-cwd', '--cwd')
args = parser.parse_args()

command_line = " ".join(argv[1:])
file = args.windows_file
with open(file) as wf:
    pw = []
    for line in wf:
        w = line.replace('\n','').split('\t')
        pw.append(w)
position_windows = pw

o = open("EM_temp","w")
o2 = open("EntropyTemp","w")
o5 = open("MethTemp","w")
o7 = open("MHLTemp","w")

for position in range(len(position_windows)):
    ID = position_windows[position][0] #id genomic region
    start_w = position_windows[position][1]
    end_w = position_windows[position][2]
    os.system('python run_MHap.py ' +command_line+ ' -ID '+ ID +' -start_w ' + start_w + ' -end_w '+ end_w)
    t = open("temp1.txt")
    for line in t:
        o.write(line)
    t2 = open("temp2.txt")
    for line2 in t2:
        o2.write(line2)
    t3 = open("temp3.txt")
    for line3 in t3:
        o5.write(line3)
    t4 = open("temp4.txt")
    for line4 in t4:
        o7.write(line4)

o.close()
o2.close()
o5.close()
o7.close()

o = open("EM_temp")
o3 = open("Estimated_Haplotype_Frequencies","w")
o3.write('Genomic_region\tWindow_start\tWindow_end\tHaplotype\tFrequency\n')
outputFile = o.readlines()
for line in outputFile:
    elem = line.replace('\n','').split('\t')
    print(line.replace('\n',''), file=o3)
o3.close()

o2 = open("EntropyTemp")
o4 = open("Estimated_Entropy","w")
o4.write('Genomic_region\tWindow_start\tWindow_end\tEstimated_Entropy\n')
outputFile2 = o2.readlines()
for line in outputFile2:
    elem = line.replace('\n','').split('\t')
    print(line.replace('\n',''), file=o4)
o4.close()
o4.close()

o5 = open("MethTemp")
o6 = open("Estimated_Methylation_Proportion","w")
o6.write('Genomic_region\tWindow_start\tWindow_end\tCpG_position\tFrom_reads\tFrom_haplotypes\n')
outputFile5 = o5.readlines()
for line in outputFile5:
    elem = line.replace('\n','').split('\t')
    print(line.replace('\n',''), file=o6)
o5.close()
o6.close()

o7 = open("MHLTemp")
o8 = open("MHL","w")
o8.write('Genomic_region\tWindow_start\tWindow_end\tMHL\n')
outputFile7 = o7.readlines()
for line in outputFile7:
    elem = line.replace('\n','').split('\t')
    print(line.replace('\n',''), file=o8)
o7.close()
o8.close()



d = str(datetime.datetime.now())
d = d.replace(' ', '_')
out_path = args.cwd

os.system('mkdir ' + out_path+'/MHap_output_'+d)
os.system('mv Estimated_Methylation_Proportion '+out_path+'/MHap_output_'+d)
os.system('mv Estimated_Entropy ' + out_path+'/MHap_output_'+d)
os.system('mv Estimated_Haplotype_Frequencies ' + out_path+'/MHap_output_'+d)
os.system('mv MHL ' + out_path+'/MHap_output_'+d)

os.system("rm temp*.txt EM_temp EntropyTemp MethTemp MHLTemp windows new_cpg_reads *.pyc")
