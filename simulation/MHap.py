# implementation of an EM algorithm to estimate the frequencies of methylation haplotypes
# python Hap_EM.py reads window_range

import argparse
from sys import *
from read_input_file import *
from possible_haplotypes import *
from seq_reads_patterns import *
from consistent import *
from initialize_EM import *
from EM import *
from entropy import *
from methylation_proportion import *
from write_output import *


#max_its = 1000  # number of iterations of EM
#epsilon = 0.00001 # convergence criteria
#min_threshold = 0.0001  # minimum threshold to display haplotype frequency
##reads_file = argv[1] # input file
#start_coordinate = int(argv[2])
#window = int(argv[3])

# argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', required=True)
parser.add_argument('-start_w', '--start_coordinate', type=int, required=True)
parser.add_argument('-len_w', '--length_window', type=int, required=True)
parser.add_argument('-max_iter', '--max_iterations', type=int, default=1000)
parser.add_argument('-conv', '--convergence_criteria', type=float, default=0.000001)
parser.add_argument('-freq_cutoff', '--threshold_min_freq', type=float, default=0.0001)
parser.add_argument('-initial_freq', '--initial_freq_em', type=str, required=False)
parser.add_argument('-id', '--id', type=str)
args = parser.parse_args()



# Reading the file:
all_reads, all_positions, meth_status, n_cpgs = read_input_file(args.infile, args.start_coordinate, args.length_window)

# array with all possible haplotypes
haplotypes = all_possible_haps(len(all_positions))

# array with all individual reads, represented as patterns
reads = read_patterns(all_reads, all_positions, meth_status)

# counting how many reads correspond to each read pattern
read_count_patterns = read_count(reads)

# find the reads consistent with each possible haplotype
consistent_hap = consistent_h(read_count_patterns, haplotypes)

# initializing dictionaries and arrays to store and update the values from each iteration of the EM
hap_freq = hap_freq(haplotypes, args.initial_freq_em)
expected_counts = expected_counts(consistent_hap)
likelihood, all_logL, logLarray = likelihood(read_count_patterns)


# EM
print 'Calculating...\t'

for its in range(args.max_iterations):
    EM(read_count_patterns, consistent_hap, hap_freq, expected_counts, all_reads)
    likelihood, logL, all_logL, logLarray = log_Likelihood(reads, consistent_hap, likelihood,
                                                           hap_freq, logLarray, all_logL, its)
    if convergence(its, args.max_iterations, logL, all_logL, args.convergence_criteria) == True:
        break


# print out the haplotypes and their frequencies
haplotypes_frequencies(hap_freq, args.threshold_min_freq)

# entropy
entropy = entropy(hap_freq, args.threshold_min_freq, n_cpgs)


# proportion of methylated CpGs
proportion_from_reads, proportion_from_haps = meth_freqs(all_positions, meth_status, hap_freq, args.threshold_min_freq, args.infile)


# write output file containing the haplotype frequencies and entropy
output_freq_entropy(hap_freq, args.threshold_min_freq, entropy, args.id)
output_meth(proportion_from_reads, proportion_from_haps, args.id)




#methylation_proportion(all_positions, meth_status, reads_file)

#print meth_status
###################################

#methfreq = open(argv[1]+"_MethFreq", "w")
#methfreq.write('ID\tPosition\tFreqA 0 \tFreqA 1 \tFreqB 0 \tFreqB 1 \n')
#
#status_proportion = {} #meth proportion in each position
#for cpg in all_positions.keys():
#    status_proportion[cpg] = []
#    m_u = []
#    m_u = meth_status[cpg].values()
#    status_proportion[cpg] = (m_u.count(['Z'])/float(len(m_u)))
##    methfreq.write(argv[1]+'\t'+str(cpg)+'\t'+str(status_proportion[cpg])+'\n')
#
#
#new_fd = {}
#for h in sorted(hap_freq.keys()):
#    len_hap = len([h][0])
#    fd = {key: 0.0 for key in range(len_hap)}
#
#for h in sorted(hap_freq.keys()):
#    if hap_freq[h] > min_threshold:
#        len_hap = len([h][0])
##        fd = {key: 0.0 for key in range(len_hap)}
#        get_hap = list(h)
#        for position in range(len_hap):
#            if get_hap[position] == '1':
#                fd[position] += float(hap_freq[h])
#        new_fd = fd
#
#
#for p in status_proportion.keys():
#    for p2 in new_fd.keys():
#        if p == p2:
#            methfreq.write(str(p)+'\t'+str(status_proportion[p])+'\t'+str(new_fd[p2])+'\n')
#



