

from numpy import *
from itertools import *


# All possible combinations of 0 and 1 with a given length n (i.e.: 2**n )
def all_possible_haps(length):
    possible_haps = list(map(''.join, product(map(str, {0,1}), repeat=length)))
    return possible_haps


# Creating patterns of reads with the same length of the haps
def read_patterns(all_reads, all_positions, meth_status):
    reads = []
    for read in all_reads.keys():
        read_binary = ""
        for cpg in all_positions.keys():
            if read not in meth_status[cpg].keys(): # read is not covering this position
                read_binary += '2'
            elif meth_status[cpg][read] in (['0'],['u'],['z']):
                read_binary += '0'
            elif meth_status[cpg][read] in (['1'],['m'],['Z']):
                read_binary += '1'
        reads.append(read_binary)
    return reads


# Counting how many reads correspond to each read pattern
def read_count(reads):
    read_count_patterns = {} # key = read pattern, value = number of reads
    for r in reads:
        if r not in read_count_patterns.keys():
            read_count_patterns[r] = 1
        else:
            read_count_patterns[r] += 1
    return read_count_patterns



# Determine the haplotypes consistent with each read by comparing the read pattern with the possible haplotype
def consistent(r,h):
    # converts each read and haplotype into list
    list_read = list(r)
    list_hap = list(h)
    counter = 0
    # Compares each position of the read to the correspondent position of the haplotype
    for i in range(len(list_read)):
        if list_read[i] == list_hap[i]:
            counter += 1
    # The haplotype is consistent with the read if the length of the read minus the number of unknown CpG status
    #   in the read is equal to the number of matched positions
    if (len(list_read) - list_read.count('2')) == counter:
        return 1
    else:
        return 0

# Find the haplotypes consistent with each read callig the function 'consistent', giving a read pattern and an haplotype as arguments
def consistent_h(read_count_patterns, haplotypes):
    consistent_hap = {}
    for r in read_count_patterns.keys():
        consistent_hap[r] = {}
        for h in haplotypes:
            consistent_hap[r][h] = consistent(r,h)
    return consistent_hap



# Initializing dictionaries and arrays to store and update the values from each iteration of the EM

def hap_freq(haplotypes):
    hap_freq = {}
    for h in haplotypes:
        hap_freq[h] = 1.0/len(haplotypes) # The initial haplotype frequency is the same for all haplotypes
    return hap_freq

def expected_counts(consistent_hap):
    expected_counts = {}
    for k in consistent_hap.keys():
        expected_counts[k] = {}
        for v in consistent_hap[k].keys():
            expected_counts[k][v] = 0
    return expected_counts

def likelihood(read_count_patterns):
    likelihood = {}
    for r in read_count_patterns.keys():
        likelihood[r] = 0
    all_logL = {}
    logLarray = []
    return likelihood, all_logL, logLarray





















