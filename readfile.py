

from numpy import *
from itertools import *

# Functions

# Reading the file:
def read_input_file(file_argv, start_coordinate, window):
    all_positions = {} # all CpGs as keys
    all_reads = {} # all reads as keys
    meth_status = {}  # CpGs as keys and the values consist of another dictionary where the reads are keys and
                      # the values are the methylation status of the CpG. > dict {position: {read: meth}, }
    with open(file_argv) as file:
        for line in file:
            fields = line.replace('\n','').split('\t')
            position = int(fields[0]) # each position corresponds to one CpG
            if (position >= start_coordinate and position < (start_coordinate + window)):
                # if the CpG is in the genomic locus of interest, those 3 dictionaries initialized above are created.
                all_positions[position] = 1
                all_reads[fields[2]] = 1
                meth_status.setdefault(int(fields[0]), {}).setdefault(fields[2], []).append(str(fields[1]))
    n_cpgs = len(all_positions.keys()) # total number of CpGs
    print 'Number of CpG sites: ', n_cpgs
    return all_reads, all_positions, meth_status, n_cpgs




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




# EM algorithm
def EM(read_count_patterns, consistent_hap, hap_freq, expected_counts, all_reads):
    # E step
    for r in read_count_patterns.keys():
        total_consistent = 0
        for h in consistent_hap[r].keys():
            if consistent_hap[r][h] == 1:
                total_consistent += hap_freq[h] # Sum of the frequencies of the haplotypes consistent with each read, the likelihood
        for h in hap_freq.keys():
            if consistent_hap[r][h] == 1:
                expected_counts[r][h] = float(read_count_patterns[r]) * hap_freq[h] / total_consistent
    # M step
    for h in hap_freq.keys():
        hap_freq[h] = 0
        for r in read_count_patterns.keys():
            hap_freq[h] += float(expected_counts[r][h]) / len(all_reads)
    return hap_freq


# log-Likelihood
def log_Likelihood(reads, consistent_hap, likelihood, hap_freq, logLarray, all_logL, its):
    logL = 0
    for r in reads: # the loop is calling each individual read for the likelihood calculation
        total_consistent = 0
        for h in consistent_hap[r].keys():
            if consistent_hap[r][h] == 1:
                total_consistent += hap_freq[h]
        # The likehood of the read is the sum of the frequencies of the haplotypes consistent with the read
        likelihood[r] = total_consistent
        logL += log(likelihood[r])
    all_logL[its] = logL  # Dictionary storing the log-Likelihood of each iteration, used to check the convergence
    #logLarray.append(logL)
    return likelihood, logL, all_logL, logLarray
    


# convergence
def convergence(its, max_its, logL, all_logL, epsilon):
    converged = 0
    counter = 0
    if its > 0:
        counter += 1
        deltaLL = logL - all_logL[its-1] # Likelihood of the current iteration minus the likelihood of the previous iteration
        if fabs(deltaLL) < epsilon: # Compares the deltaLL with the convergence criteria
            converged += 1
        if converged == counter: # If the algorithm converges before the maximun number of iterations, the code breaks and display the results
            print '\nConverged! iteration:', its
            return True
    elif its == (max_its - 1): # When the algorithm does not converge after running the given number of iterations
            print 'Convergence problem'


# Print the haplotypes and their frequencies
def haplotypes_frequencies(hap_freq, min_threshold):
    print '-'*16, "\nHaplotype\tEstimated frequency\n"
    for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        if hap_freq[h] > min_threshold:
            print h, '\t', hap_freq[h]
#        out.write(argv[1] +'\t'+ str(h) +'\t'+ str(v)+'\n')

# Shannon entropy is used to calculate the methylation entropy
def entropy(hap_freq, min_threshold, n_cpgs):
    entropy = 0
    for h, v in sorted(hap_freq.iteritems()):
        if hap_freq[h] > min_threshold:
            log2f = float(log(hap_freq[h])) / log(2)
            entropy = entropy - hap_freq[h] * log2f

    entropy = float((1.0 / n_cpgs)) * entropy
    print '\nEntropy:', entropy, '\n'
























