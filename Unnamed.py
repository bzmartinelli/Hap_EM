
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


max_its = 1000  # number of iterations of EM
epsilon = 0.00001 # convergence criteria
min_threshold = 0.0001  # minimum threshold to display haplotype frequency
reads_file = argv[1] # input file
start_coordinate = int(argv[2])
window = int(argv[3])


# reading the input file:
all_reads, all_positions, meth_status, n_cpgs = read_input_file(argv[1], start_coordinate, window)

# array with all possible haplotypes
haplotypes = all_possible_haps(len(all_positions))

# array with all individual reads, represented as patterns
reads = read_patterns(all_reads, all_positions, meth_status)

# counting how many reads correspond to each read pattern
read_count_patterns = read_count(reads)

# find the reads consistent with each possible haplotype
consistent_hap = consistent_h(read_count_patterns, haplotypes)

# initializing dictionaries and arrays to store and update the values from each iteration of the EM
hap_freq = hap_freq(haplotypes)
expected_counts = expected_counts(consistent_hap)
likelihood, all_logL, logLarray = likelihood(read_count_patterns)


# EM
print 'Calculating...\t'

for its in range(max_its):
    EM(read_count_patterns, consistent_hap, hap_freq, expected_counts, all_reads)
    likelihood, logL, all_logL, logLarray = log_Likelihood(reads, consistent_hap, likelihood,
                                                           hap_freq, logLarray, all_logL, its)
    if convergence(its, max_its, logL, all_logL, epsilon) == True:
        break


# print the haplotypes and their frequencies
haplotypes_frequencies(hap_freq, min_threshold)

# entropy
entropy = entropy(hap_freq, min_threshold, n_cpgs)

# proportion of methylated CpGs
# proportion_from_reads, proportion_from_haps = meth_freqs(all_positions, meth_status, hap_freq, min_threshold, argv[1])

# write output file containing the haplotype frequencies and entropy
output_freq_entropy(hap_freq, min_threshold, entropy, argv[1])
# output_meth(proportion_from_reads, proportion_from_haps, argv[1])





