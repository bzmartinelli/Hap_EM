
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
from mhl import *

# argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', required=True)
parser.add_argument('-wf', '--windows_file', required=True )
parser.add_argument('-gr', '--genomic_regions', required=True )
parser.add_argument('-data', '--cpg_reads_data', required=True)
parser.add_argument('-data_from', '--data_from', required=True)
parser.add_argument('-start_w', '--start_coordinate', type=int, required=True)
parser.add_argument('-end_w', '--end_window', type=int, required=True)
parser.add_argument('-w', '--window_size', type=int, default=1000)
parser.add_argument('-ncpgs', '--number_of_cpgs', type=int, default=5)
parser.add_argument('-mmin', '--meth_min', type=float, default=0)
parser.add_argument('-mmax', '--meth_max', type=float, default=100)
parser.add_argument('-max_iter', '--max_iterations', type=int, default=1000)
parser.add_argument('-conv', '--convergence_threshold', type=float, default=0.000001)
parser.add_argument('-freq_cutoff', '--min_freq_cutoff', type=float, default=0.0001)
parser.add_argument('-initial_freq', '--initial_freq_em', type=str, required=False)
parser.add_argument('-ID', '--ID_genomic_region', type=str, required=True)
parser.add_argument('-cwd', '--cwd')
args = parser.parse_args()


# Reading the file:
all_reads, all_positions, meth_status = read_input_file(args.infile, args.start_coordinate, args.end_window, args.data_from)

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
print 'Calculating...'

for its in range(args.max_iterations):
    EM(read_count_patterns, consistent_hap, hap_freq, expected_counts, all_reads)
    likelihood, logL, all_logL, logLarray = log_Likelihood(reads, consistent_hap, likelihood,
                                                           hap_freq, logLarray, all_logL, its)
    if convergence(its, args.max_iterations, logL, all_logL, args.convergence_threshold) == True:
        break


# print out the haplotypes and their frequencies
#haplotypes_frequencies(hap_freq, args.min_freq_cutoff)

# entropy
entropy = entropy(hap_freq, args.min_freq_cutoff, args.number_of_cpgs)


# proportion of methylated CpGs
proportion_from_reads, proportion_from_haps = meth_freqs(all_positions, meth_status, hap_freq, args.min_freq_cutoff, args.infile)


# write output file containing the haplotype frequencies and entropy
haplos = output_freq_entropy(hap_freq, args.min_freq_cutoff, entropy, args.ID_genomic_region, args.start_coordinate, args.end_window)
output_meth(proportion_from_reads, proportion_from_haps, args.ID_genomic_region, args.start_coordinate, args.end_window)


# MHL
MHL(haplos, args.number_of_cpgs, args.ID_genomic_region, args.start_coordinate, args.end_window)


