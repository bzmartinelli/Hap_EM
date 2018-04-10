
# Initializing dictionaries and arrays to store and update the values from each iteration of the EM
from numpy import *

def hap_freq(haplotypes, initial_freq_em):
    hap_freq = {}
    if initial_freq_em == 'random': # The initial frequencies are random, determined using Dirichlet distribution
        randfreq = random.dirichlet(ones(len(haplotypes)),size=1)
        for i, h in enumerate(haplotypes):
            hap_freq[h] = randfreq[0][i]
    else:
        for h in haplotypes:
            hap_freq[h] = 1.0/len(haplotypes) # The initial haplotype frequencies are the same for all haplotypes
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





















