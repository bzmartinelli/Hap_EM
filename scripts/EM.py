from __future__ import print_function
from numpy import *

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
def convergence(its, max_iterations, logL, all_logL, convergence_threshold):
    converged = 0
    counter = 0
    if its > 0:
        counter += 1
        deltaLL = logL - all_logL[its-1] # Likelihood of the current iteration minus the likelihood of the previous iteration
        if fabs(deltaLL) < convergence_threshold: # Compares the deltaLL with the convergence criteria
            converged += 1
        if converged == counter: # If the algorithm converges before the maximun number of iterations, the code breaks and display the results
            print('Converged!') #'iteration:', its
            return True
    elif its == (max_iterations - 1): # When the algorithm does not converge after running the given number of iterations
            print('Not converged')


# Print the haplotypes and their frequencies
#def haplotypes_frequencies(hap_freq, min_freq_cutoff):
#    print '-'*16, "\nHaplotype\tEstimated frequency\n"
#    for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
#        if hap_freq[h] > min_freq_cutoff:
#            print h, '\t', hap_freq[h]
#
#




