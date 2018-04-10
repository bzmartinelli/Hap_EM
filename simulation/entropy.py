from numpy import *


# Shannon entropy is used to calculate the methylation entropy
def entropy(hap_freq, threshold_min_freq, n_cpgs):
    entropy = 0
    for h, v in sorted(hap_freq.iteritems()):
        if hap_freq[h] > threshold_min_freq:
            log2f = float(log(hap_freq[h])) / log(2)
            entropy = entropy - hap_freq[h] * log2f

    entropy = float((1.0 / n_cpgs)) * entropy
#    print "---n_cpgs: ", n_cpgs
    print '\nEntropy:', entropy, '\n'
    return entropy
