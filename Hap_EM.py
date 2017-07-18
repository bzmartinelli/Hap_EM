# implementation of an EM algorithm to estimate the frequencies of methylation haplotypes
# python Hap_EM.py reads window_range


from sys import *
from numpy import *
from matplotlib import pyplot
from Hap_EM_Functions import *


max_its = 1000  # number of iterations of EM
epsilon = 0.0001
min_threshold = 0.0001  # minimum threshold to display haplotype frequency
start_coordinate = int(argv[2])
window = int(argv[3])

#o = open(argv[1]+"temp","w")

# Reading the file:
all_reads = {}
all_positions = {}
meth_status = {}  # dict {position: {read: meth}, }
with open(argv[1]) as file:
    for line in file:
        fields = line.replace('\n','').split('\t')
        position = int(fields[0])
        if (position >= start_coordinate and position < (start_coordinate + window)):
            all_positions[position] = 1
            all_reads[fields[2]] = 1
            meth_status.setdefault(int(fields[0]), {}).setdefault(fields[2], []).append(str(fields[1]))

n_cpgs = len(all_positions.keys())
print 'Number of CpG sites: ', n_cpgs

#


# array with all possible haplotypes
haplotypes = all_possibilities(len(all_positions))
for i, h in enumerate(haplotypes):
    haplotypes[i] = ''.join(h)


# array with all individual reads, represented as patterns
reads = read_patterns(all_reads, all_positions, meth_status)


# counting how many reads correspond to each read pattern
read_count_patterns = read_count(reads)


# find the reads consistent with each possible haplotype
consistent_hap = {}
for r in read_count_patterns.keys():
    consistent_hap[r] = {}
    for h in haplotypes:
        consistent_hap[r][h] = consistent(r,h)



# initializing dictionaries and arrays to store and update the values from each iteration of the EM

hap_freq = {}
for h in haplotypes:
    hap_freq[h] = 1.0/len(haplotypes)

expected_counts = {}
for k in consistent_hap.keys():
    expected_counts[k] = {}
    for v in consistent_hap[k].keys():
        expected_counts[k][v] = 0

likelihood = {}
for r in read_count_patterns.keys():
    likelihood[r] = 0

all_logL = {}
logLarray = []


print 'Calculating...\t'

# EM

for its in range(max_its):
    
    # E step
    for r in read_count_patterns.keys():
        total_consistent = 0
        for h in consistent_hap[r].keys():
            if consistent_hap[r][h] == 1:
                total_consistent += hap_freq[h] # sum of the frequencies of the haps consistent with the read
        for h in hap_freq.keys():
            if consistent_hap[r][h] == 1:
                expected_counts[r][h] = float(read_count_patterns[r]) * hap_freq[h] / total_consistent

    # M step
    for h in hap_freq.keys():
        hap_freq[h] = 0
        for r in read_count_patterns.keys():
            hap_freq[h] += float(expected_counts[r][h]) / len(all_reads)


    # log Likelihood
    logL = 0
    for r in reads: # the loop is calling each individual read for the likelihood calculation
        total_consistent = 0
        for h in consistent_hap[r].keys():
            if consistent_hap[r][h] == 1:
                total_consistent += hap_freq[h]
        likelihood[r] = total_consistent
        logL += log(likelihood[r])
    logLarray.append(logL)
    all_logL[its] = logL  # dict with the log likelihood in each iteration, used to check the convergence

    # checking convergence
    converged = 0
    counter = 0
    if its > 0:
        counter += 1
        deltaLL = logL - all_logL[its-1] # likelihood of the current iteration minus the likelihood of the previous iteration
        if fabs(deltaLL) < epsilon: # compares the deltaLL with the convergence criteria
            converged += 1
        if converged == counter:
            print '\nConverged! iteration:', its
            break
        elif its == (max_its - 1):
            print 'Convergence problem'

print '\nlog-Likelihood:', logL
print '-'*16, "\nHaplotype\tEstimated frequency\n"



#out = open(argv[1]+"_HapsFreqs", "w")
#out.write('ID\tHaplotypes\tFrequencies\n')


# entropy calculation
entropy = 0
for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    if hap_freq[h] > min_threshold:
        print h, '\t', hap_freq[h]
#        out.write(argv[1] +'\t'+ str(h) +'\t'+ str(v)+'\n')
        log2f = float(log(hap_freq[h])) / log(2)
        entropy = entropy - hap_freq[h] * log2f

N = n_cpgs
entropy = float((1.0 / N)) * entropy
print '\nEntropy:', entropy, '\n'



#Plotting the log-likelihood
#pyplot.title('logLikelihood')
#pyplot.xlabel('Iterations')
#pyplot.ylabel('logLikelihood')
#pyplot.plot(logLarray)
#pyplot.show()





# Methylation Frequency calculation

#methfreq = open(argv[1]+"_MethFreq", "w")
#methfreq.write('ID\tPosition\tFreqA 0 \tFreqA 1 \tFreqB 0 \tFreqB 1 \n')
#status_proportion = {} #meth proportion in each position
#for cpg in all_positions.keys():
#    status_proportion[cpg] = []
#    m_u = []
#    m_u = meth_status[cpg].values()
#    status_proportion[cpg] = [(m_u.count(['0'])/float(len(m_u))),(m_u.count(['1'])/float(len(m_u)))]
#    methfreq.write(argv[1]+'\t'+str(cpg)+'\t'+str((m_u.count(['0'])/float(len(m_u))))+'\t'+str((m_u.count(['1'])/float(len(m_u))))+'\n')
#
#
#high_freq = sorted(hap_freq.values(), reverse = True)
#high_freq = high_freq[:2]
#o.write(str(N) +'\t'+ str(high_freq[0]) +'\t'+ str(high_freq[1]) +'\t'+ str(entropy))
#o.close()
#out.close()




