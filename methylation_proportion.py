# Methylation Frequency calculation


# Calculate the proportion of the CpGs in a methylated state.
def meth_freqs(all_positions, meth_status, hap_freq, min_threshold, file_argv):
    methfreq = open(file_argv+"_Proportion_Methylated_CpGs", "w")
    methfreq.write('Position\tFrom reads\tFrom haplotypes\n')
    
    # Proportion of methylated CpGs calculated using the sequence reads
    proportion_from_reads = {} # proportion of methylated CpGs in each position
    for cpg in all_positions.keys():
        meth_states = []
        meth_states = meth_status[cpg].values() # all methylation states of each CpG
        proportion_from_reads[cpg] = (meth_states.count(['1'])/float(len(meth_states)))

    # Proportion of methylated CpGs calculated using the haplotype frequencies
    proportion_from_haps = {}
    for hap in sorted(hap_freq.keys()):
        len_hap = len([hap][0])
        meth_proportion = {key: 0.0 for key in range(len_hap)} # dictionary to store the methylation proportion for each CpG in the haplotype
    for hap in sorted(hap_freq.keys()):
        if hap_freq[hap] > min_threshold:
            len_hap = len([hap][0])
            hap_cpg = list(hap)
            for position in range(len_hap): # get the methylation state of each CpG in the haplotype
                if hap_cpg[position] == '1':
                    meth_proportion[position] += float(hap_freq[hap]) # sum of the frequencies of all haplotypes containing the given CpG in a methylated state
            proportion_from_haps = meth_proportion

    # Writing the output file containing the proportion of methylated CpGs based on the calculation using the sequence reads or the haplotype frequencies
    for cpg1 in proportion_from_reads.keys():
        for cpg2 in proportion_from_haps.keys():
            if cpg1 == cpg2:
                methfreq.write(str(cpg1)+'\t'+str(proportion_from_reads[cpg1])+'\t'+str(proportion_from_haps[cpg2])+'\n')


