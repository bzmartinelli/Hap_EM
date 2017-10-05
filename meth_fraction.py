# Methylation Frequency calculation



def meth_freqs(all_positions, meth_status, hap_freq, min_threshold, file_argv):
    methfreq = open(file_argv+"_MethFreq", "w")
    methfreq.write('Position\tFrom reads 1 \tFrom haps 1 \n')
    
    # FROM READS
    status_proportion = {} # meth proportion in each position
    for cpg in all_positions.keys():
        status_proportion[cpg] = []
        m_u = []
        m_u = meth_status[cpg].values()
        status_proportion[cpg] = (m_u.count(['1'])/float(len(m_u)))

    # FROM HAPLOTYPES
    new_fd = {}
    for h in sorted(hap_freq.keys()):
        len_hap = len([h][0])
        fd = {key: 0.0 for key in range(len_hap)}
    for h in sorted(hap_freq.keys()):
        if hap_freq[h] > min_threshold:
            len_hap = len([h][0])
            get_hap = list(h)
            for position in range(len_hap):
                if get_hap[position] == '1':
                    fd[position] += float(hap_freq[h])
            new_fd = fd

#   writing the output file containing the proportion of each CpG based on the calculation using the sequence reads 
#  or the haplotype frequencies
    for cpg1 in status_proportion.keys():
        for cpg2 in new_fd.keys():
            if cpg1 == cpg2:
                methfreq.write(str(cpg1)+'\t'+str(status_proportion[cpg1])+'\t'+str(new_fd[cpg2])+'\n')

