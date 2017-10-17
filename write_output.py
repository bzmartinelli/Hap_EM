
# output file containing the haplotype frequencies and the entropy
def output_freq_entropy(hap_freq, min_threshold, entropy, file_argv):
    out = open(file_argv+"_Frequencies_and_Entropy", "w")
    out.write('Haplotype\tFrequency\n')
    for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        if hap_freq[h] > min_threshold:
            out.write(str(h) +'\t'+ str(v)+'\n')
    out.write('\nEntropy:' +' '+ str(entropy))



# # output file containing the proportion of each CpG based on the calculation using the sequence reads or the haplotype frequencies
# def output_meth(proportion_from_reads, proportion_from_haps, file_argv):
#     methout = open(file_argv+"_Proportion_Methylated_CpGs", "w")
#     methout.write('Position\tFrom reads\tFrom haplotypes\n')
#     for cpg1 in proportion_from_reads.keys():
#         for cpg2 in proportion_from_haps.keys():
#             if cpg1 == cpg2:
#                 methout.write(str(cpg1)+'\t'+str(proportion_from_reads[cpg1])+'\t'+str(proportion_from_haps[cpg2])+'\n')


