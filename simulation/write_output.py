
# output file containing the haplotype frequencies and the entropy
def output_freq_entropy(hap_freq, threshold_min_freq, entropy, id_out):
    out = open(id_out+"_Estimated_Frequencies", "w")
    out.write('Haplotype\tFrequency\tSim_id\n')
    for h, f in sorted(hap_freq.items(), key=lambda x: x[1], reverse=True):
#            out.write(str(h) +'\t'+ str(v)+'\n')
    #for h, v in hap_freq.items():
        if hap_freq[h] > threshold_min_freq:
            out.write(str(h) +'\t'+ str(float(f)) +'\t'+ str(id_out)+'\n')
    
    out2 = open(id_out+"_Estimated_Entropy", "w")
#    out2.write('Simulated_Entropy+'\t'+ str(file_argv) +\n')
    out2.write(str(entropy) +'\t'+ str(id_out) +'\n')




# output file containing the proportion of each CpG based on the calculation using the sequence reads or the haplotype frequencies
def output_meth(proportion_from_reads, proportion_from_haps, id_out):
#    methout = open(file_argv+"_Proportion_Methylated_CpGs", "w")
#    methout.write('Position\tFrom_reads\tFrom_haplotypes\n')
    out3 = open(id_out+"_Estimated_Meth_Proportion", "w")
    for cpg1 in proportion_from_reads.keys():
        for cpg2 in proportion_from_haps.keys():
            if cpg1 == cpg2:
                out3.write(str(cpg1)+'\t'+str(proportion_from_reads[cpg1])+'\t'+str(proportion_from_haps[cpg2]) +'\t'+ str(id_out)+'\n')

