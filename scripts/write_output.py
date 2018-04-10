
# output file containing the haplotype frequencies and the entropy
def output_freq_entropy(hap_freq, threshold_min_freq, entropy, ID_gr, wstart, wend):
#    out = open("EM_temp", "w")
    out = open("temp1.txt", "w")
    for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        if hap_freq[h] > threshold_min_freq:
            out.write(str(ID_gr)+'\t'+str(wstart)+'\t'+str(wend)+'\t'+ str(h) +'\t'+ str(float(v)) +'\n')
    out2 = open("temp2.txt", "w")
    out2.write(str(ID_gr)+'\t'+str(wstart)+'\t'+str(wend)+'\t'+ str(entropy) +'\n')




# output file containing the proportion of each CpG based on the calculation using the sequence reads or the haplotype frequencies
def output_meth(proportion_from_reads, proportion_from_haps, ID_gr, wstart, wend):
    out3 = open("temp3.txt", "w")    
    kr = sorted(proportion_from_reads.keys())
    vh = []
    for k in sorted(proportion_from_haps.keys()):
        vh.append(proportion_from_haps[k])
    proportion_from_haps = {}
    for i, z in enumerate(kr):
        proportion_from_haps[z] = vh[i]
    for cpg1 in proportion_from_reads.keys():
        for cpg2 in proportion_from_haps.keys():
            if cpg1 == cpg2:
                out3.write(str(ID_gr)+'\t'+str(wstart)+'\t'+str(wend)+'\t'+ str(cpg1)+'\t'+str(proportion_from_reads[cpg1])+'\t'+str(proportion_from_haps[cpg2])+'\n')
