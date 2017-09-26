
def output(hap_freq, min_threshold, entropy, file_argv):
    out = open(file_argv+"_output", "w")
    out.write('Haplotype\tFrequency\n')
    for h, v in sorted(hap_freq.iteritems(), key=lambda (k,v): (v,k), reverse=True):
        if hap_freq[h] > min_threshold:
            out.write(str(h) +'\t'+ str(v)+'\n')
    out.write('\nEntropy:' +' '+ str(entropy))



