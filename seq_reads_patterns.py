
# Creating patterns of reads with the same length of the haplotypes
def read_patterns(all_reads, all_positions, meth_status):
    reads = []
    for read in all_reads.keys():
        read_binary = ""
        for cpg in all_positions.keys():
            if read not in meth_status[cpg].keys(): # read is not covering this position
                read_binary += '2'
            elif meth_status[cpg][read] in (['0'],['u'],['z']):
                read_binary += '0'
            elif meth_status[cpg][read] in (['1'],['m'],['Z']):
                read_binary += '1'
        reads.append(read_binary)
    return reads


# Counting how many reads correspond to each read pattern
def read_count(reads):
    read_count_patterns = {} # key = read pattern, value = number of reads
    for r in reads:
        if r not in read_count_patterns.keys():
            read_count_patterns[r] = 1
        else:
            read_count_patterns[r] += 1
    return read_count_patterns
