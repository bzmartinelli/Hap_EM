



# Functions

# function returning array containing all possible combinations of 0 and 1 with a given length n (i.e.: 2**n )
def all_possibilities(length):
    alphabet = ['0','1']
    c = [[]]
    for i in range(length):
        c = [[x]+y for x in alphabet for y in c]
    return c


# function creating patterns of reads with the same length of the haps
def read_patterns(all_reads, all_positions, meth_status):
    reads = []
    for read in all_reads.keys():
        read_binary = ""
        for cpg in all_positions.keys():
            if read not in meth_status[cpg].keys():
                read_binary += '2'
            elif meth_status[cpg][read] == ['0']:
                read_binary += '0'
            elif meth_status[cpg][read] == ['1']:
                read_binary += '1'
        reads.append(read_binary)
    return reads

# function to find the haplotypes consistent with each read
def consistent(r,h):
    list_read = list(r)
    list_hap = list(h)
    counter = 0
    for i in range(len(list_read)):
        if list_read[i] == list_hap[i]:
            counter += 1
    if (len(list_read) - list_read.count('2')) == counter:
        return 1
    else:
        return 0

# counting how many reads correspond to each read pattern
def read_count(reads):
    read_count_patterns = {} # key = read pattern, value = number of reads
    for r in reads:
        if r not in read_count_patterns.keys():
            read_count_patterns[r] = 1
        else:
            read_count_patterns[r] += 1
    return read_count_patterns


