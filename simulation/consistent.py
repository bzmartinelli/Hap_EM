
# Determine the haplotypes consistent with each read by comparing the read pattern with the possible haplotype
def consistent(r,h):
    # converts each read and haplotype into list
    list_read = list(r)
    list_hap = list(h)
    counter = 0
    # Compares each position of the read to the correspondent position of the haplotype
    for i in range(len(list_read)):
        if list_read[i] == list_hap[i]:
            counter += 1
    # The haplotype is consistent with the read if the length of the read minus the number of unknown CpG status
    #   in the read is equal to the number of matched positions
    if (len(list_read) - list_read.count('2')) == counter:
        return 1
    else:
        return 0

# Find the haplotypes consistent with each read callig the function 'consistent', giving a read pattern and an haplotype as arguments
def consistent_h(read_count_patterns, haplotypes):
    consistent_hap = {}
    for r in read_count_patterns.keys():
        consistent_hap[r] = {}
        for h in haplotypes:
            consistent_hap[r][h] = consistent(r,h)
    return consistent_hap


