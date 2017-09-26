
# All possible combinations of 0 and 1 with a given length n (i.e.: 2**n )
def all_possible_haps(length):
    possible_haps = list(map(''.join, product(map(str, {0,1}), repeat=length)))
    return possible_haps
