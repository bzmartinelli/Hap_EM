

from numpy import *
from itertools import *

# Functions

# Reading the file:
def read_input_file(file_argv, start_coordinate, window):
    all_positions = {} # all CpGs as keys
    all_reads = {} # all reads as keys
    meth_status = {}  # CpGs as keys and the values consist of another dictionary where the reads are keys and
                      # the values are the methylation status of the CpG. > dict {position: {read: meth}, }
    with open(file_argv) as file:
        for line in file:
            fields = line.replace('\n','').split('\t')
            position = int(fields[0]) # each position corresponds to one CpG
            if (position >= start_coordinate and position < (start_coordinate + window)):
                # if the CpG is in the genomic locus of interest, those 3 dictionaries initialized above are created.
                all_positions[position] = 1
                all_reads[fields[2]] = 1
                meth_status.setdefault(int(fields[0]), {}).setdefault(fields[2], []).append(str(fields[1]))
    n_cpgs = len(all_positions.keys()) # total number of CpGs
    print 'Number of CpG sites: ', n_cpgs
    return all_reads, all_positions, meth_status, n_cpgs








