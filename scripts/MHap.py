
import argparse
from sys import *
import os
from sliding_window import *
from intermediate_methylation import *
from open_gr_data import *
from write_new_input import *


# get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gr', '--genomic_regions', required=True )
parser.add_argument('-data', '--cpg_reads_data', required=True)
parser.add_argument('-data_from', '--data_from', required=True)
parser.add_argument('-w', '--window_size', type=int, default=1000)
parser.add_argument('-ncpgs', '--number_of_cpgs', type=int, default=5)
parser.add_argument('-mmin', '--meth_min', type=float, default=0)
parser.add_argument('-mmax', '--meth_max', type=float, default=100)
args = parser.parse_args()

command_line = " ".join(argv[1:])
cwd = os.getcwd()

if str(cwd) in args.cpg_reads_data:
    cpg_reads_data = args.cpg_reads_data
else:
    cpg_reads_data = cwd+'/'+args.cpg_reads_data

if str(cwd) in args.genomic_regions:
    genomic_regions = args.genomic_regions
else:
    genomic_regions = cwd+'/'+args.genomic_regions

os.chdir(path[0])

selected_cpgs, id_windows = open_genomic_regions(genomic_regions, cpg_reads_data, args.data_from, args.meth_min, args.meth_max, args.window_size, args.number_of_cpgs)

writeout(cpg_reads_data, id_windows, selected_cpgs, args.data_from, command_line)

os.system("python call.py "+ command_line + " -f new_cpg_reads " + " -wf windows " + "-cwd "+cwd)


