from sliding_window import *
from intermediate_methylation import *


def open_genomic_regions(gen_regions, cpg_reads_file, data_from, meth_min, meth_max, window_size, n_cpgs):
    with open(gen_regions) as file:
        header = file.readline()
        selected_cpgs =[]
        id_windows = {}
        for line in file:
            fields = line.replace('\n','').split('\t')
            ID = fields[0] +':'+ fields[1]+'-'+fields[2]
            chr_gr = fields[0]
            p_start = int(fields[1])
            p_end = int(fields[2])
            cpg_meth_status = read_data_file(cpg_reads_file, p_start, p_end, chr_gr, data_from)
            cpg_meth_percentage = meth_percentage(cpg_meth_status, meth_min, meth_max)
            nframes = sliding_window(cpg_meth_percentage, n_cpgs, window_size)
            for w in nframes:
                for p in w:
                    id_windows.setdefault(ID, []).append(p)
                    for i in p:
                        selected_cpgs.append(i)
    return selected_cpgs, id_windows


def read_data_file(input_file, p_start, p_end, chr_gr, data_from):
    cpg_meth_status = {}
    with open(input_file) as file:
        header = file.readline()
        for line in file:
            fields = line.replace('\n','').split('\t')
            if data_from == 'bismark':
                position = int(fields[0]) # each position corresponds to one CpG
                methylation = str(fields[1])
                chr = fields[2]
            if data_from == 'bissnp':
                position = int(fields[1]) # each position corresponds to one CpG
                methylation = str(fields[2])
                chr = fields[0]
            # methylation state can be represented in 3 ways on the input file. The code will store the information using '1' for methylated and '0' for unmethylated CpGs representation
            if methylation in ('1','m','Z'):
                methylation = '1'
            elif methylation in ('0','u','z'):
                methylation = '0'
            else:
                raise SystemExit("\nError on input file. Methylation states have to be represented by 1, m, Z (methylated CpGs) or 0, u, z (unmethylated CpGs).\n")
            if (position >= p_start and position <= p_end) and (chr == chr_gr):
                cpg_meth_status.setdefault(position,[]).append(methylation)
    return cpg_meth_status
