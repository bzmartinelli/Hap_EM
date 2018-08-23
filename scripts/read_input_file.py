
# Reading the input file:



def read_input_file(input_file, start_coordinate, end_window, data_from):
    all_positions = {} # all CpGs as keys
    all_reads = {} # all reads as keys
    meth_status = {}  # CpGs as keys and the values consist of another dictionary where the reads are keys and the values are the methylation status of the CpG. > dict {position: {read: meth}, }

    with open(input_file) as file:
        header = file.readline()
        for line in file:
            fields = line.replace('\n','').split('\t')

            if data_from == 'bismark':
                position = int(fields[3]) # each position corresponds to one CpG
                methylation = str(fields[4])
                read = fields[0]
                chr = fields[2]
            if data_from == 'bissnp':
                position = int(fields[1]) # each position corresponds to one CpG
                methylation = str(fields[2])
                read = fields[5]
                chr = fields[0]

            # methylation state can be represented in 3 ways on the input file. The code will store the information using '1' for methylated and '0' for unmethylated CpGs representation

            if methylation in ('1','m','Z'):
                methylation = '1'
            elif methylation in ('0','u','z'):
                methylation = '0'
            else:
                raise SystemExit("\nError on input file. Methylation states have to be represented by 1, m, Z ( for methylated CpGs) or 0, u, z (for unmethylated CpGs).\n")
            if (position >= start_coordinate and position <= end_window):
                # if the CpG is in the genomic locus of interest, those 3 dictionaries initialized above are created.
                all_positions[position] = 1
                all_reads[read] = 1
                meth_status.setdefault(position, {}).setdefault(read, []).append(methylation)
    return all_reads, all_positions, meth_status




