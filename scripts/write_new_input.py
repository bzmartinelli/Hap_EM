

def writeout(input_file, id_windows, selected_cpgs, data_from, command_line):
    w_out = open("windows", "w")
    data_out = open("new_cpg_reads","w")
    for id, w in id_windows.iteritems():
       for p in w:
            start = p[0]
            end = p[len(p)-1]
            w_out.write(id +'\t'+ str(start) +'\t'+ str(end) +'\n')
    with open(input_file) as file:
        header = file.readline()
        for line in file:
            fields = line.replace('\n','').split('\t')
            if data_from == 'bismark' and int(fields[0]) in selected_cpgs:
                data_out.write(line)
            if data_from == 'bissnp' and int(fields[1]) in selected_cpgs:
                data_out.write(line)
    w_out.close()
    data_out.close()
