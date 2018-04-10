

def sliding_window(cpg_meth_percentage, n_cpgs, window_size): #window_size=args.n_cpgs given by the user
    nframes = []
    line_number = sorted(cpg_meth_percentage.keys())
    frames = [ line_number[start:start + n_cpgs] for start in range(0, len(line_number)) ]
    for w in frames:
        if len(w) == n_cpgs and (int(w[n_cpgs-1]) - int(w[0]) < window_size):
            nframes.append([w])  #choose the windows
    return nframes


