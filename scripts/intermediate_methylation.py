
# calculate the percentage of methylation of individual CpGs and select only the ones with intermediate methylation
def meth_percentage(cpg_meth_status, meth_min, meth_max):
    cpg_meth_percentage = {}   # {position: % meth}
    for k in cpg_meth_status.keys():
        meth = float(cpg_meth_status[k].count('1'))*100/len(cpg_meth_status[k])
        if meth > meth_min and meth < meth_max:
            cpg_meth_percentage[k] = meth
    return cpg_meth_percentage




