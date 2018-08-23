import numpy as np

def MHL(haplos, ncpgs, ID_gr, wstart, wend):
    m = np.array(haplos)  # binary array of methylation state along a set of sequences
    mhl = 0
    weight_sum = 0
    for i in range(0,ncpgs):
        allcount = 0
        count = 0
        for j in range(0,(ncpgs-i)):  # starting position of window
            for k in range(0,len(m)):
                allmeth = 1
                for q in range(j,(j+i+1)): # determine whether this row of the matrix is all methylated in the required window
                    if m[k][q] == '0':
                        allmeth = 0
                if (allmeth == 1): #and (m[k][:] == 1)
                    count = count + 1
                allcount = allcount + 1
        prop_allmeth = float(count)/float(allcount)
#        print "\nproportion fully methylated in windows of size",i,prop_allmeth
        mhl += float(float(prop_allmeth) * (i+1))
        weight_sum += (i+1)


    mhl = float(float(mhl)/float(weight_sum))
    out = open("temp4.txt", "w")
    out.write(str(ID_gr)+'\t'+str(wstart)+'\t'+str(wend)+'\t'+ str(mhl) +'\n')










