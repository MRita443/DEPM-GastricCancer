import bct
import numpy as np
import sys
import pandas as pd

if __name__ == "__main__":
    inFile = sys.argv[1]

#print(inFile)

inp = pd.read_csv("input-matrix.csv", sep= ";", index_col = 0, decimal = ",")

mat = inp.to_numpy(copy=True)

res = bct.modularity_louvain_und_sign(mat, seed=13)

#print(res[0])

with open('output.txt', 'w') as f:
    for line in res[0]:
        f.write("%s\n" % line)
