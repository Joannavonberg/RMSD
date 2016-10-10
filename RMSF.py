#!/usr/bin/env python

import time
import sys
import pandas as pd
import numpy as np

start = time.time()

if sys.argv[1] == "mainchain":
    mainchain = True
else:
    mainchain = False

nfiles = int(sys.argv[2])

if sys.argv[3] == "cryo":
    cryo = True
else:
    cryo = False
    
# unit cell dimensions
if(cryo):
    a = 78.64
    b = 78.64
    c = 37.06
    temp = "cryo"
else:
    a = 79.3
    b = 79.3 
    c = 38.2
    temp = "RT"
    
# load crystal structure references
if not mainchain:
    refx = pd.read_csv("/work/berg/Git/ref/x_%s_protein_noH.txt" %temp, header=None, squeeze=True)
    refy = pd.read_csv("/work/berg/Git/ref/y_%s_protein_noH.txt" %temp, header=None, squeeze=True)
    refz = pd.read_csv("/work/berg/Git/ref/z_%s_protein_noH.txt" %temp, header=None, squeeze=True)

rmsd = []

def RMSD(data, ref, ucp):
    tot = 0
    tot = data.apply(MeanDifference, axis=0, ref = ref, ucp = ucp)
    tot = np.sqrt(tot/data.loc[1,:].size)
    return tot

def MeanDifference(col, ref, ucp):
    tot = 0
    for n in range(0, col.size-1):
        tot = tot + min((ucp - max(ref[n], col[n]) + min(ref[n], col[n])), abs(ref[n] - col[n])) ** 2
    tot = tot/col.size
    return tot

def load(a, n):
    tmp = pd.read_csv("%s%.0f.txt" %(a, n), header=None, squeeze=True)
    chains = {}
    begin = 0
    step = tmp.size/8
    end = step
    for letter in range(ord('A'), ord('I')):
        chains[chr(letter)] = tmp.iloc[begin:end].values
        begin += step
        end += step
    x = pd.DataFrame(chains)
    return x

def writePDB(x, y, z, filenumber, mainchain):
    #    For writing mainchain PDB files
    if mainchain:
        f = open("/work2/berg/Simulations/Unit_Cells/CorrectBox/300K_NVT_cryodim/PCA//chain_A/average.pdb", 'r')
    #For writing protein-H files
    else:
        f = open("/work2/berg/Simulations/Unit_Cells/CorrectBox/300K_NVT_cryodim/Simulation/eerste_noH.pdb")
    fout = open("%.0f.pdb" % filenumber, 'w')
    n = 0
    for line in f:
        if line[0:6] == "ATOM  ":
            newline = line[0:30] + str("%8.3f" % x[n]) +  str("%8.3f" % y[n]) + str("%8.3f" % z[n]) + line[55:-1] + "\n"
            fout.write(newline)
            n += 1
        else:
            fout.write(line)
    f.close()
    fout.close()

for n in range(1, nfiles+1):
    print("\n n is " + str(n) + "\n")
    x = load("x", n)
    y = load("y", n)
    z = load("z", n)

    x2 = x
    y2 = y
    z2 = z
""" 
The chains that have an exchange of axes seem to go wrong...
"""
    x2.loc[:,'B'] = x.loc[:,'B'] * -1
    y2.loc[:,'B'] = y.loc[:,'B'] * -1
    z2.loc[:,'B'] = z.loc[:,'B'] - 0.5 * c

    x2.loc[:,'C'] = y.loc[:,'C'] - 0.5 * b #
    y2.loc[:,'C'] = x.loc[:,'C'] + 0.5 * a
    z2.loc[:,'C'] = z.loc[:,'C'] + 0.25 * c

    x2.loc[:,'D'] = y.loc[:,'D'] * -1 + 0.5 * b #
    y2.loc[:,'D'] = x.loc[:,'D'] - 0.5 * a
    z2.loc[:,'D'] = z.loc[:,'D'] - 0.25 * c

    x2.loc[:,'E'] = x.loc[:,'E'] * -1 + 0.5 * a 
    y2.loc[:,'E'] = y.loc[:,'E'] - 0.5 * b
    z2.loc[:,'E'] = z.loc[:,'E'] * -1 + 0.75 * c

    x2.loc[:,'F'] = x.loc[:,'F'] - 0.5 * a
    y2.loc[:,'F'] = y.loc[:,'F'] * -1 + 0.5 * b
    z2.loc[:,'F'] = z.loc[:,'F'] * -1 + 0.25 * c

    x2.loc[:,'G'] = y.loc[:,'G'] #
    y2.loc[:,'G'] = x.loc[:,'G']
    z2.loc[:,'G'] = z.loc[:,'G'] * -1

    x2.loc[:,'H'] = y.loc[:,'H'] * -1 #
    y2.loc[:,'H'] = x.loc[:,'H'] * -1
    z2.loc[:,'H'] = z.loc[:,'H'] * -1 + 0.5 * c

    if n is 1:
        x_trans = []
        y_trans = []
        z_trans = []
        for letter in range(ord('A'), ord('I')):
            x_premodulo = x2.loc[505,chr(letter)]
            y_premodulo = y2.loc[505,chr(letter)]
            z_premodulo = z2.loc[505,chr(letter)]
            x_trans.append(x_premodulo%a - x_premodulo)
            y_trans.append(y_premodulo%b - y_premodulo)
            z_trans.append(z_premodulo%c - z_premodulo)
    i = 0
    for letter in range(ord('A'), ord('I')):
        x2.loc[:,chr(letter)] = x2.loc[:,chr(letter)] + x_trans[i]
        y2.loc[:,chr(letter)] = y2.loc[:,chr(letter)] + y_trans[i]
        z2.loc[:,chr(letter)] = z2.loc[:,chr(letter)] + z_trans[i]
        i+=1
        
    # for writing to PDB file    
    x3 = []
    y3 = []
    z3 = []
    for letter in range(ord('A'), ord('I')):
        for t in range(x2.loc[:,chr(letter)].size):
            x3.append(x2.loc[t,chr(letter)])
            y3.append(y2.loc[t,chr(letter)])
            z3.append(z2.loc[t,chr(letter)])

    writePDB(x3, y3, z3, n, mainchain) 
        
end = time.time()
print("This script took " + str(end - start) + " seconds to finish.")
