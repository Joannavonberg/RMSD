#!/usr/bin/env python

import time
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

start = time.time()

# unit cell dimensions
if(cryo):
    a = 78.64	# NB!!! these are cryo-dimensions!
    b = 78.64
    c = 37.06
    temp = "cryo"
else:
    a = 79.3	# NB!!! these are RT-dimensions!
    b = 79.3 
    c = 38.2
    temp = "RT"
# load references
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
        tot = tot + min((ucp - max(ref[n], col[n]) + min(ref[n], col[n])), abs(ref[n] - col[n]))^2
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

def changePDB(x, y, z, filenumber, mainchain):
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

def ChangeSign(x, col):
    x.loc[:,col] = x.loc[:,col].apply(lambda x: -x)
    
for(n in 1:nfiles):
    x = load("x", n)
    y = load("y", n)
    z = load("z", n)

    """

construct a boolean dataframe to assign certain columns for sign change 
    - initialize everything to False
    - manually change the wanted columns to True
    """
    """init = np.repeat(False, 4)
    signchangers = {}
    for letter in range(ord('A'), ord('I')):
        signchangers[chr(letter)] = init
    signchangers = pd.DataFrame(signchangers)"""

    cols = {}
    cols['A'] = np.repeat(False, 3)
    cols['B'] = [True, True, False]
    cols['C'] = [False, True, False]
    cols['D'] = [True, False, False]
    cols['E'] = [True, False, True]
    cols['F'] = [False, True, True]
    cols['G'] = [False, False, True]
    cols['H'] = np.repeat(True, 3)
    cols = pd.DataFrame(cols)
    cols.index = ['x', 'y', 'z']
    
    co = ord('w')
    for coord in x, y, z:
        co += 1
        for letter in range(ord('A'), ord('I')):
            if cols.loc[chr(co), chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: -x)

    z.loc[:,'B'] = z.loc[:,'B'].apply(lambda x: x - 0.5*c)

    

end = time.time()
print(end - start)
