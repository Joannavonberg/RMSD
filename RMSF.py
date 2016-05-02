#!/usr/bin/env python

import time
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

start = time.time()

if sys.argv[1] == "mainchain":
    mainchain = True

nfiles = int(sys.argv[2])

if sys.argv[3] == "cryo":
    cryo = True
    
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
x_trans = []
y_trans = []
z_trans = []

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

def ChangeSign(x, col):
    x.loc[:,col] = x.loc[:,col].apply(lambda x: -x)
    
for n in range(1, nfiles):
    x = load("x", n)
    y = load("y", n)
    z = load("z", n)

    """

construct a boolean dataframe for each symmetry operator and assign set value to True for monomers that need that specific operator
    - sg : sign change, x = -x
    - mh : minus half, x = x - 0.5*ucp
    - ph : plus half, x = x + 0.5*ucp
    - mq : minus quarter, x = x - 0.25*ucp
    - pq : plus quarter, x = x + 0.25*ucp
    - pt : plus three quarter, x = x + 0.75*ucp
    """

    """init = np.repeat(False, 4)
    signchangers = {}
    for letter in range(ord('A'), ord('I')):
        signchangers[chr(letter)] = init
    signchangers = pd.DataFrame(signchangers)"""

    # sign change
    sg = {}
    sg['A'] = np.repeat(False, 3)
    sg['B'] = [True, True, False]
    sg['C'] = [False, True, False]
    sg['D'] = [True, False, False]
    sg['E'] = [True, False, True]
    sg['F'] = [False, True, True]
    sg['G'] = [False, False, True]
    sg['H'] = [True, True, True]
    sg = pd.DataFrame(sg)
    #sg.index = ['x', 'y', 'z']

    # - 0.5 * ucp
    mh = {}
    mh['A'] = np.repeat(False, 3)
    mh['B'] = [False, False, True]
    mh['C'] = [True, False, False]
    mh['D'] = [False, True, False]
    mh['E'] = [False, True, False]
    mh['F'] = [True, False, False]
    mh['G'] = mh['H'] = np.repeat(False, 3)
    mh = pd.DataFrame(mh)

    # + 0.5 * ucp
    ph = {}
    ph['A'] = ph['B'] = ph['G'] = ph['H'] = np.repeat(False, 3)
    ph['C'] = [False, True, False]
    ph['D'] = [True, False, False]
    ph['E'] = [True, False, False]
    ph['F'] = [False, True, False]
    ph = pd.DataFrame(ph)

    # minus quarter
    mq = {}
    mq['A'] = mq['B'] = mq['C'] = mq['E'] = mq['F'] = mq['G'] = mq['H'] = np.repeat(False, 3)
    mq['D'] = [False, False, True]
    mq = pd.DataFrame(mq)

    # plus quarter
    pq = {}
    pq['A'] = pq['B'] = pq['D'] = pq['E'] = pq['G'] = pq['H'] = np.repeat(False, 3)
    pq['C'] = [False, False, True]
    pq['F'] = [False, False, True]
    pq = pd.DataFrame(pq)

    # plus three quarter
    pt = {}
    pt['A'] = pt['B'] = pt['C'] = pt['D'] = pt['F'] = pt['G'] = pt['H'] = np.repeat(False, 3)
    pt['E'] = [False, False, True]
    pt = pd.DataFrame(pt)

    # applying the relevant symmetry operations
    n = -1
    ucp = [a, b, c]
    for coord in x, y, z:
        n += 1
        #co += 1
        for letter in range(ord('A'), ord('I')):
            if sg.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: -x)
            if mh.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: x - 0.5*ucp[n])
            if ph.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: x + 0.5*ucp[n])
            if mq.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: x - 0.25*ucp[n])
            if pq.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: x + 0.25*ucp[n])
            if pt.loc[n, chr(letter)]:
                coord.loc[:,chr(letter)] = coord.loc[:,chr(letter)].apply(lambda x: x + 0.75*ucp[n])
    
    if n==1:
        for letter in range(ord('A'), ord('I')):
            x_premodulo = x.loc[505,chr(letter)]
            y_premodulo = y.loc[505,chr(letter)]
            z_premodulo = z.loc[505,chr(letter)]

            x_trans.append(x_premodulo.apply(lambda x: x%a - x))
            y_trans.append(y_premodulo.apply(lambda x: x%a - x))
            z_trans.append(z_premodulo.apply(lambda x: x%a - x))
    i = -1
    for letter in range(ord('A'), ord('I')):        
        x.loc[:,chr(letter)] = x.loc[:,chr(letter)].apply(lambda x: x +x_trans[i])
        y.loc[:,chr(letter)] = y.loc[:,chr(letter)].apply(lambda x: x +x_trans[i])
        z.loc[:,chr(letter)] = z.loc[:,chr(letter)].apply(lambda x: x +x_trans[i])
        i+=1

    rmsd.append(np.sqrt(np.sqrt(RMSD(x, ref = refx, ucp = a)^2 + RMSD(y, ref = refy, ucp = b)^2)+ RMSD(z2, ref = refz, ucp = c)^2))
    x2 = []
    y2 = []
    z2 = []
    for t in range(1,8):
        x2.append(x.loc[:,t])
        y2.append(y.loc[:,t])
        z2.append(z.loc[:,t])

    writePDB(x2, y2, z2, n, mainchain)
        
end = time.time()
print("This script took " + str(end - start) + " seconds to finish.")
