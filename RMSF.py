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
    # load references
else:
    a = 79.3	# NB!!! these are RT-dimensions!
    b = 79.3 
    c = 38.2
    # load references

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
    tmp =
    x = pd.DataFrame()

def changePDB():
    # paste changePDB-script here, file can be deleted afterwards
    pass
    
for(n in 1:nfiles):
    # do all the work
    
    
end = time.time()
print(end - start)
