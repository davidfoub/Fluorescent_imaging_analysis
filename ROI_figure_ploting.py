# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 10:37:20 2025

@author: David
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

def plot_masks(ops, stat, iscell, clr="rx"):
    img=ops['meanImg']
    mask=np.zeros((ops['Ly'],ops['Lx'],4))
    rng=np.random.default_rng()
    for n in range(iscell.shape[0]):
        ypix = stat[n]['ypix']
        xpix = stat[n]['xpix']
        #mask[ypix,xpix]=1
        mask[ypix,xpix]=np.append(rng.random(3),0.3)
    plt.imshow(img, vmax=np.max(img)*0.1,interpolation='none', cmap="gray")
    plt.imshow(mask)
    #plt.plot(np.array(med)[:,0],np.array(med)[:,1],"bx",ms=3.0)
    #plt.text(x,y, idx,c='1', size =4.0)
    
    return mask

ops=np.load("E:/glia projects/plasticity/data/training_A1_01apr23/A1_min15_01apr23/suite2p/plane0/ops.npy",allow_pickle=True).item()
stat=np.load("E:/glia projects/plasticity/data/training_A1_01apr23/A1_min15_01apr23/suite2p/plane0/stat.npy",allow_pickle=True)
iscell=np.load("E:/glia projects/plasticity/data/training_A1_01apr23/A1_min15_01apr23/suite2p/plane0/iscell.npy",allow_pickle=True)
plot_masks(ops, stat, iscell)