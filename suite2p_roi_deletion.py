# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 11:05:53 2023

@author: BioCraze
"""

import numpy as np

stat = np.load("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/stat.npy", allow_pickle=True)
F = np.load("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/F.npy", allow_pickle=True)
iscell = np.load("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/iscell.npy", allow_pickle=True)
spks = np.load("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/spks.npy", allow_pickle=True)
fneu = np.load("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/Fneu.npy", allow_pickle=True)

fneu0 = np.delete(fneu, [47,13,5,41,120,130,84,37,30,103,88,74,36,59,45,50,43,6,11,75,132,96],0)
F0 = np.delete(F, [47,13,5,41,120,130,84,37,30,103,88,74,36,59,45,50,43,6,11,75,132,96],0)
stat0 = np.delete(stat, [47,13,5,41,120,130,84,37,30,103,88,74,36,59,45,50,43,6,11,75,132,96],0)
iscell0 = np.delete(iscell, [47,13,5,41,120,130,84,37,30,103,88,74,36,59,45,50,43,6,11,75,132,96],0)
spks0 = np.delete(spks, [47,13,5,41,120,130,84,37,30,103,88,74,36,59,45,50,43,6,11,75,132,96],0)

np.save("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/fneu.npy", fneu0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/F.npy", F0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/iscell.npy", iscell0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/stat.npy", stat0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/APV_training_A1_20aug25/A1_APV_train_min60_20aug25/suite2p/plane0/spks.npy", spks0,allow_pickle=True)
