# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 14:03:59 2026

@author: David
"""

from tifffile import imread, imwrite
import numpy as np
from matplotlib.pyplot import imshow
import glob
import os

rootdir = 'E:/glia projects/plasticity/data/'

files = [f.path for i in glob.glob(f'{rootdir}/**/',recursive=True) for f in os.scandir(i) if f.path.endswith('.tif') and 'min' in  f.path and 'APV' in f.path]

for file in files:
    img= imread(file)
    print('read ', file)
    if img.shape[0] == 9000:
        imwrite(file[:file.rfind('.')]+'_mean'+file[file.rfind('.'):], np.mean([img[0:9000:2],img[1:9000:2]], axis=0).astype('uint16'),dtype='uint16')

    print("wrote ", file)