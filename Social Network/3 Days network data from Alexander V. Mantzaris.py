# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 20:30:45 2017

@author: Kanak
"""
import numpy as np
data = np.loadtxt('D:/UCF/STA 6908 - Alexander V. Mantzaris/Data/3 days network data from Alexander V. Mantzaris.txt')

print(data.shape)

# convert into 3D array
data = data.reshape((3,21,21))

print(data[2][9])
