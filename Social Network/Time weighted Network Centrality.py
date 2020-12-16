# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:32:15 2017

@author: Kanak
"""
import numpy as np

Acells = np.loadtxt('D:/UCF/STA 6908 - Alexander V. Mantzaris/Data/3 days network data from Alexander V. Mantzaris.txt')
node_num = 21
nday = 3
alpha = 0.1

# convert into 3D array
Acells = Acells.reshape((nday, node_num, node_num))


eyeM = np.identity(node_num, dtype=float)
S_iters = np.zeros((nday, node_num, node_num), dtype=float)
S_iters_temp = np.dot((eyeM + np.exp(-0.2)*0),(np.linalg.inv(eyeM - alpha*0))) - eyeM
S_iters_Broadcast = np.zeros((nday, node_num), dtype=float);
for ii in range(nday):
    Aday = Acells[ii,:,:];
    S_iters[ii,:,:] = np.dot((eyeM + np.exp(-0.2)*S_iters_temp),(np.linalg.inv(eyeM - alpha*Aday))) - eyeM
    S_iters_temp = S_iters[ii,:,:];
    S_iters_Broadcast[ii,:] = np.sum(S_iters[ii, :,:],axis = 1);



S_last = S_iters[ii,:,:];
S_last_Broadcast = np.sum(S_last,axis = 1);
rankIndexesS_last_Broadcast = S_last_Broadcast.argsort()[::-1]
S_last_Broadcast_sort = S_last_Broadcast[rankIndexesS_last_Broadcast]
S_last_Broadcast_rank = rankIndexesS_last_Broadcast[0:4]