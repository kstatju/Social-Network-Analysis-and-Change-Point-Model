# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:25:39 2017

@author: Kanak
"""

import numpy as np
import sys
sys.path.insert(0, 'D:/UCF/STA 6908 - Alexander V. Mantzaris/Python Code/')
import Network_Centrality as cent

#
new_data = np.loadtxt('D:/UCF/STA 6908 - Alexander V. Mantzaris/Data/enronAtensor.txt', delimiter = ',')
nnode = 151
nday = 1137

# convert into 3D array
new_data = new_data.reshape((nnode, nnode, nday))

Atensor = np.zeros((nday, nnode,nnode), dtype=int)

dayspec = np.zeros(nday, dtype=float)
for k in range(nday):
    A = new_data[:,:, k]
    #ensure symmetric, binary, zero diag
    A = np.sign(A+np.transpose(A))
    A = A - np.diag(np.diag(A))
    Atensor[k,:,:] = A
    w, v = np.linalg.eig(A)
    dayspec[k] = np.max(np.abs(w))
del(w,v,k,A)
    
#################################################################

new_data = np.loadtxt('D:/UCF/STA 6908 - Alexander V. Mantzaris/Data/3 days network data from Alexander V. Mantzaris.txt')
nnode = 21
nday = 3
new_data.shape

# convert into 3D array
new_data = new_data.reshape((nday, nnode, nnode))

Atensor = np.zeros((nday, nnode,nnode), dtype=int)

dayspec = np.zeros(nday, dtype=float)
for k in range(nday):
    A = new_data[k, :,:]
    #ensure symmetric, binary, zero diag
    A = np.sign(A+np.transpose(A))
    A = A - np.diag(np.diag(A))
    Atensor[k,:,:] = A
    w, v = np.linalg.eig(A)
    dayspec[k] = np.max(np.abs(w))
del(w,v,k,A)    
    
cent.test_symmetric(Atensor)
xx = cent.create_sym(Atensor)

B = cent.grindrod(Atensor)    
rankind = cent.dy_rank_degree(B)
print(rankind['sort_node_index'])
print(rankind['sort_node_ln_index'])


kz = cent.katz(np.sum(Atensor, 0), .015)
rkz = cent.kz_rank_degree(kz)
rkz

    
