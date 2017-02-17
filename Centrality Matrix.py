# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:25:39 2017

@author: Kanak
"""

import numpy as np
import sys
sys.path.insert(0, 'Python Code/')
import Network_Centrality as cent

#
new_data = np.loadtxt('enronAtensor.txt', delimiter = ',')
nnode = 151
nday = 1137

# convert into 3D array
new_data = new_data.reshape((nnode, nnode, nday))

Atensor = np.zeros((nday, nnode,nnode), dtype=int)

cent.test_symmetric(Atensor)
xx = cent.create_sym(Atensor)

B = cent.grindrod(Atensor)    
rankind = cent.dy_rank_degree(B)
print(rankind['sort_node_index'])
print(rankind['sort_node_ln_index'])


kz = cent.katz(np.sum(Atensor, 0), .015)
rkz = cent.kz_rank_degree(kz)
rkz


