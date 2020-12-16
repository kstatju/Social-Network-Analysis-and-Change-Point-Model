# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:28:20 2017

@author: Kanak
"""
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
nnode = 151
nday = 1137
new_data = np.loadtxt('D:/UCF/STA 6908 - Alexander V. Mantzaris/Data/enronAtensor.txt')
new_data.shape
new_data = new_data.reshape((nday, nnode, nnode))
new_data = new_data.astype(int)


Atensor = np.zeros((nday, nnode,nnode), dtype=int)

dayspec = np.zeros(nday, dtype=float)
for k in range(nday):
    A = new_data[:,:,k]
    #ensure symmetric, binary, zero diag
    A = np.sign(A+np.transpose(A))
    A = A - np.diag(np.diag(A))
    Atensor[k,:,:] = A
    w, v = np.linalg.eig(A)
    dayspec[k] = np.max(np.abs(w))

#####################################
# Dynamic Centrality 
#####################################


alpha = 0.1

Adeg = np.zeros(nday, dtype=float)

B = np.identity(nnode)
I = np.identity(nnode)
for k in reversed(range(nday)):
    Aday = Atensor[k,:,:]
    Adeg[k] = sum(sum(Aday))/2.0
    B = np.dot(np.linalg.inv(I - alpha*Aday),B)
    B = abs(B)
    B = B/np.linalg.norm(B)

Cbroad = np.sum(B, axis = 1)
Crec = np.sum(B, axis = 0)


#%!!!Alex
#%now do the ranking in terms of Hierarchy 
sortedQrankindices = (Cbroad - Crec).argsort()[::-1]
sortedQrankLNindices = (np.log(Cbroad) - np.log(Crec)).argsort()[::-1]
Cbroad = Cbroad[sortedQrankLNindices]
Crec = Crec[sortedQrankLNindices]
print(Cbroad[1:10])
print(Crec[1:10])
print(sortedQrankindices +1, '\n', sortedQrankLNindices + 1)
#%end Alex!!!

Degtemp = np.sum(Atensor,axis = 2)
Degall = sum(Degtemp)

specmax = max(dayspec)
amax = 1/specmax

#####################################
# %Do KATZ on binarized aggregate
#####################################

Astat = np.real(np.sum(Atensor,0)>0)
N = len(Astat)
w,v = np.linalg.eig(Astat)
amaxstar = 1/max(np.abs(w)) #%gives 0.0196
alph = 0.015
katz = np.linalg.inv(np.identity(N) - alph*Astat)
katz_sum = np.sum(katz,axis = 0)

katz_sum.shape
cbk = np.corrcoef(np.transpose(Cbroad), katz_sum)[0,1]
kenbk = stats.kendalltau(np.transpose(Cbroad),katz_sum)[0]
crk = np.corrcoef(np.transpose(Crec), katz_sum)[0,1]
kenrk = stats.kendalltau(np.transpose(Crec),katz_sum)[0]



sortedQrankindices_k = katz_sum.argsort()[::-1]
sortedQrankLNindices_k = np.log(katz_sum).argsort()[::-1]
Centrality_k = katz_sum[sortedQrankLNindices_k]
print(Centrality_k[1:10])
print(sortedQrankindices_k + 1,'\n', sortedQrankLNindices_k + 1)
# end
#####################################




#####################################
# Figures
#####################################


#figure(1)
#subplot(221)
Act1 = np.sum(Atensor, axis = 2)
Actall = np.sum(Act1, axis = 2)
plt.plot(Actall,'.')
plt.xlabel('Day')
plt.ylabel('Total Activity')
plt.title('gca')
plt.show()

#subplot(222)
plt.loglog(Cbroad, Crec, '.', basex=np.e, basey=np.e)
plt.xlim([1e-10,1e3])
plt.ylim([1e-15,1e3])
plt.xlabel('Broadcast')
plt.ylabel('Receive')
plt.title('gca')
plt.show()

#subplot(223)
plt.loglog(Cbroad, Degall, '.', basex=np.e, basey=np.e)
plt.xlim([1e-10,1e3])
plt.ylim([1e0,1e4])
plt.xlabel('Broadcast')
plt.ylabel('Total Degree')
plt.title('gca')
plt.show()

#subplot(224)

plt.loglog(Crec, Degall, '.', basex=np.e, basey=np.e)
plt.xlim([1e-10,1e3])
plt.ylim([1e0,1e4])
plt.xlabel('Receive')
plt.ylabel('Total Degree')
plt.title('gca')
plt.show()