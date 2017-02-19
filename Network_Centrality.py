# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 15:41:58 2017

@author: Kanak
"""
import numpy as np

######################################################
# Test symmetric and 0 diagonal matrix
######################################################

def test_symmetric(x):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return
    
    dimsize = x.shape
    if len(dimsize) == 3:
        nday = dimsize[0]
    else:
        if (x.transpose() == x).all():
            if (np.diag(x) == 0).all():
                print('Input is a symmetric matrix and diagonal elements equal 0')
                return(True)
            print('Input is a symmetric matrix but diagonal elements not equal 0')
            return(False)
        else:
            print('Input is not a symmetric matrix')
            return(False)
    i = 0
    symm = True
    diagt = True
    while i < nday and symm:
        a = x[i,:,:]
        if (a.transpose() == a).all():
            if (np.diag(a) == 0).all() and diagt:
                diagt = True
            else:
                diagt = False
            symm = True
            i += 1
        else:
            symm = False
            i += 1
    if symm and diagt:
        print('Input is a symmetric matrix and diagonal elements equal 0')
        return(True)
    elif symm and not diagt:
        print('Input is a symmetric matrix but diagonal elements not equal 0')
        return(False)
    else:
        print('Input is not a symmetric matrix')
        return(False) 
    
######################################################
# create binary symmetric and 0 diagonal matrix 
######################################################


def create_sym(x):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return
    
    dimsize = x.shape
    if len(dimsize) == 3:
        nday = dimsize[0]
        nnode = dimsize[1]
    else:
        #ensure symmetric, binary, zero diag
        x = np.sign(x + np.transpose(x))
        np.fill_diagonal(x, 0)
        return(x)

        
    Atensor = np.zeros((nday, nnode,nnode), dtype=int)
    for k in range(nday):
        A = x[k,:,:]
        A = np.sign(A+np.transpose(A))
        np.fill_diagonal(A, 0)
        Atensor[k,:,:] = A
    return(Atensor)

######################################################
# Dynamic centrality matrix 
######################################################

def grindrod(x,
             alpha = 0.1
             ):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return
    dimention = x.shape
    if dimention[1] == dimention[2]:
        nday = dimention[0]
        nnode = dimention[1]
    else:
        print ("Input data is not Symmetric")
        return

    B = np.identity(nnode)
    I = np.identity(nnode)
    for k in reversed(range(nday)):
        Aday = x[k,:,:]
        B = np.dot(np.linalg.inv(I - alpha*Aday),B)
        B = abs(B)
        B = B/np.linalg.norm(B)
    print('A. V. Mantzaris and D. J. Higham, “A model for dynamic communicators,” European Journal of Applied Mathematics, vol. 23, no. 06, pp. 659–668, 2012.')
    return(B)

    
    
######################################################
# Rank based on dynamic centrality matrix 
######################################################

def dy_rank_degree(x, ninfnode = None):
    broadcast = np.sum(x, axis = 1)
    receive = np.sum(x, axis = 0)
    lnbroadcast = np.log(broadcast)
    lnreceive = np.log(receive)
    sortedQrankindices = (broadcast - receive).argsort()[::-1]
    sortedQrankLNindices = (lnbroadcast - lnreceive).argsort()[::-1]
    s_broadcast = broadcast[sortedQrankindices]
    s_receive = receive[sortedQrankindices]
    s_lnbroadcast = lnbroadcast[sortedQrankLNindices]
    s_lnreceive = lnreceive[sortedQrankLNindices]
    if ninfnode == None:
        ninfnode = x.shape[0]
    return{'broadcase': broadcast, 'receive': receive,
           'ln_broadcase': lnbroadcast, 'ln_receive': lnreceive,
           'sort_broadcase': s_broadcast[:ninfnode],
           'sort_receive': s_receive[:ninfnode],
           'sort_ln_broadcase': s_lnbroadcast[:ninfnode],
           'sort_ln_receive': s_lnreceive[:ninfnode],
           'sort_node_index': sortedQrankindices[:ninfnode],
           'sort_node_ln_index': sortedQrankLNindices[:ninfnode]}
           
    
######################################################
# Katz centrality matrix 
######################################################
    
def katz(x,
         alpha = None
         ):
    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except:
            print ("x is not possible to convert as array")
            return
    dimention = x.shape
    if dimention[0] == dimention[1]:
        N = dimention[0]
    else:
        print ("Input data is not Symmetric")
        return
    if alpha == None:
        w,v = np.linalg.eig(x)
        alpha = np.round_(1/max(np.abs(w)), decimals = 3) #%gives 0.0196
    print('L. Katz, “A new status index derived from sociometric analysis,” Psychometrika, vol. 18, no. 1, pp. 39–43, 1953.')
    return(np.linalg.inv(np.identity(N) - alpha*x))


######################################################
# Rank based on katz centrality matrix 
######################################################


def kz_rank_degree(x, ninfnode = None):
    centrality = np.sum(x,axis = 0)
    ln_centrality = np.log(centrality)
    sortedQrankindices_k = centrality.argsort()[::-1]
    sortedQrankLNindices_k = (ln_centrality).argsort()[::-1]
    Centrality_k = centrality[sortedQrankindices_k]
    Centrality_ln_k = ln_centrality[sortedQrankLNindices_k]
    if ninfnode == None:
        ninfnode = x.shape[0]

    return{'centrality': centrality, 'ln_centrality': ln_centrality,
           'sort_centrality': Centrality_k[:ninfnode], 
           'sort_ln_centrality': Centrality_ln_k[:ninfnode],
           'sort_node_index': sortedQrankindices_k[:ninfnode],
           'sort_node_ln_index': sortedQrankLNindices_k[:ninfnode]}
   


