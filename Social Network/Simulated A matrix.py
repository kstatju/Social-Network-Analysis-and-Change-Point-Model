# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 16:28:08 2017

@author: ka746940
"""

import numpy as np

size_n = 5 # sample size of each timestamp
high_n = 3 # number of node
time_n = 4 # number of timestamp 

# Generating data contains timestamp, sender number, receiver number

dt = []
for j in range(time_n):
    time = np.repeat(j, size_n)
    sender = np.random.randint(low = 0, high = high_n, size = size_n)
    receiver = np.random.randint(low = 0, high = high_n, size = size_n)

    for i in range(size_n):
        if sender[i] != receiver[i]:
            dt.append([time[i], sender[i], receiver[i]])
data = np.array(dt)


# generating #node x #node zero matrix for A for each timestamp (time_n)
A = np.zeros((time_n, high_n, high_n))

# if sender sends a message to receiver, assign 1 otherwise 0 for each timestamp

for row in data:
    A[row[0], row[1], row[2]] = 1
    
print(data)
print(A)