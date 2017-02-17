# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 21:51:15 2017

@author: Kanak
"""
import numpy as np
import scipy.io
mat = scipy.io.loadmat('enronAtensor.mat')
data = mat.get('enronAtensor')
with open('enronAtensor.txt', 'w') as outfile:

    outfile.write('# Array shape: {0}\n'.format(data.shape))


    for data_slice in data:

        for row in data_slice:
            j = 0
            lengthr = len(row)
            for i in row:
                j +=1
                if j < lengthr:
                    outfile.write('%s,'%i)
                else:
                    outfile.write('%s\n'%i)
            #outfile.write('\n')

        # Writing out a break to indicate different slices...
        outfile.write('# New slice\n')

        
# Read the array from disk
new_data = np.loadtxt('enronAtensor.txt')

print(new_data.shape)

# convert into 3D array
new_data = new_data.reshape((151,151,1137))

# compare two data set
assert np.all(new_data == data)
