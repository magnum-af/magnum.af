#https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
#https://stackoverflow.com/questions/42855285/plotting-average-curve-for-points-in-gnuplot
#https://stackoverflow.com/questions/1796597/import-an-array-in-python

import numpy as np

from numpy import *

data = loadtxt("m.dat")
t,mz = data[:,0], data [:,3]
#print t
#print mz[:]
print "len(mz)=", len(mz)
#n = 100

import matplotlib.pyplot as plt
modes = ['full', 'same', 'valid']

for n in range (1000, 11000, 1000):
    fig, ax = plt.subplots()
    
    ax.plot(mz);
    for m in modes:
        ax.plot(np.convolve(mz, np.ones((n,))/n, mode=m));
    #plt.axis([-10, len(mz), -.1, 1.1]);
    ax.legend(['mz', 'full', 'same', 'valid'], loc='lower center');
    fig.savefig("fig%d"%n)
    #plt.show()
#plt.plot(np.convolve(np.ones((200,)), np.ones((50,))/50, mode=m));
