#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value6.dat')

plt.plot(M[:,0], M[:,1], 'r', label="signal =1")
plt.plot(M[:,0], M[:,2], 'b', label="signal =10")
plt.plot(M[:,0], M[:,3], 'g', label="signal =30")
plt.plot(M[:,0], M[:,4], 'b', label="signal =50")
plt.plot(M[:,0], M[:,5], 'g', label="signal =100")

plt.axis('equal')


plt.legend()
plt.title('Combine with some error on the signal and background in function of the background with 0.2 incertitude')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
