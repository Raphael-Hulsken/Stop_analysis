#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value.dat')

plt.plot(M[:,0], M[:,1], 'rx', label="Zbi s=2")
plt.plot(M[:,0], M[:,2], 'r', label="combine s=2")
plt.plot(M[:,0], M[:,3], 'gx', label="Zbi s=5")
plt.plot(M[:,0], M[:,4], 'g', label="combine s=5")
plt.plot(M[:,0], M[:,5], 'bx', label="Zbi s=10")
plt.plot(M[:,0], M[:,6], 'b', label="combine s=10")


plt.axis('equal')


plt.legend()
plt.title('Significance in function of the background with 0.2 error')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
