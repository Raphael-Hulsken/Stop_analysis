#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value3.dat')

plt.plot(M[:,0], M[:,1], 'rx', label="lnN S=5")
plt.plot(M[:,0], M[:,4], 'bx', label="gmN S=5")
plt.plot(M[:,0], M[:,7], 'gx', label="lnU S=5")
plt.plot(M[:,0], M[:,2], 'ro', label="lnN S=10")
plt.plot(M[:,0], M[:,5], 'bo', label="gmN S=10")
plt.plot(M[:,0], M[:,8], 'go', label="lnU S=10")
plt.plot(M[:,0], M[:,3], 'r', label="lnN S=20")
plt.plot(M[:,0], M[:,6], 'b', label="gmN S=20")
plt.plot(M[:,0], M[:,9], 'g', label="lnU S=20")

plt.axis('equal')


plt.legend()
plt.title('Combine with different function in function of the background with 0.2 error')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
