#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value4.dat')

plt.plot(M[:,0], M[:,1], 'rx', label="S=10")
plt.plot(M[:,0], M[:,3], 'bx', label="S=30")
plt.plot(M[:,0], M[:,5], 'gx', label="S=50")

plt.plot(M[:,0], M[:,2], 'ro', label="Zbi S=10")
plt.plot(M[:,0], M[:,4], 'bo', label="Zbi S=30")
plt.plot(M[:,0], M[:,6], 'go', label="Zbi S=50")


plt.axis('equal')


plt.legend()
plt.title('Combine with different Signal in function of the background with 0.2 error')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
