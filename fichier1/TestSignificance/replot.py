#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value2.dat')

plt.plot(M[:,0], M[:,1], 'rx', label="combine s=2")
plt.plot(M[:,0], M[:,2], 'bx', label="combine s=5")
plt.plot(M[:,0], M[:,3], 'gx', label="combine s=10")
plt.plot(M[:,0], M[:,5], 'ro', label="combine s=2 error 0.2")
plt.plot(M[:,0], M[:,6], 'bo', label="combine s=5 error 0.2")
plt.plot(M[:,0], M[:,7], 'go', label="combine s=10 error 0.2")

plt.plot(M[:,0], M[:,9], 'rs', label="combine s=2 error 0.4")
plt.plot(M[:,0], M[:,10], 'bs', label="combine s=5 error 0.4")
plt.plot(M[:,0], M[:,11], 'gs', label="combine s=10 error 0.4")

plt.plot(M[:,0], M[:,13], 'r', label="combine s=2 error 0.6")
plt.plot(M[:,0], M[:,14], 'b', label="combine s=5 error 0.6")
plt.plot(M[:,0], M[:,15], 'g', label="combine s=10 error 0.6")


plt.axis('equal')


plt.legend()
plt.title('Significance in function of the background with 0.2 error')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
