#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

M = np.loadtxt('value5.dat')

plt.plot(M[:,0], M[:,1], 'r', label="no incertitude signal =5")
plt.plot(M[:,0], M[:,4], 'rx', label="Signal error = 1.1 et bkg error =1.1 signal=5")
plt.plot(M[:,0], M[:,7], 'ro', label="Signal error = 0.9 et bkg error =0.9 signal=5")
plt.plot(M[:,0], M[:,10], 'r^', label="Signal error = 1.1 et bkg error =0.9 signal=5")
plt.plot(M[:,0], M[:,13], 'rs', label="Signal error = 0.9 et bkg error =1.1 signal=5")
plt.plot(M[:,0], M[:,2], 'g', label="no incertitude signal=10")
plt.plot(M[:,0], M[:,5], 'gx', label="Signal error = 1.1 et bkg error =1.1 signal=10")
plt.plot(M[:,0], M[:,8], 'go', label="Signal error = 0.9 et bkg error =0.9 signal=10")
plt.plot(M[:,0], M[:,11], 'g^', label="Signal error = 1.1 et bkg error =0.9 signal=10")
plt.plot(M[:,0], M[:,14], 'gs', label="Signal error = 0.9 et bkg error =1.1 signal=10")
plt.plot(M[:,0], M[:,3], 'b', label="no incertitude  signal=15")
plt.plot(M[:,0], M[:,6], 'bx', label="Signal error = 1.1 et bkg error =1.1 signal=15")
plt.plot(M[:,0], M[:,9], 'bo', label="Signal error = 0.9 et bkg error =0.9 signal=15")
plt.plot(M[:,0], M[:,12], 'b^', label="Signal error = 1.1 et bkg error =0.9 signal=15")
plt.plot(M[:,0], M[:,15], 'bs', label="Signal error = 0.9 et bkg error =1.1 signal=15")

plt.axis('equal')


plt.legend()
plt.title('Combine with some error on the signal and background in function of the background with 0.2 incertitude')
plt.xlabel('b')
plt.ylabel('significance')
plt.xscale('log')
#plt.savefig("plot_data_file.png")

plt.show()
