


import matplotlib.pyplot as plt

nWorkers = [ 1,  2,  3,  4,  5,  6,  8,  10,  14, 20]

realTime = [ 178, 80, 62, 50, 42, 36, 26, 29, 24, 25 ]

yieldRaw = [ 89.3056, 90.2440, 90.2091, 90.1743, 90.1394, 90.1045, 90.0347, 89.9650, 89.9281, 89.9420 ]

estimatedTime = [realTime[0] * 1 / n for n in nWorkers]

trueYield = yieldRaw[0]
yieldError = [ (y - trueYield)/trueYield * 100 for y in yieldRaw]

plt.subplot(2, 1, 1)
plt.plot(nWorkers, realTime, 'bo', label="Measured time")
plt.plot(nWorkers, estimatedTime, '-', label="Time(1 worker) / nWorkers")
plt.ylabel('Running time')

legend = plt.legend(loc='upper right', shadow=False)

plt.subplot(2, 1, 2)
plt.plot(nWorkers, yieldError, 'r.-', label="Relative error on yield (in %)")
plt.xlabel('Number of workers')
plt.ylabel('Rel. error on yield')

legend2 = plt.legend(loc='upper right', shadow=False)

plt.show()
