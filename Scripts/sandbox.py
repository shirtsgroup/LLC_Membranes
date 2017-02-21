#!/home/ben/anaconda2/bin/python2
import time
import numpy as np

times = []

for j in range(100):
    start = time.time()
    count = 0
    for i in range(1000000):
        count += 1
    end = time.time()
    times.append(end - start)

print "Average: %s, Standard Deviation: %s" %(np.mean(times), np.std(times))