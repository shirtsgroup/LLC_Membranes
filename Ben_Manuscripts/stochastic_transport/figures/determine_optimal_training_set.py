#!/usr/bin/env python

"""
Brute force find the best way to partition the trajectories
"""

import numpy as np

msds = np.array([2.5597, 0.0868, 1.3591, 0.0318, 4.7476, 8.4055, 0.0374, 1.5502, 0.086, 0.341, 0.0198, 0.1932, 4.759, 39.1757, 0.0521, 4.979, 1.0219, 9.886, 1.8209, 0.0786, 1.8244, 0.0869, 0.3532, 8.8082])

target = np.mean(msds)
ntrials = 1000
means = np.zeros(ntrials)
indices = np.zeros([ntrials, 12])
variances = np.zeros(ntrials)  # difference in variance between two halves

for t in range(ntrials):
	ndx = np.random.choice(24, size=12, replace=False)
	not_ndx = np.array([x for x in range(24) if x not in ndx])
	means[t] = np.mean(msds[ndx])
	indices[t, :] = ndx
	variances[t] = np.abs(msds[ndx].var() - msds[not_ndx].var()) 

winner = np.argmin(variances)
#print(means[winner], variances[winner])
#exit()


#winner = np.argmin(np.abs(means - target))
first_set = indices[winner, :].astype(int)
second_set = np.array([x for x in range(24) if x not in indices[winner, :]])
print('First Set:\nmean=%.3f, var=%.2f, indices=%s' % (means[winner], msds[first_set].var(), first_set))
print('Second Set:\nmean=%.3f, var=%.2f, indices=' % (msds[second_set].mean(), msds[second_set].var()), second_set)
