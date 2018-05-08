#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

raw_data = {'nmon': ['4', '5', '6', '7', '8'], 'layered': [3.71, 4.18, 4.69, 4.66, 5.00], 'layered_err' : [0.09, 0.13, 0.11, 0.05, 0.10], 'offset': [3.84, 4.17, 4.72, 4.82, 5.33], 'offset_err': [0.07, 0.12, 0.07, 0.09, 0.05]}
layered_disordered = [4.05]
layered_disordered_err = [0.14]
offset_disordered = [4.032]
offset_disordered_err = [0.07]

#raw_data = {'nmon': ['4', '5', '6', '7', '8'], 'four': [3.72, 3.84, 0.04, 0.02], 'five': [4.20, 4.23, 0.04, 0.04], 'six': [4.69, 4.72, 0.04, 0.03], 'seven' : [ 4.66, 4.82, 0.02, 0.03], 'eight': [5.00, 5.33, 0.04, 0.02]}

df = pd.DataFrame(raw_data, columns = ['nmon', 'layered', 'offset', 'layered_err', 'offset_err'])

width = 0.2
pos = list(range(len(df['layered'])))
pos = np.array([p/2 for p in pos])
pos[2:] += 2*width

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(pos, df['layered'], width, label='Sandwiched (d=3.7)', yerr=df['layered_err'])
plt.bar([p + width for p in pos], df['offset'], width, label='Parallel Displaced (d=3.7)', yerr=df['offset_err'])
plt.bar([p + 2*width for p in pos][1], layered_disordered, width, label='Sandwiched (d=5)', yerr=layered_disordered_err)
plt.bar([p + 3*width for p in pos][1], offset_disordered, width, label='Parallel Displaced (d=5)', yerr=offset_disordered_err)

ax.set_xticks([p + .5*width for p in pos])
ax.set_xticklabels(df['nmon'], fontsize=14)
ax.set_ylabel('Pore-to-Pore Distance (nm)', fontsize=14)
ax.set_xlabel('Number of columns per pore', fontsize=14)
plt.tick_params(axis='both', labelsize=14)
ax.plot([-0.20, 2.8], [4.12, 4.12], "k--", label='Experiment')
#plt.legend(['Sandwiched', 'Parallel Displaced', 'Experiment'],fontsize=14)
plt.legend(fontsize=14)
plt.ylim(3.5, 5.5)
plt.xlim(-0.2, 2.8)
plt.tight_layout()
plt.savefig('p2p.png')
plt.show()
