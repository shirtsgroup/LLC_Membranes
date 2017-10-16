#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

raw_data = {'nmon': ['4', '5', '6', '7', '8'], 'layered': [3.72, 4.20, 4.69, 4.66, 5.00], 'layered_err' : [0.04, 0.04, 0.04, 0.02, 0.04], 'offset': [3.84, 4.23, 4.72, 4.82, 5.33], 'offset_err': [0.02, 0.04, 0.03, 0.03, 0.02]}

#raw_data = {'nmon': ['4', '5', '6', '7', '8'], 'four': [3.72, 3.84, 0.04, 0.02], 'five': [4.20, 4.23, 0.04, 0.04], 'six': [4.69, 4.72, 0.04, 0.03], 'seven' : [ 4.66, 4.82, 0.02, 0.03], 'eight': [5.00, 5.33, 0.04, 0.02]}

df = pd.DataFrame(raw_data, columns = ['nmon', 'layered', 'offset', 'layered_err', 'offset_err'])

pos = list(range(len(df['layered'])))
pos = [p/2 for p in pos]
width = 0.2

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(pos, df['layered'], width, label=df['nmon'][0], yerr=df['layered_err'])
plt.bar([p + width for p in pos], df['offset'], width, label=df['nmon'][1], yerr=df['offset_err'])

ax.set_xticks([p + .5*width for p in pos])
ax.set_xticklabels(df['nmon'], fontsize=14)
ax.set_ylabel('Pore-to-Pore Distance', fontsize=14)
ax.set_xlabel('Number of monomers per layer', fontsize=14)
plt.legend(['Sandwiched', 'Parallel Displaced'],fontsize=14)
plt.tick_params(axis='both', labelsize=14)
plt.ylim(3.5, 5.5)
plt.tight_layout()
plt.savefig('p2p.png')
plt.show()
