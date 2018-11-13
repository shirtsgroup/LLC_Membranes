#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#raw_data = {'nmon': ['Sandwiched', 'Parallel\n Displaced'], 'ordered': [4.38, 4.39], 'disordered': [5.11, 4.85] }
raw_data = {'nmon': ['Sandwiched', 'Parallel\n Displaced'], 'ordered': [4.38, 4.39], 'disordered': [5.16, 4.84] }

df = pd.DataFrame(raw_data, columns = ['nmon', 'ordered', 'disordered'])

pos = list(range(len(df['ordered'])))
pos = [p/2 for p in pos]
width = 0.2

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(pos, df['ordered'], width, label=df['nmon'][0])
plt.bar([p + width for p in pos], df['disordered'], width, label=df['nmon'][1])

ax.set_xticks([p + .5*width for p in pos])
ax.set_xticklabels(df['nmon'], fontsize=14)
ax.set_ylabel('Distance between layers', fontsize=14)
ax.set_xlabel('Head group placement', fontsize=14)
plt.legend(['3.7 $\AA$ initial spacing', '5.0 $\AA$ initial spacing'],fontsize=18, prop={'size':10})
plt.tick_params(axis='both', labelsize=14)
#ax.plot([-0.20, 2.4], [4.12, 4.12], "k--")
plt.ylim(3.5, 5.5)
#plt.xlim(-0.2, 2.4)
ax.set_aspect('1')
plt.tight_layout()
plt.savefig('dbwl.png')
plt.show()
