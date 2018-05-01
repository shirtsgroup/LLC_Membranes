#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

raw_data = {'nmon': ['Sandwiched', 'Parallel\n Displaced'], 'ordered': [4.20, 4.23], 'ordered_error':[0.04, 0.04],'disordered': [4.061, 4.115], 'disordered_error':[0.065, 0.034] }

df = pd.DataFrame(raw_data, columns = ['nmon', 'ordered', 'disordered', 'ordered_error', 'disordered_error'])

pos = list(range(len(df['ordered'])))
pos = [p/2 for p in pos]
width = 0.2

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(pos, df['ordered'], width, label=df['nmon'][0], yerr=df['ordered_error'])
plt.bar([p + width for p in pos], df['disordered'], width, label=df['nmon'][1], yerr=df['disordered_error'])

ax.set_xticks([p + .5*width for p in pos])
ax.set_xticklabels(df['nmon'], fontsize=14)
ax.set_ylabel('Pore-to-pore distance', fontsize=14)
ax.set_xlabel('Head group placement', fontsize=14)
plt.legend(['3.7 $\AA$ initial spacing', '5.0 $\AA$ initial spacing'],fontsize=18, prop={'size':10})
plt.tick_params(axis='both', labelsize=14)
#ax.plot([-0.20, 2.4], [4.12, 4.12], "k--")
plt.ylim(3.5, 4.5)
ax.set_aspect(1)
#plt.xlim(-0.2, 2.4)
plt.tight_layout()
plt.savefig('p2p2.png')
plt.show()
