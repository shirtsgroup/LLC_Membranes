#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd

# method comparison
datas =	[{'label':'Aligned\n Experiment', 'color':'b', 'height':1.3*10**-5, 'error':0.1*10**-5}, 
        {'label':'Nernst\n Einstein', 'color':'g', 'height':1.33*10**-4, 'error':0.04444*10**-4},
        {'label':'Collective\n Diffusion','color':'r', 'height': 1.37*10**-4, 'error':0.45*10**-4}]

# crosslink comparison
#datas =	[{'label':'Aligned\n Experiment', 'color':'b', 'height':1.3*10**-5, 'error':0.1*10**-5}, 
#        {'label':'Simulated \n Crosslinked', 'color':'g', 'height':0.89*10**-4, 'error':0.03*10**-4},
#        {'label':'Simulated \n Uncrosslinked','color':'r', 'height':1.33*10**-4, 'error':0.0444*10**-4}]

width=0.9

i = 0
for data in datas:
	plt.bar(i, data['height'], width, align='center', alpha=0.5, color=data['color'], yerr=data['error'])
	i += 1

labels = [data['label'] for data in datas]
pos = [i for i in range(len(datas))]
plt.tick_params(axis='both', labelsize=14)
plt.xticks(pos, labels)
plt.ylabel('Ionic Conductivity (S/m)', fontsize=14)
plt.yscale('log')
plt.ylim(10**-7, 10**-3)
plt.tight_layout()
plt.savefig('IC.png')
plt.show()
