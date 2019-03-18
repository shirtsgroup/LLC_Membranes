#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

opacity = 0.8

custom = [Patch(facecolor='xkcd:blue', alpha=opacity),
          Patch(facecolor='xkcd:orange', alpha=opacity),
          Patch(facecolor='xkcd:green', alpha=opacity),
          Patch(facecolor='xkcd:red', alpha=opacity)]

fig, ax = plt.subplots(figsize=(6.5, 1.2))

ax.legend(custom, ['n = 1', 'n = 1', 'n = 3', 'n = 4'], ncol=4, title='Number of Simultaneous Hydrogen Bonds', fontsize=14, title_fontsize=14)
plt.gca().set_axis_off()
plt.tight_layout()
plt.savefig('hbond_sensitivity_legend.pdf')
plt.show()
