#!/usr/bin/env python

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(20, 20))
colors = ['xkcd:blue', 'xkcd:olive', 'xkcd:orangered', 'xkcd:magenta', 'xkcd:gold']
f = lambda m,c: ax.plot([],[], color=c, linewidth=4)[0]
handles = [f("s", colors[i]) for i in range(5)]
print(handles)
labels = ["Parallel Displaced ($d$ = 3.7 $\AA$)", "Sandwiched ($d$ = 3.7 $\AA$)", "Parallel Displaced ($d$ = 5 $\AA$)", "Sandwiched ($d$ = 5 $\AA$)", "Solvated Parallel Displaced"]
legend = plt.legend(handles, labels, loc=3, framealpha=1, frameon=False, ncol=3, fontsize=14)

def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

export_legend(legend)
plt.show()
