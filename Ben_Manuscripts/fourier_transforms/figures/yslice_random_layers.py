#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

#perfect = np.load('perfect.npz')
#random_layers = np.load('random_layers.npz')
perfect = np.load('perfect_100pores.npz')
random_layers = np.load('random_layers_100pores.npz')

plt.plot(perfect['freq_y'], perfect['slice'], label='Perfect Crystal', linewidth=2)
plt.plot(random_layers['freq_y'], random_layers['slice'], label='Random Layers')
plt.xlabel('$q_y (\AA^{-1})$')
plt.ylabel('Intensity')
plt.legend()
plt.tight_layout()
#plt.savefig('perfect_v_random_layers.png')
plt.savefig('perfect_v_random_100pores.png')
plt.show()
