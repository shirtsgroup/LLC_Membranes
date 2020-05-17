#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from hdphmm.generate_timeseries import GenARData
from LLC_Membranes.llclib import file_rw

com, dt = file_rw.load_object('/home/ben/github/hdphmm/notebooks/comr.pl')

plt.plot(com[:, 2, 1])
plt.show()

