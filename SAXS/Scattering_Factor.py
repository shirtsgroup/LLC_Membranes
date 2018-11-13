# calculate atomic scattering factor

import math
import numpy as np
import matplotlib.pyplot as plt

# fitted parameters for sodium (International Crystallographic Tables, Volume C, Table 6.1.1.3)
a = [3.25650, 3.93620, 1.39980, 1.00320]
b = [2.66710, 6.11530, 0.200100, 14.0390]
c = 0.404000


def scattering_factor(a, b, c, d):
    f = 0
    for i in range(0, 4):
        f += a[i]*math.exp(-b[i]*d**2)
    f += c
    return f


braggs = np.linspace(0, 2, 200)  # values of 1/2d which is equal to sin(theta)/lamdba by the braggs relationship

f = np.zeros((len(braggs), 1))
for i in range(0, len(f)):
    f[i] = scattering_factor(a, b, c, braggs[i])

plt.figure(1)
plt.plot(braggs, f)
plt.show()
