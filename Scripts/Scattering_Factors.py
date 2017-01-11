#!/usr/bin/bash


class ScatteringFactor(object):
    """
    Scattering factors for individual atoms are characterized by the following parameters:

    An equation used to calculate the scattering factor as a function of q is given as follows:
        f(sin(theta)/lambda) = f(q/4*pi) = (sum from i= 1 to i = 4 of ai*exp(-bi*sin^2(theta)/lambda^2)) + c

    Attributes:
        name: A string representing the atom's name.
        (The following are parameters that are found in International Crystallographic Tables, Volume C, Table 6.1.1.3)
        a: A list of the a parameters used to fit the above equation
        b: A list of b parameters used to fit in the above equation
        c: A single parameter c used to fit the above equation
    """

    def __init__(self, name, a, b, c):
        self.name = name
        self.a = a
        self.b = b
        self.c = c

NA = ScatteringFactor('NA', [3.25650, 3.93620, 1.39980, 1.00320], [2.66710, 6.11530, 0.200100, 14.0390], 0.404000)
