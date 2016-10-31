import numpy as np
import pdb
import matplotlib.pyplot as plt

# start with a number of particles
Nparticles = 100
Nbootstraps = 2000

# Time step in ps
dt = 0.1

# what fraction of the MSD graph do I show
fracshow = .3  # 30% of the graph

# Total time in ps simulated
T = 200

# number of points
nT = int(T/dt)

# Dimensions: do it in 3D
Dim = 3

# State array: how many states are there? i.e. bound and unbound
Nstates = 2

# State rate matrix: what is the rate of going back and forth between states?
# Diagonals are 0 (can't really define a rate going from self to self)

rmatrix = np.zeros([Nstates, Nstates], dtype=float)
rmatrix[0, 1] = 1  # timescale in ps from going from 1->2   # you can adjust these
rmatrix[1, 0] = 100  # timescale in ps from going from 2->1   # you can adjust these

# State transition probability: the odds of being in state i at one dt and state j at the next dt.
# Off diagonals : time step / rmatrix[i, j] --> If it takes 1 ps to go from state 1 to 2, and the time step is 0.1 ps
# then the probability of it making that transition at a given time, t, is 0.1 ps / 1 ps = 0.1
# The diagonals are 1 - sum of the other rates. unitless (depends on the time step)
tmatrix = np.zeros([Nstates, Nstates], dtype=float)
for i in range(0, Nstates):
    rowsum = 0
    for j in range(0, Nstates):
        if i != j:
            tmatrix[i, j] = dt/rmatrix[i, j]
            rowsum += tmatrix[i, j]
    tmatrix[i, i] = 1.0 - rowsum

equilibrium = np.zeros([Nstates, Nstates], dtype=float)
# Equilibrium ratios of the different states i to j - unitless
# assume detailed balance such that rate_1/rate_2 = probability_1/probability_2, determined from rmatrix
for i in range(0, Nstates):
    equilibrium[i, i] = 1  # normalization in each row is 1, other states are in proportion to the amount of i.
    for j in range(0, Nstates):
        if i != j:
            equilibrium[i, j] = rmatrix[j, i]/rmatrix[i, j]

# Define the state array for each particle
states = np.zeros([Nparticles, nT], dtype=int)

# define the type of motion.  For now, we will assume it is diffusive,
# indicated by 'D' We can use 'B' for bound, and 'P' for some
# arbitrary distance distribution. Only 'D' and 'B' are implemented for now.

#types = np.array(['D']*Nstates)  # adjust these
types = np.array(['D','B'])  # the first is pure diffusive, the 2nd is bound

# constants for transport.
# Used both types = 'D' and 'B', units are per state in nm^2 / ps.
# diffusion constants
dconsts = np.zeros([Nstates], dtype = float)
dconsts[0] = 0.002                   # fiddle with these!  NOTE: even if the diffusion is zero, switching into a diffusive state resets the center for the bound state.  I think this is probably the correct behavior.
dconsts[1] = 0.005                    # you can fiddle with these!

# if types = 'B', we need both a center location (per particle,
# defined when the state is turned on)), a diffusion constant (from above) and
# a force constant kconsts, in units  beta * k =  nm^-2
centers = np.zeros([Dim, Nparticles], dtype = float)  # start as all zeros
kconsts = np.zeros([Nstates], dtype = float)
kconsts[0] = 0
kconsts[1] = 50  # how tightly it is bound in the bound state

# Define the location array for each particle at each time point
x = np.zeros([Dim, Nparticles, nT], dtype=float)

# we can leave these as zeros because of now we really only care about
# dx.  For more complex models, that are concentration dependent, we
# would need to figure out what particle was at each point.

# first, initialize the state matrix. Figure out the total weight. Can
# be done with any state equilibrium

totalw = 0
initial_prob = np.zeros(Nstates, dtype=float)  # [0, 0]
for i in range(0, Nstates):
   totalw += equilibrium[0, i]
for i in range(0, Nstates):
    initial_prob[i] = equilibrium[0, i]/totalw

# now, initialize the state according to the equilibrium probability distribution
probstate = np.random.rand(Nparticles)  # fills an array of length Nparticles, with random numbers drawn between 0 and 1

# for each particles, see
for n in range(0, Nparticles):
    for i in range(0, Nstates):
        if probstate[n] < initial_prob[i]:
            states[n, 0] = i  # states has a column for each time step. Each row is a different particle at each step.
            # The number that ends up in states[n, 0] refers to what state its in
            break
        else:
            probstate[n] -= initial_prob[i]  # not sure why this is here?

for t in range(nT-1):  # now fill in the rest of the time steps

    # we need a random number for each state to determine if it switches state
    probstate = np.random.rand(Nparticles)

    for n in range(0, Nparticles):
        state = states[n, t]  # state of particle n at time t starting from initialized state distribution

        # Determine how far it moves.
        if types[states[n, t]] == 'D':
            # the displacement will be drawn from a random distribution with mean zero and
            # standard deviation sqrt(2Ddt) in each direction --> i.e. Gaussian
            # see for example http://iopscience.iop.org/article/10.1088/1478-3967/1/3/001/fulltext/#pb183249bib26
            s = np.sqrt(2*dconsts[state]*dt)  # standard deviation
            dvecs = s*np.random.randn(3)  # https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.randn.html
            # np.random.randn gives three random numbers chosen from the standard normal distribution
            # for random sampling, use sigma*np.random.randn() + mu (mu = mean, sigma = standard deviation)
            x[:, n, t + 1] = x[:, n, t] + dvecs  # Set the particle position in the next time step
        elif types[states[n, t]] == 'B':  # If the particle is in the bound state
            if t > 0 and types[states[n, t-1]] != 'B':  # if it just switched into the bound state, it's new location is the center.
                centers[:, n] = x[:, n, t]
            # We want a Gaussian process with some correlation time, since we don't want it to jump,
            # See http://www.ks.uiuc.edu/Services/Class/PHYS498/LectureNotes/chp4.pdf, equation 4.120, or https://softphys.files.wordpress.com/2014/03/lang_eqn_harmonic_field.pdf
            # should be ornstein-furth formula? <(dx)^2> = 2kT/mb^2 ( bt - 1 + exp(-betat)
            # possibly https://ocw.mit.edu/courses/mathematics/18-366-random-walks-and-diffusion-fall-2006/study-materials/lec14.pdf
            # I think this is right: a random force, and a drift in the direction determined by the force constant.
            s = np.sqrt(2*dconsts[state]*dt)
            dvecs = s*np.random.randn(3) - dconsts[state]*kconsts[state]*dt*(x[:,n,t]-centers[:,n])
            x[:, n, t+1] = x[:, n,t] + dvecs
        else:
            print "This non-diffusive motion not implemented!"
            exit()

        # we now need to determine if the state moves to state i. It
        # will move with probability dt/rate, as stored in
        # tmatrix[i,j]

        for j in range(0, Nstates):
            if probstate[n] < tmatrix[state, j]:
                states[n, t+1] = j
                break
            else:
                probstate[n] -= tmatrix[state, j]  # again -- why?

exit()
# Now calculate the MSD of the particles as a function of time
MSD = np.zeros([nT], dtype=float)  # find msd at each time step
itau = 0

#### straightforward - slow?
if 0:
    while itau < nT:
        ncount = Nparticles*(nT-itau)
        for n in range(0, Nparticles):
            xn = x[:, n, :]  # a trajectory of a single particle - x,y,z position of particle n for all time steps
            for t in range(nT-itau):  # run for total number of time points from itau to nT
                xo = xn[:, t + itau] - xn[:, t]
                MSD[itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
        MSD[itau] /= ncount
        itau += 1

##### using Fast Fourier transforms - from
##### http://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft


def autocorrFFT(x):
    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n= N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)
    return res/n  # this is the autocorrelation in convention A

def msd_fft(r):
    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = Q/(N-m)
    return S1-2*S2

MSD = np.zeros([nT], dtype=float)
MSDs = np.zeros([nT, Nparticles], dtype=float)  # a set of MSDs per particle
for n in range(0, Nparticles):
    MSDs[:, n] = msd_fft(x[:, n, :].T)
    MSD += MSDs[:, n]
MSD /= Nparticles


eMSDs = np.zeros([nT, Nbootstraps], dtype=float)  # a set of MSDs per particle (time step?)
# Now, bootstrap over number of particles, assuming that the particles are sufficiently independent
for b in range(0, Nbootstraps):
    indices = np.random.randint(0, Nparticles, Nparticles)  # randomly select N of the particles with replacement
    # ^ makes a list of length Nparticles each with a random number from 0 to  Nparticles
    for n in range(0, Nparticles):
        eMSDs[:, b] += MSDs[:, indices[n]]  # for this bootstrap trial b, add the MSDs of a randomly selected particle
        # to the particle at each time step
    eMSDs[:, b] /= Nparticles  # Divide every timestep by Nparticles -- average the MSDs

# The result is Nbootstraps trajectories

limits = np.zeros([2, nT], dtype=float)  # a set of MSDs per particle
# now, let's determine a 95\% error bound for each tau (out of the
# Nbootstrapped MSD's, use that for the error bars
for t in range(0, nT):
    limits[0, t] = np.abs(np.percentile(eMSDs[t, :], 2.5) - MSD[t])
    limits[1, t] = np.abs(np.percentile(eMSDs[t, :], 97.5) - MSD[t])

# plot just the first half of the MSD (avoid the noisy part)
endMSD = int(np.floor(nT*fracshow))
# plot only 100 bars total
errorevery = int(np.ceil(fracshow*nT/100.0))
plt.errorbar(dt*np.array(range(0,endMSD)),MSD[:endMSD],yerr=[limits[0,:endMSD],limits[1,:endMSD]],errorevery=errorevery)
plt.ylabel('MSD')
plt.xlabel('time (ps)')
# I'm not sure why the error appears when plt.show is done.
plt.show()