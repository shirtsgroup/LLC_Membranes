#!/usr/bin/env python

import numpy as np
from LLC_Membranes.timeseries.ctrwsim import CTRW
from LLC_Membranes.timeseries.fbmsim import FractionalBrownianMotion
from LLC_Membranes.machine_learning.features import Features
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import tqdm


labels_dict = {"CTRW": 0, "FBM": 1, "sFBM": 2}
models = ["CTRW", "FBM", "sFBM"]
# create n trajectories where a fraction, split, are CTRWs and the rest follow FBM
splits = [.5, .5, 0]  # percent that are CTRW

# fration of data to use to train model
train_frac = 0.75  # fraction of data used to train model

# Trajectory details
steps = 2000  # number of steps in each trajectory
ntraj = 1000
hop_sigma = 1
alpha = 0.6
noise = 0

n = [int(s*ntraj) for s in splits]
dividers = np.cumsum(n)

labels = np.zeros(ntraj)
for i, d in enumerate(dividers):
    if i == 0:
        labels[:dividers[i]] = labels_dict[models[i]]
    else:
        labels[dividers[i - 1]:dividers[i]] = labels_dict[models[i]]

all_features = np.zeros([ntraj, 4])
all_features[:, 0] = labels

# Generate continuous time random walk trajectories
ctrw = CTRW(steps, n[0], hop_dist='Gaussian', dwell_dist='power', hop_sigma=hop_sigma, alpha=alpha,
            dt=1, H=0.5, padding=1)
ctrw.generate_trajectories(fixed_time=True, noise=noise)

ctrw_features = Features(ctrw.z_interpolated.T[:, :, np.newaxis])

ctrw_features.compute_msd(axis=0)
ctrw_features.determine_stationarity()
ctrw_features.compute_autocorrelation()


all_features[:dividers[0], 1:] = np.concatenate((ctrw_features.msd[:, np.newaxis],
                                                     ctrw_features.stationarity[:, np.newaxis],
                                                     ctrw_features.autocorrelation[:, np.newaxis]), axis=1)

# Generate fractional brownian motion trajectories
hurst_distribution = np.random.uniform(size=n[1])*0.5
scale = 1
fbm_trajectories = np.zeros([steps, n[1]])
print('Generating FBM trajectories...')
for i, h in enumerate(tqdm.tqdm(hurst_distribution)):
    fbm_trajectories[:, i] = FractionalBrownianMotion(steps, h, length=steps, scale=scale, progress=False).trajectories[:-1, 0]

fbm_features = Features(fbm_trajectories[:, :, np.newaxis])
fbm_features.compute_msd(axis=0)
fbm_features.determine_stationarity()
fbm_features.compute_autocorrelation()

all_features[dividers[0]:dividers[1], 1:] = np.concatenate((fbm_features.msd[:, np.newaxis],
                                                            fbm_features.stationarity[:, np.newaxis],
                                                            fbm_features.autocorrelation[:, np.newaxis]), axis=1)

# Generate sFBM trajectories
# sfbm = CTRW(steps, int(split*ntraj), hop_dist='Gaussian', dwell_dist='power', hop_sigma=hop_sigma, alpha=alpha,
#             dt=1, H=H, padding=1, )

np.random.shuffle(all_features)  # shuffle all of the rows of the feature array
RF = RandomForestClassifier(n_estimators=100)
RF.fit(all_features[:int(train_frac*ntraj), 1:], all_features[:int(train_frac*ntraj), 0])
predictions = RF.predict(all_features[int(train_frac*ntraj):, 1:])

# Accuracy of predictions
print(sum(predictions == all_features[int(train_frac*ntraj):, 0]) / len(predictions))
print(RF.feature_importances_)