#! /usr/bin/env python

import numpy as np
import pandas as pd
from mcfly import modelgen, find_architecture, storage
from LLC_Membranes.timeseries.ctrwsim import CTRW
from LLC_Membranes.timeseries.fbmsim import FractionalBrownianMotion
from keras.models import load_model
import tqdm

labels_dict = {"CTRW": 0, "FBM": 1, "sFBM": 2}
diffusion_models = ["CTRW", "FBM", "sFBM"]
# create n trajectories where a fraction, split, are CTRWs and the rest follow FBM
splits = [.5, .5, 0]  # percent that are CTRW

# fration of data to use to train model
fracs = [0.75, 0.05, 0.2]  # fraction of trajectories for [training, validation, testing]

# Trajectory details
steps = 2000  # number of steps in each trajectory
ntraj = 1000  # number of trajectories (i.e. number of samples)

# mcfly parameters
n_models = 20  # generate 'n_models' different deep learning models with random parameters. This will help decide the
# best parameters to use for training the full data set.

trajectories = np.zeros([ntraj, steps])  # initialize array for all trajectories

n = [int(s*ntraj) for s in splits if s != 0]  # number of trajectories to generate for each diffusion model
dividers = np.cumsum(n)

# labels organized according to one-hot-encoding
labels = np.zeros([ntraj, len(n)])
for i, d in enumerate(dividers):
    if i == 0:
        labels[:dividers[i], labels_dict[diffusion_models[i]]] = 1
    else:
        labels[dividers[i - 1]:dividers[i], labels_dict[diffusion_models[i]]] = 1

# CTRW trajectories. Should modify to make a distribution of all of these parameters
hop_sigma = 1
alpha = 0.6
noise = 0.1

print('Generating CTRW trajectories...')
ctrw = CTRW(steps, n[0], hop_dist='Gaussian', dwell_dist='power', hop_sigma=hop_sigma, alpha=alpha,
            dt=1, H=0.5, padding=1)
ctrw.generate_trajectories(fixed_time=True, noise=noise)

trajectories[:dividers[0]] = ctrw.z_interpolated

# Generate fractional brownian motion trajectories
hurst_distribution = np.random.uniform(size=n[1])*0.5
scale = 1
fbm_trajectories = np.zeros([steps, n[1]])
print('Generating FBM trajectories...')
for i, h in enumerate(tqdm.tqdm(hurst_distribution)):
    fbm_trajectories[:, i] = FractionalBrownianMotion(steps, h, length=steps, scale=scale, progress=False).trajectories[:-1, 0]

trajectories[dividers[0]:dividers[1]] = fbm_trajectories.T

# shuffle the trajectories so we can ensure randomized data for training, validating and testing
shuffle_indices = np.random.permutation(ntraj)
np.random.shuffle(shuffle_indices)  # shuffle all of the rows of the feature array
trajectories = trajectories[shuffle_indices, :]
labels = labels[shuffle_indices, :]

# generate and test accuracy of models
models = modelgen.generate_models((ntraj, steps, 1), number_of_classes=len(n), number_of_models=n_models)

# models_to_print = range(len(models))
# for i, item in enumerate(models):
#     if i in models_to_print:
#         model, params, model_types = item
#         print("-------------------------------------------------------------------------------------------------------")
#         print("Model " + str(i))
#         print(" ")
#         print("Hyperparameters:")
#         print(params)
#         print(" ")
#         print("Model description:")
#         model.summary()
#         print(" ")
#         print("Model type:")
#         print(model_types)
#         print(" ")

set_dividers = np.cumsum([int(n*ntraj) for n in fracs])

# Partition Data into training, validation and test sets
X_train = trajectories[:set_dividers[0], :, np.newaxis]  # need to add an axis for X values.
X_val = trajectories[set_dividers[0]:set_dividers[1], :, np.newaxis]
X_test = trajectories[set_dividers[1]:, :, np.newaxis]
y_train_binary = labels[:set_dividers[0], :]
y_val_binary = labels[set_dividers[0]:set_dividers[1], :]  # labels of validation data
y_test_binary = labels[set_dividers[1]:, :]

outputfile = "test.json"

# test all of the models with a subset of the data
epochs = 5  # number of times to train
histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(X_train, y_train_binary,
                                                                                  X_val, y_val_binary, models,
                                                                                  nr_epochs=epochs, subset_size=500,
                                                                                  verbose=True, outputfile=outputfile)

print('Value accuracies:\n', val_accuracies)  # honestly not sure what this is

best_model_index = np.argmax(val_accuracies)  # choose best model based on highest accuracy
best_model, best_params, best_model_types = models[best_model_index]

print('Best Model Type: %s' % best_model_types)
print('Best Parameters: %s' % best_params)

# now train using the best model with all of the data
history = best_model.fit(X_train, y_train_binary, epochs=1, validation_data=(X_val, y_val_binary))

# Construct confusion matrix : columns are predicted, rows are truth
probs = best_model.predict_proba(X_val)
predicted = probs.argmax(axis=1)
y_index = y_val_binary.argmax(axis=1)  # basically find where the nonzero values are in each row
confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
confusion_matrix.index = [diffusion_models[i] for i in confusion_matrix.index]
confusion_matrix.columns = [diffusion_models[i] for i in confusion_matrix.columns]
confusion_matrix.reindex(columns=[l for l in diffusion_models], fill_value=0)
print('Confusion Matrix:\n', confusion_matrix)

score_test = best_model.evaluate(X_test, y_test_binary, verbose=True)
print('Score of best model: ' + str(score_test))
modelname = 'best_model.h5'
best_model.save(modelname)
