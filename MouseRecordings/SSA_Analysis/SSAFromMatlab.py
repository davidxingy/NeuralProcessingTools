import numpy as np
import numpy.random as npr
from matplotlib import pyplot as plt
from scipy.linalg import orth
import time
import os
from os.path import dirname, join
import matlab.engine

from ssa.models import fit_ssa, weighted_pca, weighted_rrr
from ssa.util import get_sample_weights

import scipy.io as sio

eng = matlab.engine.start_matlab()
baseDir = os.getcwd()

if multiPop:
    model, latent, y_pred, losses = \
        fit_ssa(X=sourceData, Y=targetData, R=nDims)

    matSaveDict = {"latents":latent.detach().numpy(), "weightsSource":model.fc1.weight.detach().numpy(), \
        "weightsTarget":model.fc2.weight.detach().numpy(), "loss":losses}

else:
    print(nDims)

    model, latent, y_pred, losses = fit_ssa(X=sourceData, R=nDims)

    matSaveDict = {"latents":latent.detach().numpy(), "weightsSource":model.fc2.weight.detach().numpy(), \
        "weightsTarget":[], "loss":losses}

sio.savemat(join(baseDir,'TmpSaveData.mat'),matSaveDict)    