import numpy as np
import numpy.random as npr
from matplotlib import pyplot as plt
from scipy.linalg import orth
import time
from os.path import dirname, join

from contrastive import CPCA

import scipy.io as sio

cpca = CPCA(standardize=False, n_components=nDims)
projected_data = cpca.fit_transform(foregroundMat, backgroundMat, alpha_selection='manual', alpha_value=alpha)
