import numpy as np
import numpy.random as npr
from matplotlib import pyplot as plt
from scipy.linalg import orth
import time
from os.path import dirname, join

from ssa.models import fit_ssa, weighted_pca, weighted_rrr
from ssa.util import get_sample_weights

import scipy.io as sio
import h5py


baseDIR = 'Z:\\David\\ArenaRecordings\\NeuropixelsTest\\D024-111022-ArenaRecording\\ProcessedData\\SSA\\LabeledTransitions\\'
regionNames = ["Striatum", "Cortex", "Cortex_and_Striatum"] 
nDims = 20

#next, do single population SSA analysis for each of the populations
for iRegion in regionNames:
    loadedMat = sio.loadmat(join(baseDIR,iRegion))
    transData = loadedMat['dataCatCent']
    controlData = loadedMat['controlDataCatCent']
    model, latent, y_pred, losses = fit_ssa(X=transData, R=nDims)
    model_control, latent_control, y_pred_control, losses_control = fit_ssa(X=controlData, R=nDims)

    matSaveDict = {"latents":latent.detach().numpy(), "weights":model.fc2.weight.detach().numpy(), \
        "loss":losses, "latentsControl":latent_control.detach().numpy(), \
        "weightsControl":model_control.fc2.weight.detach().numpy(), "loss_control":losses_control}

    sio.savemat(join(baseDIR,iRegion+'_SSA.mat'),matSaveDict)


#first, do multi-population SSA analysis for cortex and striatum, at varying time lags
shiftMatVars = {}
shiftFile = h5py.File(join(baseDIR,'StrCtxShift.mat'))
for varName, varData in shiftFile.items():
    shiftMatVars[varName] = np.transpose(np.array(varData))

ctxShiftData = shiftMatVars['ctxShiftCatCent']
strShiftData = shiftMatVars['strShiftCatCent']
    
strDataSize = strShiftData.shape
nShifts = strDataSize[2]

shiftLatents_strSource = {}
shiftWeights_strSource = {}
shiftWeights2_strSource = {}
shiftLosses_strSource = {}
shiftLatents_strSourceControl = {}
shiftWeights_strSourceControl = {}
shiftWeights2_strSourceControl = {}
shiftLosses_strSourceControl = {}

shiftLatents_ctxSource = {}
shiftWeights_ctxSource = {}
shiftWeights2_ctxSource = {}
shiftLosses_ctxSource = {}
shiftLatents_ctxSourceControl = {}
shiftWeights_ctxSourceControl = {}
shiftWeights2_ctxSourceControl = {}
shiftLosses_ctxSourceControl = {}

for iShift in range(nShifts):
    thisStrShift = strShiftData[:,:,iShift]

    #first do striatum as the source
    model_strSource, latent_strSource, y_pred_strSource, losses_strSource = \
        fit_ssa(X=thisStrShift, Y=ctxShiftData, R=nDims)

    shiftLatents_strSource["shift"+str(iShift+1)] = latent_strSource.detach().numpy()
    shiftWeights_strSource["shift"+str(iShift+1)] = model_strSource.fc1.weight.detach().numpy()
    shiftWeights2_strSource["shift"+str(iShift+1)] = model_strSource.fc2.weight.detach().numpy()
    shiftLosses_strSource["shift"+str(iShift+1)] = losses_strSource

    #do control trials
    loadedMat = sio.loadmat(join(baseDIR,'StrCtxControlShift_'+str(iShift+1)+'.mat'))
    strControlData = loadedMat['strControlShiftCatCent']
    ctxControlData = loadedMat['ctxControlShiftCatCent']
    model_strSourceControl, latent_strSourceControl, y_pred_strSourceControl, losses_strSourceControl = \
        fit_ssa(X=strControlData, Y=ctxControlData, R=nDims)

    shiftLatents_strSourceControl["shift"+str(iShift+1)] = latent_strSourceControl.detach().numpy()
    shiftWeights_strSourceControl["shift"+str(iShift+1)] = model_strSourceControl.fc1.weight.detach().numpy()
    shiftWeights2_strSourceControl["shift"+str(iShift+1)] = model_strSourceControl.fc2.weight.detach().numpy()
    shiftLosses_strSourceControl["shift"+str(iShift+1)] = losses_strSourceControl

    #next do cortex as the source
    model_ctxSource, latent_ctxSource, y_pred_ctxSource, losses_ctxSource = \
        fit_ssa(X=ctxShiftData, Y=thisStrShift, R=nDims)

    shiftLatents_ctxSource["shift"+str(iShift)] = latent_ctxSource.detach().numpy()
    shiftWeights_ctxSource["shift"+str(iShift)] = model_ctxSource.fc1.weight.detach().numpy()
    shiftWeights2_ctxSource["shift"+str(iShift)] = model_ctxSource.fc2.weight.detach().numpy()
    shiftLosses_ctxSource["shift"+str(iShift)] = losses_ctxSource

    #do control trials with cortex as source
    model_ctxSourceControl, latent_ctxSourceControl, y_pred_ctxSourceControl, losses_ctxSourceControl = \
    fit_ssa(X=ctxControlData, Y=ctxControlData, R=nDims)

    shiftLatents_ctxSourceControl["shift"+str(iShift+1)] = latent_ctxSourceControl.detach().numpy()
    shiftWeights_ctxSourceControl["shift"+str(iShift+1)] = model_ctxSourceControl.fc1.weight.detach().numpy()
    shiftWeights2_ctxSourceControl["shift"+str(iShift+1)] = model_ctxSourceControl.fc2.weight.detach().numpy()
    shiftLosses_ctxSourceControl["shift"+str(iShift+1)] = losses_ctxSourceControl

    #save this shift (just in case)
    matSaveDict = {"strSourceLatents":latent_strSource.detach().numpy(), \
        "strSourceWeights":model_strSource.fc1.weight.detach().numpy(), \
        "strSourceWeights2":model_strSource.fc2.weight.detach().numpy(), \
        "strSourceLoss":losses_strSource, "strSourceLatentsControl":latent_strSourceControl.detach().numpy(), \
        "strSourceWeightsControl":model_strSourceControl.fc1.weight.detach().numpy(), \
        "strSourceWeights2Control":model_strSourceControl.fc2.weight.detach().numpy(), \
        "strSourceLoss_control":losses_strSourceControl, \
        "ctxSourceLatents":latent_ctxSource.detach().numpy(), \
        "ctxSourceWeights":model_ctxSource.fc1.weight.detach().numpy(), \
        "ctxSourceWeights2":model_ctxSource.fc2.weight.detach().numpy(), \
        "ctxSourceLoss":losses_ctxSource, "ctxSourceLatentsControl":latent_ctxSourceControl.detach().numpy(), \
        "ctxSourceWeightsControl":model_ctxSourceControl.fc1.weight.detach().numpy(), \
        "ctxSourceWeights2Control":model_ctxSourceControl.fc2.weight.detach().numpy(), \
        "ctxSourceLoss_control":losses_ctxSourceControl}

    sio.savemat(join(baseDIR,'shiftSSA'+str(iShift+1)+'.mat'),matSaveDict)

#save all shifts
matSaveDict = {"shiftLatents_strSource":shiftLatents_strSource, \
        "shiftWeights_strSource":shiftWeights_strSource, \
        "shiftWeights2_strSource":shiftWeights2_strSource, \
        "shiftLosses_strSource":shiftLosses_strSource, \
        "shiftLatents_strSourceControl":shiftLatents_strSourceControl, \
        "shiftWeights_strSourceControl":shiftWeights_strSourceControl, \
        "shiftWeights2_strSourceControl":shiftWeights2_strSourceControl, \
        "shiftLosses_strSourceControl":shiftLosses_strSourceControl, \
        "shiftLatents_ctxSource":shiftLatents_ctxSource, \
        "shiftWeights_ctxSource":shiftWeights_ctxSource, \
        "shiftWeights2_ctxSource":shiftWeights2_ctxSource, \
        "shiftLosses_ctxSource":shiftLosses_ctxSource, \
        "shiftLatents_ctxSourceControl":shiftLatents_ctxSourceControl, \
        "shiftWeights_ctxSourceControl":shiftWeights_ctxSourceControl, \
        "shiftWeights2_ctxSourceControl":shiftWeights2_ctxSourceControl, \
        "shiftLosses_ctxSourceControl":shiftLosses_ctxSourceControl}

sio.savemat(join(baseDIR,'shiftSSA.mat'),matSaveDict)

print('Done')
