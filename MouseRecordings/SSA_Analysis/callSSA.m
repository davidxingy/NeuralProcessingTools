function ssaResults = callSSA(sourceData, targetData, nDims, pathToEnvExe)
% set the python enviornment to SSA
currentPyEnv = pyenv;

if ~strcmpi(currentPyEnv.Executable,pathToEnvExe)
    terminate(pyenv)
    pyenv('Version',pathToEnvExe,"ExecutionMode","OutOfProcess")
end

thisFilePath = mfilename('fullpath');
thisFileFolder = fileparts(thisFilePath);

if isempty(targetData)
%     [trajs, weightsSource, modelLoss] = pyrunfile(fullfile(thisFileFolder,"SSAFromMatlab.py"),string(["trajs","weightsSource","weightsTarget","modelLoss"]),...
%         sourceData = py.numpy.array(sourceData), nDims = int32(nDims), multiPop = false);
    matSaveDict = pyrunfile(fullfile(thisFileFolder,"SSAFromMatlab.py"), string("matSaveDict"),...
        sourceData = py.numpy.array(sourceData), nDims = int32(nDims), multiPop = false);

else
%     [trajs, weightsSource, weightsTarget, modelLoss] = pyrunfile(fullfile(thisFileFolder,"SSAFromMatlab.py"),string(["trajs","weightsSource","weightsTarget","modelLoss"]),...
%         sourceData = py.numpy.array(sourceData),targetData = py.numpy.array(targetData), nDims = int32(nDims), multiPop = true);

    matSaveDict = pyrunfile(fullfile(thisFileFolder,"SSAFromMatlab.py"),string(["matSaveDict"]),...
        sourceData = py.numpy.array(sourceData),targetData = py.numpy.array(targetData), nDims = int32(nDims), multiPop = true);
end

load('TmpSaveData')
ssaResults.trajs = double(latents);
ssaResults.weightsSource = double(weightsSource);
ssaResults.weightsTarget = double(weightsTarget);
ssaResults.loss = double(loss);

% outputs = struct(matSaveDict);
% ssaResults.trajs = double(outputs.latents);
% ssaResults.weightsSource = single(outputs.weightsSource);
% ssaResults.weightsTarget = single(outputs.weightsTarget);
% ssaResults.losses = double(outputs.losses);


% 