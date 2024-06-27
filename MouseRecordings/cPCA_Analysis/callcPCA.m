function outputs = callcPCA(foregroundMat, backgroundMat, nDims, alpha, pathToEnvExe)
% set the python enviornment to cPCA
currentPyEnv = pyenv;

if ~strcmpi(currentPyEnv.Executable,pathToEnvExe)
    pyenv('Version',pathToEnvExe)
end

thisFilePath = mfilename('fullpath');
thisFileFolder = fileparts(thisFilePath);
outputs = pyrunfile(fullfile(thisFileFolder,"cPCAFromMatlab.py"),"projected_data",...
    foregroundMat = foregroundMat,backgroundMat = backgroundMat, nDims = int32(nDims), alpha = alpha);

outputs = double(outputs);


% 