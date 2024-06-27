function [cPCATrials, cPCAControl]  = transitionsCPCA(baseDir,pythonEnvDir,nDims,alphaValue)

load(fullfile(baseDir,'ProcessedData','LabeledTransitions.mat'))
load(fullfile(baseDir,'ProcessedData','EpochedData10ms.mat'))

% get background data
behvFRs = {};
nSamples = [];
backgroundBehaviors = behaviors([7]);

for iBehv = 1:length(backgroundBehaviors)

    %go through each epoch, don't use the first and last 200ms (make sure
    %its fully within the movement period)
    boutFRs = {};
    for iBout = 1:length(behavioralData.(backgroundBehaviors{iBehv}).boutFRs)

        boutDuration = size(behavioralData.(backgroundBehaviors{iBehv}).boutFRs{iBout},2);

        if boutDuration < 41
            continue
        end

        boutFRs{iBout} = behavioralData.(backgroundBehaviors{iBehv}).boutFRs{iBout}(1:254,20:end-20);
    end

    behvFRs{iBehv} = cat(2,boutFRs{:});

    %remove nans
    behvFRs{iBehv}(:,any(isnan(behvFRs{iBehv}))) = [];

    nSamples(iBehv) = size(behvFRs{iBehv},2);

end

% subsample so all behaviors have the same number of time points
for iBehv = 1:length(behvFRs)
    behvFRs{iBehv} = behvFRs{iBehv}(:,randperm(nSamples(iBehv),min(nSamples)));
end

backgroundData = cat(2,behvFRs{:});
backgroundDataCent = backgroundData' - mean(backgroundData');


% get forground data and run for all behaviors
outputs = {};
for iBehv = 1:length(periEventFR)

    %do it for striatum and cortex
    foregroundData = permute(periEventFR{iBehv}{3}(:,170:320,:), [3 2 1]);
    foregroundDataCent = foregroundData(:,:)' - mean(foregroundData(:,:)');

    if isempty(foregroundDataCent)
        continue
    end

    outputs = callcPCA(foregroundDataCent, backgroundDataCent, nDims, alphaValue, pythonEnvDir);

%     [V, ~] = eig(cov(foregroundDataCent) - alphaValue * cov(backgroundDataCent));

    cPCATrials(iBehv).ssaProjs = reshape(outputs,size(foregroundData,2),...
        size(outputs,1)/size(foregroundData,2),size(outputs,2));
    cPCATrials(iBehv).ssaWeights = foregroundDataCent'/outputs';

end

% also do control
controlData = permute(controlFR{3}(:,170:320,:), [3 2 1]);
controlDataCent = controlData(:,:)' - mean(controlData(:,:)');

outputs = callcPCA(controlDataCent, backgroundDataCent, nDims, alphaValue, pythonEnvDir);

cPCAControl.ssaProjs = reshape(outputs,size(controlData,2),...
    size(outputs,1)/size(controlData,2),size(outputs,2));
cPCAControl.ssaWeights = controlDataCent'/outputs';



% 
