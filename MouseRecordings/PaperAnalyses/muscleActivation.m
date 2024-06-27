clear 
close all

use10msBins = false;

% if just want to do single region decoding without testing across regions
% (to save time)
doCrossRegion = true;

% baseDir = 'X:\David\ArenaRecordings\D026-032923-ArenaRecording';
% recordingName = 'D026-032923-ArenaRecording';

baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording';
recordingName = 'D024-111022-ArenaRecording';

% baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording';
% recordingName = 'D020-062922-ArenaRecording';

load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'regionAssignmentsFiltered','reduction','origDownsampEMGInd',...
    'behvLabelsDown','analyzedBehaviors','regionWatershedLabels','regionAssignmentsNoBound','regionBehvAssignments')

% load(fullfile(baseDir,'ProcessedData','RegionEpochedData10ms'),'regionFRsOrig','regionEMGs','origDownsampEMGInd',...
%     'behvLabelsDown','analyzedBehaviors','regionWatershedLabels')

cd(fullfile(baseDir,'ProcessedData'))

% get EMG
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
load(fullfile(baseDir,'ProcessedData','neuronDataStruct.mat'))
load(fullfile(baseDir,'Neuropixels','artifactTimestamps.mat'))
% load(fullfile(baseDir,'ProcessedData','EpochedData1ms.mat'))
load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))
load(fullfile(baseDir,'ProcessedData',[recordingName '_ProcessedEMG_MetaData.mat']))

badEMGChans = [];

% load neural data
if use10msBins
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins10msGauss.mat'), 'cortexInds', 'striatumInds','allFRs')
else
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds','allFRs')
end
nNeurons = size(allFRs,1);

% downsample EMG further to 10ms if using 10ms binned FRs
if use10msBins
    for iChan = 1:size(downsampEMG,1)
        avedSig = conv(repmat(0.1,1,10),downsampEMG(iChan,:));
        emg10ms(iChan,:) = avedSig(6:10:end);
    end
    downsampEMG = emg10ms;
end

nonNanNeurons = find(~all(isnan(allFRs),2));

% go through biceps to find muscle activations
fs = 1000;
threshChan = 4;
thresh = 70;
baselineFlucThresh = 50;
baselineChanThresh = 50;
neurBinSize = 1; %in ms

preThreshDuration = 500; %in ms
preThreshbuffer = 50;

neurWindow = [200 100]; %in ms

threshCrossings = find(downsampEMG(threshChan,2:end) > thresh & downsampEMG(threshChan,1:end-1) <= thresh);





% 

