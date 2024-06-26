function calcFiringRates(baseRecordingFolder, binSize, gaussStd, depthCtxStr, depthStr, minSpikes, minAmp)
% function calcFiringRates(baseRecordingFolder, binSize, gaussStd, depthCtxStr, depthStr)
% 
% function to generate neural firing rate time series from a set of
% timestamps. Will bin acorrding to binSize and convolve with a guassian
% with standard deviation size gaussStd. For the mouse-arena project, will
% also segment the neurons into cortex and striatum based on depth. If no
% input is provided for those, will just use default values 1000mm and
% 2600mm as end depths for M1 and Str respectively. Will also remove
% neurons with really low spikes (set by minSpikes, default 100) or really
% low waveform amplitudes (set by minAmp, default 0.3)

load(fullfile(baseRecordingFolder,'ProcessedData','VideoSyncFrames.mat'))
load(fullfile(baseRecordingFolder,'ProcessedData','neuronDataStruct.mat'))
load(fullfile(baseRecordingFolder,'Neuropixels','artifactTimestamps.mat'))

maxSamples = frameNeuropixelSamples{1}{end}(end);

if isempty(depthCtxStr)
    depthCtxStr = 2600; %in mm
end
if isempty(depthStr)
    depthStr = 1000; %in mm
end
if isempty(minSpikes)
    minSpikes = 100; %in number of spikes
end
if isempty(minAmp)
    minAmp = 0.3;
end


% bin and smooth FRs
for iNeuron = 1:length(neuronDataStruct)
    
    spikesPerBin = histcounts(neuronDataStruct(iNeuron).timeStamps,1:(binSize*30):maxSamples);
    smoothedFRs(iNeuron,:) = convGauss(spikesPerBin, binSize, gaussStd,0);
    
end

% remove neurons with too small amplitude or too few spikes
rejectedNeurons = unique([find(cellfun(@length,{neuronDataStruct.timeStamps})<100)...
    find([neuronDataStruct.amplitude]<0.3)]);

% segment into cortex and striatum
cortexInds = setdiff(find([neuronDataStruct.depth] > depthCtxStr), rejectedNeurons);
striatumInds = setdiff(find([neuronDataStruct.depth] < depthCtxStr & [neuronDataStruct.depth] > depthStr), rejectedNeurons);
cortexFRs = smoothedFRs(cortexInds,:);
striatumFRs = smoothedFRs(striatumInds,:);
smoothedFRs(rejectedNeurons,:)=[];
allFRs = smoothedFRs;
noNanFRs = smoothedFRs;

% set time bins that have artifact to NaN
artifacts = histcounts(artifactTS,1:(30*binSize):maxSamples);
artifactBins = find(convGauss(artifacts, binSize, gaussStd,0));

cortexFRs(:,artifactBins) = nan;
striatumFRs(:,artifactBins) = nan;
allFRs(:,artifactBins) = nan;

% save
save(fullfile(baseRecordingFolder,'ProcessedData',['NeuralFiringRates' num2str(binSize) 'msBins' num2str(gaussStd) 'msGauss']),...
    'cortexFRs','striatumFRs','allFRs','noNanFRs','cortexInds','striatumInds','rejectedNeurons',"-v7.3")



% 
