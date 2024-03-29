% function allChansFiltData = optoUMAPOverlay(baseDir)

baseDir = 'X:\David\ArenaRecordings\D036-101123-ArenaRecording';

% load UMAP projection
load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','watershedRegions','annotatedBehvLabels')

% load in sync data
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))


% load EMG data
load(fullfile(baseDir,'ProcessedData', 'EMG1ms.mat'))

%get meta data as well
sessionFiles = string(ls(fullfile(baseDir,'ProcessedData')));
metaDataFile = sessionFiles(contains(sessionFiles,'ProcessedEMG_MetaData'));

if isempty(metaDataFile)
    error('Unable to find EMG Meta data file!')
end
load(fullfile(baseDir,'ProcessedData', metaDataFile))

% get stim windows
windowSize = 200;
pulseWindowStarts = [];
pulseWindowEnds= [];

controlPulseWindowStarts = [];
controlPulseWindowEnds = [];
for iPulse = 1:length(downsampleLaserOnsetInds)
    stimWindow{iPulse} = downsampleLaserOnsetInds(iPulse):downsampleLaserOnsetInds(iPulse)+windowSize;
    pulseWindowStarts(iPulse) = downsampleLaserOnsetInds(iPulse);
    pulseWindowEnds(iPulse) = downsampleLaserOnsetInds(iPulse)+windowSize;
    
    controlPulseWindowStarts(iPulse) = downsampleLaserOnsetInds(iPulse)-windowSize;
    controlPulseWindowEnds(iPulse) = downsampleLaserOnsetInds(iPulse);
end
allStimWindowInds = unique(cat(2,stimWindow{:}));

% get emg corresponding to the UMAP
% origDownsampEMGInd = origDownsampEMGInd(annotatedBehvLabels(1:50:end));
reducEMGInds = origDownsampEMGInd;

% % % get the subset of timepoints which correspond to stim windows and no stim
% % % windows
% % stimUseInds = [];
% % stimNoUseInds = [];
% % for iReducInd = 1:length(origDownsampEMGInd)
% %     if any(origDownsampEMGInd(iReducInd)==allStimWindowInds)
% %         stimUseInds = [stimUseInds iReducInd];
% %     else
% %         stimNoUseInds = [stimNoUseInds iReducInd];
% %     end
% % end

stimUseInds = {}; 
for iPulse = 1:length(pulseWindowStarts)

    stimUseInds{iPulse} = find(origDownsampEMGInd>pulseWindowStarts(iPulse) & origDownsampEMGInd<pulseWindowEnds(iPulse));
    stimControlInds{iPulse} = find(origDownsampEMGInd>controlPulseWindowStarts(iPulse) & origDownsampEMGInd<controlPulseWindowEnds(iPulse));

end

stimUseInds = cat(2,stimUseInds{:});
stimNoUseInds = setdiff(1:length(origDownsampEMGInd),stimUseInds);

stimControlInds = cat(2,stimControlInds{:});
stimControlNoUseInds = setdiff(1:length(origDownsampEMGInd),stimControlInds);

origDownsampEMGIndStim = origDownsampEMGInd(stimUseInds);
origDownsampEMGIndNoStim = origDownsampEMGInd(stimNoUseInds);

origDownsampEMGIndStimControl = origDownsampEMGInd(stimControlInds);
origDownsampEMGIndNoStimControl = origDownsampEMGInd(stimControlNoUseInds);

emgStim = downsampEMG(:,origDownsampEMGIndStim);
emgNoStim = downsampEMG(:,origDownsampEMGIndNoStim);

emgStimControl = downsampEMG(:,origDownsampEMGIndStimControl);
emgNoStimControl = downsampEMG(:,origDownsampEMGIndNoStimControl);

% get subsample to build null distribution
nShuffs = 100;
for iShuff = 1:nShuffs
    shuffInds = randperm(length(stimNoUseInds),length(stimUseInds));
    emgNoStimShuff{iShuff} = downsampEMG(:,origDownsampEMGIndNoStim(shuffInds));
    stimNoUseShuffInds{iShuff} = stimNoUseInds(shuffInds);
    
    shuffInds = randperm(length(stimControlNoUseInds),length(stimControlInds));
    emgNoStimShuffControl{iShuff} = downsampEMG(:,origDownsampEMGIndNoStimControl(shuffInds));
    stimNoUseShuffIndsControl{iShuff} = stimControlNoUseInds(shuffInds);
end

densityGaussStd = 0.2;

[heatMapStim, heatMapFiltStim, gaussKernal] = calcActivityMap(baseDir,emgStim,stimUseInds,densityGaussStd);
[heatMapNoStim, heatMapFiltNoStim] = calcActivityMap(baseDir,emgNoStim,stimNoUseInds,densityGaussStd);

[heatMapStimControl, heatMapFiltStimControl, gaussKernal] = calcActivityMap(baseDir,emgStimControl,stimControlInds,densityGaussStd);
[heatMapNoStimControl, heatMapFiltNoStimControl] = calcActivityMap(baseDir,emgNoStimControl,stimControlNoUseInds,densityGaussStd);

heatMapDiff = heatMapNoStim - heatMapStim;
heatMapDiffControl = heatMapNoStimControl - heatMapStimControl;

for iChan = 1:size(heatMapDiff,3)
    heatMapDiffFilt(:,:,iChan) = fftshift(real(ifft2(fft2(gaussKernal).*fft2(heatMapDiff(:,:,iChan)))));
end

for iShuff = 1:nShuffs
    [heatMapShuff(:,:,:,iShuff), heatMapFiltShuff(:,:,:,iShuff)] = ...
        calcActivityMap(baseDir,emgNoStimShuff{iShuff},stimNoUseShuffInds{iShuff},densityGaussStd);
    heatMapShuffDiff(:,:,:,iShuff) = heatMapNoStim - heatMapShuff(:,:,:,iShuff);
    
    for iChan = 1:size(heatMapShuffDiff,3)
        heatMapShuffDiffFilt(:,:,iChan,iShuff) = fftshift(real(ifft2(fft2(gaussKernal).*fft2(heatMapShuffDiff(:,:,iChan,iShuff)))));
    end
    
    [heatMapShuffControl(:,:,:,iShuff), heatMapFiltShuffControl(:,:,:,iShuff)] = ...
        calcActivityMap(baseDir,emgNoStimShuffControl{iShuff},stimNoUseShuffIndsControl{iShuff},densityGaussStd);
    heatMapShuffDiffControl(:,:,:,iShuff) = heatMapNoStimControl - heatMapShuffControl(:,:,:,iShuff);
    
end

figure
tiledlayout(2,4)
for iChan = 1:size(heatMapDiff,3)
    nexttile
    imagesc(heatMapFiltNoStim(:,:,iChan))
    title(channelNames{iChan})
    colorbar
end

figure
tiledlayout(2,4)
for iChan = 1:size(heatMapDiff,3)
    nexttile
    imagesc(heatMapDiffFilt(:,:,iChan))
    title(channelNames{iChan})
    colorbar
end


figure
tiledlayout(2,4)
for iChan = 1:size(heatMapDiff,3)
    nexttile
    diffMap = heatMapDiff(:,:,iChan) - nanmean(heatMapShuffDiff(:,:,iChan,:),4);
    zScoreMap = diffMap./(nanstd(heatMapShuffDiff(:,:,iChan,:),[],4)+100);
%     zScoreMap(std(heatMapShuffDiffFilt(:,:,iChan,:),[],4)<30) = 0;
    zScoreMap(isnan(zScoreMap)) = 0;
    zScoreMapFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(zScoreMap))));
%     imagesc(zScoreMap)
    imagesc(zScoreMapFilt)
    title(channelNames{iChan})
    colorbar
end
              
figure
tiledlayout(2,4)
for iChan = 1:size(heatMapDiff,3)
    nexttile
    diffMap = heatMapDiffControl(:,:,iChan) - nanmean(heatMapShuffDiffControl(:,:,iChan,:),4);
    zScoreMap = diffMap./(nanstd(heatMapShuffDiffControl(:,:,iChan,:),[],4)+100);
%     zScoreMap(std(heatMapShuffDiffFilt(:,:,iChan,:),[],4)<30) = 0;
    zScoreMap(isnan(zScoreMap)) = 0;
    zScoreMapFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(zScoreMap))));
    %     imagesc(zScoreMap)
    imagesc(zScoreMapFilt)
    title(channelNames{iChan})
    colorbar
end

x=1;

% 



function [heatMap, heatMapFilt, gaussKernal] = calcActivityMap(baseDir,reducSigs,reducIndsToUse,densityGaussStd)
load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','watershedRegions','annotatedBehvLabels')

% get watershedding region splitting
[regionBoundaryIndsX, regionBoundaryIndsY] = find(watershedRegions==0);

% Generate grid
nGridPoints = length(gridXInds);
rangeVals = [min(reduction(reducIndsToUse,1)) max(reduction(reducIndsToUse,1)); ...
    min(reduction(reducIndsToUse,2)) max(reduction(reducIndsToUse,2))]*1.5;

gridIndsX = gridXInds;
gridIndsY = gridYInds;
[meshGridX,meshGridY] = meshgrid(gridIndsX,gridIndsY);

% get the gaussian kernal to convolve with for better visualization
gaussKernal = exp(-.5.*((meshGridX-mean([max(max(meshGridX)) min(min(meshGridX))])).^2 + (meshGridY-mean([max(max(meshGridY)) min(min(meshGridY))])).^2)./densityGaussStd^2) ./ (2*pi*densityGaussStd^2);
gaussKernal = gaussKernal./(sum(sum(gaussKernal)));

% for each point in the UMAP reduction, get the bin which it belongs to when doing 2D binning using the grid points
[~,~,~,binIndX, binIndY] = histcounts2(reduction(reducIndsToUse,1),reduction(reducIndsToUse,2),gridIndsX,gridIndsY);
binLabels = sub2ind([nGridPoints, nGridPoints], binIndX, binIndY);

% for each 2D bin, get how many time points are inside it (i.e. how many of
% the umap reduction points are in each bin)
binNPoints = accumarray(binLabels, ones(1,length(binLabels)));

% also, sometimes the bin has no reduction points, in which case, don't use
% it for firing rate calculations (e.g. it's a on the corner of the grid)
accumarrayRealOutputs = find(binNPoints~=0);

for iChan = 1:size(reducSigs,1)

    binSigSums = accumarray(binLabels,reducSigs(iChan,:));
    binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*1000;

    spaceAveZScores = nan(nGridPoints,nGridPoints);
    binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
    binAveZScores = binAveSig;

    %map the firing rates into the actual bins within the grid now (using
    %the bin labels as the index). Only the bins with anything in them are
    %mapped
    spaceAveZScores(unique(binLabels)) = binAveZScores;
    spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));

    heatMap(:,:,iChan) = spaceAveZScores;
    heatMapFilt(:,:,iChan) = spaceAveZScoresFilt;

end

end


% 

