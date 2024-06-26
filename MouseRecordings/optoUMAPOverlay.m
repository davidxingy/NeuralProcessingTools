% function allChansFiltData = optoUMAPOverlay(baseDir)
clearvars -except D040Sessions
% baseDir = 'X:\David\ArenaRecordings\D040-110223-ArenaRecording';

baseDirs = D040Sessions;

for iSession = 1:length(baseDirs)

    % load UMAP projection
    if iSession == 1
        load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'analyzedBehaviors','behvLabelsNoArt','regionBehvAssignments','regionWatershedLabels',...
            'gridXInds','gridYInds','watershedRegions','annotatedBehvLabels');
        % load in sync data
        load(fullfile(baseDirs{iSession},'ProcessedData','VideoSyncFrames.mat'))
    end

    load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'reduction','origDownsampEMGInd','regionAssignmentsNoBound')

    % load EMG data
    load(fullfile(baseDirs{iSession},'ProcessedData','EMG1ms'))

    %get meta data as well
    sessionFiles = string(ls(fullfile(baseDirs{iSession},'ProcessedData')));
    metaDataFile = sessionFiles(contains(sessionFiles,'ProcessedEMG_MetaData'));
    load(fullfile(baseDirs{iSession},'ProcessedData', metaDataFile))

    sessTotalPulses = length(downsampleLaserOnsetInds);

    if iSession == 2
        %         usedPulseInds = downsampleLaserOnsetInds(1:3900);
    end

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

    nShuffs = 100;
    % nullPulseMinOffset = 700;
    % nullPulseMaxOffset = 1000;
    % for iShuff = 1:nShuffs
    %     %get null distribution pulses
    %     for iPulse = 1:length(downsampleLaserOnsetInds)
    %
    %         gotPulse = false;
    %         if iPulse == 1
    %             tryPulse = 10000 + randi(nullPulseMaxOffset-nullPulseMinOffset) + nullPulseMinOffset;
    %         else
    %             tryPulse = nullPulses(iShuff,iPulse-1) + randi(nullPulseMaxOffset-nullPulseMinOffset) + nullPulseMinOffset;
    %         end
    %         while ~gotPulse
    %             if any(tryPulse >= (controlPulseWindowStarts-windowSize - 50) & tryPulse <= (pulseWindowEnds + 50))
    %                 tryPulse = pulseWindowEnds(find(tryPulse <= (pulseWindowEnds + 50),1)) + 50 + randi(nullPulseMaxOffset-nullPulseMinOffset);
    %             elseif tryPulse > max(origDownsampEMGInd)
    %                 error('Unable to get enough null pulses')
    %             else
    %                 nullPulses(iShuff,iPulse) = tryPulse;
    %                 gotPulse = true;
    %             end
    %         end
    %
    %         nullPulseWindowStarts(iShuff,iPulse) = nullPulses(iShuff,iPulse);
    %         nullPulseWindowEnds(iShuff,iPulse) = nullPulses(iShuff,iPulse) + windowSize;
    %
    %     end
    % end



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

    stimUseIndsCell = {};
    stimControlIndsCell = {};
    for iPulse = 1:length(pulseWindowStarts)

        stimUseIndsCell{iPulse} = find(origDownsampEMGInd>pulseWindowStarts(iPulse) & origDownsampEMGInd<pulseWindowEnds(iPulse));
        stimControlIndsCell{iPulse} = find(origDownsampEMGInd>controlPulseWindowStarts(iPulse) & origDownsampEMGInd<controlPulseWindowEnds(iPulse));

        %     for iShuff = 1:nShuffs
        %         stimNullInds{iShuff}{iPulse} = find(origDownsampEMGInd>nullPulseWindowStarts(iShuff,iPulse) & origDownsampEMGInd<nullPulseWindowEnds(iShuff,iPulse));
        %     end

    end

    stimUseInds = cat(2,stimUseIndsCell{:});
    stimNoUseInds = setdiff(1:length(origDownsampEMGInd),stimUseInds);

    stimControlInds = cat(2,stimControlIndsCell{:});
    stimControlNoUseInds = setdiff(1:length(origDownsampEMGInd),[stimControlInds stimUseInds]);

    % for iShuff = 1:nShuffs
    %     stimNullShuffInds{iShuff} = cat(2,stimNullInds{iShuff}{:});
    %     origDownsampEMGIndStimNull{iShuff} = origDownsampEMGInd(stimNullShuffInds{iShuff});
    %     emgStimNull{iShuff} = downsampEMG(:,origDownsampEMGIndStimNull{iShuff});
    % end

    origDownsampEMGIndStim = origDownsampEMGInd(stimUseInds);
    origDownsampEMGIndNoStim = origDownsampEMGInd(stimNoUseInds);

    origDownsampEMGIndStimControl = origDownsampEMGInd(stimControlInds);
    origDownsampEMGIndNoStimControl = origDownsampEMGInd(stimControlNoUseInds);


    emgStim = downsampEMG(:,origDownsampEMGIndStim);
    emgNoStim = downsampEMG(:,origDownsampEMGIndNoStim);

    emgStimControl = downsampEMG(:,origDownsampEMGIndStimControl);
    emgNoStimControl = downsampEMG(:,origDownsampEMGIndNoStimControl);

    % get subsample to build null distribution
    for iShuff = 1:nShuffs

        shuffInds = randperm(length(stimNoUseInds),length(stimUseInds));
        emgNoStimShuff{iShuff} = downsampEMG(:,origDownsampEMGIndNoStim(shuffInds));
        stimNoUseShuffInds{iShuff} = stimNoUseInds(shuffInds);

        shuffInds = randperm(length(stimControlNoUseInds),length(stimControlInds));
        emgNoStimShuffControl{iShuff} = downsampEMG(:,origDownsampEMGIndNoStimControl(shuffInds));
        stimNoUseShuffIndsControl{iShuff} = stimControlNoUseInds(shuffInds);

    end

    sessionsReductions{iSession} = reduction;
    if iSession == 1
        reducIndsOffset = 0;
    else
        reducIndsOffset = sum(cellfun(@(x) size(x,1),sessionsReductions(1:iSession-1)));
    end

    sessionsStimUseInds{iSession} = stimUseInds + reducIndsOffset;
    sessionsStimNoUseInds{iSession} = stimNoUseInds + reducIndsOffset;
    sessionsStimControlInds{iSession} = stimControlInds + reducIndsOffset;
    sessionsStimControlNoUseInds{iSession} = stimControlNoUseInds + reducIndsOffset;
    for iShuff = 1:nShuffs
        sessionsNoUseShuffInds{iSession,iShuff} = stimNoUseShuffInds{iShuff} + reducIndsOffset;
        sessionsNoUseControlShuffInds{iSession,iShuff} = stimNoUseShuffIndsControl{iShuff} + reducIndsOffset;
    end

    sessionsEmgStim{iSession} = emgStim;
    sessionsEmgNoStim{iSession} = emgNoStim;
    sessionsEmgStimControl{iSession} = emgStimControl;
    sessionsEmgNoStimControl{iSession} = emgNoStimControl;
    for iShuff = 1:nShuffs
        sessionsEmgNoStimShuff{iSession,iShuff} = emgNoStimShuff{iShuff};
        sessionsEmgNoStimControlShuff{iSession,iShuff} = emgNoStimShuffControl{iShuff};
    end

end

% consolicate all the sessions together for calculating heatmaps
allReductions = cat(1,sessionsReductions{:});
allStimUseInds = cat(2,sessionsStimUseInds{:});
allStimNoUseInds = cat(2,sessionsStimNoUseInds{:});
allStimControlInds = cat(2,sessionsStimControlInds{:});
allStimControlNoUseInds = cat(2,sessionsStimControlNoUseInds{:});
for iShuff = 1:nShuffs
    allNoUseShuffInds{iShuff} = cat(2,sessionsNoUseShuffInds{:,iShuff});
    allNoUseControlShuffInds{iShuff} = cat(2,sessionsNoUseControlShuffInds{:,iShuff});
end
clear sessionsNoUseShuffInds sessionsNoUseControlShuffInds
allEmgStim = cat(2,sessionsEmgStim{:});
allEmgNoStim = cat(2,sessionsEmgNoStim{:});
allEmgStimControl = cat(2,sessionsEmgStimControl{:});
allEmgNoStimControl = cat(2,sessionsEmgNoStimControl{:});
for iShuff = 1:nShuffs
    allEmgNoStimShuff{iShuff} = cat(2,sessionsEmgNoStimShuff{:,iShuff});
end
clear sessionsEmgNoStimShuff
for iShuff = 1:nShuffs
    allEmgNoStimControlShuff{iShuff} = cat(2,sessionsEmgNoStimControlShuff{:,iShuff});
end
clear sessionsEmgNoStimControlShuff

% calculate heat maps for both the stim and no stim time points
densityGaussStd = 0.1;

% make input struct for the heat map functions
stimHeatMapInput.reduction = allReductions;
stimHeatMapInput.origDownsampEMGInd = origDownsampEMGInd;
stimHeatMapInput.gridXInds = gridXInds;
stimHeatMapInput.gridYInds = gridYInds;
stimHeatMapInput.watershedRegions = watershedRegions;
stimHeatMapInput.annotatedBehvLabels = [];
stimHeatMapInput.densityGaussStd = densityGaussStd;

stimHeatMapInput.reducSigs = allEmgStim;
stimHeatMapInput.reducIndsToUse = allStimUseInds;

noStimHeatMapInput = stimHeatMapInput;
noStimHeatMapInput.reducSigs = allEmgNoStim;
noStimHeatMapInput.reducIndsToUse = allStimNoUseInds;

controlStimHeatMapInput = stimHeatMapInput;
controlStimHeatMapInput.reducSigs = allEmgStimControl;
controlStimHeatMapInput.reducIndsToUse = allStimControlInds;

controlNoStimHeatMapInput = stimHeatMapInput;
controlNoStimHeatMapInput.reducSigs = allEmgNoStimControl;
controlNoStimHeatMapInput.reducIndsToUse = allStimControlNoUseInds;

for iShuff = 1:nShuffs
    shuffHeatMapInput(iShuff) = stimHeatMapInput;
    shuffHeatMapInput(iShuff).reducSigs = allEmgNoStimShuff{iShuff};
    shuffHeatMapInput(iShuff).reducIndsToUse = allNoUseShuffInds{iShuff};

    shuffControlHeatMapInput(iShuff) = stimHeatMapInput;
    shuffControlHeatMapInput(iShuff).reducSigs = allEmgNoStimControlShuff{iShuff};
    shuffControlHeatMapInput(iShuff).reducIndsToUse = allNoUseControlShuffInds{iShuff};
end

% do heat map generation
[heatMapStim, heatMapFiltStim, gaussKernal, stimNPoints] = calcActivityMap(stimHeatMapInput);
[heatMapNoStim, heatMapFiltNoStim, noStimNPoints] = calcActivityMap(noStimHeatMapInput);

% do it for the control as well
[heatMapStimControl, heatMapFiltStimControl, gaussKernal] = calcActivityMap(controlStimHeatMapInput);
[heatMapNoStimControl, heatMapFiltNoStimControl] = calcActivityMap(controlNoStimHeatMapInput);

% the effect size is the difference between the stim and no stim
heatMapDiff = heatMapNoStim - heatMapStim;
heatMapDiffControl = heatMapNoStimControl - heatMapStimControl;

% apply gaussian filter for visualization
for iChan = 1:size(heatMapDiff,3)
    chanDiff = heatMapDiff(:,:,iChan);
    chanDiff(isnan(chanDiff)) = 0;
    heatMapDiffFilt(:,:,iChan) = fftshift(real(ifft2(fft2(gaussKernal).*fft2(chanDiff))));
end

% now get the null distribution effect size heat maps
for iShuff = 1:nShuffs

    shuffHeatMapInput = stimHeatMapInput;
    shuffHeatMapInput.reducSigs = allEmgNoStimShuff{iShuff};
    shuffHeatMapInput.reducIndsToUse = allNoUseShuffInds{iShuff};

    [heatMapShuff(:,:,:,iShuff), heatMapFiltShuff(:,:,:,iShuff), ~ , shuffNPoints] = ...
        calcActivityMap(shuffHeatMapInput);
    heatMapShuffDiff(:,:,:,iShuff) = heatMapNoStim - heatMapShuff(:,:,:,iShuff);

%     [heatMapShuff(:,:,:,iShuff), heatMapFiltShuff(:,:,:,iShuff)] = ...
%         calcActivityMap(baseDir,emgStimNull{iShuff},stimNullShuffInds{iShuff},densityGaussStd);
%     heatMapShuffDiff(:,:,:,iShuff) = heatMapNoStim - heatMapShuff(:,:,:,iShuff);
%     
%     for iChan = 1:size(heatMapShuffDiff,3)
%         heatMapShuffDiffFilt(:,:,iChan,iShuff) = fftshift(real(ifft2(fft2(gaussKernal).*fft2(heatMapShuffDiff(:,:,iChan,iShuff)))));
%     end
    
    shuffControlHeatMapInput = stimHeatMapInput;
    shuffControlHeatMapInput.reducSigs = allEmgNoStimControlShuff{iShuff};
    shuffControlHeatMapInput.reducIndsToUse = allNoUseControlShuffInds{iShuff};

    [heatMapShuffControl(:,:,:,iShuff), heatMapFiltShuffControl(:,:,:,iShuff)] = ...
        calcActivityMap(shuffControlHeatMapInput);
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



function [heatMap, heatMapFilt, gaussKernal, nPointsMap] = calcActivityMap(inputs)

reduction = inputs.reduction;
origDownsampEMGInd = inputs.origDownsampEMGInd;
gridXInds = inputs.gridXInds;
gridYInds = inputs.gridYInds;
watershedRegions = inputs.watershedRegions;
annotatedBehvLabels = inputs.annotatedBehvLabels;
reducSigs = inputs.reducSigs;
reducIndsToUse = inputs.reducIndsToUse;
densityGaussStd = inputs.densityGaussStd;

reducIndsToUse = reducIndsToUse-200;
badInds = find(reducIndsToUse <= 0);
reducIndsToUse(badInds) = [];
reducSigs(:,badInds) = [];

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
    binNPointsSig = binNPoints(accumarrayRealOutputs);

    spaceAveZScores = nan(nGridPoints,nGridPoints);
    spaceAveNPoints = nan(nGridPoints,nGridPoints);
    binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
    binAveZScores = binAveSig;

    %map the firing rates into the actual bins within the grid now (using
    %the bin labels as the index). Only the bins with anything in them are
    %mapped
    spaceAveZScores(unique(binLabels)) = binAveZScores;
    spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));
    spaceAveNPoints(unique(binLabels)) = binNPointsSig;

    heatMap(:,:,iChan) = spaceAveZScores;
    nPointsMap(:,:,iChan) = spaceAveNPoints;
    heatMapFilt(:,:,iChan) = spaceAveZScoresFilt;

end

end


% 

