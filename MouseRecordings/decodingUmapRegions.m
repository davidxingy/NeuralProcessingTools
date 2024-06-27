clear 
close all

use10msBins = true;

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

% normalize EMG and firing rate
normalizedEMG = (downsampEMG-nanmean(downsampEMG,2))./nanstd(downsampEMG,[],2);
normalizedFRs = (allFRs-nanmean(allFRs,2))./nanstd(allFRs,[],2);

nonNanNeurons = find(~all(isnan(normalizedFRs),2));

% regionWatershedLabels = regionWatershedLabels(behvAlignPerm);

% go through each region, get the time inds for the region
for iRegion = 1:length(regionWatershedLabels)

    regionTimeInds = find(regionAssignmentsNoBound == regionWatershedLabels(iRegion));  
    regionTimeInds(regionTimeInds > floor(frameEMGSamples{1}{end}(end)/20)) = [];

    if use10msBins
        origEmgRegionInds = unique(round(origDownsampEMGInd(regionTimeInds)/10));
        origEmgRegionInds(origEmgRegionInds==0) = [];
    else
        origEmgRegionInds = origDownsampEMGInd(regionTimeInds);
    end
    regionEMGsInds{iRegion} = origEmgRegionInds;

    if use10msBins
        currentDir = pwd;
        cd(fullfile(baseDir,'ProcessedData'))
        regionNeurInds{iRegion} = round(NeurEMGSync(origEmgRegionInds*200,...
            frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
        cd(currentDir)

        outOfBoundInds = find(regionNeurInds{iRegion} > floor(frameNeuropixelSamples{1}{end}(end)/300));
        nanInds = find(isnan(regionNeurInds{iRegion}));
        zeroInds = find(regionNeurInds{iRegion}==0);
        
        regionNeurInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
        regionEMGsInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
    else
        regionNeurInds{iRegion} = round(NeurEMGSync(origDownsampEMGInd(regionTimeInds)*20,...
            frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
        outOfBoundInds = regionNeurInds{iRegion} > floor(frameNeuropixelSamples{1}{end}(end)/30);
        nanInds = find(isnan(regionNeurInds{iRegion}));
        zeroInds = find(regionNeurInds{iRegion}==0);
        regionNeurInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
        regionEMGsInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
    end
    regionEMGs{iRegion} = normalizedEMG(:,regionEMGsInds{iRegion});
    regionFRs{iRegion} = normalizedFRs(:,regionNeurInds{iRegion});
    regionFRsOrig{iRegion} = allFRs(:,regionNeurInds{iRegion});

    %remove Nans
    nanEMGInds = find(any(isnan(regionEMGs{iRegion})));
    nanFRInds = find(any(isnan(regionFRs{iRegion}(nonNanNeurons,:))));
    regionEMGs{iRegion}(:,unique([nanEMGInds nanFRInds])) = [];
    regionFRs{iRegion}(:,unique([nanEMGInds nanFRInds])) = [];
    regionFRsOrig{iRegion}(:,unique([nanEMGInds nanFRInds])) = [];
    regionNeurInds{iRegion}(unique([nanEMGInds nanFRInds])) = [];
    regionEMGsInds{iRegion}(unique([nanEMGInds nanFRInds])) = [];

    neurIndsToUse = 1:size(regionFRs{iRegion},1);

    regionBehvNames{iRegion} = string(join(analyzedBehaviors(regionBehvAssignments{iRegion}),'/'));
end

% % don't use those artifact points
% downsampEMG(:,artifactInEMGInds) = nan;

% go through each region and find neurons with low firing rates
leastPoints = min(cellfun(@length,regionFRs));
regionMeanFRsCell = cellfun(@(x) nanmean(x,2),regionFRsOrig,'UniformOutput',false);
regionMeanFRs = cat(2,regionMeanFRsCell{:});

goodNeuronsAll = find(any(regionMeanFRs*(1000/(10*use10msBins))>0.2,2));
reallyGoodNeuronsAll = find(all(regionMeanFRs*(1000/(10*use10msBins))>0.2,2));

goodNeuronsStr = intersect(goodNeuronsAll,1:length(striatumInds));
goodNeuronsCtx = intersect(goodNeuronsAll,length(striatumInds)+1:size(normalizedFRs,1));

reallyGoodNeuronsStr = intersect(reallyGoodNeuronsAll,1:length(striatumInds));
reallyGoodNeuronsCtx = intersect(reallyGoodNeuronsAll,length(striatumInds)+1:size(normalizedFRs,1));

nRegions = length(unique(regionAssignmentsFiltered));

%decoding parameters
nHist = 5;
histSkip = 1;
nFolds = 10;

% go through each region and get chunks of continuous data
minChunkDuration = 20;

regionBothHist = {};
regionStrHist = {};
regionCtxHist = {};
regionEMGHist = {};
for iRegion = 1:nRegions

    regionInds = regionNeurInds{iRegion};
    regionChunkStarts = [1, find(diff(regionInds)~=1)+1];
    regionChunkStops = [find(diff(regionInds)~=1), length(regionInds)];

    goodChunks = find((regionChunkStops-regionChunkStarts)>minChunkDuration);

    regionChunkStarts = regionChunkStarts(goodChunks);
    regionChunkStops = regionChunkStops(goodChunks);

    % add history
    bothHist = {};
    strHist = {};
    ctxHist = {};
    emgHist = {};
    for iChunk = 1:length(regionChunkStarts)

        emgInds = regionEMGsInds{iRegion}(regionChunkStarts(iChunk):regionChunkStops(iChunk));
        neurInds = regionNeurInds{iRegion}(regionChunkStarts(iChunk):regionChunkStops(iChunk));

        if any(isnan(neurInds))
            neurHist{iChunk} = [];
            emgHist{iChunk} = [];
            continue
        end

        bothHist{iChunk} = addHistory(allFRs(goodNeuronsAll,neurInds), nHist, 2, histSkip, false);
        strHist{iChunk} = addHistory(allFRs(goodNeuronsStr,neurInds), nHist, 2, histSkip, false);
        ctxHist{iChunk} = addHistory(allFRs(goodNeuronsCtx,neurInds), nHist, 2, histSkip, false);
        emgHist{iChunk} = downsampEMG(:,emgInds);
        
        bothHist{iChunk}(:,1:nHist*histSkip) = [];
        strHist{iChunk}(:,1:nHist*histSkip) = [];
        ctxHist{iChunk}(:,1:nHist*histSkip) = [];
        emgHist{iChunk}(:,1:nHist*histSkip,:) = [];
    end

    %don't use any points which have artifacts
    
    regionBothHist{iRegion} = cat(2,bothHist{:});
    regionStrHist{iRegion} = cat(2,strHist{:});
    regionCtxHist{iRegion} = cat(2,ctxHist{:});
    regionEMGHist{iRegion} = cat(2,emgHist{:});

    clear bothHist strHist ctxHist emgHist
    
    nanInds = any(isnan(regionBothHist{iRegion}),1) | any(isnan(regionEMGHist{iRegion}),1);

    regionBothHist{iRegion}(:,nanInds) = [];
    regionStrHist{iRegion}(:,nanInds) = [];
    regionCtxHist{iRegion}(:,nanInds) = [];
    regionEMGHist{iRegion}(:,nanInds) = [];

    %z-score the EMG
    regionEMGHistNormalized{iRegion} = (regionEMGHist{iRegion}-mean(regionEMGHist{iRegion},2))./std(regionEMGHist{iRegion},[],2);

end

clear allFRs normalizedFRs
nTrainPoints = floor(min([cellfun(@(x) size(x,2),regionBothHist)])/nFolds * (nFolds-1));
% allBadNeurons = unique(cat(1,regionBadNeurons{:}));
% goodNeurons = setdiff(1:size(allFRs,1),allBadNeurons);
% goodNeuronsHist = {};
% for iHist = 1:nHist+1
%     goodNeuronsHist{iHist} = goodNeurons+(iHist-1)*nNeurons;
% end
% allGoodNeuronsHist = cat(2,goodNeuronsHist{:});

% train models
for iRegion = 1:nRegions 

    %get indices for testing and training
    foldTestInds{iRegion} = divideBlocks(1:size(regionBothHist{iRegion},2), nFolds);

    %do cross validation
    for iFold = 1:nFolds

        trainInds = setdiff(1:size(regionBothHist{iRegion},2),foldTestInds{iRegion}{iFold});
        
        %downsample so same number of training points across all behavior
        %regions, if doing cross-region decoding
        if doCrossRegion
            trainInds = trainInds(randperm(length(trainInds),nTrainPoints));
        end

        emgCent = regionEMGHistNormalized{iRegion}(:,trainInds)-mean(regionEMGHistNormalized{iRegion}(:,trainInds),2);
        bothCent = regionBothHist{iRegion}(:,trainInds)-mean(regionBothHist{iRegion}(:,trainInds),2);
        strCent = regionStrHist{iRegion}(:,trainInds)-mean(regionStrHist{iRegion}(:,trainInds),2);
        ctxCent = regionCtxHist{iRegion}(:,trainInds)-mean(regionCtxHist{iRegion}(:,trainInds),2);

        %     %get optimal ridge regression parameter usgin CV
        %     optimalParam = findBestRidgeParam(neurCent,emgCent, [0.001 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5]);

        for iChan = 1:size(emgCent,1)
            decoderWeightsBoth{iRegion,iFold}(iChan,:) = ridge(emgCent(iChan,:)', bothCent', 0.1);
            decoderWeightsStr{iRegion,iFold}(iChan,:) = ridge(emgCent(iChan,:)', strCent', 0.1);
            decoderWeightsCtx{iRegion,iFold}(iChan,:) = ridge(emgCent(iChan,:)', ctxCent', 0.1);
        end

        %get polynomial fit
        thisEMGEstBoth = bothCent'*decoderWeightsBoth{iRegion,iFold}';
        thisEMGEstStr = strCent'*decoderWeightsStr{iRegion,iFold}';
        thisEMGEstCtx = ctxCent'*decoderWeightsCtx{iRegion,iFold}';

        for iChan = 1:size(thisEMGEstBoth,2)
            polyFitWeightsBoth{iRegion,iFold}(iChan,:) = polyfit(thisEMGEstBoth(:,iChan),emgCent(iChan,:),3);
            polyFitWeightsStr{iRegion,iFold}(iChan,:) = polyfit(thisEMGEstStr(:,iChan),emgCent(iChan,:),3);
            polyFitWeightsCtx{iRegion,iFold}(iChan,:) = polyfit(thisEMGEstCtx(:,iChan),emgCent(iChan,:),3);
        end

        if doCrossRegion
            % get data for training model based on all regions combined
            combTrainInds{iRegion,iFold} = trainInds(randperm(nTrainPoints,floor(nTrainPoints/nRegions)));
            combTrainEMG{iRegion,iFold} = regionEMGHistNormalized{iRegion}(:,combTrainInds{iRegion,iFold});
            combTrainBoth{iRegion,iFold} = regionBothHist{iRegion}(:,combTrainInds{iRegion,iFold});
            combTrainStr{iRegion,iFold} = regionStrHist{iRegion}(:,combTrainInds{iRegion,iFold});
            combTrainCtx{iRegion,iFold} = regionCtxHist{iRegion}(:,combTrainInds{iRegion,iFold});

            %get data for training model for just 2 regions
            twoTrainInds{iRegion,iFold} = trainInds(randperm(nTrainPoints,floor(nTrainPoints/2)));
            twoTrainEMG{iRegion,iFold} = regionEMGHistNormalized{iRegion}(:,twoTrainInds{iRegion,iFold});
            twoTrainBoth{iRegion,iFold} = regionBothHist{iRegion}(:,twoTrainInds{iRegion,iFold});
            twoTrainStr{iRegion,iFold} = regionStrHist{iRegion}(:,twoTrainInds{iRegion,iFold});
            twoTrainCtx{iRegion,iFold} = regionCtxHist{iRegion}(:,twoTrainInds{iRegion,iFold});
        end

    end
end

if doCrossRegion
    % train model on 2 regions combined
    for iRegion1 = 1:nRegions
        for iRegion2 = 1:nRegions
            if iRegion1 > iRegion2

                for iFold = 1:nFolds
                    %combine the data from the two regions
                    thisTwoTrainEMG = cat(2,twoTrainEMG{iRegion1,iFold},twoTrainEMG{iRegion2,iFold});
                    thisTwoTrainBoth = cat(2,twoTrainBoth{iRegion1,iFold},twoTrainBoth{iRegion2,iFold});
                    thisTwoTrainStr = cat(2,twoTrainStr{iRegion1,iFold},twoTrainStr{iRegion2,iFold});
                    thisTwoTrainCtx = cat(2,twoTrainCtx{iRegion1,iFold},twoTrainCtx{iRegion2,iFold});

                    twoTrainEMGCent = thisTwoTrainEMG - nanmean(thisTwoTrainEMG,2);
                    twoTrainBothCent = thisTwoTrainBoth - nanmean(thisTwoTrainBoth,2);
                    twoTrainStrCent = thisTwoTrainStr - nanmean(thisTwoTrainStr,2);
                    twoTrainCtxCent = thisTwoTrainCtx - nanmean(thisTwoTrainCtx,2);

                    % linear weights
                    for iChan = 1:size(emgCent,1)
                        decoderWeightsTwoBoth{iRegion1,iRegion2,iFold}(iChan,:) = ridge(twoTrainEMGCent(iChan,:)', twoTrainBothCent', 0.1);
                        decoderWeightsTwoStr{iRegion1,iRegion2,iFold}(iChan,:) = ridge(twoTrainEMGCent(iChan,:)', twoTrainStrCent', 0.1);
                        decoderWeightsTwoCtx{iRegion1,iRegion2,iFold}(iChan,:) = ridge(twoTrainEMGCent(iChan,:)', twoTrainCtxCent', 0.1);
                    end

                    %get polynomial fit
                    twoEMGEstBoth = twoTrainBothCent'*decoderWeightsTwoBoth{iRegion1,iRegion2,iFold}';
                    twoEMGEstStr = twoTrainStrCent'*decoderWeightsTwoStr{iRegion1,iRegion2,iFold}';
                    twoEMGEstCtx = twoTrainCtxCent'*decoderWeightsTwoCtx{iRegion1,iRegion2,iFold}';

                    for iChan = 1:size(twoEMGEstBoth,2)
                        polyFitWeightsTwoBoth{iRegion1,iRegion2,iFold}(iChan,:) = polyfit(twoEMGEstBoth(:,iChan),twoTrainEMGCent(iChan,:),3);
                        polyFitWeightsTwoCtx{iRegion1,iRegion2,iFold}(iChan,:) = polyfit(twoEMGEstStr(:,iChan),twoTrainEMGCent(iChan,:),3);
                        polyFitWeightsTwoStr{iRegion1,iRegion2,iFold}(iChan,:) = polyfit(twoEMGEstCtx(:,iChan),twoTrainEMGCent(iChan,:),3);
                    end

                end

            end
        end
    end

    % train model based on all regions combined
    for iFold = 1:nFolds
        combTrainEMGCat = cat(2,combTrainEMG{:,iFold});
        combTrainBothCat = cat(2,combTrainBoth{:,iFold});
        combTrainStrCat = cat(2,combTrainStr{:,iFold});
        combTrainCtxCat = cat(2,combTrainCtx{:,iFold});

        combTrainEMGCent = combTrainEMGCat - nanmean(combTrainEMGCat,2);
        combTrainBothCent = combTrainBothCat - nanmean(combTrainBothCat,2);
        combTrainStrCent = combTrainStrCat - nanmean(combTrainStrCat,2);
        combTrainCtxCent = combTrainCtxCat - nanmean(combTrainCtxCat,2);

        % linear weights
        for iChan = 1:size(combTrainEMGCent,1)
            decoderWeightsCombBoth{iFold}(iChan,:) = ridge(combTrainEMGCent(iChan,:)', combTrainBothCent', 0.1);
            decoderWeightsCombStr{iFold}(iChan,:) = ridge(combTrainEMGCent(iChan,:)', combTrainStrCent', 0.1);
            decoderWeightsCombCtx{iFold}(iChan,:) = ridge(combTrainEMGCent(iChan,:)', combTrainCtxCent', 0.1);
        end

        %get polynomial fit
        combEMGEstBoth = combTrainBothCent'*decoderWeightsCombBoth{iFold}';
        combEMGEstStr = combTrainStrCent'*decoderWeightsCombStr{iFold}';
        combEMGEstCtx = combTrainCtxCent'*decoderWeightsCombCtx{iFold}';

        for iChan = 1:size(combEMGEstBoth,2)
            polyFitWeightsCombBoth{iFold}(iChan,:) = polyfit(combEMGEstBoth(:,iChan),combTrainEMGCent(iChan,:),3);
            polyFitWeightsCombStr{iFold}(iChan,:) = polyfit(combEMGEstStr(:,iChan),combTrainEMGCent(iChan,:),3);
            polyFitWeightsCombCtx{iFold}(iChan,:) = polyfit(combEMGEstCtx(:,iChan),combTrainEMGCent(iChan,:),3);
        end
    end

end
clear allNeurHist neurCent

% now do decoding

% go through each region
for iRegion = 1:nRegions

    %go through each CV fold
    estSingleFoldsBoth = {};
    estSingleFoldsStr = {};
    estSingleFoldsCtx = {};
    estTwoFoldsBoth = {};
    estTwoFoldsStr = {};
    estTwoFoldsCtx = {};
    estCombFoldsBoth = {};
    estCombFoldsStr = {};
    estCombFoldsCtx = {};
    for iFold = 1:nFolds

        %get the testing inds
        thisTestInds = foldTestInds{iRegion}{iFold};

        %get testing input and output
        testEMG = regionEMGHistNormalized{iRegion}(:,thisTestInds);
        testNeurBoth = regionBothHist{iRegion}(:,thisTestInds);
        testNeurStr = regionStrHist{iRegion}(:,thisTestInds);
        testNeurCtx = regionCtxHist{iRegion}(:,thisTestInds);

        testEMGCent = testEMG - mean(testEMG,2);
        testNeurBothCent = testNeurBoth - mean(testNeurBoth,2);
        testNeurStrCent = testNeurStr - mean(testNeurStr,2);
        testNeurCtxCent = testNeurCtx - mean(testNeurCtx,2);
        
        %decode
        %first use single trained decoder
        for iChan = 1:size(testEMG,1)
            estSingleFoldsBoth{iFold}(iChan,:) = polyval(polyFitWeightsBoth{iRegion,iFold}(iChan,:),testNeurBothCent'*decoderWeightsBoth{iRegion,iFold}(iChan,:)');
            estSingleFoldsStr{iFold}(iChan,:) = polyval(polyFitWeightsStr{iRegion,iFold}(iChan,:),testNeurStrCent'*decoderWeightsStr{iRegion,iFold}(iChan,:)');
            estSingleFoldsCtx{iFold}(iChan,:) = polyval(polyFitWeightsCtx{iRegion,iFold}(iChan,:),testNeurCtxCent'*decoderWeightsCtx{iRegion,iFold}(iChan,:)');
        end

        if doCrossRegion

            %next use two region trained decoders
            for iRegion2 = 1:nRegions

                if iRegion2 == iRegion
                    continue
                end

                if iRegion2 < iRegion
                    thisPolyWeightsBoth = polyFitWeightsTwoBoth{iRegion,iRegion2,iFold};
                    thisPolyWeightsStr = polyFitWeightsTwoStr{iRegion,iRegion2,iFold};
                    thisPolyWeightsCtx = polyFitWeightsTwoCtx{iRegion,iRegion2,iFold};

                    thisDecoderWeightsBoth = decoderWeightsTwoBoth{iRegion, iRegion2,iFold};
                    thisDecoderWeightsStr = decoderWeightsTwoStr{iRegion, iRegion2,iFold};
                    thisDecoderWeightsCtx = decoderWeightsTwoCtx{iRegion, iRegion2,iFold};

                else
                    thisPolyWeightsBoth = polyFitWeightsTwoBoth{iRegion2,iRegion,iFold};
                    thisPolyWeightsStr = polyFitWeightsTwoStr{iRegion2,iRegion,iFold};
                    thisPolyWeightsCtx = polyFitWeightsTwoCtx{iRegion2,iRegion,iFold};

                    thisDecoderWeightsBoth = decoderWeightsTwoBoth{iRegion2, iRegion,iFold};
                    thisDecoderWeightsStr = decoderWeightsTwoStr{iRegion2, iRegion,iFold};
                    thisDecoderWeightsCtx = decoderWeightsTwoCtx{iRegion2, iRegion,iFold};
                end

                for iChan = 1:size(testEMG,1)
                    estTwoFoldsBoth{iRegion2,iFold}(iChan,:) = polyval(thisPolyWeightsBoth(iChan,:),testNeurBothCent'*thisDecoderWeightsBoth(iChan,:)');
                    estTwoFoldsStr{iRegion2,iFold}(iChan,:) = polyval(thisPolyWeightsStr(iChan,:),testNeurStrCent'*thisDecoderWeightsStr(iChan,:)');
                    estTwoFoldsCtx{iRegion2,iFold}(iChan,:) = polyval(thisPolyWeightsCtx(iChan,:),testNeurCtxCent'*thisDecoderWeightsCtx(iChan,:)');
                end
            end

            %finally, all regions trained decoder
            for iChan = 1:size(testEMG,1)
                estCombFoldsBoth{iFold}(iChan,:) = polyval(polyFitWeightsCombBoth{iFold}(iChan,:),testNeurBothCent'*decoderWeightsCombBoth{iFold}(iChan,:)');
                estCombFoldsStr{iFold}(iChan,:) = polyval(polyFitWeightsCombStr{iFold}(iChan,:),testNeurStrCent'*decoderWeightsCombStr{iFold}(iChan,:)');
                estCombFoldsCtx{iFold}(iChan,:) = polyval(polyFitWeightsCombCtx{iFold}(iChan,:),testNeurCtxCent'*decoderWeightsCombCtx{iFold}(iChan,:)');
            end

        end

    end

    %concatenate all the folds and calc performance metric
    estEMGSingleBoth{iRegion} = cat(2,estSingleFoldsBoth{:});
    estEMGSingleStr{iRegion} = cat(2,estSingleFoldsStr{:});
    estEMGSingleCtx{iRegion} = cat(2,estSingleFoldsCtx{:});

    if doCrossRegion

        for iRegion2 = 1:nRegions
            if iRegion2 == iRegion
                continue
            end
            estEMGTwoBoth{iRegion,iRegion2} = cat(2,estTwoFoldsBoth{iRegion2,:});
            estEMGTwoStr{iRegion,iRegion2} = cat(2,estTwoFoldsStr{iRegion2,:});
            estEMGTwoCtx{iRegion,iRegion2} = cat(2,estTwoFoldsCtx{iRegion2,:});
        end

        estEMGCombBoth{iRegion} = cat(2,estCombFoldsBoth{:});
        estEMGCombStr{iRegion} = cat(2,estCombFoldsStr{:});
        estEMGCombCtx{iRegion} = cat(2,estCombFoldsCtx{:});

    end

    %calc decoding performance
    singlePerformanceBoth{iRegion} = calcPerformanceMetrics(estEMGSingleBoth{iRegion},regionEMGHistNormalized{iRegion});
    singlePerformanceStr{iRegion} = calcPerformanceMetrics(estEMGSingleStr{iRegion},regionEMGHistNormalized{iRegion});
    singlePerformanceCtx{iRegion} = calcPerformanceMetrics(estEMGSingleCtx{iRegion},regionEMGHistNormalized{iRegion});

    if doCrossRegion

        for iRegion2 = 1:nRegions
            if iRegion2 == iRegion
                continue
            end
            twoPerformanceBoth{iRegion,iRegion2} = calcPerformanceMetrics(estEMGTwoBoth{iRegion,iRegion2},regionEMGHistNormalized{iRegion});
            twoPerformanceStr{iRegion,iRegion2} = calcPerformanceMetrics(estEMGTwoStr{iRegion,iRegion2},regionEMGHistNormalized{iRegion});
            twoPerformanceCtx{iRegion,iRegion2} = calcPerformanceMetrics(estEMGTwoCtx{iRegion,iRegion2},regionEMGHistNormalized{iRegion});
        end

        combPerformanceBoth{iRegion} = calcPerformanceMetrics(estEMGCombBoth{iRegion},regionEMGHistNormalized{iRegion});
        combPerformanceStr{iRegion} = calcPerformanceMetrics(estEMGCombStr{iRegion},regionEMGHistNormalized{iRegion});
        combPerformanceCtx{iRegion} = calcPerformanceMetrics(estEMGCombCtx{iRegion},regionEMGHistNormalized{iRegion});

    end

end


%%
% get global performance metrics (i.e. not for individual time points)
figure('Color','w')
tiledlayout(nRegions,1)
for iRegion = 1:nRegions

    nexttile
    %find nan points and remove
    regionNans = unique([regionDecodedNanInds{iRegion} testEMGNanInds{iRegion}]);

    realEMGRegion = testEMGNormalizedRegions{iRegion};
    realEMGRegion(:,regionNans) = [];

    %do for each decoder trained on one region
    for iRegionTrain = 1:nRegions

        estEMGRegion = estEMGRegions{iRegion}{iRegionTrain};
        estEMGRegion(:,regionNans) = [];
        decodingPerformance{iRegion}{iRegionTrain} = calcPerformanceMetrics(estEMGRegion,realEMGRegion);

    end

    % decoding with decoder trained on two regions
    for iTwoRegions = 1:length(twoRegionLabels)

        estEMGRegion = estEMGRegionsTwo{iRegion}{iTwoRegions};
        estEMGRegion(:,regionNans) = [];
        decodingPerformanceTwo{iRegion}{iTwoRegions} = calcPerformanceMetrics(estEMGRegion,realEMGRegion);

    end

    % decoding with decoder trained on all regions
    estEMGRegion = estEMGRegionsComb{iRegion};
    estEMGRegion(:,regionNans) = [];
    decodingPerformanceComb{iRegion} = calcPerformanceMetrics(estEMGRegion,realEMGRegion);

    % plot
    emgChans = 2:7;

    hold on
    for iRegionTrain = 1:nRegions
        plot(repmat(iRegionTrain,length(decodingPerformance{iRegion}{iRegionTrain}.CC(emgChans))),...
            decodingPerformance{iRegion}{iRegionTrain}.CC(emgChans),'.','MarkerSize',10,'Color',lines(1),'MarkerFaceColor','k')
    end
%     for iTwoRegions = 1:length(twoRegionLabels)
%         plot(repmat(iTwoRegions+3,length(decodingPerformanceTwo{iRegion}{iTwoRegions}.CC(emgChans))),...
%             decodingPerformanceTwo{iRegion}{iTwoRegions}.CC(emgChans),'o','MarkerSize',4,'Color','k','MarkerFaceColor','k')
%     end
%     plot(repmat(7,length(decodingPerformanceComb{iRegion}.CC(emgChans))),...
%             decodingPerformanceComb{iRegion}.CC(emgChans),'o','MarkerSize',4,'Color','k','MarkerFaceColor','k')

    plot(1:7,[cellfun(@(x) mean(x.CC),decodingPerformance{iRegion})],'.-','MarkerSize',20,'Color',[0.3 0.3 0.3])
    title(['Decoding Region ' num2str(iRegion)])
    xlim([0 4])
    set(gca,'XTick',1:3)
    set(gca,'fontsize',15)
    set(gca,'linewidth',2)
    set(gca,'fontsize',15)
    set(gca,'XTickLabelRotation',45)
    if iRegion == nRegions
        set(gca,'XTickLabel',{'Climbing Region Trained','Walking Region Trained','Eat/Grooming Region Trained','Region 1+2 Trained',...
            'Region 1+3 Trained','Region 2+3 Trained','All region trained'})
    else
        set(gca,'XTick',[])
    end
    box off
    ylabel('Decoder CC')

end

% calc absolute error at each time point


% calc performance metrics using sliding windows
windowSize = 100;

% save(fullfile(baseDir,'ProcessedData','uMAPDecoding.mat'),'decodingPerformance','decodingPerformanceTwo','decodingPerformanceComb',...
%     'testEMGNormalizedRegions','estEMGRegions','estEMGRegionsTwo','estEMGRegionsComb','regionDecodedNanInds','testEMGNanInds')


function dataHist = addHistorySkip(data, nHist, skip)

    for iHist = 1:nHist
        
        shiftAmount = iHist*skip;
        shiftData{iHist} = [zeros(size(data,1), shiftAmount) data(:, 1:end-shiftAmount)];
        
    end

    dataHist = [data; cat(1,shiftData{:})];

end



function optimalParam = findBestRidgeParam(input,output,parameterVals)
% do 5-fold cross validation to determine the optimal ridge regression
% parameter

nFolds = 5;
foldInds = divideBlocks(1:size(input,2),nFolds);

for iFold = 1:nFolds
    
    %get training set and mean center
    trainInds = setdiff(1:size(input,2),foldInds{iFold});
    foldInputTrain = input(:,trainInds);
    foldInputTrainCent = foldInputTrain - mean(foldInputTrain,2);
    foldOutputTrain = output(:,trainInds);
    foldOutputTrainCent = foldOutputTrain - mean(foldOutputTrain,2);

    %get testing set and mean center
    foldInputTest = input(:,foldInds{iFold});
    foldInputTestCent = foldInputTest - mean(foldInputTest,2);
    foldOutputTest = output(:,foldInds{iFold});
    foldOutputTestCent{iFold} = foldOutputTest - mean(foldOutputTest,2);

    %test all ridge values
    for iChan = 1:size(foldOutputTest,1)
        
        ridgeParams = ridge(foldOutputTrainCent(iChan,:)', foldInputTrainCent', parameterVals);
        ridgeEst = foldInputTestCent'*ridgeParams;

        for iParam = 1:length(parameterVals)
            outputEst{iParam}{iFold}(iChan,:) = ridgeEst(:,iParam);
        end

    end

end

% concatenate all folds and then calculate performance metrics
for iParam = 1:length(parameterVals)

    allFoldOutputEst{iParam} = cat(2,outputEst{iParam}{:});
    allFoldOutputTest = cat(2,foldOutputTestCent{:});
    decodeMetrics{iParam} = calcPerformanceMetrics(allFoldOutputEst{iParam},allFoldOutputTest);

end


end


% 
