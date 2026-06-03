function allChansFiltData = activityUMAPOverlay(sessionName, sessionProbe, frFile, activityType, varargin)
% activityUMAPOverlay(baseDir, frFile, activityType, saveResultsDir, nShuff)
% 
% Overlay recorded time-series into UMAP space by binning into grids and
% taking average within each grid bin. Also does some smoothing by
% convolving with a gaussian. For firing rate time-series, assumes it's in
% 10ms time bins.
% 
% Inputs:
% baseDir           - Folder containing the data for the recording session
% 
% frFile            - File contanin the firing rate data
% 
% activityType      - String specifying the type of activity to overlay on
%                     the UMAP space. Can be either 'FiringRate','SSA', or
%                     'EMG'
% 
% saveResultsDir    - Optional input, string specifying the folder to save
%                     the overlay plots and results to. If empty, will not 
%                     save any plots. Default empty.
% 
% nShuff            - Optional input, number of time-shift shuffles to do
%                     for permutation tests. Default 200.
% 
% doSubsamp         - Optional input, boolean specifying whether we should
%                     subsample number of time points in each behavior
%                     regions so that the same number of time points are
%                     used across all regions. Default false.
% 
% nResamp           - Optional input, number of bootstrap resamples to do
%                     on the time points to get variace for average 
%                     in behavior clusters. Default 500.
% 
% David Xing, last updated 6/10/2025

% parse optional input parameters
narginchk(4, 8)

if length(varargin)>=1
    if isempty(varargin{1})
        saveResultsDir = [];
    else
        saveResultsDir = varargin{1};
    end
else
    saveResultsDir = [];
end

if length(varargin)>=2
    if isempty(varargin{2})
        nShuff = 200;
    else
        nShuff = varargin{2};
    end
else
    nShuff = 200;
end

if length(varargin)>=3
    if isempty(varargin{3})
        doSubsamp = false;
    else
        doSubsamp = varargin{3};
    end
else
    doSubsamp = false;
end


if length(varargin)>=4
    if isempty(varargin{4})
        nResamp = 500;
    else
        nResamp = varargin{4};
    end
else
    nResamp = 500;
end

rng(2025)
useHumanAnnoBouts = false;

filePaths = getMouseDataNames(sessionName(1:4),sessionName,sessionProbe);

% load UMAP projection
load(filePaths.UMAPFile,'reduction','origDownsampEMGInd','gridXInds','gridYInds','umapArtifactInds','classifierLabels',...
    'musclePercentileGroupLabels','watershedRegions','regionWatershedLabels','regionAssignmentsFiltered','behvLabelsNoArt','regionBehvAssignments')

if ~exist('umapArtifactInds')
    umapArtifactInds = [];
end

% for sparsity calculations, use 100x100 grid rather than 1000x1000 grid
gridXIndsDownsamp = linspace(gridXInds(1),gridXInds(end),100);
gridYIndsDownsamp = linspace(gridYInds(1),gridYInds(end),100);

% load in sync data
load(filePaths.VideoSyncFrames)

% just do it for every time point in reduction
% sync UMAP to neural timepoints
currentDir = pwd;
cd(filePaths.processedDataFolder)
reducNeurInds = round(NeurEMGSync(origDownsampEMGInd*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
cd(currentDir);
maxNeurSamples = floor(frameNeuropixelSamples{1}{end}(end)/300);

% don't use any points past the end of video
pastEndPoints = find(reducNeurInds > maxNeurSamples);
reducNeurInds(pastEndPoints) = [];

% don't use artifact points 
reducIndsToUse = 1:length(reducNeurInds);

sessionArtifacts = consolidateArtifactInds(sessionName, sessionProbe);
artInds = sessionArtifacts.allArtNeurInds10ms30msSmooth;
artInds(artInds > maxNeurSamples) = [];

removePoints = find(ismember(reducNeurInds,artInds));
reducIndsToUse(removePoints) = [];

if ~isdir(fullfile(filePaths.processedDataFolder,saveResultsDir))
    mkdir(fullfile(filePaths.processedDataFolder,saveResultsDir))
end

switch lower(activityType)

    case 'firingrate'
        % load firing rate data
        load(fullfile(filePaths.processedDataFolder,frFile),'allFRs','striatumInds','cortexInds')
        
        % get neural data corresponding to the UMAP time points
        reducFRs = allFRs(:,reducNeurInds(reducIndsToUse));
        clear allFRs

        % don't use any points where there's nans in the FR data
        nanReducPoints = any(isnan(reducFRs),1);
        reducIndsToUse(nanReducPoints) = [];

        reducSigs = reducFRs(:,~nanReducPoints);
        clear reducFRs

    case {'pcaall','pcacortex','pcastriatum'}
        % load firing rate data
        load(fullfile(filePaths.processedDataFolder,frFile),'allFRs','striatumInds','cortexInds')

        %run PCA
        
        switch lower(activityType)
            case 'pcaall'
                pcaGoodNeurons = goodNeuronsAll;
                [pcaWeights, pcaTrajs] = pca(allFRs(pcaGoodNeurons,:)');

            case 'pcacortex'
                pcaGoodNeurons = intersect(goodNeuronsAll,cortexInds);
                [pcaWeights, pcaTrajs] = pca(allFRs(pcaGoodNeurons,:)');

            case 'pcastriatum'
                pcaGoodNeurons = intersect(goodNeuronsAll,striatumInds);
                [pcaWeights, pcaTrajs] = pca(allFRs(pcaGoodNeurons,:)');
        end

        % get the data corresponding to the UMAP time points
        reducPCATrajs = pcaTrajs(reducNeurInds(reducIndsToUse),:)';
        clear pcaTrajs

        % don't use any points where there's nans in the data
        nanReducPoints = any(isnan(reducPCATrajs),1);
        reducIndsToUse(nanReducPoints) = [];

        reducSigs = reducPCATrajs(:,~nanReducPoints);
        clear reducPCATrajs allFRs

    case 'emg'
        % load EMG data
        load(filePaths.EMG1ms)

        %get meta data as well
        load(filePaths.emgMetaData)

%         % Z-score
%         zScoredEMG = (downsampEMG-mean(downsampEMG,2))./std(downsampEMG,[],2);

        % get emg corresponding to the UMAP
        reducEMGInds = origDownsampEMGInd;
        reducSigs = downsampEMG(:,reducEMGInds(reducIndsToUse));

    case {'ssaall','ssacortex','ssastriatum'}
        %load SSA projections
        load(fullfile(filePaths.processedDataFolder, 'SSA','allTimePointsSSA.mat'))
        load(fullfile(filePaths.processedDataFolder,'NeuralFiringRates10msBins30msGauss.mat'),'allFRs');
        nanInds = find(any(isnan(allFRs)));
        
        %just put in the nans again, then index into the SSA trajs, then
        %remove nans, probably easiest way to make it consistent with the
        %neural firing rate inputs
        ssaTrajWithNans = zeros(nSSADims,size(allFRs,2));
        ssaTrajWithNans(:,nanInds) = nan;

        switch lower(activityType)
            case 'ssastriatum'
                ssaTrajWithNans(:,setdiff(1:size(allFRs,2),nanInds)) = ssaResults{1}.trajs';
            case 'ssacortex'
                ssaTrajWithNans(:,setdiff(1:size(allFRs,2),nanInds)) = ssaResults{2}.trajs';
            case 'ssaall'
                ssaTrajWithNans(:,setdiff(1:size(allFRs,2),nanInds)) = ssaResults{3}.trajs';
        end

        reducSSAs = ssaTrajWithNans(:,reducNeurInds(reducIndsToUse));
        nanReducPoints = any(isnan(reducSSAs),1);
        reducIndsToUse(nanReducPoints) = [];
        reducSigs = reducSSAs(:,~nanReducPoints);
        clear reducSSAs
        clear allFRs
end

% get watershedding region splitting
[regionBoundaryIndsX, regionBoundaryIndsY] = find(watershedRegions==0);

% Generate grid
nGridPoints = length(gridXInds);
nGridPointsDownsamp = length(gridXIndsDownsamp);
[meshGridX,meshGridY] = meshgrid(gridXInds,gridYInds);

% get the gaussian kernal to convolve with for better visualization
densityGaussStd = 0.3;
gaussKernal = exp(-.5.*((meshGridX-mean([max(max(meshGridX)) min(min(meshGridX))])).^2 + (meshGridY-mean([max(max(meshGridY)) min(min(meshGridY))])).^2)./densityGaussStd^2) ./ (2*pi*densityGaussStd^2);
gaussKernal = gaussKernal./(sum(sum(gaussKernal)));

% for each point in the UMAP reduction, get the bin which it belongs to when doing 2D binning using the grid points 
[~,~,~,binIndX, binIndY] = histcounts2(reduction(reducIndsToUse,1),reduction(reducIndsToUse,2),gridXInds,gridYInds);
binLabels = sub2ind([nGridPoints, nGridPoints], binIndX, binIndY);

% also do it for the coarser grid points for sparsity calculation
[~,~,~,binIndXDownsamp, binIndYDownsamp] = histcounts2(reduction(reducIndsToUse,1),reduction(reducIndsToUse,2),gridXIndsDownsamp,gridYIndsDownsamp);
binLabelsDownsamp = sub2ind([nGridPointsDownsamp, nGridPointsDownsamp], binIndXDownsamp, binIndYDownsamp);


% for each 2D bin, get how many time points are inside it (i.e. how many of
% the umap reduction points are in each bin)
binNPoints = accumarray(binLabels, ones(1,length(binLabels)));
binNPointsDownsamp = accumarray(binLabelsDownsamp, ones(1,length(binLabelsDownsamp)));

% also, sometimes the bin has no reduction points, in which case, don't use
% it for firing rate calculations (e.g. it's a on the corner of the grid)
accumarrayRealOutputs = find(binNPoints~=0);
accumarrayRealOutputsDownsamp = find(binNPointsDownsamp~=0);

% Changed code to do signal calculations not just on UMAP watershed
% regions, but also on the 10 behaviors based on the supervised classifier
allUsedLabels = {regionAssignmentsFiltered, classifierLabels, musclePercentileGroupLabels};

if iscell(regionWatershedLabels)
    regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
else
    regionWatershedLabelsCat = regionWatershedLabels;
end

allUsedLabelIds = {regionWatershedLabelsCat, 1:10, 1:8};

% loop through each behavior label type
behvAveSigs = [];
behvPrctSigs = [];
behvAveSigsResamp = [];
behvPrctSigsResamp = [];
behvAveSigsShuff = [];
behvPrctSigsShuff = [];

for iLabelType = 1:length(allUsedLabels)

    usedLabels = allUsedLabels{iLabelType};

    % now, if we want to do permutation test shuffles, loop through all the
    % subsequent calculations here
    nTotalInds = size(reducSigs,2);
    for iShuff = 1:nShuff + 1

        tic
        if iShuff == 1
            %first loop, just do non-shuffled
            itrReducSigs = reducSigs;
        else
            %do perm shuffles
            randShift(iShuff-1) = randperm(size(reducSigs,2),1);

            %shift should be at least 30 seconds away from original time point
            while randShift(iShuff-1) < 30000 && randShift(iShuff-1) > size(reducSigs,2)-30000
                randShift(iShuff-1) = randperm(size(reducSigs,2),1);
            end

            itrReducSigs = circshift(reducSigs,randShift(iShuff-1),2);
            reducSigsShuff = circshift(reducSigs,randShift(iShuff-1),2);

            % instead of doing time shifting, actually permute the indices for
            % calculating controls for sparsity (since shifting relative to
            % behavior doesn't serve as a control as behavior isn't used in
            % the sparsity calculation)
            reducSigsPerm = reducSigs(:,randperm(size(reducSigs,2)));
        end

        % go through each neuron and plot the average firing rates in the space
        for iChan = 1:size(itrReducSigs,1)

            % For the UMAP overlay stuff, only do it for the UMAP labels,
            % don't need to re-run this part for the classifier labels
            if iLabelType == 1
                % get average firing rate in each spatial bin
                binSigSums = accumarray(binLabels,itrReducSigs(iChan,:));
                binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*100;

                % also use coarser bins for sparsity calculation
                binSigSums = accumarray(binLabelsDownsamp,itrReducSigs(iChan,:));
                binAveSigDownsamp = binSigSums(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;

                spaceAveZScores = zeros(nGridPoints,nGridPoints);
                binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
                binAveZScores = binAveSig;

                %map the firing rates into the actual bins within the grid now (using
                %the bin labels as the index). Only the bins with anything in them are
                %mapped
                spaceAveZScores(unique(binLabels)) = binAveZScores;

                spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));
                allChansFiltData(:,:,iChan) = spaceAveZScoresFilt;

                if iShuff == 1
                    %don't need to plot perm shuffles
                    figH = figure('Visible','on','Units','pixels','OuterPosition',[200 50 500 1000]);
                    tiledlayout(2,1,"TileSpacing","tight")
                    nexttile
                    imagesc(gridXInds, gridYInds, spaceAveZScoresFilt')
                    hold on
                    plot((gridXInds(regionBoundaryIndsY)),gridYInds(regionBoundaryIndsX),'k.')
                    colorbar
                    set(gca,'YDir','normal')

                end

            end

            %now the average activity within each of the divided regions
            %     watershedRegionsAdj = watershedRegions';
            %     watershedRegionsBins = watershedRegionsAdj(:);
            %     reducUsedLabels = watershedRegionsBins(binLabels);
            reducUsedLabels = usedLabels(reducIndsToUse);
            usedLabelIds = allUsedLabelIds{iLabelType};

            for iBehv = 1:length(usedLabelIds)
                regionNPoints(iBehv) = sum(reducUsedLabels == usedLabelIds(iBehv));
            end
            regionMinPoints = min(regionNPoints);

            % use human annotated labels to get bouts for estimating
            % intra-behavior variance
            if useHumanAnnoBouts

                for iBehv = 1:length(regionBehvAssignments)

                    regionBehvs = regionBehvAssignments{iBehv};
                    thisBehvBoutMeans = {};
                    shedRegionNBouts = [];
                    for iBehv = 1:length(regionBehvs)

                        behvInds = find(behvLabelsNoArt(reducIndsToUse)==regionBehvs(iBehv));
                        boutStarts = [1 find(diff(behvInds) > 1)+1];
                        boutEnds = [find(diff(behvInds) > 1) length(behvInds)];
                        shedRegionNBouts(iBehv) = length(boutStarts);

                        for iBout = 1:length(boutStarts)
                            boutSigs = itrReducSigs(iChan,behvInds(boutStarts(iBout)):behvInds(boutEnds(iBout)));
                            thisBehvBoutMeans{iBehv}(iBout) = nanmean(boutSigs)*1000;
                        end


                    end

                    regionAnnoBoutAveSigs{iChan,iBehv} = cat(2,thisBehvBoutMeans{:});
                    regionAnnoBoutLabels{iBehv} = repmat(iBehv,1,sum(shedRegionNBouts));

                end

            end

            for iBehv = 1:length(usedLabelIds)
                
                itrReducInds = find(reducUsedLabels == usedLabelIds(iBehv));

                % get bouts
                boutStarts = itrReducInds([1 find(diff(itrReducInds) > 1)+1]);
                boutEnds = itrReducInds([find(diff(itrReducInds) > 1) length(itrReducInds)]);
                boutDurations = boutEnds - boutStarts;

                % for resample do resampling of continuous blocks
                bootStrapBlockDur = 500;
                doBlockResamp = true;
                behvBoutRegionMeans = {};
                for iBout = 1:length(boutStarts)
                    if boutDurations(iBout) < bootStrapBlockDur
                        continue
                    end
                    blockStartInds = boutStarts(iBout):bootStrapBlockDur:boutEnds(iBout);

                    for iBlock = 1:length(blockStartInds)
                        if blockStartInds(iBlock)+bootStrapBlockDur-1 > size(itrReducSigs,2)
                            continue
                        end
                        behvBoutRegionMeans{iBout}(iBlock) = nanmean(itrReducSigs(iChan,blockStartInds(iBlock):blockStartInds(iBlock)+bootStrapBlockDur-1));
                    end
                end

                thisBehvRegionBlockMeans{iBehv} = cat(2,behvBoutRegionMeans{:});

                % next do averages within bouts, but only use "good" bouts that
                % satisfy a set of criterion
                boutStarts = boutStarts(boutDurations > 500);
                boutEnds = boutEnds(boutDurations > 500);

                boutEnds = boutEnds(diff(boutStarts)>1000);
                boutStarts = boutStarts(diff(boutStarts)>1000);
                behvBoutDurs(iBehv) = length(boutStarts);
                behvBoutLabels{iBehv} = repmat(iBehv,1,length(boutStarts));

                for iBout = 1:length(boutStarts)

                    boutSigs = itrReducSigs(iChan,boutStarts(iBout):boutEnds(iBout));
                    behvBoutAveSigs{iChan,iBehv}(iBout) = nanmean(boutSigs)*1000;

                end

            end % of behaviors loop for calculating bouts

            regionMinBlocks = min(cellfun(@length,thisBehvRegionBlockMeans));

            for iBehv = 1:length(usedLabelIds)
                %if we want to downsample to same number of points across behavior
                %regions, do it now
                itrReducInds = find(reducUsedLabels == usedLabelIds(iBehv));
                if doSubsamp
                    randDownsampInds = randperm(length(itrReducInds),regionMinPoints);
                    itrReducInds = itrReducInds(randDownsampInds);
                    thisBehvRegionBlockMeans{iBehv} = thisBehvRegionBlockMeans{iBehv}(...
                        randperm(length(thisBehvRegionBlockMeans{iBehv}),regionMinBlocks));
                end

                if iShuff == 1
                    %save the indices that we ended up using
                    if iLabelType == 1
                        umapRegionReducInds{iBehv} = itrReducInds;
                    elseif iLabelType == 2
                        classifierReducInds{iBehv} = itrReducInds;
                    elseif iLabelType == 3
                        emgPercentileReducInds{iBehv} = itrReducInds;
                    end
                end

                % now loop through and do bootstrap resamples if we want to
                % (only do it for the actual data, not any permutation
                % shuffles)
                if iShuff == 1
                    nExtraLoops = nResamp;
                else
                    nExtraLoops = 0;
                end

                for iResamp = 1:nExtraLoops+1

                    if iResamp == 1
                        %first loop is just non-resampled indices
                        itrUsedInds = itrReducInds;
                    elseif ~doBlockResamp
                        %otherwise do resampling same number of points with replacement
                        itrUsedInds = itrReducInds(randsample(length(itrReducInds),length(itrReducInds),true));
                    else
                        %do resampling of the blocks
                        resampMeans = thisBehvRegionBlockMeans{iBehv}(randsample(length(thisBehvRegionBlockMeans{iBehv}),length(thisBehvRegionBlockMeans{iBehv}),true));
                        behvAveSigsResamp(iChan,iBehv,iResamp-1) = mean(resampMeans)*1000;
                        continue
                    end

                    behvReducSigs = itrReducSigs(iChan,itrUsedInds);
                    behvReducSigsNonZeros = behvReducSigs(behvReducSigs ~= 0);
                    aveSigs = mean(behvReducSigs)*1000;
                    prctSigs = prctile(behvReducSigsNonZeros,90)*1000;

                    %also compute rates not by averging just the raw firing rates,
                    %but by taking the filtered values from the map
                    % this part only do for umap watershed labels
                    if iLabelType == 1
                        regionFiltValues = spaceAveZScoresFilt(binLabels(itrUsedInds));
                        regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
                        aveSigsFilt = mean(regionFiltValues);
                        prctSigsFilt = prctile(regionFiltValuesNonZeros,99);
                    end

                    % save values to arrays across channels and regions
                    if iResamp == 1 && iShuff == 1
                        behvAveSigs(iChan,iBehv) = aveSigs;
                        behvPrctSigs(iChan,iBehv) = prctSigs;

                        if iLabelType == 1
                            regionAveSigsFilt(iChan,iBehv) = aveSigsFilt;
                            regionPrctSigsFilt(iChan,iBehv) = prctSigsFilt;
                        end
                    elseif iShuff == 1
                        behvAveSigsResamp(iChan,iBehv,iResamp-1) = aveSigs;
                        behvPrctSigsResamp(iChan,iBehv,iResamp-1) = prctSigs;

                        if iLabelType == 1
                            regionAveSigsFiltResamp(iChan,iBehv,iResamp-1) = aveSigsFilt;
                            regionPrctSigsFiltResamp(iChan,iBehv,iResamp-1) = prctSigsFilt;
                        end
                    else
                        behvAveSigsShuff(iChan,iBehv,iShuff-1) = aveSigs;
                        behvPrctSigsShuff(iChan,iBehv,iShuff-1) = prctSigs;

                        if iLabelType == 1
                            regionAveSigsFiltShuff(iChan,iBehv,iShuff-1) = aveSigsFilt;
                            regionPrctSigsFiltShuff(iChan,iBehv,iShuff-1) = prctSigsFilt;
                        end
                    end

                end %resample loop

            end %behavior region loop

            %calculate anova eta^2
            %         [pVal(iChan),tbl] = anova1(cat(2,behvBoutAveSigs{iChan,:}),cat(2,regionBoutLabels{:}),'off');
            %         etaSquared(iChan) = tbl{2,2}/tbl{4,2};
            %
            %         [~,tbl] = anova1(cat(2,regionAnnoBoutAveSigs{iChan,:}),cat(2,regionAnnoBoutLabels{:}),'off');
            %         etaSquaredAnno(iChan) = tbl{2,2}/tbl{4,2};
            %
            %         [~,tbl] = anova1(squeeze(behvAveSigsResamp(iChan,:,:))',[],'off');
            %         etaSquaredResamp(iChan) = tbl{2,2}/tbl{4,2};

            %calculate sparsity and spike information metrics from Jung et al 1994,
            %Skaggs et al, 1996
            %only use bins that have at least 100 points
            if iLabelType == 1
                hasMinPointsBins = find(binNPointsDownsamp(accumarrayRealOutputsDownsamp) >= 100);
                probVisitBin = binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins))/...
                    sum(binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins)));

                if iShuff == 1
                    sparsity(iChan) = sum(probVisitBin.*binAveSigDownsamp(hasMinPointsBins))^2 / ...
                        sum(probVisitBin.*binAveSigDownsamp(hasMinPointsBins).^2);
                else
                    % instead of doing time shifting, actually permute the indices for
                    % calculating controls for sparsity (since shifting relative to
                    % behavior doesn't serve as a control as behavior isn't used in
                    % the sparsity calculation)
                    binSigSumsPerm = accumarray(binLabelsDownsamp,reducSigsPerm(iChan,:));
                    binAveSigPerm = binSigSumsPerm(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;
                    sparsityShuff(iChan,iShuff-1) = sum(probVisitBin.*binAveSigPerm(hasMinPointsBins))^2 / ...
                        sum(probVisitBin.*binAveSigPerm(hasMinPointsBins).^2);
                end

                % add to plots
                if iShuff == 1

                    % naming scheme
                    switch lower(activityType)

                        case 'firingrate'
                            if iChan > length(striatumInds)
                                region = 'Cortex';
                                neuronNum = iChan - length(striatumInds);
                            else
                                region = 'Striatum';
                                neuronNum = iChan;
                            end
                            title([region ', Neuron ' num2str(neuronNum)])
                            figFilename = [region '_Neuron' num2str(neuronNum) '.png'];

                        case 'emg'
                            title(channelNames{iChan})
                            figFilename = ['Muscle_' replace(channelNames{iChan},' ','-') '.png'];

                        case 'ssaall'
                            title(['SSA (all neurons) dim ' num2str(iChan)])
                            figFilename = ['SSAAll_Dim' num2str(iChan) '.png'];
                        case 'ssastriatum'
                            title(['SSA (striatum neurons) dim ' num2str(iChan)])
                            figFilename = ['SSAStr_Dim' num2str(iChan) '.png'];
                        case 'ssacortex'
                            title(['SSA (cortex neurons) dim ' num2str(iChan)])
                            figFilename = ['SSACtx_Dim' num2str(iChan) '.png'];
                        case 'pcaall'
                            title(['PCA (all neurons) dim ' num2str(iChan)])
                            figFilename = ['PCAAll_Dim' num2str(iChan) '.png'];
                        case 'pcastriatum'
                            title(['PCA (striatum neurons) dim ' num2str(iChan)])
                            figFilename = ['PCAStr_Dim' num2str(iChan) '.png'];
                        case 'pcacortex'
                            title(['PCA (cortex neurons) dim ' num2str(iChan)])
                            figFilename = ['PCACtx_Dim' num2str(iChan) '.png'];

                    end

                    nexttile
                    plot(regionAveSigsFilt(iChan,:))
                    hold on
                    plot(behvAveSigs(iChan,:))
                    plot(regionPrctSigsFilt(iChan,:))
                    %             title(['Sparsity = ' num2str(sparsity(iChan)) ', Bout EtaSqr = ' num2str(etaSquaredAnno(iChan)) ', Bootstrap EtaSqr = ' num2str(etaSquaredResamp(iChan))])
                    legend('Filtered Average','Value Averaged','Filtered 99 Percentile','box','off')

                    % save plots if desired
                    if ~isempty(saveResultsDir)
                        saveas(figH,fullfile(filePaths.processedDataFolder,saveResultsDir,figFilename))
                    end
                    close(figH)
                end % of adding to plots if block

            end %of label type if block for sparsity calc and adding to plot

        end %activity signal channel loop

        disp(toc)

    end %perm shuffle loop

    % save umap region vs classifier behvs to different variable names
    if iLabelType == 1
        regionAveSigs = behvAveSigs;
        regionPrctSigs = behvPrctSigs;
        regionAveSigsResamp = behvAveSigsResamp;
        regionPrctSigsResamp = behvPrctSigsResamp;
        regionAveSigsShuff = behvAveSigsShuff;
        regionPrctSigsShuff = behvPrctSigsShuff;
    elseif iLabelType == 2
        classifierAveSigs = behvAveSigs;
        classifierPrctSigs = behvPrctSigs;
        classifierAveSigsResamp = behvAveSigsResamp;
        classifierPrctSigsResamp = behvPrctSigsResamp;
        classifierAveSigsShuff = behvAveSigsShuff;
        classifierPrctSigsShuff = behvPrctSigsShuff;
    elseif iLabelType == 3
        emgPercentileAveSigs = behvAveSigs;
        emgPercentilePrctSigs = behvPrctSigs;
        emgPercentileAveSigsResamp = behvAveSigsResamp;
        emgPercentilePrctSigsResamp = behvPrctSigsResamp;
        emgPercentileAveSigsShuff = behvAveSigsShuff;
        emgPercentilePrctSigsShuff = behvPrctSigsShuff;
    end

end % of loop across behavior label types

% % % also do permutation shuffles and do the same calculations
% % for iShuff = 1:nShuff
% %     
% % %     indShuff = randperm(nTotalInds);
% % %     reducFRsNoNansShuff = reducFRsNoNans(:,indShuff);
% % 
% %     tic
% %     for iChan = 1:size(reducSigs,1)
% % 
% %         binSigSums = accumarray(binLabels,reducSigsShuff(iChan,:));
% %         binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*100;
% % 
% %         binSigSums = accumarray(binLabelsDownsamp,reducSigsShuff(iChan,:));
% %         binAveSigDownsamp = binSigSums(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;
% % 
% %         binSigSumsPerm = accumarray(binLabelsDownsamp,reducSigsPerm(iChan,:));
% %         binAveSigPerm = binSigSumsPerm(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;
% % 
% %         spaceAveZScores = zeros(nGridPoints,nGridPoints);
% %         binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
% %         binAveZScores = binAveSig;
% % 
% %         spaceAveZScores(unique(binLabels)) = binAveZScores;
% % 
% %         spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));
% % 
% %         regionWatershedLabelsCat = regionWatershedLabels;
% %         for iShedRegion = 1:length(regionWatershedLabelsCat)
% %             shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
% %             shedReducSigs = reducSigsShuff(iChan,shedReducInds);
% %             regionAveSigsShuff(iChan,iShedRegion,iShuff) = mean(shedReducSigs)*1000;
% %             shedReducSigsNonZeros = shedReducSigs(shedReducSigs ~= 0);
% %             regionPrctSigsShuff(iChan,iShedRegion,iShuff) = prctile(shedReducSigsNonZeros,90)*1000;
% % 
% %             %get the bins at which the points are in
% %             regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
% % 
% %             %and average across that instead
% %             regionAveSigsFiltShuff(iChan,iShedRegion,iShuff) = mean(regionFiltValues);
% %             regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
% %             regionPrctSigsFiltShuff(iChan,iShedRegion,iShuff) = prctile(regionFiltValuesNonZeros,99);
% %         end
% % 
% %         %sparsity
% %         hasMinPointsBins = find(binNPointsDownsamp(accumarrayRealOutputsDownsamp) >= 100);
% %         probVisitBin = binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins))/...
% %             sum(binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins)));
% %         sparsityShuff(iChan,iShuff) = sum(probVisitBin.*binAveSigPerm(hasMinPointsBins))^2 / ...
% %         sum(probVisitBin.*binAveSigPerm(hasMinPointsBins).^2);
% % 
% %     end
% %     disp(toc)
% % 
% % end

commonSaveVars = {'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt','umapRegionReducInds','classifierReducInds'...
    'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff',...
    'classifierAveSigs','classifierPrctSigs','classifierAveSigsShuff','classifierPrctSigsShuff','emgPercentileReducInds',...
    'emgPercentileAveSigs','emgPercentilePrctSigs','emgPercentileAveSigsShuff','emgPercentilePrctSigsShuff'};

switch lower(activityType)
    case 'firingrate'
        outputFileName = 'NeuronRegionProps';
        extraSaveVars = {};
        
    case 'emg'
        outputFileName = 'MuscleRegionProps';
        extraSaveVars = {};

    case 'ssaall'
        outputFileName = 'SSAAllRegionProps';
        extraSaveVars = {};

    case 'ssastriatum'
        outputFileName = 'SSAStrRegionProps';
        extraSaveVars = {};

    case 'ssacortex'
        outputFileName = 'SSACtxRegionProps';
        extraSaveVars = {};

    case 'pcaall'
        outputFileName = 'PCAAllRegionProps';
        extraSaveVars = {'pcaWeights','pcaGoodNeurons'};

    case 'pcastriatum'
        outputFileName = 'PCAStrRegionProps';
        extraSaveVars = {'pcaWeights','pcaGoodNeurons'};

    case 'pcacortex'
        outputFileName = 'PCACtxRegionProps';
        extraSaveVars = {'pcaWeights','pcaGoodNeurons'};

end


if ~isempty(saveResultsDir)
    save(fullfile(filePaths.processedDataFolder,saveResultsDir,outputFileName),commonSaveVars{:},extraSaveVars{:});
end


% 
