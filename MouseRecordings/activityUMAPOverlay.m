function allChansFiltData = activityUMAPOverlay(baseDir, frFile, activityType, varargin)
% activityUMAPOverlay(baseDir, frFile, activityType, saveResultsDir, nShuff)
% 
% Function to cluster UMap embedding space into distinct regions using
% watershed and manual assignment of clusters. Will also do some cleaning
% up of the time-series of the assigned regions. Will also save the region
% assignments to the same file that contains the umap results.
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
% David Xing, last updated 2/16/2024

% parse optional input parameters
narginchk(3, 5)

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

% load UMAP projection
load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','watershedRegions','regionWatershedLabels')

% for sparsity calculations, use 100x100 grid rather than 1000x1000 grid
gridXIndsDownsamp = linspace(gridXInds(1),gridXInds(end),100);
gridYIndsDownsamp = linspace(gridYInds(1),gridYInds(end),100);

% load in sync data
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

% just do it for every time point in reduction
% sync UMAP to neural timepoints
currentDir = pwd;
cd(fullfile(baseDir,'ProcessedData'))
reducNeurInds = round(NeurEMGSync(origDownsampEMGInd*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
cd(currentDir);
maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/300);
reducIndsToUse = find(reducNeurInds < maxNeurSamples);

switch lower(activityType)

    case 'firingrate'
        % load firing rate data
        load(fullfile(baseDir,'ProcessedData',frFile),'allFRs','striatumInds','cortexInds')
        
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
        load(fullfile(baseDir,'ProcessedData',frFile),'allFRs','striatumInds','cortexInds')

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
        load(fullfile(baseDir,'ProcessedData', 'EMG1ms.mat'))

        %get meta data as well
        sessionFiles = string(ls(fullfile(baseDir,'ProcessedData')));
        metaDataFile = sessionFiles(contains(sessionFiles,'ProcessedEMG_MetaData'));

        if isempty(metaDataFile)
            error('Unable to find EMG Meta data file!')
        end
        load(fullfile(baseDir,'ProcessedData', metaDataFile))

%         % Z-score
%         zScoredEMG = (downsampEMG-mean(downsampEMG,2))./std(downsampEMG,[],2);

        % get emg corresponding to the UMAP
        reducEMGInds = origDownsampEMGInd;
        reducSigs = downsampEMG(:,reducEMGInds(reducIndsToUse));

    case {'ssaall','ssacortex','ssastriatum'}
        %load SSA projections
        load(fullfile(baseDir,'ProcessedData', 'SSA','allTimePointsSSA.mat'))
        load(fullfile(baseDir, 'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs');
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

% go through each neuron and plot the average firing rates in the space
for iChan = 1:size(reducSigs,1)

    % get average firing rate in each spatial bin
    binSigSums = accumarray(binLabels,reducSigs(iChan,:));
    binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*100;

    % also use coarser bins for sparsity calculation
    binSigSums = accumarray(binLabelsDownsamp,reducSigs(iChan,:));
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

    figH = figure('Visible','on','Units','pixels','OuterPosition',[200 50 500 1000]);
    tiledlayout(2,1,"TileSpacing","tight")
    nexttile
    imagesc(gridXInds, gridYInds, spaceAveZScoresFilt')
    hold on
    plot((gridXInds(regionBoundaryIndsY)),gridYInds(regionBoundaryIndsX),'k.')
    colorbar
    set(gca,'YDir','normal')

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

    %now the average activity within each of the divided regions
    watershedRegionsAdj = watershedRegions';
    watershedRegionsBins = watershedRegionsAdj(:);
    reducWatershedRegions = watershedRegionsBins(binLabels);

    if iscell(regionWatershedLabels)
        regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
    else
        regionWatershedLabelsCat = regionWatershedLabels;
    end

    for iShedRegion = 1:length(regionWatershedLabelsCat)
        shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
        shedReducSigs = reducSigs(iChan,shedReducInds);
        regionAveSigs(iChan,iShedRegion) = mean(shedReducSigs)*1000;
        shedReducSigsNonZeros = shedReducSigs(shedReducSigs ~= 0);
        regionPrctSigs(iChan,iShedRegion) = prctile(shedReducSigsNonZeros,90)*1000;

        %get the bins at which the points are in
        regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
        %and average across that instead
        regionAveSigsFilt(iChan,iShedRegion) = mean(regionFiltValues);
        regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
        regionPrctSigsFilt(iChan,iShedRegion) = prctile(regionFiltValuesNonZeros,99);
    end

    %calculate sparsity and spike information metrics from Jung et al 1994,
    %Skaggs et al, 1996
    %only use bins that have at least 100 points
    hasMinPointsBins = find(binNPointsDownsamp(accumarrayRealOutputsDownsamp) >= 100);
    probVisitBin = binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins))/...
        sum(binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins)));
    sparsity(iChan) = sum(probVisitBin.*binAveSigDownsamp(hasMinPointsBins))^2 / ...
        sum(probVisitBin.*binAveSigDownsamp(hasMinPointsBins).^2);

    nexttile
    plot(regionAveSigsFilt(iChan,:))
    hold on
    plot(regionAveSigs(iChan,:))
    plot(regionPrctSigsFilt(iChan,:))
    legend('Filtered Average','Value Averaged','Filtered 99 Percentile','box','off')

    if ~isempty(saveResultsDir)
        saveas(figH,fullfile(baseDir,'ProcessedData',saveResultsDir,figFilename))
    end
    close(figH)

end

% also do permutation test and do the same calculations
% do permutation shuffles for ssa and neurons
nTotalInds = size(reducSigs,2);
for iShuff = 1:nShuff
    
%     indShuff = randperm(nTotalInds);
%     reducFRsNoNansShuff = reducFRsNoNans(:,indShuff);

    randShift(iShuff) = randperm(size(reducSigs,2),1);

    %shift should be at least 30 seconds away from original time point
    while randShift(iShuff) < 30000 && randShift(iShuff) > size(reducSigs,2)-30000
        randShift(iShuff) = randperm(size(reducSigs,2),1);
    end

    reducSigsShuff = circshift(reducSigs,randShift(iShuff),2);

    % instead of doing time shifting, actually permute the indices for
    % calculating controls for sparsity (since shifting relative to
    % behavior doesn't serve as a control as behavior isn't used in
    % the sparsity calculation)
    reducSigsPerm = reducSigs(:,randperm(size(reducSigs,2)));

    tic
    for iChan = 1:size(reducSigs,1)

        binSigSums = accumarray(binLabels,reducSigsShuff(iChan,:));
        binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*100;

        binSigSums = accumarray(binLabelsDownsamp,reducSigsShuff(iChan,:));
        binAveSigDownsamp = binSigSums(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;

        binSigSumsPerm = accumarray(binLabelsDownsamp,reducSigsPerm(iChan,:));
        binAveSigPerm = binSigSumsPerm(accumarrayRealOutputsDownsamp)./binNPointsDownsamp(accumarrayRealOutputsDownsamp)*100;

        spaceAveZScores = zeros(nGridPoints,nGridPoints);
        binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
        binAveZScores = binAveSig;

        spaceAveZScores(unique(binLabels)) = binAveZScores;

        spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));

        regionWatershedLabelsCat = regionWatershedLabels;
        for iShedRegion = 1:length(regionWatershedLabelsCat)
            shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
            shedReducSigs = reducSigsShuff(iChan,shedReducInds);
            regionAveSigsShuff(iChan,iShedRegion,iShuff) = mean(shedReducSigs)*1000;
            shedReducSigsNonZeros = shedReducSigs(shedReducSigs ~= 0);
            regionPrctSigsShuff(iChan,iShedRegion,iShuff) = prctile(shedReducSigsNonZeros,90)*1000;

            %get the bins at which the points are in
            regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));

            %and average across that instead
            regionAveSigsFiltShuff(iChan,iShedRegion,iShuff) = mean(regionFiltValues);
            regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
            regionPrctSigsFiltShuff(iChan,iShedRegion,iShuff) = prctile(regionFiltValuesNonZeros,99);
        end

        %sparsity
        hasMinPointsBins = find(binNPointsDownsamp(accumarrayRealOutputsDownsamp) >= 100);
        probVisitBin = binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins))/...
            sum(binNPointsDownsamp(accumarrayRealOutputsDownsamp(hasMinPointsBins)));
        sparsityShuff(iChan,iShuff) = sum(probVisitBin.*binAveSigPerm(hasMinPointsBins))^2 / ...
        sum(probVisitBin.*binAveSigPerm(hasMinPointsBins).^2);

    end
    disp(toc)

end

if ~isempty(saveResultsDir)
    switch activityType
        case 'firingrate'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/NeuronRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff');

        case 'emg'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/MuscleRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff');

        case 'ssaall'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/SSAAllRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff');

        case 'ssastriatum'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/SSAStrRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff');
        
        case 'ssacortex'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/SSACtxRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','sparsity','sparsityShuff');

        case 'pcaall'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/PCAAllRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','pcaWeights','pcaGoodNeurons','sparsity','sparsityShuff');

        case 'pcastriatum'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/PCAStrRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','pcaWeights','pcaGoodNeurons','sparsity','sparsityShuff');

        case 'pcacortex'

            save(fullfile(baseDir,'ProcessedData/UMAPFRs/PCACtxRegionProps'),'regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
                ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff','pcaWeights','pcaGoodNeurons','sparsity','sparsityShuff');

    end
end


% 
