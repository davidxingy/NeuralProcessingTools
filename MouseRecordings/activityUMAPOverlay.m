function activityUMAPOverlay(baseDir, frFile, activityType, varargin)
% activityUMAPOverlay(baseDir, frFile, activityType, savePlotDir, nShuff)
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
% savePlotDir       - Optional input, string specifying the folder to save
%                     the overlay plots to. If empty, will not save any
%                     plots. Default empty.
% 
% nShuff            - Optional input, number of time-shift shuffles to do
%                     for permutation tests. Default 200.
% 
% David Xing, last updated 3/20/2023

% parse optional input parameters
narginchk(3, 5)

if length(varargin)>=1
    if isempty(varargin{1})
        savePlotDir = [];
    else
        savePlotDir = varargin{1};
    end
else
    savePlotDir = [];
end

if length(varargin)>=2
    if isempty(varargin{1})
        nShuff = 200;
    else
        nShuff = varargin{1};
    end
else
    nShuff = 200;
end

% load UMAP projection
load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridInds','watershedRegions','regionWatershedLabels')

% load in sync data
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

% get mapping from EMG to neural time points
% just do linear interpolation
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30)-round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
    (round(frameEMGSamples{1}{end}(end)/20)-round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

% just do it for every time point in reduction
reducNeurInds = round(NeurEMGSync(origDownsampEMGInd*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/30);
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

    case 'ssa'
        %load SSA projections
        % load(fullfile(baseDir,'ProcessedData', 'SSA','allTimePointsSSA.mat'))
        % ssaInput = load(fullfile(baseDir, 'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs');

end

% get watershedding region splitting
[regionBoundaryIndsX, regionBoundaryIndsY] = find(watershedRegions==0);

% Generate grid
nGridPoints = length(gridInds);
rangeVals = [min(reduction(reducIndsToUse,1)) max(reduction(reducIndsToUse,1)); ...
    min(reduction(reducIndsToUse,2)) max(reduction(reducIndsToUse,2))]*1.5;

gridIndsX = gridInds;
gridIndsY = gridInds;
[meshGridX,meshGridY] = meshgrid(gridIndsX,gridIndsY);

% get the gaussian kernal to convolve with for better visualization
densityGaussStd = 0.3;
gaussKernal = exp(-.5.*(meshGridX.^2 + meshGridY.^2)./densityGaussStd^2) ./ (2*pi*densityGaussStd^2);
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

% go through each SSA dimension

% % also get SSA low dimensional projections
% for iRegion = 1:3
%     ssaData = ssaResults{iRegion}.trajs;
%     origSSAInds = 1:size(ssaInput.allFRs,2);
%     ssaNanInds = any(isnan(ssaInput.allFRs));
%     origSSAInds(ssaNanInds) = [];
% 
%     withNanSSA = nan(max(origSSAInds),size(ssaData,2));
%     withNanSSA(origSSAInds,:) = ssaData;
% 
%     reducSSAInds = floor(reducNeurInds(reducIndsToUse)/10);
%     reducSSA = withNanSSA(reducSSAInds,:)';
% 
%     for iSSA = 1:size(reducSSA,1)
% 
%         nonNanInds = find(~isnan(reducSSA(iSSA,:)));
%         binFRSums = accumarray(binLabels(nonNanInds),reducSSA(iSSA,nonNanInds));
%         binAveFRs = binFRSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs);
% 
%         spaceAveZScores = zeros(nGridPoints,nGridPoints);
%         binAveZScores = (binAveFRs - mean(binAveFRs))/std(binAveFRs);
%         binAveZScores = binAveFRs;
% 
%         spaceAveZScores(unique(binLabels)) = binAveZScores;
% 
%         spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));
% 
%         figH = figure('Visible','off','Units','pixels','OuterPosition',[200 50 500 1000]);
%         tiledlayout(2,1,"TileSpacing","tight")
%         nexttile
%         imagesc(gridIndsX, gridIndsY, fliplr(spaceAveZScoresFilt'))
%         hold on
%         plot(-1*(gridIndsX(regionBoundaryIndsY)),gridIndsY(regionBoundaryIndsX),'k.')
%         colorbar
% 
% %             xlim([-12.5,9])
% %             ylim([-9,7])
% %         xlim([-14,8])
% %         ylim([-8,6])
%     xlim([-12,8])
%     ylim([-8,12])
% 
%         %now the average latent value within each of the divided regions
%         watershedRegionsAdj = watershedRegions';
%         watershedRegionsBins = watershedRegionsAdj(:);
%         reducWatershedRegions = watershedRegionsBins(binLabels);
% 
%         regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
%         for iShedRegion = 1:length(regionWatershedLabelsCat)
%             shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
%             ssaReduc = reducSSA(iSSA,shedReducInds);
%             regionAveSSA{iRegion}(iSSA,iShedRegion) = nanmean(abs(ssaReduc));
%             ssaReducNonZeros = ssaReduc(ssaReduc ~= 0);
%             regionPrctSSA{iRegion}(iSSA,iShedRegion) = prctile(abs(ssaReducNonZeros),99)*1000;
% 
%             %get the bins at which the points are in
%             regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
%             %and average across that instead
%             regionAveSSAFilt{iRegion}(iSSA,iShedRegion) = mean(abs(regionFiltValues));
%             regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
%             regionPrctSSAFilt{iRegion}(iSSA,iShedRegion) = prctile(abs(regionFiltValuesNonZeros),99);
% 
%         end
% 
%         nexttile
%         plot(regionAveSSAFilt{iRegion}(iSSA,:))
%         hold on
%         plot(regionAveSSA{iRegion}(iSSA,:))
%         plot(regionPrctSSAFilt{iRegion}(iSSA,:))
% 
%         legend('Filtered Average','Value Averaged','Filtered 99 Percentile','box','off')
%         title([regionNames{iRegion} ' SSA Latent Dim ' num2str(iSSA)])
%         saveas(figH,['./UMAPFRs/', regionNames{iRegion} '_SSALatent' num2str(iSSA) '.png'])
%         close(figH)
% 
% 
%     end
% 
% end


% % go through each muscle
% for iMusc = 1:size(reducEMGs,1)
%     
%     binFRSums = accumarray(binLabels,reducEMGs(iMusc,:));
%     binAveFRs = binFRSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs);
% 
%     spaceAveZScores = zeros(nGridPoints,nGridPoints);
%     binAveZScores = (binAveFRs - mean(binAveFRs))/std(binAveFRs);
%     binAveZScores = binAveFRs;
% 
%     spaceAveZScores(unique(binLabels)) = binAveZScores;
% 
%     spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));
% 
%     figH = figure('Visible','on');
%     imagesc(gridIndsX, gridIndsY, fliplr(spaceAveZScoresFilt'))
%     hold on
%     plot(-1*(gridIndsX(regionBoundaryIndsY)),gridIndsY(regionBoundaryIndsX),'k.')
%     colorbar
% 
%     xlim([-12.5,9])
%     ylim([-9,7])
% 
%     title(['Muscle ' channelNames{iMusc}])
% 
%     saveas(figH,['./UMAPFRs/', 'Muscle_' channelNames{iMusc} '.png'])
% 
% end


% go through each neuron and plot the average firing rates in the space
for iChan = 1:size(reducSigs,1)

    binSigSums = accumarray(binLabels,reducSigs(iChan,:));
    binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*1000;

    spaceAveZScores = zeros(nGridPoints,nGridPoints);
    binAveZScores = (binAveSig - mean(binAveSig))/std(binAveSig);
    binAveZScores = binAveSig;

    %map the firing rates into the actual bins within the grid now (using
    %the bin labels as the index). Only the bins with anything in them are
    %mapped
    spaceAveZScores(unique(binLabels)) = binAveZScores;

    spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));

    figH = figure('Visible','off','Units','pixels','OuterPosition',[200 50 500 1000]);
    tiledlayout(2,1,"TileSpacing","tight")
    nexttile
    imagesc(gridIndsX, gridIndsY, fliplr(spaceAveZScoresFilt'))
    hold on
    plot(-1*(gridIndsX(regionBoundaryIndsY)),gridIndsY(regionBoundaryIndsX),'k.')
    colorbar

%     xlim([-12.5,9])
%     ylim([-9,7])
%     xlim([-12,8])
%     ylim([-8,12])
%     xlim([-4,8])
%     ylim([-10,6])
%     xlim([-10,10])
%     ylim([-6,9])

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
            figFilename = ['./UMAPFRs/' region '_Neuron' num2str(neuronNum) '.png'];

        case 'emg'
            title(channelNames{iChan})
            figFilename = ['./UMAPFRs/Muscle_' replace(channelNames{iChan},' ','-') '.png'];
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

    nexttile
    plot(regionAveSigsFilt(iChan,:))
    hold on
    plot(regionAveSigs(iChan,:))
    plot(regionPrctSigsFilt(iChan,:))
    legend('Filtered Average','Value Averaged','Filtered 99 Percentile','box','off')


    saveas(figH,figFilename)
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

    tic
    for iChan = 1:size(reducSigs,1)

        binSigSums = accumarray(binLabels,reducSigsShuff(iChan,:));
        binAveSig = binSigSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*1000;

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

    end
    disp(toc)

end

switch activityType
    case 'firingrate'

        save('./UMAPFRs/NeuronRegionProps','regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
            ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff');

    case 'emg'

        save('./UMAPFRs/MuscleRegionProps','regionAveSigs','regionPrctSigs','regionAveSigsFilt','regionPrctSigsFilt'...
            ,'regionAveSigsShuff','regionPrctSigsShuff','regionAveSigsFiltShuff','regionPrctSigsFiltShuff');

end



% 
