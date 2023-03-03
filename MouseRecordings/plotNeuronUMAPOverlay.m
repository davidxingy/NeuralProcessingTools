clear
close all

baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording';
frFile = 'NeuralFiringRates1msBins10msGauss.mat';
umapDownSamp = 1;

% load firing rate data
load(fullfile(baseDir,'ProcessedData',frFile))
clear striatumFRs cortexFRs noNanFRs
% load UMAP projection
load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridInds','watershedRegions','regionWatershedLabels')
% load in sync data
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

% also load EMG data (for looking at how each muscle is distributed)    
load(fullfile(baseDir,'ProcessedData', 'EMG1ms.mat'))
load(fullfile(baseDir,'ProcessedData', 'D024-111022-ArenaRecording_ProcessedEMG_MetaData.mat'))

% % also load SSA Data (and corresponding input data) 
% load(fullfile(baseDir,'ProcessedData', 'SSA','allTimePointsSSA.mat'))
% ssaInput = load(fullfile(baseDir, 'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs');

% get watershedding region splitting
[regionBoundaryIndsX, regionBoundaryIndsY] = find(watershedRegions==0);

downSampleAmount = 1;

% get mapping from EMG to neural time points
% just do linear interpolation
emgNeurSlope = (round(frameNeuropixelSamples{1}{end}(end)/30)-round(frameNeuropixelSamples{1}{1}(1)/30)) / ...
    (round(frameEMGSamples{1}{end}(end)/20)-round(frameEMGSamples{1}{1}(1)/20));
emgNeurOffset = round(frameNeuropixelSamples{1}{1}(1)/30) - emgNeurSlope*round(frameEMGSamples{1}{1}(1)/20);

% just do it for every time point in reduction
reducNeurInds = round(round(origDownsampEMGInd(1:umapDownSamp:end)*emgNeurSlope + emgNeurOffset)/1);
maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/30);
reducIndsToUse = find(reducNeurInds < maxNeurSamples);

% get neural data corresponding to the UMAP time points
reducFRs = allFRs(:,reducNeurInds(reducIndsToUse));
clear allFRs

% don't use any points where there's nans in the FR data
nanReducPoints = any(isnan(reducFRs),1);
reducIndsToUse(nanReducPoints) = [];
reducFRsNoNans = reducFRs(:,~nanReducPoints);
clear reducFRs

% get emg
reducEMGInds = origDownsampEMGInd(1:umapDownSamp:end);
reducEMGs = downsampEMG(:,reducEMGInds(reducIndsToUse));

% Generate grid
nGridPoints = length(gridInds);
rangeVals = [min(reduction(reducIndsToUse,1)) max(reduction(reducIndsToUse,1)); ...
    min(reduction(reducIndsToUse,2)) max(reduction(reducIndsToUse,2))]*1.5;

% rangeVals = [-15 15; -15 15];
% xx = linspace(rangeVals(1,1),rangeVals(1,2),nGridPoints);
% yy = linspace(rangeVals(2,1),rangeVals(2,2),nGridPoints);
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

% figure;
% imagesc(xx,yy,fliplr(density))
% colorbar
% hold on
% plot(-xx(jj),xx(ii),'k.')
% xlim([-12.5,9])
% ylim([-9,7])


% % go through each SSA dimension
% 
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
%         figH = figure('Visible','on');
%         imagesc(gridIndsX, gridIndsY, fliplr(spaceAveZScoresFilt'))
%         hold on
%         plot(-1*(gridIndsX(regionBoundaryIndsY)),gridIndsY(regionBoundaryIndsX),'k.')
%         colorbar
% 
%         %     xlim([-12.5,9])
%         %     ylim([-9,7])
% %         xlim([-14,8])
% %         ylim([-8,6])
%     xlim([-12,8])
%     ylim([-8,12])
% 
%         title(['SSA Latent Dim ' num2str(iSSA)])
% 
%         saveas(figH,['./UMAPFRs/', regionNames{iRegion} '_SSALatent' num2str(iSSA) '.png'])
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
%             regionAveSSA{iRegion}(iSSA,iShedRegion) = nanmean(ssaReduc);
%             ssaReducNonZeros = ssaReduc(ssaReduc ~= 0);
%             regionPrctSSA{iRegion}(iSSA,iShedRegion) = prctile(ssaReducNonZeros,99.99)*1000;
% 
%             %get the bins at which the points are in
%             regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
%             %and average across that instead
%             regionAveSSAFilt{iRegion}(iSSA,iShedRegion) = mean(regionFiltValues);
%             regionFiltValuesNonZeros = regionFiltValues(regionFiltValues ~= 0);
%             regionPrctSSAFilt(iSSA,iShedRegion) = prctile(regionFiltValuesNonZeros,99);
% 
%         end
% 
%         %     %also do permutation
%         %     for iShuff = 1:nShuff
%         %         shuffSSA = reducSSA(iSSA,indShuff(iShuff,:));
%         %
%         %         for iShedRegion = 1:length(regionWatershedLabelsCat)
%         %             shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
%         %             ssaReduc = shuffSSA(shedReducInds);
%         %             regionAveSSAShuff(iSSA,iShedRegion,iShuff) = nanmean(ssaReduc);
%         %             ssaReducNonZeros = ssaReduc(ssaReduc ~= 0);
%         %             regionPrctSSAShuff(iSSA,iShedRegion,iShuff) = prctile(ssaReducNonZeros,99.99)*1000;
%         %
%         %             %get the bins at which the points are in
%         %             regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
%         %             %and average across that instead
%         %             regionAveSSAFiltShuff(iSSA,iShedRegion,iShuff) = mean(regionFiltValues);
%         %
%         %         end
%         %
%         %     end
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
for iNeuron = 1:size(reducFRsNoNans,1)

    binFRSums = accumarray(binLabels,reducFRsNoNans(iNeuron,:));
    binAveFRs = binFRSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*1000;

    spaceAveZScores = zeros(nGridPoints,nGridPoints);
    binAveZScores = (binAveFRs - mean(binAveFRs))/std(binAveFRs);
    binAveZScores = binAveFRs;

    %map the firing rates into the actual bins within the grid now (using
    %the bin labels as the index). Only the bins with anything in them are
    %mapped
    spaceAveZScores(unique(binLabels)) = binAveZScores;

    spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));

    figH = figure('Visible','off','Units','pixels','OuterPosition',[200 50 500 1000]);
    tiledlayout(2,1,"TileSpacing","tight")
%     nexttile
%     imagesc(xx, yy, spaceAveZScores)
    nexttile
    imagesc(gridIndsX, gridIndsY, fliplr(spaceAveZScoresFilt'))
    hold on
    plot(-1*(gridIndsX(regionBoundaryIndsY)),gridIndsY(regionBoundaryIndsX),'k.')
    colorbar

%     xlim([-12.5,9])
%     ylim([-9,7])
%     xlim([-12,8])
%     ylim([-8,12])
    xlim([-4,8])
    ylim([-10,6])

    if iNeuron > length(striatumInds)
        region = 'Cortex';
        neuronNum = iNeuron - length(striatumInds);
    else
        region = 'Striatum';
        neuronNum = iNeuron;
    end
    title([region ', Neuron ' num2str(neuronNum)])

    %now the average firing rate within each of the divided regions
    watershedRegionsAdj = watershedRegions';
    watershedRegionsBins = watershedRegionsAdj(:);
    reducWatershedRegions = watershedRegionsBins(binLabels);

    regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
    for iShedRegion = 1:length(regionWatershedLabelsCat)
        shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
        neuronReducFR = reducFRsNoNans(iNeuron,shedReducInds);
        regionAveFRs(iNeuron,iShedRegion) = mean(neuronReducFR)*1000;
        neuronReducFRNonZeros = neuronReducFR(neuronReducFR ~= 0);
        regionPrctFRs(iNeuron,iShedRegion) = prctile(neuronReducFRNonZeros,90)*1000;

        %get the bins at which the points are in
        regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));
        %and average across that instead
        regionAveFRFilt(iNeuron,iShedRegion) = mean(regionFiltValues);
        neuronReducFRNonZeros = regionFiltValues(regionFiltValues ~= 0);
        regionPrctFRFilt(iNeuron,iShedRegion) = prctile(regionFiltValues,99);
    end

    nexttile
    plot(regionAveFRFilt(iNeuron,:))
    hold on
    plot(regionAveFRs(iNeuron,:))

    saveas(figH,['./UMAPFRs/' region '_Neuron' num2str(neuronNum) '.png'])
    close(figH)

end

% also do permutation test and do the same calculations
% do permutation shuffles for ssa and neurons
nShuff = 200;
nTotalInds = size(reducFRsNoNans,2);
for iShuff = 1:nShuff   
    
%     indShuff = randperm(nTotalInds);
%     reducFRsNoNansShuff = reducFRsNoNans(:,indShuff);

    randShift(iShuff) = randperm(size(reducFRsNoNans,2),1);
    reducFRsNoNansShuff = circshift(reducFRsNoNans,randShift(iShuff),2);

    tic
    for iNeuron = 1:size(reducFRsNoNans,1)

        binFRSums = accumarray(binLabels,reducFRsNoNansShuff(iNeuron,:));
        binAveFRs = binFRSums(accumarrayRealOutputs)./binNPoints(accumarrayRealOutputs)*1000;

        spaceAveZScores = zeros(nGridPoints,nGridPoints);
        binAveZScores = (binAveFRs - mean(binAveFRs))/std(binAveFRs);
        binAveZScores = binAveFRs;

        spaceAveZScores(unique(binLabels)) = binAveZScores;

        spaceAveZScoresFilt = fftshift(real(ifft2(fft2(gaussKernal).*fft2(spaceAveZScores))));

        regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
        for iShedRegion = 1:length(regionWatershedLabelsCat)
            shedReducInds = find(reducWatershedRegions == regionWatershedLabelsCat(iShedRegion));
            neuronReducFR = reducFRsNoNansShuff(iNeuron,shedReducInds);
            regionAveFRsShuff(iNeuron,iShedRegion,iShuff) = mean(neuronReducFR)*1000;
            neuronReducFRNonZeros = neuronReducFR(neuronReducFR ~= 0);
            regionPrctFRsShuff(iNeuron,iShedRegion,iShuff) = prctile(neuronReducFRNonZeros,90)*1000;

            %get the bins at which the points are in
            regionFiltValues = spaceAveZScoresFilt(binLabels(shedReducInds));

            %and average across that instead
            regionAveFRFiltShuff(iNeuron,iShedRegion,iShuff) = mean(regionFiltValues);
            neuronReducFRNonZeros = regionFiltValues(regionFiltValues ~= 0);
            regionPrctFRFiltShuff(iNeuron,iShedRegion,iShuff) = prctile(regionFiltValues,99);
        end

    end
    disp(toc)

end



% 
