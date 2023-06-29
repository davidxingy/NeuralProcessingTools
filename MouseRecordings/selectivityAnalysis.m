% function selectivityAnalysis(baseDir)
% selectivityAnalysis(baseDir)
% function to calculate selectivity metrics from the UMAP regions activity overlay


load(fullfile(baseDir,'ProcessedData','UMAPFRs','NeuronRegionProps.mat'))
load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds')
load(fullfile(baseDir,'ProcessedData','UMAP.mat'), 'regionBehvAssignments','analyzedBehaviors')
load(fullfile(baseDir,'ProcessedData','neuronDataStruct.mat'))

allInds = [striatumInds cortexInds];
depths = [neuronDataStruct(allInds).depth];
ts = {neuronDataStruct(allInds).timeStamps};

if exist("regionAveSigs")
    regionAveFRs = regionAveSigs;
    regionAveFRsShuff = regionAveSigsShuff;
end

nRegions = size(regionAveFRs,2);
for iRegion = 1:nRegions
    regionBehvNames{iRegion} = string(join(analyzedBehaviors(regionBehvAssignments{iRegion}),'/'));
end

multiBehvSpec = calcMultiBehvSpec(regionAveFRs);
multiBehvSpecShuff = calcMultiBehvSpec(regionAveFRsShuff);
[singleBehvSpec, allSingleBehvSpec, singleBehvSpecRegion] = calcSingleBehvSpec2(regionAveFRs);
[singleBehvSpecShuff, allSingleBehvSpecShuff, singleBehvSpecRegionShuff] = calcSingleBehvSpec2(regionAveFRsShuff);

allSingleBehvSpecCorrected = allSingleBehvSpec - squeeze(nanmean(allSingleBehvSpecShuff));
singleBehvSpecCorrected = singleBehvSpec - nanmean(singleBehvSpecShuff);
allMultiBehvSpecCorrected = multiBehvSpec - nanmean(multiBehvSpecShuff,2);

strSpec = sortrows([allMultiBehvSpecCorrected(1:length(striatumInds))...
    singleBehvSpecCorrected(1:length(striatumInds))'...
    (1:length(striatumInds))' singleBehvSpecRegion(1:length(striatumInds))']);
strSpec(isnan(strSpec(:,1)),:) = [];
[cdfFreqStr, cdfSpecStr] = ecdf(strSpec(:,1));

ctxSpec = sortrows([allMultiBehvSpecCorrected(length(striatumInds)+1:end)...
    singleBehvSpecCorrected(length(striatumInds)+1:end)'...
    (length(striatumInds)+1:length(multiBehvSpec))' singleBehvSpecRegion(length(striatumInds)+1:end)']);
ctxSpec(isnan(ctxSpec(:,1)),:) = [];
[cdfFreqCtx, cdfSpecCtx] = ecdf(ctxSpec(:,1));

figure;
cdfplot(strSpec(:,1))
set(get(gca,'Children'),'Color','k')
hold on
cdfplot(ctxSpec(:,1))
set(get(gca,'Children'),'Color','k')

plotColors = lines(7);
for iRegion = 1:7
    
    regionNeurons = find(strSpec(:,4)==iRegion);
    plotH(iRegion) = plot(cdfSpecStr(regionNeurons),cdfFreqStr(regionNeurons),'.','MarkerSize',10,'color',plotColors(iRegion,:));

    regionNeurons = find(ctxSpec(:,4)==iRegion);
    plot(cdfSpecCtx(regionNeurons),cdfFreqCtx(regionNeurons),'.','MarkerSize',10,'color',plotColors(iRegion,:));

end

legend(plotH,regionBehvNames,'box', 'off')

% look correlation across behaviors

regionAveFRsCtx = regionAveFRs(length(striatumInds)+1:end,:);
regionAveFRsCtx(any(isnan(regionAveFRsCtx),2) | any(isinf(regionAveFRsCtx),2),:) = [];
regionAveFRsStr = regionAveFRs(1:length(striatumInds),:);
regionAveFRsStr(any(isnan(regionAveFRsStr),2) | any(isinf(regionAveFRsStr),2),:) = [];

singleBehvSpecCtx = allSingleBehvSpecCorrected(length(striatumInds)+1:end,:);
singleBehvSpecCtx(any(isnan(singleBehvSpecCtx),2) | any(isinf(singleBehvSpecCtx),2),:) = [];
singleBehvSpecStr = allSingleBehvSpecCorrected(1:length(striatumInds),:);
singleBehvSpecStr(any(isnan(singleBehvSpecStr),2) | any(isinf(singleBehvSpecStr),2),:) = [];

for iRegion1 = 1:nRegions
    for iRegion2 = 1:nRegions



        % first do all neurons
        linFit = fitlm([regionAveFRsStr(:,iRegion1); regionAveFRsCtx(:,iRegion1)],...
            [regionAveFRsStr(:,iRegion2); regionAveFRsCtx(:,iRegion2)]);
        frsRegionFitR2{1}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        frsRegionFitSlope{1}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

        linFit = fitlm([singleBehvSpecStr(:,iRegion1); singleBehvSpecCtx(:,iRegion1)],...
            [singleBehvSpecStr(:,iRegion2); singleBehvSpecCtx(:,iRegion2)]);
        specsRegionFitR2{1}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        specsRegionFitSlope{1}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

        %next, just striatum
        linFit = fitlm(regionAveFRsStr(:,iRegion1), regionAveFRsStr(:,iRegion2));
        frsRegionFitR2{2}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        frsRegionFitSlope{2}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

        linFit = fitlm(singleBehvSpecStr(:,iRegion1), singleBehvSpecStr(:,iRegion2));
        specsRegionFitR2{2}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        specsRegionFitSlope{2}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

        %next, just cortex
        linFit = fitlm(regionAveFRsCtx(:,iRegion1), regionAveFRsCtx(:,iRegion2));
        frsRegionFitR2{3}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        frsRegionFitSlope{3}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

        linFit = fitlm(singleBehvSpecCtx(:,iRegion1), singleBehvSpecCtx(:,iRegion2));
        specsRegionFitR2{3}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        specsRegionFitSlope{3}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);

    end
end

cgObj = clustergram(frsRegionFitR2{1});
clustOrdering = cellfun(@(x) str2num(x),cgObj.RowLabels);


% singleBehvSpecCorrected = singleBehvSpec - mean

colorMap = jet(nRegions);

figure;
randShifts = randn(1,length(depths))*20;
rasterplot(cellfun(@(x) x/30000, ts,'un',0),'times','.',gca,[],depths+randShifts)

for iRegion = 1:nRegions

    hold on

    regionNeurons{iRegion} = find(singleBehvSpecRegion==iRegion & singleBehvSpec-nanmean(singleBehvSpecShuff) > 0.5);
    
    legH(iRegion) = plot([0 0],[1000 1001],'Color',colorMap(iRegion,:));

    if isempty(regionNeurons{iRegion})
        continue
    end

    rasterplot(cellfun(@(x) x/30000, ts(regionNeurons{iRegion}),'un',0),...
        'times','.',gca,[],depths(regionNeurons{iRegion})+randShifts(regionNeurons{iRegion}),{'cData', colorMap(iRegion,:),'sizedata',30})
    

end

legNames = {'Eating','Grooming','Still/Rear','Walk/Jump','Climb Down','ClimbUp'};
legend(legH,legNames)
ylim([0 4000])


figure;
tiledlayout(nRegions,1)
for iRegion = 1:nRegions
    nexttile
    histH(iRegion) = histogram(depths(regionNeurons{iRegion}),0:100:4000,'FaceColor',colorMap(iRegion,:));
    title(legNames{iRegion})
    xlim([0 4000])
    box off
    ylabel('# Neurons')
end

% legend(histH,legNames)




function [singleBehvSpec singleBehvSpecRegion] = calcSingleBehvSpec1(regionMetrics)

if ndims(regionMetrics) == 3
    inputIsShuff = true;
    nShuff = size(regionMetrics,3);
else
    inputIsShuff = false;
    nShuff = 1;
end

for iShuff = 1:nShuff

    if inputIsShuff
        shuffMetrics = squeeze(regionMetrics(:,:,iShuff));
    else
        shuffMetrics = regionMetrics;
    end

    for iRegion = 1:size(shuffMetrics,2)

        if inputIsShuff
            thisRegionMetrics = shuffMetrics(:,iRegion);
        else
            thisRegionMetrics = squeeze(shuffMetrics(:,iRegion,:));
        end
        otherRegionMetrics = shuffMetrics;
        otherRegionMetrics(:,iRegion) = [];
        [otherRegionSorted, otherRegionSortedInds] = sort(otherRegionMetrics,2);
        nextHighestMetrics = otherRegionSorted(:,end);

        specs(:,iRegion) = (thisRegionMetrics - nextHighestMetrics)./(thisRegionMetrics + nextHighestMetrics);
    end
    [singleBehvSpec(iShuff,:), singleBehvSpecRegion(iShuff,:)] = max(specs,[],2);

end

end



function [singleBehvSpec, allSingleBehvSpecs, singleBehvSpecRegion] = calcSingleBehvSpec2(regionMetrics)
if ndims(regionMetrics) == 3
    inputIsShuff = true;
    nShuff = size(regionMetrics,3);
else
    inputIsShuff = false;
    nShuff = 1;
end

for iShuff = 1:nShuff

    if inputIsShuff
        shuffMetrics = squeeze(regionMetrics(:,:,iShuff));
    else
        shuffMetrics = regionMetrics;
    end

    for iRegion = 1:size(shuffMetrics,2)

        if inputIsShuff
            thisRegionMetrics = shuffMetrics(:,iRegion);
        else
            thisRegionMetrics = squeeze(shuffMetrics(:,iRegion,:));
        end
        otherRegionMetrics = shuffMetrics;
        otherRegionMetrics(:,iRegion) = [];

        specs(:,iRegion) = (thisRegionMetrics - mean(otherRegionMetrics,2))./(thisRegionMetrics);% + mean(otherRegionMetrics,2));
    end
    [singleBehvSpec(iShuff,:), singleBehvSpecRegion(iShuff,:)] = max(specs,[],2);
    allSingleBehvSpecs(iShuff,:,:) = specs;

end

allSingleBehvSpecs = squeeze(allSingleBehvSpecs);


end


function multiBehvSpec = calcMultiBehvSpec(regionMetrics)
regionMetricsNormalized = (regionMetrics./sum(regionMetrics,2));
multiBehvSpec = squeeze(sum(regionMetricsNormalized.^2,2));

end


% 
