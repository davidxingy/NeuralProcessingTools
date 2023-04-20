function selectivityAnalysis(baseDir)
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

multiBehvSpec = calcMultiBehvSpec(regionAveFRs);
multiBehvSpecShuff = calcMultiBehvSpec(regionAveFRsShuff);
[singleBehvSpec, singleBehvSpecRegion] = calcSingleBehvSpec2(regionAveFRs);
[singleBehvSpecShuff, singleBehvSpecRegionShuff] = calcSingleBehvSpec2(regionAveFRsShuff);


nRegions = size(regionAveFRs,2);
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



function [singleBehvSpec singleBehvSpecRegion] = calcSingleBehvSpec2(regionMetrics)
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

end



function multiBehvSpec = calcMultiBehvSpec(regionMetrics)
regionMetricsNormalized = (regionMetrics./sum(regionMetrics,2));
multiBehvSpec = squeeze(sum(regionMetricsNormalized.^2,2));



% 
