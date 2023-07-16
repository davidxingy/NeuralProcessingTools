% function selectivityAnalysis(baseDir)
% selectivityAnalysis(baseDir)
% function to calculate selectivity metrics from the UMAP regions activity overlay

clear

% do analysis for each of the datasets
recordingSessions = {
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording', ...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording' ...
};

behvAlignPerm = [
    1 2 3 4 5 6 7; ...
    7 6 3 5 4 2 1; ...
    1 2 4 5 3 6 7 ...
    ];

for iSess = 1:length(recordingSessions)

    baseDir = recordingSessions{iSess};
    load(fullfile(baseDir,'ProcessedData','UMAPFRs','NeuronRegionProps.mat'))
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds')
    load(fullfile(baseDir,'ProcessedData','UMAP.mat'), 'regionBehvAssignments','analyzedBehaviors')
    load(fullfile(baseDir,'ProcessedData','neuronDataStruct.mat'))

    allInds = [striatumInds cortexInds];
    depths = [neuronDataStruct(allInds).depth];
    ts = {neuronDataStruct(allInds).timeStamps};

    % load in both neural and EMG selectivity metrics
    regionAveFRs = regionAveSigs;
    regionAveFRsShuff = regionAveSigsShuff;

    load(fullfile(baseDir,'ProcessedData','UMAPFRs','MuscleRegionProps.mat'))

    regionAveEMGs = regionAveSigs;
    regionAveEMGsShuff = regionAveSigsShuff;

    % get the names associated with each region
    nRegions = size(regionAveFRs,2);
    for iRegion = 1:nRegions
        regionBehvNames{iSess,iRegion} = string(join(analyzedBehaviors(regionBehvAssignments{iRegion}),'/'));
    end

    % get both the multibehavior specificity and single behavior specificty for
    % each neuron
    fieldNames = {'str','ctx','emg'};
    % do it for strialtal and cortical neurons as well as muscles
    for iType = 1:length(fieldNames)

        if iType == 1
            regionAveVals.(fieldNames{iType}){iSess} = regionAveFRs(1:length(striatumInds),:);
            regionAveValsShuff = regionAveFRsShuff(1:length(striatumInds),:,:);
        elseif iType == 2
            regionAveVals.(fieldNames{iType}){iSess} = regionAveFRs(length(striatumInds)+1:end,:);
            regionAveValsShuff = regionAveFRsShuff(length(striatumInds)+1:end,:,:);
        elseif iType == 3
            regionAveVals.(fieldNames{iType}){iSess} = regionAveEMGs;
            regionAveValsShuff = regionAveEMGsShuff;
        end

        multiBehvSpec.(fieldNames{iType}){iSess} = calcMultiBehvSpec(regionAveVals.(fieldNames{iType}){iSess});
        multiBehvSpecShuff.(fieldNames{iType}){iSess} = calcMultiBehvSpec(regionAveValsShuff);
        [singleBehvSpec.(fieldNames{iType}){iSess}, allSingleBehvSpec.(fieldNames{iType}){iSess}, singleBehvSpecRegion.(fieldNames{iType}){iSess}] = calcSingleBehvSpec2(regionAveFRs);
        [singleBehvSpecShuff.(fieldNames{iType}){iSess}, allSingleBehvSpecShuff.(fieldNames{iType}){iSess}, singleBehvSpecRegionShuff.(fieldNames{iType}){iSess}] = calcSingleBehvSpec2(regionAveFRsShuff);

        % correct for low firing rate neurons
        allSingleBehvSpecCorrected.(fieldNames{iType}){iSess} = allSingleBehvSpec.(fieldNames{iType}){iSess} - squeeze(nanmean(allSingleBehvSpecShuff.(fieldNames{iType}){iSess}));
        singleBehvSpecCorrected.(fieldNames{iType}){iSess} = singleBehvSpec.(fieldNames{iType}){iSess} - nanmean(singleBehvSpecShuff.(fieldNames{iType}){iSess});
        allMultiBehvSpecCorrected.(fieldNames{iType}){iSess} = multiBehvSpec.(fieldNames{iType}){iSess} - nanmean(multiBehvSpecShuff.(fieldNames{iType}){iSess},2);

        %calc CDF for the multibehavioral selectivity
        multiSpecCdfVals = -0.4:0.02:1;
        multiSpecCdfFreq.(fieldNames{iType}){iSess} = calcCDF(allMultiBehvSpecCorrected.(fieldNames{iType}){iSess},multiSpecCdfVals);

    end

end

% make summary plot for multibehavior selectivity

allSessStrSpec = cat(1,multiSpecCdfFreq.str{:});
allSessCtxSpec = cat(1,multiSpecCdfFreq.ctx{:});
allSessEmgSpec = cat(1,multiSpecCdfFreq.emg{:});

plotColors = lines(3);
figure;
hold on;
for iSess = 1:length(recordingSessions)
    firstHitMax = find(allSessStrSpec(iSess,:)==1,1);
    plot(multiSpecCdfVals(1:firstHitMax),allSessStrSpec(iSess,1:firstHitMax),'Color',[plotColors(2,:) 0.2],'LineWidth',1.5);

    firstHitMax = find(allSessCtxSpec(iSess,:)==1,1);
    plot(multiSpecCdfVals(1:firstHitMax),allSessCtxSpec(iSess,1:firstHitMax),'Color',[plotColors(1,:) 0.2],'LineWidth',1.5);

    firstHitMax = find(allSessEmgSpec(iSess,:)==1,1);
    plot(multiSpecCdfVals(1:firstHitMax),allSessEmgSpec(iSess,1:firstHitMax),'Color',[plotColors(3,:) 0.2],'LineWidth',1.5);
end

firstHitMax = find(mean(allSessStrSpec)==1,1);
legH(1) = plot(multiSpecCdfVals(1:firstHitMax),mean(allSessStrSpec(:,1:firstHitMax)),'Color',plotColors(2,:),'LineWidth',3);
firstHitMax = find(mean(allSessCtxSpec)==1,1);
legH(2) = plot(multiSpecCdfVals(1:firstHitMax),mean(allSessCtxSpec(:,1:firstHitMax)),'Color',plotColors(1,:),'LineWidth',3);
firstHitMax = find(mean(allSessEmgSpec)==1,1);
legH(3) = plot(multiSpecCdfVals(1:firstHitMax),mean(allSessEmgSpec(:,1:firstHitMax)),'Color',plotColors(3,:),'LineWidth',3);

legend(legH,'Striatum','Cortex','Muscles','box','off','fontsize',13)
xlabel('Adjusted Skew Metric')
ylabel('Frequency')
box off
set(gcf,'color','w')
set(gca,'LineWidth',1.5)
set(gca,'fontsize',13)
set(gca,'TickDir','out')

% % plotting region of highest specificity
% plotColors = lines(7);
% for iRegion = 1:7
%     
%     regionNeurons = find(strSpec(:,4)==iRegion);
%     plotH(iRegion) = plot(cdfSpecStr(regionNeurons),cdfFreqStr(regionNeurons),'.','MarkerSize',10,'color',plotColors(iRegion,:));
% 
%     regionNeurons = find(ctxSpec(:,4)==iRegion);
%     plot(cdfSpecCtx(regionNeurons),cdfFreqCtx(regionNeurons),'.','MarkerSize',10,'color',plotColors(iRegion,:));
% 
% end

% make example plot of all neurons sorted by region
for iSess = 1:length(recordingSessions)
    strAligned{iSess} = regionAveVals.str{iSess}(:,behvAlignPerm(iSess,:));
    ctxAligned{iSess} = regionAveVals.ctx{iSess}(:,behvAlignPerm(iSess,:));
end

strAll = cat(1,strAligned{:});
maxSortPermStr = sortByLevelRecursive(strAll);
plotDataStr = strAll(maxSortPermStr,:)./sum(strAll(maxSortPermStr,:),2);
plotDataStr(any(isnan(plotDataStr),2),:) = [];

ctxAll = cat(1,ctxAligned{:});
maxSortPermCtx = sortByLevelRecursive(ctxAll);
plotDataCtx = ctxAll(maxSortPermCtx,:)./sum(ctxAll(maxSortPermCtx,:),2);
plotDataCtx(any(isnan(plotDataCtx),2),:) = [];

figure
imagesc([plotDataStr; plotDataCtx])
colorbar
hold on
line([0 7.5],[size(plotDataStr,1) size(plotDataStr,1)],'linewidth',2,'color','r','linestyle','--')

title('Behavioral Selectivities')
ylabel('Neuron')
set(gca,'XTickLabel',regionBehvNames(1,:))
set(gca,'XTickLabelRotation',90)
set(gcf,'Color','w')



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


function freq = calcCDF(inputData,range)

% remove nans
inputData(isnan(inputData)) = [];

for i = 1:length(range)

    freq(i) = sum(inputData <= range(i)) / length(inputData);

end

end



function sortPerm = sortByLevelRecursive(data)

[~, maxRegion] = max(data,[],2);
maxSortTable = sortrows([maxRegion,(1:length(maxRegion))']);
sortPerm = maxSortTable(:,2);

if size(data,2) == 1
    return
end

for iSub = 1:size(data,2)

    subInds = find(maxSortTable(:,1)==iSub);
    subLevelData = data(sortPerm(subInds),:);
    subLevelData(:,iSub) = [];

    subPerm = sortByLevelRecursive(subLevelData);
    sortPerm(subInds) = sortPerm(subInds(subPerm));

end

end


% 
