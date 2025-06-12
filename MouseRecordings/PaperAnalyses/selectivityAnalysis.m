% function selectivityAnalysis(baseDir)
% selectivityAnalysis(baseDir)
% function to calculate selectivity metrics from the UMAP regions activity overlay

clear
% close all

% do analysis for each of the datasets
recordingSessions = {
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording', ...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording' ...
};

behvAlignPerm = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump Down','Walk','Misc/Rearing/Still','Groom','Eat'};

for iSess = 1:length(recordingSessions)

    baseDir = recordingSessions{iSess};
    load(fullfile(baseDir,'ProcessedData','UMAPFRs','NeuronRegionProps.mat'))
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds')
    load(fullfile(baseDir,'ProcessedData','UMAP.mat'), 'regionBehvAssignments','analyzedBehaviors')
    load(fullfile(baseDir,'ProcessedData','neuronDataStruct.mat'))
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'),'cortexFRs','striatumFRs')
    
    maxStrFRs{iSess} = max(striatumFRs,[],2);
    stdStrFRs{iSess} = nanstd(striatumFRs,[],2);
    maxCtxFRs{iSess} = max(cortexFRs,[],2);
    stdCtxFRs{iSess} = nanstd(cortexFRs,[],2);

    allInds = [striatumInds cortexInds];
    depths = [neuronDataStruct(allInds).depth];
    ts = {neuronDataStruct(allInds).timeStamps};

    % load in both neural and EMG selectivity metrics
    regionAveFRs = regionAveSigs;
    regionAveFRsShuff = regionAveSigsShuff(:,:,1:100);
    sparsityFRs = sparsity;
    sparsityFRsShuff = sparsityShuff;

    load(fullfile(baseDir,'ProcessedData','UMAPFRs','MuscleRegionProps.mat'))

    regionAveEMGs = regionAveSigs;
    regionAveEMGsShuff = regionAveSigsShuff(:,:,1:100);
    sparsityEMGs = sparsity;
    sparsityEMGsShuff = sparsityShuff;

    % get the names associated with each region
    nRegions = size(regionAveFRs,2);
    for iRegion = 1:nRegions
        regionBehvNames{iSess,iRegion} = string(join(analyzedBehaviors(regionBehvAssignments{iRegion}),'/'));
    end

    % keep track of which neurons are interneurons
    strPeakToValleyCutoff = 12.45;
    ctxPeakToValleyCutoff = 11.25;
    isInterneuron.str{iSess} = [neuronDataStruct(striatumInds).peakToValley] < strPeakToValleyCutoff;
    isInterneuron.ctx{iSess} = [neuronDataStruct(cortexInds).peakToValley] < ctxPeakToValleyCutoff;

    % get both the multibehavior specificity and single behavior specificty for
    % each neuron
    fieldNames = {'str','ctx','emg'};
    % do it for strialtal and cortical neurons as well as muscles
    for iType = 1:length(fieldNames)

        if iType == 1
            regionAveVals.(fieldNames{iType}){iSess} = regionAveFRs(1:length(striatumInds),:);
            regionAveValsShuff.(fieldNames{iType}){iSess} = regionAveFRsShuff(1:length(striatumInds),:,:);

            sparsityVals.(fieldNames{iType}){iSess} = sparsityFRs(1:length(striatumInds));
            sparsityValsShuff.(fieldNames{iType}){iSess} = sparsityFRsShuff(1:length(striatumInds),:);

            [singleBehvSpec.(fieldNames{iType}){iSess}, allSingleBehvSpec.(fieldNames{iType}){iSess}, singleBehvSpecRegion.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveFRs(1:length(striatumInds),:));
            [singleBehvSpecShuff.(fieldNames{iType}){iSess}, allSingleBehvSpecShuff.(fieldNames{iType}){iSess}, singleBehvSpecRegionShuff.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveFRsShuff(1:length(striatumInds),:,:));

        elseif iType == 2
            regionAveVals.(fieldNames{iType}){iSess} = regionAveFRs(length(striatumInds)+1:end,:);
            regionAveValsShuff.(fieldNames{iType}){iSess} = regionAveFRsShuff(length(striatumInds)+1:end,:,:);

            sparsityVals.(fieldNames{iType}){iSess} = sparsityFRs(length(striatumInds)+1:end);
            sparsityValsShuff.(fieldNames{iType}){iSess} = sparsityFRsShuff(length(striatumInds)+1:end,:);

            [singleBehvSpec.(fieldNames{iType}){iSess}, allSingleBehvSpec.(fieldNames{iType}){iSess}, singleBehvSpecRegion.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveFRs(length(striatumInds)+1:end,:));
            [singleBehvSpecShuff.(fieldNames{iType}){iSess}, allSingleBehvSpecShuff.(fieldNames{iType}){iSess}, singleBehvSpecRegionShuff.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveFRsShuff(length(striatumInds)+1:end,:,:));

        elseif iType == 3
            regionAveVals.(fieldNames{iType}){iSess} = regionAveEMGs;
            regionAveValsShuff.(fieldNames{iType}){iSess} = regionAveEMGsShuff;

            sparsityVals.(fieldNames{iType}){iSess} = sparsityEMGs;
            sparsityValsShuff.(fieldNames{iType}){iSess} = sparsityEMGsShuff;

            [singleBehvSpec.(fieldNames{iType}){iSess}, allSingleBehvSpec.(fieldNames{iType}){iSess}, singleBehvSpecRegion.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveEMGs);
            [singleBehvSpecShuff.(fieldNames{iType}){iSess}, allSingleBehvSpecShuff.(fieldNames{iType}){iSess}, singleBehvSpecRegionShuff.(fieldNames{iType}){iSess}] = ...
                calcSingleBehvSpec2(regionAveEMGsShuff);
        end

        regionMaxFRs.(fieldNames{iType}){iSess} = max(regionAveVals.(fieldNames{iType}){iSess},[],2)/10;

        multiBehvSpec.(fieldNames{iType}){iSess} = calcMultiBehvSpec(regionAveVals.(fieldNames{iType}){iSess});
        multiBehvSpecShuff.(fieldNames{iType}){iSess} = calcMultiBehvSpec(regionAveValsShuff.(fieldNames{iType}){iSess});

        % correct for low firing rate neurons
        allSingleBehvSpecCorrected.(fieldNames{iType}){iSess} = allSingleBehvSpec.(fieldNames{iType}){iSess} - squeeze(nanmean(allSingleBehvSpecShuff.(fieldNames{iType}){iSess}));
        singleBehvSpecCorrected.(fieldNames{iType}){iSess} = singleBehvSpec.(fieldNames{iType}){iSess} - nanmean(singleBehvSpecShuff.(fieldNames{iType}){iSess});
        allMultiBehvSpecCorrected.(fieldNames{iType}){iSess} = multiBehvSpec.(fieldNames{iType}){iSess} - nanmean(multiBehvSpecShuff.(fieldNames{iType}){iSess},2);

        %calc CDF for the multibehavioral selectivity
        multiSpecCdfVals = 0:0.01:1;
        multispecToUse.(fieldNames{iType}){iSess} = multiBehvSpec.(fieldNames{iType}){iSess}(regionMaxFRs.(fieldNames{iType}){iSess}>0.2);
%         multispecToUse.(fieldNames{iType}){iSess} = allMultiBehvSpecCorrected.(fieldNames{iType}){iSess};
        multiSpecCdfFreq.(fieldNames{iType}){iSess} = calcCDF(multispecToUse.(fieldNames{iType}){iSess},multiSpecCdfVals);

        % get multibehavior spec shuffle values
        catMultiBehvSpecShuff.(fieldNames{iType}){iSess} = reshape(multiBehvSpecShuff.(fieldNames{iType}){iSess}(regionMaxFRs.(fieldNames{iType}){iSess}>0.2,:), 1, []);

        % also get CDF for sparsity metric
        sparsityCdfVals = 0:0.01:1;
        sparsityToUse.(fieldNames{iType}){iSess} = sparsityVals.(fieldNames{iType}){iSess}(regionMaxFRs.(fieldNames{iType}){iSess}>0.2);
        sparsityCdfFreq.(fieldNames{iType}){iSess} = calcCDF(sparsityToUse.(fieldNames{iType}){iSess},sparsityCdfVals);

        % get values of the shuffle/permutation for sparsity too
        catSparsityShuff.(fieldNames{iType}){iSess} = reshape(sparsityValsShuff.(fieldNames{iType}){iSess}(regionMaxFRs.(fieldNames{iType}){iSess}>0.2,:), 1, []);

        %as an alternative for single behavior spec, use statistical
        %approach to determine modulation of behavior
        shuffMean = mean(regionAveValsShuff.(fieldNames{iType}){iSess},3);
        shuffStd = std(regionAveValsShuff.(fieldNames{iType}){iSess},[],3);
        modulation.(fieldNames{iType}){iSess} = regionAveVals.(fieldNames{iType}){iSess} >= shuffMean + 3*shuffStd;
        for iRegion = 1:7
            regionSpec.(fieldNames{iType}){iSess}(:,iRegion) = modulation.(fieldNames{iType}){iSess}(:,iRegion) & ~any(modulation.(fieldNames{iType}){iSess}(:,setdiff(1:7,iRegion)),2);
        end
        
    end

end


% % make example plots for single neurons
% exampleSession = 3;
% exampleNeurons = [27 160 89 83];
% exampleNeuronLabels = {'Striatum Neuron 27','Cortex Neuron 35', 'Striatum Neuron 89','Striatum Neuron 83'};
% 
% baseDir = recordingSessions{exampleSession};
% load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','watershedRegions','regionWatershedLabels','behvLabelsNoArt')
% 
% frFile = 'NeuralFiringRates1msBins10msGauss.mat';
% activityType = 'firingrate';
% plotData = activityUMAPOverlay(baseDir, frFile, activityType,[],1);
% 
% nGridPoints = length(gridXInds);
% densityGaussStd = 0.2;
% [meshGridX,meshGridY] = meshgrid(gridXInds,gridYInds);
% 
% reducLimsX1 = min(min(reduction(:,1)))*1.5-2;
% reducLimsX2 = max(max(reduction(:,1)))*1.5+2;
% reducLimsY1 = min(min(reduction(:,2)))*1.5-2;
% reducLimsY2 = max(max(reduction(:,2)))*1.5+2;
% 
% [~,~,density] = findPointDensity(reduction(behvLabelsNoArt ~= 0,:),densityGaussStd,nGridPoints,[reducLimsX1 reducLimsX2 reducLimsY1 reducLimsY2]);
% 
% tileH = tiledlayout(2,2);
% tileH.TileSpacing = 'tight';
% tileH.Padding = 'tight';
% 
% for iEx = 1:length(exampleNeurons)
%     nexttile
%     imH = imagesc(gridXInds,gridYInds,plotData(:,:,exampleNeurons(iEx))');
%     alphaData = logisticTransparency(density,0.0005,0.0001);
%     imH.AlphaData = alphaData;
%     [boundaryYs,boundaryXs] = find(watershedRegions==0);
%     hold on
%     plot(gridXInds(boundaryXs),gridYInds(boundaryYs),'k.')
%     set(gca,'ydir','normal')
%     axis off
%     colorH = colorbar;
%     if iEx == 2 || iEx == 4
%         ylabel(colorH,'Firing Rate (Hz)');
%     end
%     xlim([-2.75   10.60])
%     ylim([-7.90    5.00])
%     title(exampleNeuronLabels{iEx})
% end
% 
% set(gcf,'color','w')

% make single behavior modulation and speicificty plots
ctxFigH = figure;
hold on;
set(gca,'XTick',1.5:3:21)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'LineWidth',1)
set(gca,'tickdir','out')
set(gca,'fontsize',12)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'color','w')
ylabel('Fraction of Population')
title('Cortex')
strFigH = figure;
hold on;
set(gca,'XTick',1.5:3:21)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'LineWidth',1)
set(gca,'tickdir','out')
set(gca,'fontsize',12)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'color','w')
ylabel('Fraction of Population')
title('Striatum')

sessPlotJitter = [0 0 0];%[-0.1 0 0.1];
plotColor = lines(4);
for iSess = 1:length(recordingSessions)

    strSpecFraction(iSess,:) = sum(regionSpec.str{iSess}(:,behvAlignPerm(iSess,:)))/size(modulation.str{iSess},1);
    strModFraction(iSess,:) = sum(modulation.str{iSess}(:,behvAlignPerm(iSess,:)))/size(modulation.str{iSess},1);

    ctxSpecFraction(iSess,:) = sum(regionSpec.ctx{iSess}(:,behvAlignPerm(iSess,:)))/size(modulation.ctx{iSess},1);
    ctxModFraction(iSess,:) = sum(modulation.ctx{iSess}(:,behvAlignPerm(iSess,:)))/size(modulation.ctx{iSess},1);

    figure(strFigH)
    plot((1:3:21)+sessPlotJitter(iSess),strModFraction(iSess,:),'.','MarkerSize',10,'color',plotColor(4,:))
    plot((2:3:21)+sessPlotJitter(iSess),strSpecFraction(iSess,:),'.','MarkerSize',10,'color',plotColor(3,:)*0.9)

    figure(ctxFigH)
    plot((1:3:21)+sessPlotJitter(iSess),ctxModFraction(iSess,:),'.','MarkerSize',10,'color',plotColor(4,:))
    plot((2:3:21)+sessPlotJitter(iSess),ctxSpecFraction(iSess,:),'.','MarkerSize',10,'color',plotColor(3,:)*0.9)
end

figure(strFigH)
barH = bar(1:3:21,mean(strModFraction),0.333,'EdgeColor','none','FaceColor',plotColor(4,:)*1.2);
barH.FaceAlpha = 0.5;
barH = bar(2:3:21,mean(strSpecFraction),0.333,'EdgeColor','none','FaceColor',plotColor(3,:)*1.05);
barH.FaceAlpha = 0.5;

legend('Behavior Modulated','Single Behavior Specific','box','off')

figure(ctxFigH)
barH = bar(1:3:21,mean(ctxModFraction),0.333,'EdgeColor','none','FaceColor',plotColor(4,:)*1.2);
barH.FaceAlpha = 0.5;
barH = bar(2:3:21,mean(ctxSpecFraction),0.333,'EdgeColor','none','FaceColor',plotColor(3,:)*1.05);
barH.FaceAlpha = 0.5;

legend('Behavior Modulated','Single Behavior Specific','box','off')

% do t-ttests and show stats about modulation and specificty for ctx vs
% striatum
[~, pMod] = ttest(strModFraction(:),ctxModFraction(:));
[~, pSpec] = ttest(sum(strSpecFraction,2),sum(ctxSpecFraction,2));

disp(['Mean fraction mod M1 neurons: ' num2str(mean(ctxModFraction(:))) ', str neurons: ' num2str(mean(strModFraction(:))),...
    ', t-test p = ' num2str(pMod)])

disp(['Total spec M1 neurons: ' num2str(mean(sum(ctxSpecFraction,2))) ', str neurons: ' num2str(mean(sum(strSpecFraction,2))),...
    ', t-test p = ' num2str(pSpec)])

% also look at the number of behaviors each neuron was modulated by, test for
% difference between cortex and striatum
aveStrModBehvs = cellfun(@(x) mean(sum(x,2)), modulation.str);
aveCtxModBehvs = cellfun(@(x) mean(sum(x,2)), modulation.ctx);

[~, pNModPerCell] = ttest(aveStrModBehvs,aveCtxModBehvs);

disp(['Mean # behaviors mod per cell for M1: ' num2str(mean(aveCtxModBehvs)) ', str: ' num2str(mean(aveStrModBehvs)),...
    ', t-test p = ' num2str(pNModPerCell)])

% make summary plot for multibehavior selectivity
allSessStrSpec = cat(1,multiSpecCdfFreq.str{:});
allSessCtxSpec = cat(1,multiSpecCdfFreq.ctx{:});
allSessEmgSpec = cat(1,multiSpecCdfFreq.emg{:});

plotColors = lines(4);
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


catSessStrSpec = calcCDF(cat(1,multispecToUse.str{:}),multiSpecCdfVals);
firstHitMax = find(catSessStrSpec==1,1);
legH(1) = plot(multiSpecCdfVals(1:firstHitMax),catSessStrSpec(1:firstHitMax),'Color',plotColors(2,:),'LineWidth',3);

catSessCtxSpec = calcCDF(cat(1,multispecToUse.ctx{:}),multiSpecCdfVals);
firstHitMax = find(catSessCtxSpec==1,1);
legH(2) = plot(multiSpecCdfVals(1:firstHitMax),catSessCtxSpec(1:firstHitMax),'Color',plotColors(1,:),'LineWidth',3);

catSessEmgSpec = calcCDF(cat(1,multispecToUse.emg{:}),multiSpecCdfVals);
firstHitMax = find(catSessEmgSpec==1,1);
legH(3) = plot(multiSpecCdfVals(1:firstHitMax),catSessEmgSpec(1:firstHitMax),'Color',plotColors(3,:),'LineWidth',3);

%get shuff mean, std, and 95% confidence interval
strShuffMean = nanmean(cat(2,catMultiBehvSpecShuff.str{:}));
strShuffStd = nanstd(cat(2,catMultiBehvSpecShuff.str{:}));
strShuffLowerConf = prctile(cat(2,catMultiBehvSpecShuff.str{:}),2.5);
strShuffUpperConf = prctile(cat(2,catMultiBehvSpecShuff.str{:}),97.5);

%use confidence interval rather than mean + std
patch([strShuffLowerConf, strShuffLowerConf, strShuffUpperConf, strShuffUpperConf],...
    [0.9, 1, 1, 0.9], [1 1 1 1],'facecolor',plotColors(2,:),'EdgeColor','none','facealpha',1)

ctxShuffMean = nanmean(cat(2,catMultiBehvSpecShuff.ctx{:}));
ctxShuffStd = nanstd(cat(2,catMultiBehvSpecShuff.ctx{:}));
ctxShuffLowerConf = prctile(cat(2,catMultiBehvSpecShuff.ctx{:}),2.5);
ctxShuffUpperConf = prctile(cat(2,catMultiBehvSpecShuff.ctx{:}),97.5);

patch([ctxShuffLowerConf, ctxShuffLowerConf, ctxShuffUpperConf, ctxShuffUpperConf],...
    [0.8, 0.9, 0.9, 0.8], [1 1 1 1],'facecolor',plotColors(1,:),'EdgeColor','none','facealpha',1)

emgShuffMean = nanmean(cat(2,catMultiBehvSpecShuff.emg{:}));
emgShuffStd = nanstd(cat(2,catMultiBehvSpecShuff.emg{:}));
emgShuffLowerConf = prctile(cat(2,catMultiBehvSpecShuff.emg{:}),2.5);
emgShuffUpperConf = prctile(cat(2,catMultiBehvSpecShuff.emg{:}),97.5);

patch([emgShuffLowerConf, emgShuffLowerConf, emgShuffUpperConf, emgShuffUpperConf],...
    [0.7, 0.8, 0.8, 0.7], [1 1 1 1],'facecolor',plotColors(3,:),'EdgeColor','none','facealpha',1)


legend(legH,'Striatum','Cortex','Muscles','box','off','fontsize',13)
xlabel('Adjusted Skew Metric')
ylabel('Frequency')
box off
set(gcf,'color','w')
set(gca,'LineWidth',1)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'TickLength',[0.02 0.05])
set(gca,'fontsize',13)
set(gca,'TickDir','out')
xlim([0.3 1]);

% do t-ttests for bias metric
[~, pBias] = kstest2(cat(1,multispecToUse.str{:}),cat(1,multispecToUse.ctx{:}));

disp(['Mean bias metric M1 neurons: ' num2str(mean(cat(1,multispecToUse.ctx{:}))) ', str neurons: ' num2str(mean(cat(1,multispecToUse.str{:}))),...
    ', ks-test p = ' num2str(pBias)])

% also do it for sparsity metric
% make summary plot for multibehavior selectivity
allSessStrSpar = cat(1,sparsityCdfFreq.str{:});
allSessCtxSpar = cat(1,sparsityCdfFreq.ctx{:});
allSessEmgSpar = cat(1,sparsityCdfFreq.emg{:});

plotColors = lines(4);
figure;
hold on;
for iSess = 1:length(recordingSessions)
    firstHitMax = find(allSessStrSpar(iSess,:)==1,1);
    plot(sparsityCdfVals(1:firstHitMax),allSessStrSpar(iSess,1:firstHitMax),'Color',[plotColors(2,:) 0.2],'LineWidth',1.5);

    firstHitMax = find(allSessCtxSpar(iSess,:)==1,1);
    plot(sparsityCdfVals(1:firstHitMax),allSessCtxSpar(iSess,1:firstHitMax),'Color',[plotColors(1,:) 0.2],'LineWidth',1.5);

    firstHitMax = find(allSessEmgSpar(iSess,:)==1,1);
    plot(sparsityCdfVals(1:firstHitMax),allSessEmgSpar(iSess,1:firstHitMax),'Color',[plotColors(3,:) 0.2],'LineWidth',1.5);
end


catSessStrSpar = calcCDF(cat(2,sparsityToUse.str{:}),sparsityCdfVals);
firstHitMax = find(catSessStrSpar==1,1);
legH(1) = plot(sparsityCdfVals(1:firstHitMax),catSessStrSpar(1:firstHitMax),'Color',plotColors(2,:),'LineWidth',3);

catSessCtxSpar = calcCDF(cat(2,sparsityToUse.ctx{:}),sparsityCdfVals);
firstHitMax = find(catSessCtxSpar==1,1);
legH(2) = plot(sparsityCdfVals(1:firstHitMax),catSessCtxSpar(1:firstHitMax),'Color',plotColors(1,:),'LineWidth',3);

catSessEmgSpar = calcCDF(cat(2,sparsityToUse.emg{:}),sparsityCdfVals);
firstHitMax = find(catSessEmgSpar==1,1);
legH(3) = plot(sparsityCdfVals(1:firstHitMax),catSessEmgSpar(1:firstHitMax),'Color',plotColors(3,:),'LineWidth',3);

%get shuff mean, std, and 95% confidence interval
strShuffMean = nanmean(cat(2,catSparsityShuff.str{:}));
strShuffStd = nanstd(cat(2,catSparsityShuff.str{:}));
strShuffLowerConf = prctile(cat(2,catSparsityShuff.str{:}),2.5);
strShuffUpperConf = prctile(cat(2,catSparsityShuff.str{:}),97.5);

%use confidence interval rather than mean + std
patch([strShuffLowerConf, strShuffLowerConf, strShuffUpperConf, strShuffUpperConf],...
    [0, 0.1, 0.1, 0], [1 1 1 1],'facecolor',plotColors(2,:),'EdgeColor','none','facealpha',1)

ctxShuffMean = nanmean(cat(2,catSparsityShuff.ctx{:}));
ctxShuffStd = nanstd(cat(2,catSparsityShuff.ctx{:}));
ctxShuffLowerConf = prctile(cat(2,catSparsityShuff.ctx{:}),2.5);
ctxShuffUpperConf = prctile(cat(2,catSparsityShuff.ctx{:}),97.5);

patch([ctxShuffLowerConf, ctxShuffLowerConf, ctxShuffUpperConf, ctxShuffUpperConf],...
    [0.1, 0.2, 0.2, 0.1], [1 1 1 1],'facecolor',plotColors(1,:),'EdgeColor','none','facealpha',1)

emgShuffMean = nanmean(cat(2,catSparsityShuff.emg{:}));
emgShuffStd = nanstd(cat(2,catSparsityShuff.emg{:}));
emgShuffLowerConf = prctile(cat(2,catSparsityShuff.emg{:}),2.5);
emgShuffUpperConf = prctile(cat(2,catSparsityShuff.emg{:}),97.5);

patch([emgShuffLowerConf, emgShuffLowerConf, emgShuffUpperConf, emgShuffUpperConf],...
    [0.2, 0.3, 0.3, 0.2], [1 1 1 1],'facecolor',plotColors(3,:),'EdgeColor','none','facealpha',1)

legend(legH,'Striatum','Cortex','Muscles','box','off','fontsize',13)
xlabel('Sparsity')
ylabel('Frequency')
box off
set(gcf,'color','w')
set(gca,'LineWidth',1)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'TickLength',[0.02 0.05])
set(gca,'fontsize',13)
set(gca,'TickDir','out')
xlim([0 1]);

% do t-ttests for sparsity
[~, pSpars] = kstest2(cat(2,sparsityVals.str{:}),cat(2,sparsityVals.ctx{:}));

disp(['Mean sparsity M1 neurons: ' num2str(nanmean(cat(2,sparsityVals.ctx{:}))) ', str neurons: ' num2str(mean(cat(2,sparsityVals.str{:}))),...
    ', ks-test p = ' num2str(pSpars)])

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

plotXOffsets = [0 -0.06 0.06];

% align behaviors across sessions
singleBehvSpecRegionAligned = singleBehvSpecRegion;
for iSess = 1:length(recordingSessions)
    strAligned{iSess} = regionAveVals.str{iSess}(:,behvAlignPerm(iSess,:));
    ctxAligned{iSess} = regionAveVals.ctx{iSess}(:,behvAlignPerm(iSess,:));
    strShuffAligned{iSess} = regionAveValsShuff.str{iSess}(:,behvAlignPerm(iSess,:),:);
    ctxShuffAligned{iSess} = regionAveValsShuff.ctx{iSess}(:,behvAlignPerm(iSess,:),:);
    plotXVals(iSess,:) = (1:2)+plotXOffsets(iSess);

    %align for single behavior regions
    for iBehvRegion = 1:size(behvAlignPerm,2)
        singleBehvSpecRegionAligned.str{iSess}(singleBehvSpecRegion.str{iSess} == iBehvRegion) = behvAlignPerm(iSess,iBehvRegion);
        singleBehvSpecRegionAligned.ctx{iSess}(singleBehvSpecRegion.ctx{iSess} == iBehvRegion) = behvAlignPerm(iSess,iBehvRegion);
        singleBehvSpecRegionAligned.emg{iSess}(singleBehvSpecRegion.emg{iSess} == iBehvRegion) = behvAlignPerm(iSess,iBehvRegion);
    end
end


strAllLowNeurons = cat(1,strAligned{:});
strAll = strAllLowNeurons;
strNormAll = strAll./(cat(1,stdStrFRs{:})*1000);
strShuff = cat(1,strShuffAligned{:});
strNormShuff = strShuff./(cat(1,stdStrFRs{:})*1000);
strMultiSpecs = cat(1,multiBehvSpec.str{:});
strMultiSpecs = strMultiSpecs(max(strAll,[],2)>0.2);
strAll = strAll(max(strAll,[],2)>0.2,:);
strNormAll = strNormAll(max(strAll,[],2)>0.2,:);
strNormAll(any(isnan(strNormAll),2),:) = [];
strNormShuff = strNormShuff(max(strAll,[],2)>0.2,:,:);
strNormShuff(any(any(isnan(strNormShuff),2),3),:,:) = [];

maxSortPermStr = sortByLevelRecursive(strAll);
maxSortPermStr = sortByLevelAndSpec(strAll,strMultiSpecs);
plotDataStr = strAll(maxSortPermStr,:)./sum(strAll(maxSortPermStr,:),2);
plotDataStr(any(isnan(plotDataStr),2),:) = [];

ctxAllLowNeurons = cat(1,ctxAligned{:});
ctxAll = ctxAllLowNeurons;
ctxNormAll = ctxAll./(cat(1,stdCtxFRs{:})*1000);
ctxShuff = cat(1,ctxShuffAligned{:});
ctxNormShuff = ctxShuff./(cat(1,stdCtxFRs{:})*1000);
ctxMultiSpecs = cat(1,multiBehvSpec.ctx{:});
ctxMultiSpecs = ctxMultiSpecs(max(ctxAll,[],2)>0.2);
ctxAll = ctxAll(max(ctxAll,[],2)>0.2,:);
ctxNormAll = ctxNormAll(max(ctxAll,[],2)>0.2,:);
ctxNormAll(any(isnan(ctxNormAll),2),:) = [];
ctxNormShuff = ctxNormShuff(max(ctxAll,[],2)>0.2,:,:);
ctxNormShuff(any(any(isnan(ctxNormShuff),2),3),:,:) = [];

% maxSortPermCtx = sortByLevelRecursive(ctxAll);
maxSortPermCtx = sortByLevelAndSpec(ctxAll,ctxMultiSpecs);
plotDataCtx = ctxAll(maxSortPermCtx,:)./sum(ctxAll(maxSortPermCtx,:),2);
plotDataCtx(any(isnan(plotDataCtx),2),:) = [];

%Get single behavior specificity
singleSpecStr = cat(2,singleBehvSpec.str{:});
singleSpecRegionStr = cat(2,singleBehvSpecRegionAligned.str{:});
singleSpecStr = singleSpecStr(max(strAllLowNeurons,[],2)>0.2);
singleSpecRegionStr = singleSpecRegionStr(max(strAllLowNeurons,[],2)>0.2);
highSpecInds = find(singleSpecStr > 0.67);
strHighSpecRegions = singleSpecRegionStr(highSpecInds);

singleSpecCtx = cat(2,singleBehvSpec.ctx{:});
singleSpecRegionCtx = cat(2,singleBehvSpecRegionAligned.ctx{:});
singleSpecCtx = singleSpecCtx(max(ctxAllLowNeurons,[],2)>0.2);
singleSpecRegionCtx = singleSpecRegionCtx(max(ctxAllLowNeurons,[],2)>0.2);
highSpecInds = find(singleSpecCtx > 0.67);
ctxHighSpecRegions = singleSpecRegionCtx(highSpecInds);

% make histograms of single behaviors specificty distribution
figure;
yyaxis("left")
histogram(strHighSpecRegions*3,'EdgeColor','none','LineWidth',2)
ylabel('# of striatum neurons')
yyaxis("right")
histogram(ctxHighSpecRegions*3+1,'EdgeColor','none','LineWidth',2)
ylabel('# of cortex neurons')

set(gca,'XTick',(1:7)*3+1)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'XTickLabelRotation',45)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'fontsize',13)
set(gca,'TickDir','out')
box off


% make example plot of all neurons sorted by region
figure
imagesc([plotDataStr; plotDataCtx])
colorH = colorbar;
ylabel(colorH,'Normalized Average Firing Rate');
hold on
line([0 7.5],[size(plotDataStr,1) size(plotDataStr,1)],'linewidth',2,'color','r','linestyle','--')

title('Behavioral Selectivities')
ylabel('Neuron')
set(gca,'XTick',1:7)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'XTickLabelRotation',45)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'fontsize',13)
set(gca,'TickDir','out')
box off


% look correlation across behaviors

% regionAveFRsCtx = strAll(length(striatumInds)+1:end,:);
% regionAveFRsCtx(any(isnan(regionAveFRsCtx),2) | any(isinf(regionAveFRsCtx),2),:) = [];
regionAveFRsCtx = ctxNormAll;%./max(ctxAll,[],2);
regionAveFRsCtxShuff = ctxNormShuff;
% regionAveFRsStr = regionAveFRs(1:length(striatumInds),:);
% regionAveFRsStr(any(isnan(regionAveFRsStr),2) | any(isinf(regionAveFRsStr),2),:) = [];
regionAveFRsStr = strNormAll;%./max(strAll,[],2);
regionAveFRsStrShuff = strNormShuff;


for iRegion1 = 1:nRegions
    for iRegion2 = 1:nRegions


%         % first do all neurons
%         linFit = fitlm([regionAveFRsStr(:,iRegion1); regionAveFRsCtx(:,iRegion1)],...
%             [regionAveFRsStr(:,iRegion2); regionAveFRsCtx(:,iRegion2)]);
%         frsRegionFitR2{1}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
%         frsRegionFitSlope{1}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);
%         frsRegionCC{1}(iRegion1,iRegion2) = corr([regionAveFRsStr(:,iRegion1); regionAveFRsCtx(:,iRegion1)],...
%             [regionAveFRsStr(:,iRegion2); regionAveFRsCtx(:,iRegion2)]);


        %next, just striatum
        linFit = fitlm(regionAveFRsStr(:,iRegion1), regionAveFRsStr(:,iRegion2));
        frsRegionFitR2{2}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        frsRegionFitSlope{2}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);
        frsRegionCC{2}(iRegion1,iRegion2) = corr(regionAveFRsStr(:,iRegion1), regionAveFRsStr(:,iRegion2));


        %next, just cortex
        linFit = fitlm(regionAveFRsCtx(:,iRegion1), regionAveFRsCtx(:,iRegion2));
        frsRegionFitR2{3}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
        frsRegionFitSlope{3}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);
        frsRegionCC{3}(iRegion1,iRegion2) = corr(regionAveFRsCtx(:,iRegion1), regionAveFRsCtx(:,iRegion2));

        %do for sessions individually
        for iSess = 1:length(recordingSessions)

            sessStrNorm = strAligned{iSess}./(stdStrFRs{iSess}*1000);
            sessStrNorm(any(isnan(sessStrNorm),2),:) = [];
            sessCtxNorm = ctxAligned{iSess}./(stdCtxFRs{iSess}*1000);
            sessCtxNorm(any(isnan(sessCtxNorm),2),:) = [];

            %striatum
            linFit = fitlm(sessStrNorm(:,iRegion1), sessStrNorm(:,iRegion2));
            frsRegionStrSessFitR2{iSess}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
            frsRegionStrSessFitSlope{iSess}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);
            frsRegionStrSessCC{iSess}(iRegion1,iRegion2) = corr(sessStrNorm(:,iRegion1), sessStrNorm(:,iRegion2));

            %cortex
            linFit = fitlm(sessCtxNorm(:,iRegion1), sessCtxNorm(:,iRegion2));
            frsRegionCtxSessFitR2{iSess}(iRegion1,iRegion2) = linFit.Rsquared.Ordinary;
            frsRegionCtxSessFitSlope{iSess}(iRegion1,iRegion2) = linFit.Coefficients.Estimate(2);
            frsRegionCtxSessCC{iSess}(iRegion1,iRegion2) = corr(sessCtxNorm(:,iRegion1), sessCtxNorm(:,iRegion2));

        end

        %do shuffs
        for iShuff = 1:size(strNormShuff,3)

            linFit = fitlm([regionAveFRsStrShuff(:,iRegion1,iShuff); regionAveFRsCtxShuff(:,iRegion1,iShuff)],...
                [regionAveFRsStrShuff(:,iRegion2,iShuff); regionAveFRsCtxShuff(:,iRegion2,iShuff)]);
            frsRegionFitShuffR2{1}(iRegion1,iRegion2,iShuff) = linFit.Rsquared.Ordinary;
            frsRegionFitShuffSlope{1}(iRegion1,iRegion2,iShuff) = linFit.Coefficients.Estimate(2);
            frsRegionShuffCC{1}(iRegion1,iRegion2,iShuff) = corr([regionAveFRsStrShuff(:,iRegion1,iShuff); regionAveFRsCtxShuff(:,iRegion1,iShuff)],...
                [regionAveFRsStrShuff(:,iRegion2,iShuff); regionAveFRsCtxShuff(:,iRegion2,iShuff)]);

            %next, just striatum
            linFit = fitlm(regionAveFRsStrShuff(:,iRegion1,iShuff), regionAveFRsStrShuff(:,iRegion2,iShuff));
            frsRegionFitShuffR2{2}(iRegion1,iRegion2,iShuff) = linFit.Rsquared.Ordinary;
            frsRegionFitShuffSlope{2}(iRegion1,iRegion2,iShuff) = linFit.Coefficients.Estimate(2);
            frsRegionShuffCC{2}(iRegion1,iRegion2,iShuff) = corr(regionAveFRsStrShuff(:,iRegion1,iShuff), regionAveFRsStrShuff(:,iRegion2,iShuff));

            %next, just cortex
            linFit = fitlm(regionAveFRsCtxShuff(:,iRegion1,iShuff), regionAveFRsCtxShuff(:,iRegion2,iShuff));
            frsRegionFitShuffR2{3}(iRegion1,iRegion2,iShuff) = linFit.Rsquared.Ordinary;
            frsRegionFitShuffSlope{3}(iRegion1,iRegion2,iShuff) = linFit.Coefficients.Estimate(2);
            frsRegionShuffCC{3}(iRegion1,iRegion2,iShuff) = corr(regionAveFRsCtxShuff(:,iRegion1,iShuff), regionAveFRsCtxShuff(:,iRegion2,iShuff));

            %do each session individually
            for iSess = 1:length(recordingSessions)

                sessStrNorm = strShuffAligned{iSess}(:,:,iShuff)./(stdStrFRs{iSess}*1000);
                sessStrNorm(any(isnan(sessStrNorm),2),:) = [];
                sessCtxNorm = ctxShuffAligned{iSess}(:,:,iShuff)./(stdCtxFRs{iSess}*1000);
                sessCtxNorm(any(isnan(sessCtxNorm),2),:) = [];

                %striatum
                linFit = fitlm(sessStrNorm(:,iRegion1), sessStrNorm(:,iRegion2));
                frsRegionStrSessShuffFitR2{iSess}(iRegion1,iRegion2,iShuff) = linFit.Rsquared.Ordinary;
                frsRegionStrSessShuffFitSlope{iSess}(iRegion1,iRegion2,iShuff) = linFit.Coefficients.Estimate(2);
                frsRegionStrSessShuffCC{iSess}(iRegion1,iRegion2,iShuff) = corr(sessStrNorm(:,iRegion1), sessStrNorm(:,iRegion2));

                %cortex
                linFit = fitlm(sessCtxNorm(:,iRegion1), sessCtxNorm(:,iRegion2));
                frsRegionCtxSessShuffFitR2{iSess}(iRegion1,iRegion2,iShuff) = linFit.Rsquared.Ordinary;
                frsRegionCtxSessShuffFitSlope{iSess}(iRegion1,iRegion2,iShuff) = linFit.Coefficients.Estimate(2);
                frsRegionCtxSessShuffCC{iSess}(iRegion1,iRegion2,iShuff) = corr(sessCtxNorm(:,iRegion1), sessCtxNorm(:,iRegion2));

            end

        end

    end
end

cgObj = clustergram(frsRegionFitR2{2});
clustOrdering = cellfun(@(x) str2num(x),cgObj.RowLabels);

exampleSess = 3;
figure;
tileH = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
nexttile 
plot(strAligned{exampleSess}(:,1)./(stdStrFRs{exampleSess}*1000),strAligned{exampleSess}(:,2)./(stdStrFRs{exampleSess}*1000),...
    'o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(1,:))
linFit = fitlm(strAligned{exampleSess}(:,1)./(stdStrFRs{exampleSess}*1000), strAligned{exampleSess}(:,2)./(stdStrFRs{exampleSess}*1000));
xRange = get(gca,'XLim');
line([0 xRange(2)],[linFit.Coefficients.Estimate(1) linFit.Coefficients.Estimate(1)+linFit.Coefficients.Estimate(2)*xRange(2)],'linestyle','--','linewidth',1.5,'color','r')
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'fontsize',13)
set(gca,'TickDir','out')
box off
xlabel('Climb Up z-scored firing rate')
ylabel('Climb Down z-scored firing rate')
text(0.1, 1.2,['Correlation = ' num2str(frsRegionStrSessCC{exampleSess}(1,2))],'fontsize',14);
ylim([-0.02 1.3])
xlim([-0.02 1.3])

nexttile 
plot(strAligned{exampleSess}(:,1)./(stdStrFRs{exampleSess}*1000),strAligned{exampleSess}(:,7)./(stdStrFRs{exampleSess}*1000),...
    'o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(1,:))
linFit = fitlm(strAligned{exampleSess}(:,1)./(stdStrFRs{exampleSess}*1000), strAligned{exampleSess}(:,7)./(stdStrFRs{exampleSess}*1000));
xRange = get(gca,'XLim');
line([0 xRange(2)],[linFit.Coefficients.Estimate(1) linFit.Coefficients.Estimate(1)+linFit.Coefficients.Estimate(2)*xRange(2)],'linestyle','--','linewidth',1.5,'color','r')
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'fontsize',13)
set(gca,'TickDir','out')
box off
xlabel('Climb Up z-scored firing rate')
ylabel('Eating z-scored firing rate')
text(0.1, 1.2,['Correlation = ' num2str(frsRegionStrSessCC{exampleSess}(1,7))],'fontsize',14);
ylim([-0.02 1.3])
xlim([-0.02 1.3])

figure
tileH = tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
nexttile
imagesc(1-frsRegionStrSessCC{exampleSess})
box off
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
cH = colorbar;
cH.Label.String = 'Linear Fit R^2';
% caxis([0.4 1])
title('Striatum')

nexttile
imagesc(1-frsRegionCtxSessCC{exampleSess})
box off
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
cH = colorbar;
cH.Label.String = 'Linear Fit R^2';
% caxis([0.4 1])
title('Cortex')

set(gcf,'color','w')

figure('Units','normalized','Color','w','OuterPosition',[0.1 0.1 0.7 0.4]);

strDendAx = axes('Units','normalize','OuterPosition',[0.00 0.05 0.35 0.9]);
plotMatrix = 1-frsRegionStrSessCC{exampleSess};
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],strDendAx,{'Linewidth',2,'Color','k'});
% plot shuff
plotMatrix = 1-frsRegionStrSessShuffCC{exampleSess}(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,strDendAx,{'Linewidth',2,'color',[0.2 0.2 0.2 0.5]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
ylabel('1 - Correlation')
ylim([0 0.5])


ctxPlotsAx = axes('Units','normalize','OuterPosition',[0.30 0.05 0.35 0.9]);
plotMatrix = 1-frsRegionCtxSessCC{exampleSess};
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],ctxPlotsAx,{'Linewidth',2,'Color','k'});
% plot shuff
plotMatrix = 1-frsRegionCtxSessShuffCC{exampleSess}(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,ctxPlotsAx,{'Linewidth',2,'color',[0.2 0.2 0.2 0.5]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
ylabel('1 - Correlation')
ylim([0 0.5])

% Now do plot showing hiearchy
strShuffPlotsAx = axes('Units','normalize','OuterPosition',[0.65 0.05 0.27 0.9]);
controlPlotColors = lines(3);
plotXOffsets = [0 -0.06 0.06];
for iSess = 1:length(recordingSessions)

    strMat = frsRegionStrSessCC{iSess};
    strMat(logical(diag(ones(length(strMat),1),0))) = 0;
    strSpread(iSess) = std(squareform(strMat));
    ctxMat = frsRegionCtxSessCC{iSess};
    ctxMat(logical(diag(ones(length(ctxMat),1),0))) = 0;
    ctxSpread(iSess) = std(squareform(ctxMat));
    for iShuff = 1:size(strNormShuff,3)
        strMat = frsRegionStrSessShuffCC{iSess}(:,:,iShuff);
        strMat(logical(diag(ones(length(strMat),1),0))) = 0;
        strSpreadShuff(iSess,iShuff) = std(squareform(strMat));
        ctxMat = frsRegionCtxSessShuffCC{iSess}(:,:,iShuff);
        ctxMat(logical(diag(ones(length(ctxMat),1),0))) = 0;
        ctxSpreadShuff(iSess,iShuff) = std(squareform(ctxMat));
    end

    plot(plotXVals(iSess,:),[ctxSpread(iSess) strSpread(iSess)],'o-','MarkerSize',5,'Color',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3])
    hold on

    %plot controls
    plotH = scatter(repmat(plotXVals(iSess,1),1,size(ctxSpreadShuff,2))+0.005,ctxSpreadShuff(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXVals(iSess,2),1,size(strSpreadShuff,2))+0.005,strSpreadShuff(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)

%     % add legend
%     legendH = legend('Random Rotation Control','Label Shuffle Control','Box','off','FontSize',11);
%     for iLabel = 1:length(legendH.String)
%         legendH.String{iLabel} = ['\color[rgb]{' num2str(controlPlotColors(iLabel,:)) '} ' legendH.String{iLabel}];
%     end

end

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Std dev. of correlation across behaviors')
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Striatum','Cortex'})
ylim([0 0.3])
xlim([0.5 2.5])
box off

% singleBehvSpecCorrected = singleBehvSpec - mean

% colorMap = jet(nRegions);
% 
% figure;
% randShifts = randn(1,length(depths))*20;
% rasterplot(cellfun(@(x) x/30000, ts,'un',0),'times','.',gca,[],depths+randShifts)
% 
% for iRegion = 1:nRegions
% 
%     hold on
% 
%     regionNeurons{iRegion} = find(singleBehvSpecRegion==iRegion & singleBehvSpec-nanmean(singleBehvSpecShuff) > 0.5);
%     
%     legH(iRegion) = plot([0 0],[1000 1001],'Color',colorMap(iRegion,:));
% 
%     if isempty(regionNeurons{iRegion})
%         continue
%     end
% 
%     rasterplot(cellfun(@(x) x/30000, ts(regionNeurons{iRegion}),'un',0),...
%         'times','.',gca,[],depths(regionNeurons{iRegion})+randShifts(regionNeurons{iRegion}),{'cData', colorMap(iRegion,:),'sizedata',30})
%     
% 
% end
% 
% legNames = {'Eating','Grooming','Still/Rear','Walk/Jump','Climb Down','ClimbUp'};
% legend(legH,legNames)
% ylim([0 4000])
% 
% 
% figure;
% tiledlayout(nRegions,1)
% for iRegion = 1:nRegions
%     nexttile
%     histH(iRegion) = histogram(depths(regionNeurons{iRegion}),0:100:4000,'FaceColor',colorMap(iRegion,:));
%     title(legNames{iRegion})
%     xlim([0 4000])
%     box off
%     ylabel('# Neurons')
% end
% 
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
multiBehvSpec = sqrt(squeeze(sum(regionMetricsNormalized.^2,2)));

end


function freq = calcCDF(inputData,range)

% remove nans
inputData(isnan(inputData)) = [];

for i = 1:length(range)

    freq(i) = sum(inputData <= range(i)) / length(inputData);

end

end



function sortPerm = sortByLevelRecursive(data)
% sort neural population by grouping each of the behavioral regions with highest firing
% rates together, and then subggrouping by the second highest firing rate, then thrird, ect
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



function sortPerm = sortByLevelAndSpec(data,spec)
% sort neural population by first grouping each of the behavioral regions
% with highest firing rate, then sort within each group by the
% multi-behavioral specificity
[~, maxRegion] = max(data,[],2);
maxSortTable = sortrows([maxRegion,(1:length(maxRegion))']);
sortPerm = maxSortTable(:,2);

for iRegion = 1:size(data,2)

    regionInds = find(maxSortTable(:,1)==iRegion);
    regionSpecs = spec(sortPerm(regionInds));
    [~, specSortPerm] = sort(regionSpecs);

    sortPermRegion = sortPerm(regionInds);
    sortPerm(regionInds) = sortPermRegion(specSortPerm);

end



end



% 
