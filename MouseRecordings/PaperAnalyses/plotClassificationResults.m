clear; close all

load('X:\David\AnalysesData\ClassificationAnalysis.mat')

behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rearing/Still','Groom','Eat'};

% first make confusion matrices
exampleSess = 1;
realClasses = cat(2,testLabels{exampleSess,:});
classLabels = unique(realClasses);
strPredClasses = cat(1,strRFPredLabels{exampleSess,:})';
ctxPredClasses = cat(1,ctxRFPredLabels{exampleSess,:})';

for iClass = 1:length(classLabels)

    realClassInds = find(realClasses==classLabels(iClass));
    thisClassStrPreds = strPredClasses(realClassInds);
    thisClassCtxPreds = ctxPredClasses(realClassInds);
    classNPoints(iClass) = length(realClassInds);

    for iClassCross = 1:length(classLabels)

        strConfMat(iClass,iClassCross) = sum(thisClassStrPreds==classLabels(iClassCross))/length(realClassInds);
        ctxConfMat(iClass,iClassCross) = sum(thisClassCtxPreds==classLabels(iClassCross))/length(realClassInds);

    end

end

figure('color','w')
tiledlayout(1,2,"TileSpacing",'tight','Padding','tight')

confMats{1} = strConfMat;
confMats{2} = ctxConfMat;

for iRegion = 1:length(confMats)

    nexttile
    imagesc(confMats{iRegion})
    colormap hot
    hold on
    for iClass = 1:length(classLabels)
        accuracy = num2str(round(confMats{iRegion}(iClass,iClass)*100)/100);
        if iClass==3
            textColor = 'w';
        else
            textColor = 'k';
        end
        text(iClass-0.45,iClass,accuracy,'Color',textColor)
    end
    box off
    clim([0 1])
    cH = colorbar;
    cH.Label.String = 'Accuracy';

    set(gca,'XTick',1:length(classLabels))
    set(gca,'YTick',1:length(classLabels))
    set(gca,'XTickLabel',behvRegionLabels)
    set(gca,'YTickLabel',behvRegionLabels)
    set(gca,'TickDir','out')
    set(gca,'FontSize',13)

    xlabel('Predicted Class')
    ylabel('Real Class')

end

% plot overall classification results
scatterData{1} = strAccuracy;
scatterData{2} = ctxAccuracy;

[plotH, barH] = barScatterPlot(scatterData,'none',ones(1,2),repmat({[-0.1 0 0.1]},1,2),[]);
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'Striatum','Cortex'})
xlim([0.5 2.5])

ylabel('Classifier Accuracy')


% do PCA histograms
tmp = cellfun(@(x) x(:),pcaProjCombRegionsStr,'UniformOutput',0);
strCombWeights = cat(1,tmp{:});

tmp = cellfun(@(x) x(:),pcaProjCombRegionsCtx,'UniformOutput',0);
ctxCombWeights = cat(1,tmp{:});

tmp = cellfun(@(x) x(:),pcaProjStr,'UniformOutput',0);
strIndividualWeights = cat(1,tmp{:});
        
tmp = cellfun(@(x) x(:),pcaProjCtx,'UniformOutput',0);
ctxIndividualWeights = cat(1,tmp{:});

strHistCounts = histcounts(strIndividualWeights,-1:0.01:1)/length(strIndividualWeights);
ctxHistCounts = histcounts(ctxIndividualWeights,-1:0.01:1)/length(ctxIndividualWeights);

plotColors = lines(2);

bH = bar(-1:0.01:0.99,ctxHistCounts,1);
bH.FaceAlpha = 0.6;
bH.FaceColor = plotColors(1,:);
hold on
bH = bar(-1:0.01:0.99,strHistCounts,1);
bH.FaceAlpha = 0.6;
bH.FaceColor = plotColors(2,:);

xlim([-0.3 0.3])
ylabel('Fraction of Weights')
xlabel('PCA Weight')
box off
set(gca,'LineWidth',1.5)
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'Color','w')

% calculate sparsity of the weights as the L1 norm of the weights
for iSess = 1:size(pcaProjStr)
    for iRegion = 1:size(pcaProjStr,2)

        nDimsStr(iSess,iRegion) = length(varExpStr{iSess,iRegion});% find(cumsum(varExpStr{iSess,iRegion})/sum(varExpStr{iSess,iRegion})>0.9,1);
        nDimsCtx(iSess,iRegion) = length(varExpCtx{iSess,iRegion});%find(cumsum(varExpCtx{iSess,iRegion})/sum(varExpCtx{iSess,iRegion})>0.9,1);

        directSparse = sum(abs(pcaProjStr{iSess,iRegion}(:,1:nDimsStr(iSess,iRegion))));
        strWeightsSpars{iSess,iRegion} = (directSparse-1)/(sqrt(length(varExpStr{iSess,iRegion}))-1);
        directSparse = sum(abs(pcaProjCtx{iSess,iRegion}(:,1:nDimsCtx(iSess,iRegion))));
        ctxWeightsSpars{iSess,iRegion} = (directSparse-1)/(sqrt(length(varExpCtx{iSess,iRegion}))-1);

    end
end

strWeightsSparsAll = cat(2,strWeightsSpars{:});
ctxWeightsSparsAll = cat(2,ctxWeightsSpars{:});

bar([mean(strWeightsSparsAll) mean(ctxWeightsSparsAll)])
hold on
errorbar([mean(strWeightsSparsAll) mean(ctxWeightsSparsAll)],[std(strWeightsSparsAll)/sqrt(length(strWeightsSparsAll)) std(ctxWeightsSparsAll)/sqrt(length(ctxWeightsSparsAll))],'.','Marker','none')
xlim([0 3])
ylim([0 1])

% 
