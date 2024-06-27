clear
close all

sessionNames = {'pcaAlignment_D020.mat','pcaAlignment_D024.mat','pcaAlignment_D026.mat'};

plotXOffsets = [0 -0.08 0.08];

spreadFigH = figure;
hold on
meanFigH = figure;
hold on

exampleEMGSess = 3;
exampleNeurSess = 3;


inputData = 'umapregions';
metric = 'angle';

switch lower(inputData)
    case 'umapregions'
        load('X:\David\AnalysesData\PCASubspaces.mat')
        behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rearing/Still','Groom','Eat'};

        allBehvAlignPerms = [
            1 2 3 4 5 6 7; ...
            1 2 3 4 5 6 7; ...
            1 2 4 5 3 6 7 ...
            ];

    case 'humanannotated'
        load('X:\David\AnalysesData\PCASubspacesAnnotated.mat')
        behvRegionLabels = {'Climb Up','Climb Down','Jump Down','Walk Flat','Walk Grid','Rear','Groom','Eat'};

        allBehvAlignPerms = repmat(1:8,3,1);

end

behvAlignPerm = allBehvAlignPerms(exampleEMGSess,:);

switch lower(metric)
    case 'alignment'
        pcaMetricCtx = 1-pcaAlignCtx;
        pcaMetricStr = 1-pcaAlignStr;
        pcaMetricEmg = 1-pcaAlignEmg;
        pcaRotCtx = 1-pcaAlignRotCtx;
        pcaRotStr = 1-pcaAlignRotStr;
        pcaRotEmg = 1-pcaAlignRotEmg;
        pcaShiftCtx = 1-pcaAlignShiftCtx;
        pcaShiftStr = 1-pcaAlignShiftStr;
        pcaShiftEmg = 1-pcaAlignShiftEmg;
        metricLabel = 'Alignment Metric';
        plotMatrixInv = 1;

    case 'angle'
        pcaMetricCtx = pcaAngleCtx;
        pcaMetricStr = pcaAngleStr;
        pcaMetricEmg = pcaAngleEmg;
        pcaRotCtx = pcaAngleRotCtx;
        pcaRotStr = pcaAngleRotStr;
        pcaRotEmg = pcaAngleRotEmg;
        pcaShiftCtx = pcaAngleShiftCtx;
        pcaShiftStr = pcaAngleShiftStr;
        pcaShiftEmg = pcaAngleShiftEmg;
        metricLabel = 'Principal Angle';
        plotMatrixInv = 0;

end

for iSess = 1:length(sessionNames)

    %calc summary metrics for the random rotations controls
    for iShuff = 1:size(pcaRotCtx,4)

        ctxSpreadRot(iSess,iShuff) = std(squareform(squeeze(pcaRotCtx(iSess,:,:,iShuff))));
        ctxMeanRot(iSess,iShuff) = mean(squareform(squeeze(pcaRotCtx(iSess,:,:,iShuff))));
        strSpreadRot(iSess,iShuff) = std(squareform(squeeze(pcaRotStr(iSess,:,:,iShuff))));
        strMeanRot(iSess,iShuff) = mean(squareform(squeeze(pcaRotStr(iSess,:,:,iShuff))));
        emgSpreadRot(iSess,iShuff) = std(squareform(squeeze(pcaRotEmg(iSess,:,:,iShuff))));
        emgMeanRot(iSess,iShuff) = mean(squareform(squeeze(pcaRotEmg(iSess,:,:,iShuff))));

    end

    %calc summary metrics for the behavior shift controls
    for iShuff = 1:size(pcaShiftCtx,4)

        ctxSpreadShuff(iSess,iShuff) = std(squareform(squeeze(pcaShiftCtx(iSess,:,:,iShuff))));
        ctxMeanShuff(iSess,iShuff) = mean(squareform(squeeze(pcaShiftCtx(iSess,:,:,iShuff))));
        strSpreadShuff(iSess,iShuff) = std(squareform(squeeze(pcaShiftStr(iSess,:,:,iShuff))));
        strMeanShuff(iSess,iShuff) = mean(squareform(squeeze(pcaShiftStr(iSess,:,:,iShuff))));
        emgSpreadShuff(iSess,iShuff) = std(squareform(squeeze(pcaShiftEmg(iSess,:,:,iShuff))));
        emgMeanShuff(iSess,iShuff) = mean(squareform(squeeze(pcaShiftEmg(iSess,:,:,iShuff))));

    end

    ctxSpread(iSess) = std(squareform(squeeze(pcaMetricCtx(iSess,:,:))));
    ctxMean(iSess) = mean(squareform(squeeze(pcaMetricCtx(iSess,:,:))));

    strSpread(iSess) = std(squareform(squeeze(pcaMetricStr(iSess,:,:))));
    strMean(iSess) = mean(squareform(squeeze(pcaMetricStr(iSess,:,:))));

    emgSpread(iSess) = std(squareform(squeeze(pcaMetricEmg(iSess,:,:))));
    emgMean(iSess) = mean(squareform(squeeze(pcaMetricEmg(iSess,:,:))));

    allSessMetricsCtx{iSess} = squareform(squeeze(pcaMetricCtx(iSess,:,:)));
    allSessMetricsStr{iSess} = squareform(squeeze(pcaMetricStr(iSess,:,:)));
    allSessMetricsEmg{iSess} = squareform(squeeze(pcaMetricEmg(iSess,:,:)));

    plotXVals(iSess,:) = (1:2)+plotXOffsets(iSess);

end


% plot EMG Spread separately
exampleEMG = squeeze(pcaMetricEmg(exampleEMGSess,:,:));
exampleEMGShuff = squeeze(pcaShiftEmg(exampleEMGSess,:,:,1));


% plot example hiearchy
plotColors = lines(7);
% first do cortex
figure('Units','normalized','OuterPosition',[0.1 0.1 0.445 0.7])
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
plotTitle = 'Cortex';
if plotMatrixInv
    plotMatrix = 1-pcaMetricCtx(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:));
else
    plotMatrix = pcaMetricCtx(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:));
end
imagesc(squeeze(plotMatrix))
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabelRotation',30)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
cH = colorbar;
% caxis([0 0.75])
cH.Label.String = metricLabel;
title(plotTitle)

% next do striatum
nexttile
plotTitle = 'Striatum';
if plotMatrixInv
    plotMatrix = 1-pcaMetricStr(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:));
else
    plotMatrix = pcaMetricStr(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:));
end
imagesc(squeeze(plotMatrix))
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabelRotation',30)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
cH = colorbar;
% caxis([0 0.75])
cH.Label.String = metricLabel;
title(plotTitle)

% plot dendrograms
% first cortex
nexttile
plotMatrix = squeeze(pcaMetricCtx(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:)));
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});
% add an example shift control dendrogram
plotMatrix = squeeze(pcaShiftCtx(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:),1));
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(3,:) 0.7]})
plotMatrix = squeeze(pcaRotCtx(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:),1));
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(4,:) 0.7]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
ylabel(metricLabel)
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')
ylim([0 0.3])
axH = gca;

% next striatum
nexttile
plotMatrix = squeeze(pcaMetricStr(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:)));
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});
% add an example shift control dendrogram
plotMatrix = squeeze(pcaShiftStr(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:),1));
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(3,:) 0.7]})
plotMatrix = squeeze(pcaRotStr(exampleNeurSess,allBehvAlignPerms(exampleNeurSess,:),allBehvAlignPerms(exampleNeurSess,:),1));
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(4,:) 0.7]})


set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
ylabel(metricLabel)
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')
ylim([0 0.6])
axH = gca;

% make summary plots
plotJitter = [-0.3 0 0.3];
figure;
barScatterPlot({ctxSpread;strSpread},'none',[1; 1],repmat({plotJitter},2,1),[],true);
hold on

% also plot controls
nShifts = size(strSpreadShuff,2);
for iAnimal = 1:length(sessionNames)

    controlJitter = randn(1,nShifts)/100;
    
    % plot shuffle controls
    % first do cortex
    plotH = scatter(controlJitter+plotJitter(iAnimal)+1-0.05,ctxSpreadShuff(iAnimal,:),10,plotColors(3,:),'filled');
    alpha(plotH,0.2)
    % then striatum
    plotH = scatter(controlJitter+plotJitter(iAnimal)+3-0.05,strSpreadShuff(iAnimal,:),10,plotColors(3,:),'filled');
    alpha(plotH,0.2)

    % plot rotation controls
    % first do cortex
    plotH = scatter(controlJitter+plotJitter(iAnimal)+1+0.05,ctxSpreadRot(iAnimal,:),10,plotColors(4,:),'filled');
    alpha(plotH,0.2)
    % then striatum
    plotH = scatter(controlJitter+plotJitter(iAnimal)+3+0.05,strSpreadRot(iAnimal,:),10,plotColors(4,:),'filled');
    alpha(plotH,0.2)

end

xlim([0 4])
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylabel('Hierarchy spread (s.d.)')


controlPlotColors = lines(3);
for iSess = 1:length(sessionNames)

    %first plot spread
    figure(spreadFigH);
    plot(plotXVals(iSess,:),[ctxSpread(iSess) strSpread(iSess)],'o-','MarkerSize',5,'Color',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3])
    
    %plot controls
    plotH = scatter(repmat(plotXVals(iSess,1),1,size(ctxSpreadRot,2))+0.005,ctxSpreadRot(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXVals(iSess,2),1,size(strSpreadRot,2))+0.005,strSpreadRot(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
%     plotH = scatter(repmat(plotXVals(iSess,3),1,size(emgSpreadRot,2)),emgSpreadRot(iSess,:),10,controlPlotColors(1,:),'filled');
%     alpha(plotH,0.2)

    plotH = scatter(repmat(plotXVals(iSess,1),1,size(ctxSpreadShuff,2))-0.005,ctxSpreadShuff(iSess,:),10,controlPlotColors(2,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXVals(iSess,2),1,size(strSpreadShuff,2))-0.005,strSpreadShuff(iSess,:),10,controlPlotColors(2,:),'filled');
    alpha(plotH,0.2)
%     plotH = scatter(repmat(plotXVals(iSess,3),1,size(emgSpreadShuff,2)),emgSpreadShuff(iSess,:),10,controlPlotColors(2,:),'filled');
%     alpha(plotH,0.2)

    % add legend
    legendH = legend('Random Rotation Control','Label Shuffle Control','Box','off','FontSize',11);
    for iLabel = 1:length(legendH.String)
        legendH.String{iLabel} = ['\color[rgb]{' num2str(controlPlotColors(iLabel,:)) '} ' legendH.String{iLabel}];
    end


    %next plot means
    figure(meanFigH);
    plot(plotXVals(iSess,:),[ctxMean(iSess) strMean(iSess)],'o-','MarkerSize',5,'Color',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3])
    
    %plot controls
    plotH = scatter(repmat(plotXVals(iSess,1),1,size(ctxMeanRot,2))+0.005,ctxMeanRot(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXVals(iSess,2),1,size(strMeanRot,2))+0.005,strMeanRot(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
%     plotH = scatter(repmat(plotXVals(iSess,3),1,size(emgMeanRot,2)),emgMeanRot(iSess,:),10,controlPlotColors(1,:),'filled');
%     alpha(plotH,0.2)

    plotH = scatter(repmat(plotXVals(iSess,1),1,size(ctxMeanShuff,2))-0.005,ctxMeanShuff(iSess,:),10,controlPlotColors(2,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXVals(iSess,2),1,size(strMeanShuff,2))-0.005,strMeanShuff(iSess,:),10,controlPlotColors(2,:),'filled');
    alpha(plotH,0.2)
%     plotH = scatter(repmat(plotXVals(iSess,3),1,size(emgMeanShuff,2)),emgMeanShuff(iSess,:),10,controlPlotColors(2,:),'filled');
%     alpha(plotH,0.2)

    % add legend
    legendH = legend('Random Rotation Control','Label Shuffle Control','Box','off','FontSize',11);
    for iLabel = 1:length(legendH.String)
        legendH.String{iLabel} = ['\color[rgb]{' num2str(controlPlotColors(iLabel,:)) '} ' legendH.String{iLabel}];
    end

end


% make pretty
figure(spreadFigH)
set(gcf,'Color','w')
set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
ylabel(['Std dev. of ' metricLabel])
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylim([0 0.2])
xlim([0.5 2.5])

% make pretty
figure(meanFigH)
set(gcf,'Color','w')
set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'FontSize',14)
ylabel(['Mean. of ' metricLabel])
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylim([0 0.5])
xlim([0.5 2.5])



% plot EMG Spread separately
emgFigH = figure;
%plot spread and controls
plot(plotXOffsets, emgSpread,'o','MarkerSize',5,'Color',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3])
hold on
for iSess = 1:length(sessionNames)
    plotH = scatter(repmat(plotXOffsets(iSess),1,size(emgSpreadRot,2))+0.005,emgSpreadRot(iSess,:),10,controlPlotColors(1,:),'filled');
    alpha(plotH,0.2)
    plotH = scatter(repmat(plotXOffsets(iSess),1,size(emgSpreadRot,2))+0.005,emgSpreadShuff(iSess,:),10,controlPlotColors(2,:),'filled');
end

% add legend
legendH = legend('Random Rotation Control','Label Shuffle Control','Box','off','FontSize',11);
for iLabel = 1:length(legendH.String)
    legendH.String{iLabel} = ['\color[rgb]{' num2str(controlPlotColors(iLabel,:)) '} ' legendH.String{iLabel}];
end

% make pretty

figure(spreadFigH);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylabel(['Std dev of ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.4)
set(gcf,'Color','w')
xlim([0.5 2.5])

figure(meanFigH);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylabel(['Mean of ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.4)
set(gcf,'Color','w')
xlim([0.5 2.5])

figure(emgFigH);
box off
set(gca,'XTick',[0])
set(gca,'XTickLabel',{'EMG'})
ylabel(['Std dev of ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.4)
set(gcf,'Color','w')
xlim([-1 1])

% plot EMG hiearchy 
figure('Units','normalized','OuterPosition',[0.1 0.1 0.445 0.7])
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile
plotTitle = 'D026 EMG';
if plotMatrixInv
    plotMatrix = 1-exampleEMG(behvAlignPerm,behvAlignPerm,1);
else
    plotMatrix = exampleEMG(behvAlignPerm,behvAlignPerm,1);
end
imagesc(plotMatrix)
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabelRotation',30)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
cH = colorbar;
% caxis([0 0.55])
cH.Label.String = metricLabel;
title(plotTitle)

% plot dendrograms
nexttile
plotMatrix = exampleEMG(behvAlignPerm,behvAlignPerm,1);
% plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});
% plotMatrix = exampleEMGShuff(behvAlignPerm,behvAlignPerm,1);
% plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
% customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[lines(1) 0.3]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
ylabel(metricLabel)
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')
ylim([0 0.45])
axH = gca;



% plot EMG correlations
figure('Units','normalized','OuterPosition',[0.1 0.1 0.445 0.7])
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')

allEMGs = [allSessMetricsEmg{:}];
allStrs = [allSessMetricsStr{:}];
allCtxs = [allSessMetricsCtx{:}];

% plot emg vs cortex
nexttile
scatter(allEMGs,allCtxs,15,[0.3 0.3 0.3],'filled')
linReg = fitlm(allEMGs,allCtxs);
hold on
plot(get(gca,'XLim'),get(gca,'XLim')*linReg.Coefficients.Estimate(2)+linReg.Coefficients.Estimate(1),'Color',lines(1),'LineWidth',2,'LineStyle','--')
text(0.05,0.38,['R^2 = ' num2str(corr(allEMGs',allCtxs'))],'fontSize',14)
xlabel(['EMG ' metricLabel])
ylabel(['Cortex ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.5)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
xlim([0 0.8])

% plot emg vs striatum
nexttile
scatter(allEMGs,allStrs,15,[0.3 0.3 0.3],'filled')
linReg = fitlm(allEMGs,allStrs);
hold on
plot(get(gca,'XLim'),get(gca,'XLim')*linReg.Coefficients.Estimate(2)+linReg.Coefficients.Estimate(1),'Color',lines(1),'LineWidth',2,'LineStyle','--')
text(0.05,0.75,['R^2 = ' num2str(corr(allEMGs',allStrs'))],'fontSize',14)
xlabel(['EMG ' metricLabel])
ylabel(['Striatum ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.5)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
xlim([0 0.8])

% plot cortex vs striatum
nexttile
scatter(allCtxs,allStrs,15,[0.3 0.3 0.3],'filled')
linReg = fitlm(allCtxs,allStrs);
hold on
plot(get(gca,'XLim'),get(gca,'XLim')*linReg.Coefficients.Estimate(2)+linReg.Coefficients.Estimate(1),'Color',lines(1),'LineWidth',2,'LineStyle','--')
text(0.03,0.75,['R^2 = ' num2str(corr(allCtxs',allStrs'))],'fontSize',14)
xlabel(['Cortex ' metricLabel])
ylabel(['Striatum ' metricLabel])
set(gca,'TickDir','out')
set(gca,'FontSize',13)
set(gca,'LineWidth',1.5)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gcf,'Color','w')


% 
