clear; close all

load('X:\David\AnalysesData\CCACorrelations.mat')

% plot example time series
plotColors = lines(3);
tileH = tiledlayout(2,1);
tileH.Padding = "tight"; tileH.TileSpacing = "tight";
plotH(1) = nexttile;
plot(ccaEMGTrajCtx{1}(100:1300,1),'color',[0.5 0.5 0.5],'linewidth',1.5);
hold on;
plot(ccaNeurTrajCtx{1}(100:1300,1),'color',plotColors(1,:),'linewidth',1.5);
line([0 200],[-3 -3],'color','k','linewidth',3)
axis off


plotH(2) = nexttile;
plot(ccaEMGTrajStr{1}(100:1300,1),'color',[0.5 0.5 0.5],'linewidth',2);
hold on;
plot(ccaNeurTrajStr{1}(100:1300,1),'color',plotColors(2,:),'linewidth',1.5);
axis off
linkaxes(plotH,'x');

set(gcf,'color','w')

% get the 95% interval for all the shift trials
ctxShiftsAllBehvCat = squeeze(permute(cannonCorrsCtxAllBehvShift,[4 3 2 1]));
ctxShiftsAllBehvCat = ctxShiftsAllBehvCat(:,:);

ctxShiftsAllBehvLowerPrct = prctile(ctxShiftsAllBehvCat,2.5,2);
ctxShiftsAllBehvUpperPrct = prctile(ctxShiftsAllBehvCat,97.5,2);

strShiftsAllBehvCat = squeeze(permute(cannonCorrsStrAllBehvShift,[4 3 2 1]));
strShiftsAllBehvCat = strShiftsAllBehvCat(:,:);

strShiftsAllBehvLowerPrct = prctile(strShiftsAllBehvCat,2.5,2);
strShiftsAllBehvUpperPrct = prctile(strShiftsAllBehvCat,97.5,2);

% plot the summary CCA corrlations for all animals for all behaviors, as a
% summery value, take the average across all 4 CC's
plotColors = lines(3);

% first do all behaviors combined
scatterData{1,1} = squeeze(mean(cannonCorrsCtxAllBehv,3));
scatterData{1,2} = squeeze(mean(cannonCorrsStrAllBehv,3));

plot(squeeze(cannonCorrsCtxAllBehv)','Color',[plotColors(1,:) 0.5],'LineWidth',1)
hold on
plot(squeeze(cannonCorrsStrAllBehv)','Color',[plotColors(2,:) 0.5],'LineWidth',1)
plot(mean(squeeze(cannonCorrsCtxAllBehv)),'.-','Color',[plotColors(1,:)],'LineWidth',2,'MarkerSize',15)
plot(mean(squeeze(cannonCorrsStrAllBehv)),'.-','Color',[plotColors(2,:)],'LineWidth',2,'MarkerSize',15)
box off
plot(1:4,ctxShiftsAllBehvLowerPrct,'--','Color',plotColors(1,:),'LineWidth',1.5)
plot(1:4,ctxShiftsAllBehvUpperPrct,'--','Color',plotColors(1,:),'LineWidth',1.5)
plot(1:4,strShiftsAllBehvLowerPrct,'--','Color',plotColors(2,:),'LineWidth',1.5)
plot(1:4,strShiftsAllBehvUpperPrct,'--','Color',plotColors(2,:),'LineWidth',1.5)
% shadedErrorBar(1:4,ctxShiftMean,ctxShiftRange,'transparent',1,'lineProps',{'color',plotColors(1,:),'linestyle','none'})
% shadedErrorBar(1:4,strShiftMean,strShiftRange,'transparent',1,'lineProps',{'color',plotColors(2,:),'linestyle','none'})
set(gca,'XTick',1:4)
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
set(gca,'TickDir','out')
help canoncorr
xlabel('Canonical variable')
ylabel('Canonical correlation')
set(gcf,'Color','w')
set(gca,'Color','k')
set(gca,'Color','w')
set(gca,'XColor','k')
set(gca,'YColor','k')

% make bar plots using averages
% next break out for each behavior individually
scatterData(2:size(cannonCorrsCtx,1)+1,1) = ...
    mat2cell(mean(cannonCorrsCtx,3),ones(size(cannonCorrsCtx,1),1),size(cannonCorrsCtx,2));
scatterData(2:size(cannonCorrsStr,1)+1,2) = ...
    mat2cell(mean(cannonCorrsStr,3),ones(size(cannonCorrsStr,1),1),size(cannonCorrsStr,2));

[plotH, barH] = barScatterPlot(scatterData,'none',ones(size(cannonCorrsCtx,1)+1,2),[],[1 2]);
ylabel('Mean Canonical Corr')
set(gca,'XTickLabel',{'All','Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'})
legendH = legend(barH, 'Cortex','Striatum','box','off');
for iLabel = 1:length(legendH.String)
    legendH.String{iLabel} = ['\color[rgb]{' num2str(barH(iLabel).FaceColor) '} ' legendH.String{iLabel}];
end




% 

