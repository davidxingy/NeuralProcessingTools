function lookDecoding

close all

sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData'};


allBehvAlignPerms = [...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

allMuscleAlignPerms = [...
    1 2 3 4 5 6 7 8; ...
    1 2 3 4 5 6 7 8; ...
    1 2 3 4 6 5 8 7 ...
    ];

behvRegionLabels = {'Climb Up','Climb Down','Jump Down/Across','Walk Flat/Grid','Rearing/Still','Grooming','Eating'};
muscleLabels = {'Right Biceps','Right Triceps','Right ECR','Right PL','Right TA','Right Gastr','Left Triceps','Left Biceps'};

for iSess = 1:length(sessionDirs)

    behvAlignPerm = allBehvAlignPerms(iSess,:);
    muscleAlignPerm = allMuscleAlignPerms(iSess,:);

    load(fullfile(sessionDirs{iSess},'singleRegionDecoding'))
    

    %get R2, CC, and MSE
    ctxR2Cell = cellfun(@(x) x.R2, singlePerformanceCtx,'UniformOutput',0);
    ctxR2(:,:,iSess) = [ctxR2Cell{behvAlignPerm}];
    ctxR2(:,:,iSess) = ctxR2(muscleAlignPerm,:,iSess);
    ctxCCCell = cellfun(@(x) x.CC, singlePerformanceCtx,'UniformOutput',0);
    ctxCC(:,:,iSess) = [ctxCCCell{behvAlignPerm}];
    ctxCC(:,:,iSess) = ctxCC(muscleAlignPerm,:,iSess);
    ctxMSECell = cellfun(@(x) x.MSE, singlePerformanceCtx,'UniformOutput',0);
    ctxMSE(:,:,iSess) = [ctxMSECell{behvAlignPerm}];
    ctxMSE(:,:,iSess) = ctxMSE(muscleAlignPerm,:,iSess);

    strR2Cell = cellfun(@(x) x.R2, singlePerformanceStr,'UniformOutput',0);
    strR2(:,:,iSess) = [strR2Cell{behvAlignPerm}];
    strR2(:,:,iSess) = strR2(muscleAlignPerm,:,iSess);
    strCCCell = cellfun(@(x) x.CC, singlePerformanceStr,'UniformOutput',0);
    strCC(:,:,iSess) = [strCCCell{behvAlignPerm}];
    strCC(:,:,iSess) = strCC(muscleAlignPerm,:,iSess);
    strMSECell = cellfun(@(x) x.MSE, singlePerformanceStr,'UniformOutput',0);
    strMSE(:,:,iSess) = [strMSECell{behvAlignPerm}];
    strMSE(:,:,iSess) = strMSE(muscleAlignPerm,:,iSess);
    
end

% make cortex vs stiratum figures
scatterData(:,1) = mat2cell(squeeze(mean(ctxR2(1:4,:,:),1)),ones(1,size(ctxR2,2)),size(ctxR2,3));
scatterData(:,2) = mat2cell(squeeze(mean(strR2(1:4,:,:),1)),ones(1,size(strR2,2)),size(strR2,3));

[plotH, barH] = barScatterPlot(scatterData,'none',ones(size(scatterData,1),2),[],[1 2]);
ylabel('Decoding R^2')
set(gca,'XTickLabel',{'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'})
legendH = legend(barH, 'Cortex','Striatum','box','off');
for iLabel = 1:length(legendH.String)
    legendH.String{iLabel} = ['\color[rgb]{' num2str(barH(iLabel).FaceColor) '} ' legendH.String{iLabel}];
end
ylim([-0.25 0.5])
set(gca,'YTick',-0.25:0.25:0.5)

muscleFigH = figure;
hold on

ctxR2Muscles = ctxR2(:,:);
strR2Muscles = strR2(:,:);

barH(1) = bar((1:8)*3,mean(ctxR2Muscles,2),'FaceColor',plotColors(1,:),'BarWidth',0.3);
errorbar((1:8)*3,mean(ctxR2Muscles,2),std(ctxR2Muscles,[],2)/sqrt(size(ctxR2Muscles,2)),'.','Color',[0.3 0.3 0.3])
barH(2) = bar((1:8)*3+1,mean(strR2Muscles,2),'FaceColor',plotColors(2,:),'BarWidth',0.3);
errorbar((1:8)*3+1,mean(strR2Muscles,2),std(strR2Muscles,[],2)/sqrt(size(strR2Muscles,2)),'.','Color',[0.3 0.3 0.3])

set(gca,'XTick',(1:8)*3+0.5)
set(gca,'XTickLabels',muscleLabels)
ylabel('R2')
legend(barH,'Cortex','Striatum')

%Make example figure
exampleSession = 2;
exampleRegion = 1;
exampleMuscle = 2;

exampleTimes = 28725:29175;

exampleFigureH = figure;
load(fullfile(sessionDirs{exampleSession},'singleRegionDecoding'))
hold on
samplePlotH(1) = plot(regionEMGHistNormalized{allBehvAlignPerms(exampleSession,exampleRegion)}(exampleMuscle,exampleTimes),'color',[0.2 0.2 0.2 0.5],'LineWidth',2.5);
samplePlotH(2) = plot(estEMGSingleCtx{allBehvAlignPerms(exampleSession,exampleRegion)}(exampleMuscle,exampleTimes),'color',plotColors(1,:),'LineWidth',2);
samplePlotH(3) = plot(estEMGSingleStr{allBehvAlignPerms(exampleSession,exampleRegion)}(exampleMuscle,exampleTimes),'color',plotColors(2,:),'LineWidth',2);
legend(samplePlotH,['Real ' muscleLabels{exampleMuscle}],'Cortex Decoded','Striatum Decoded','box','off','fontsize',12)

line([0 50],[-1 -1],'color','k','linewidth',2)
% ylabel('Normalized EMG')
title([behvRegionLabels{exampleRegion} ', Ctx R2 = ' num2str(singlePerformanceCtx{allBehvAlignPerms(exampleSession,exampleRegion)}.R2(exampleMuscle))...
    ', Str R2 = ' num2str(singlePerformanceStr{allBehvAlignPerms(exampleSession,exampleRegion)}.R2(exampleMuscle))])

axis off
set(gcf,'color','w')



% 
