clear

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

use10msBins = true;
nBehvShifts = 5;

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};

for iSess = 1:length(allDirs)
    
    load(fullfile(allDirs{iSess},'ProcessedData','EMG1ms.mat'))
    load(fullfile(allDirs{iSess},'ProcessedData','UMAP.mat'),'regionAssignmentsFiltered','reduction','origDownsampEMGInd',...
        'behvLabelsDown','analyzedBehaviors','regionWatershedLabels','regionAssignmentsFiltered','regionBehvAssignments')
    
    if use10msBins
        load(fullfile(allDirs{iSess},'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs','cortexInds','striatumInds')
    else
        load(fullfile(allDirs{iSess},'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'),'allFRs','cortexInds','striatumInds')
    end

    load(fullfile(allDirs{iSess},'ProcessedData','VideoSyncFrames.mat'))

    emg10ms = [];
    if use10msBins
        for iChan = 1:size(downsampEMG,1)
            avedSig = conv(repmat(0.1,1,10),downsampEMG(iChan,:));
            emg10ms(iChan,:) = avedSig(6:10:end);
        end
        downsampEMG = emg10ms;
    end

    usedEmgInds = (1:length(downsampEMG)-1000);
    currentDir = pwd;
    cd(fullfile(allDirs{iSess},'ProcessedData'))
    matchedNeurInds = round(NeurEMGSync(usedEmgInds*200,...
        frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
    cd(currentDir)

    % some emg time points don't have a corresponding neural sample point
    nanNeurInds = find(isnan(matchedNeurInds));
    matchedNeurInds(nanNeurInds) = [];
    usedEmgInds(nanNeurInds) = [];
    
    strData = allFRs(1:length(striatumInds),matchedNeurInds);
    ctxData = allFRs(length(striatumInds)+1:end,matchedNeurInds);
    emgData = downsampEMG(:,usedEmgInds);

    %remove nans
    nanInds = find(any(isnan(allFRs(:,matchedNeurInds))));
    strData(:,nanInds) = [];
    ctxData(:,nanInds) = [];
    emgData(:,nanInds) = [];

    %look at correlations
    strCorrs{iSess} = corr(strData',emgData');
    ctxCorrs{iSess} = corr(ctxData',emgData');

    %get controls and also lags
%     lags = -99:100;
%     minShiftAmount = 1000;
%     for iShift = 1:1200
% 
%         if iShift <=1000
%             shiftAmount = randperm(size(downsampEMG,2)-2*minShiftAmount,1)+minShiftAmount;
%         else
%             shiftAmount = lags(iShift-1000);
%         end
% 
%         shiftFRs = circshift(allFRs,shiftAmount,2);
%         nanInds = find(any(isnan(shiftFRs(:,matchedNeurInds))));
%         strDataShift = shiftFRs(1:length(striatumInds),matchedNeurInds);
%         ctxDataShift = shiftFRs(length(striatumInds)+1:end,matchedNeurInds);
%         emgDataShift = downsampEMG(:,usedEmgInds);
% 
%         strDataShift(:,nanInds) = [];
%         ctxDataShift(:,nanInds) = [];
%         emgDataShift(:,nanInds) = [];
% 
%         if iShift <=1000
%             strCorrsShift{iSess}(:,:,iShift) = corr(strDataShift',emgDataShift');
%             ctxCorrsShift{iSess}(:,:,iShift) = corr(ctxDataShift',emgDataShift');
%         else
%             strCorrsLag{iSess}(:,:,iShift-1000) = corr(strDataShift',emgDataShift');
%             ctxCorrsLag{iSess}(:,:,iShift-1000) = corr(ctxDataShift',emgDataShift');
%         end
% 
%     end
% 
%     strCutoff = mean(abs(strCorrsShift{iSess}),3)+3*std(abs(strCorrsShift{iSess}),[],3);
%     ctxCutoff = mean(abs(ctxCorrsShift{iSess}),3)+3*std(abs(ctxCorrsShift{iSess}),[],3);
% 
%     strModNeurs{iSess} = any(abs(strCorrs{iSess}) > strCutoff,2);
%     ctxModNeurs{iSess} = any(abs(ctxCorrs{iSess}) > ctxCutoff,2);
% 
%     strCorrZscore{iSess} = (abs(strCorrs{iSess})-mean(abs(strCorrsShift{iSess}),3)) ./ std(abs(strCorrsShift{iSess}),[],3);
%     ctxCorrZscore{iSess} = (abs(ctxCorrs{iSess})-mean(abs(ctxCorrsShift{iSess}),3)) ./ std(abs(ctxCorrsShift{iSess}),[],3);
% 
%     strCorrPVals{iSess} = (sum(repmat(abs(strCorrs{iSess}),1,1,1000) <= abs(strCorrsShift{iSess}),3)+1)/1001;
%     ctxCorrPVals{iSess} = (sum(repmat(abs(ctxCorrs{iSess}),1,1,1000) <= abs(ctxCorrsShift{iSess}),3)+1)/1001;


    % go through each region, get the time inds for the region
    % do it for a bunch of time shifts as well
    for iShift = 1:nBehvShifts+1

        %first loop don't do any shifting, otherwise, shift by at least 10s
        if iShift > 1
            minShiftAmount = 10000;
            shiftAmount = randperm(length(regionAssignmentsFiltered)-2*minShiftAmount,1)+minShiftAmount;
            regionLabelsToUse = circshift(regionAssignmentsFiltered,shiftAmount);
        else
            regionLabelsToUse = regionAssignmentsFiltered;
        end

        behvAlignPerms = allBehvAlignPerms(iSess,:);
        regionWatershedLabels = regionWatershedLabels(behvAlignPerms);

        for iRegion = 1:length(regionWatershedLabels)

            regionTimeInds = find(regionLabelsToUse == regionWatershedLabels(iRegion));
            regionTimeInds(regionTimeInds > floor(frameEMGSamples{1}{end}(end)/20)) = [];

            % first get indices for the EMG timeseries
            if use10msBins
                origEmgRegionInds = unique(round(origDownsampEMGInd(regionTimeInds)/10));
                origEmgRegionInds(origEmgRegionInds==0) = [];
            else
                origEmgRegionInds = origDownsampEMGInd(regionTimeInds);
            end
            regionEMGsInds{iRegion} = origEmgRegionInds;

            % next get it for the neural data
            if use10msBins
                currentDir = pwd;
                cd(fullfile(allDirs{iSess},'ProcessedData'))
                regionNeurInds{iRegion} = round(NeurEMGSync(origEmgRegionInds*200,...
                    frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
                cd(currentDir)

                outOfBoundInds = find(regionNeurInds{iRegion} > floor(frameNeuropixelSamples{1}{end}(end)/300));
                nanInds = find(isnan(regionNeurInds{iRegion}));
                zeroInds = find(regionNeurInds{iRegion}==0);

                %remove bad time points (nan inds,out of bound inds, ect)
                regionNeurInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
                regionEMGsInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
            else
                regionNeurInds{iRegion} = round(NeurEMGSync(origDownsampEMGInd(regionTimeInds)*20,...
                    frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
                outOfBoundInds = regionNeurInds{iRegion} > floor(frameNeuropixelSamples{1}{end}(end)/30);
                nanInds = find(isnan(regionNeurInds{iRegion}));
                zeroInds = find(regionNeurInds{iRegion}==0);
                regionNeurInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
                regionEMGsInds{iRegion}([outOfBoundInds nanInds zeroInds]) = [];
            end

            regionEMGs{iRegion} = downsampEMG(:,regionEMGsInds{iRegion});
            regionFRs{iRegion} = allFRs(:,regionNeurInds{iRegion});

            % don't use any nans in the data
            nanEMGInds = find(any(isnan(regionEMGs{iRegion})));
            nanFRInds = find(any(isnan(regionFRs{iRegion})));
            regionEMGs{iRegion}(:,unique([nanEMGInds nanFRInds])) = [];
            regionFRs{iRegion}(:,unique([nanEMGInds nanFRInds])) = [];
            regionNeurInds{iRegion}(unique([nanEMGInds nanFRInds])) = [];
            regionEMGsInds{iRegion}(unique([nanEMGInds nanFRInds])) = [];

            %save the data for later making example plots
            if iShift == 1
                savedRegionEMGs{iSess,iRegion} = regionEMGs{iRegion};
                savedRegionFRs{iSess,iRegion} = regionFRs{iRegion};
            end

        end

        % to keep number of time points consistent across regions,
        % downsample to the minimum number of points across regions
        regionMinPoints = min(cellfun(@(x) size(x,2),regionFRs));

        % calculate correlations
        for iRegion = 1:length(regionWatershedLabels)
            usedTimePoints = randperm(size(regionFRs{iRegion},2),regionMinPoints);
            corrs = corr(regionFRs{iRegion}(length(striatumInds)+1:end,usedTimePoints)',regionEMGs{iRegion}(:,usedTimePoints)');
            if iShift == 1
                regionCorrs{iSess,iRegion} =  corrs;
            else
                regionCorrsShift{iSess,iRegion,iShift-1} = corrs;
            end
        end

        % now go through and summerize cross region stability using
        % different metrics
        for iRegion1 = 1:length(regionWatershedLabels)
            for iRegion2 = 1:length(regionWatershedLabels)

                if iRegion1 == iRegion2
                    if iShift == 1
                        meanCorrDiff(iRegion1,iRegion2,iSess) = 0;
                        spearmanCorr(iRegion1,iRegion2,iSess) = 1;
                        r2Fit(iRegion1,iRegion2,iSess) = 1;
                    else
                        meanCorrDiffShift(iRegion1,iRegion2,iSess,iShift-1) = 0;
                        spearmanCorrShift(iRegion1,iRegion2,iSess,iShift-1) = 1;
                        r2FitShift(iRegion1,iRegion2,iSess,iShift-1) = 1;
                    end
                    continue
                end

                if iShift == 1
                    corrs1 = regionCorrs{iSess,iRegion1}(:,1:4);
                    corrs2 = regionCorrs{iSess,iRegion2}(:,1:4);
                else
                    corrs1 = regionCorrsShift{iSess,iRegion1,iShift-1}(:,1:4);
                    corrs2 = regionCorrsShift{iSess,iRegion2,iShift-1}(:,1:4);
                end

                %remove nans
                nanNeurs = unique([find(any(isnan(corrs1),2)); find(any(isnan(corrs2),2))]);
                corrs1(nanNeurs,:) = [];
                corrs2(nanNeurs,:) = [];

                %get average change in correlation
                thisCorrDiffs = mean(abs(corrs1(:)-corrs2(:)));

                %get spearman correlation across all neuron-muscle combos
                thisSpearCorrs = corr(corrs1(:),corrs2(:),'type','spearman');

                %get linear regression R2
                linReg = fitlm(corrs1(:),corrs2(:));
                thisR2 = linReg.Rsquared.Ordinary;

                if iShift == 1
                    meanCorrDiff(iRegion1,iRegion2,iSess) = thisCorrDiffs;
                    spearmanCorr(iRegion1,iRegion2,iSess) = thisSpearCorrs;
                    r2Fit(iRegion1,iRegion2,iSess) = thisR2;
                else
                    meanCorrDiffShift(iRegion1,iRegion2,iSess,iShift-1) = thisCorrDiffs;
                    spearmanCorrShift(iRegion1,iRegion2,iSess,iShift-1) = thisSpearCorrs;
                    r2FitShift(iRegion1,iRegion2,iSess,iShift-1) = thisR2;
                end

            end
        end

        %calculate hierarchy
        if iShift == 1
            corrDiffSpread(iSess) = std(squareform(meanCorrDiff(:,:,iSess)));
            spearmanCorrSpread(iSess) = std(squareform(1-spearmanCorr(:,:,iSess)));
            r2FitSpread(iSess) = std(squareform(1-r2Fit(:,:,iSess)));
        else
            corrDiffSpreadShift(iSess,iShift-1) = std(squareform(squeeze(meanCorrDiffShift(:,:,iSess,iShift-1))));
            spearmanCorrSpreadShift(iSess,iShift-1) = std(squareform(1-squeeze(spearmanCorrShift(:,:,iSess,iShift-1))));
            r2FitSpreadShift(iSess,iShift-1) = std(squareform(1-squeeze(r2FitShift(:,:,iSess,iShift-1))));
        end


    end

end

exampleAnimal = 3;
exampleRegion = 1;
exampleCloseRegion = 2;
exampleFarRegion = 7;
exampleNeuron = 148;
exampleMusc = 3;
exampleTime = 14000:17000;
exampleCloseTime = 6000:9000;
exampleFarTime = 27000:30000;

% make plots
plotColors = lines(7);

% first do example time series of neurons and muscles across behaviors
figure
tiledlayout(2,3)
nexttile
plot((0:3000)/100,savedRegionFRs{exampleAnimal,exampleRegion}(exampleNeuron,exampleTime))
nexttile
plot((0:3000)/100,savedRegionFRs{exampleAnimal,exampleCloseRegion}(exampleNeuron,exampleCloseTime))
nexttile
plot((0:3000)/100,savedRegionFRs{exampleAnimal,exampleFarRegion}(exampleNeuron,exampleFarTime))
nexttile
plot((0:3000)/100,savedRegionEMGs{exampleAnimal,exampleRegion}(exampleMusc,exampleTime))
nexttile
plot((0:3000)/100,savedRegionEMGs{exampleAnimal,exampleCloseRegion}(exampleMusc,exampleCloseTime))
nexttile
plot((0:3000)/100,savedRegionEMGs{exampleAnimal,exampleFarRegion}(exampleMusc,exampleFarTime))


% example scatter plot of all neur-musc corrs
figure('Color','w')
tiledlayout(2,1)
nexttile
exampleCorrs1 = regionCorrs{exampleAnimal,exampleRegion}(:,1:4);
exampleCorrs2 = regionCorrs{exampleAnimal,exampleCloseRegion}(:,1:4);
exampleCorrs3 = regionCorrs{exampleAnimal,exampleFarRegion}(:,1:4);
plot(exampleCorrs1(:),exampleCorrs2(:),'.')
xlims = get(gca,'XLim');
line(xlims,xlims,'color','k','linestyle','--')
box off
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'LineWidth',0.5)
set(gca,'FontSize',7)
set(gca,'TickDir','out')
xlabel('Climb up corr')
ylabel('Climb down corr')

nexttile
plot(exampleCorrs1(:),exampleCorrs3(:),'.')
xlims = get(gca,'XLim');
line(xlims,xlims,'color','k','linestyle','--')
box off
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'LineWidth',0.5)
set(gca,'FontSize',7)
set(gca,'TickDir','out')
xlabel('Climb up corr')
ylabel('Eating corr')


% next do example distance matrix
plotMat = meanCorrDiff(:,:,iSess);
figure
tiledlayout(2,1)
nexttile
imagesc(plotMat)
set(gca,'XTick',1:length(behvRegionLabels))
set(gca,'YTick',1:length(behvRegionLabels))
set(gca,'XTickLabelRotation',30)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gcf,'Color','w')
set(gca,'LineWidth',0.5)
set(gca,'FontSize',7)
set(gca,'TickLength',[0 0])
cH = colorbar;
caxis([0 0.12])
cH.Label.String = 'Mean corr difference';

% plot dendrogram
nexttile
[~, dendPerm] = customDendrogram(linkage(squareform(plotMat)),[],gca,{'Linewidth',2,'Color','k'});
% add an example shift control dendrogram
shiftExampleInd = find(corrDiffSpreadShift(iSess,1:99)==median(corrDiffSpreadShift(iSess,1:99)));
plotMatrix = meanCorrDiffShift(:,:,iSess,shiftExampleInd);
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(3,:) 0.7]})

set(gca,'TickDir','out')
set(gca,'LineWidth',0.5)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',7)
ylabel('Mean corr difference')
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')

% plot overall hierarchy spread across all animals
figure;
plotJitter = [-0.2 0 0.2];
barScatterPlot({corrDiffSpread},'none',[1; 1],repmat({plotJitter},2,1),[],true);
hold on

% also plot controls
for iSess = 1:length(allDirs)
    controlJitter = randn(1,size(corrDiffSpreadShift,2))/100;
    plotH = scatter(controlJitter+plotJitter(iSess)+1,corrDiffSpreadShift(iSess,:),10,plotColors(3,:),'filled');
    alpha(plotH,0.2)
end

ylabel('Hierarchy spread (s.d.)')


save('X:\David\AnalysesData\singleNeuronMuscCorrs','strCorrs','ctxCorrs','strCorrsShift','ctxCorrsShift','strCorrZscore','ctxCorrZscore','strCorrPVals','ctxCorrPVals','strModNeurs','ctxModNeurs','lags')
save('X:\David\AnalysesData\singleNeuronMuscCorrsRegions','regionCorrs','regionCorrsShift','savedRegionEMGs','savedRegionFRs','meanCorrDiff','meanCorrDiffShift','spearmanCorr','spearmanCorrShift','r2Fit','r2FitShift')

% 
