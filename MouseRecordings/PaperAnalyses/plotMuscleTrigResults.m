clear
load('X:\David\AnalysesData\MovementInitData.mat')

xTimes = -1*periTransTimes(1):periTransTimes(2);
plotColors = lines(3);
behvLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

useNonINs = false;

regionCrossingCtxControlNorm = cellfun(@(x) x - mean(x(:,1:50),2),regionCrossingCtxControlMean,'un',0);
regionCrossingStrControlNorm = cellfun(@(x) x - mean(x(:,1:50),2),regionCrossingStrControlMean,'un',0);

clear regionCrossingCtxControlMean regionCrossingStrControlMean

% first combine across all behavioral regions, average across neurons
avePopResponse = {};
for iAnimal = 1:size(regionCrossingCtx,1)

    load(fullfile(allDirs{iAnimal},'ProcessedData','neuronDataStruct.mat'))
    load(fullfile(allDirs{iAnimal},'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'striatumInds','cortexInds')
    
    spikeWidths = [neuronDataStruct.peakToValley];
    nonINCortexInds = find(spikeWidths(cortexInds) > 11.25);
    nonINStriatunmInds = find(spikeWidths(striatumInds) > 12.45);

    if useNonINs
        usedCortexInds = nonINCortexInds;
        usedStriatumInds = nonINStriatunmInds;
    else
        usedCortexInds = 1:length(cortexInds);
        usedStriatumInds = 1:length(striatumInds);
    end

    for iMusc = 1:size(regionCrossingCtx,3)

        % first concatenate across all regions for making example plots
        ctxAllRegions = cat(1,regionCrossingCtx{iAnimal,:,iMusc});
        ctxAllRegions = mean(ctxAllRegions(:,:,usedCortexInds),3);
        strAllRegions = cat(1,regionCrossingStr{iAnimal,:,iMusc});
        strAllRegions = mean(strAllRegions(:,:,nonINStriatunmInds),3);
        emgAllRegions = cat(1,regionCrossingEmg{iAnimal,:,iMusc});
        emgAllRegions = squeeze(emgAllRegions(:,:,iMusc))/50;

        % also concatenate across all regions and all muscles for making
        % summary plots
        ctxAllRegionsAllMusc = cat(1,regionCrossingCtx{iAnimal,:,:});
        ctxAllRegionsAllMusc = mean(ctxAllRegionsAllMusc(:,:,usedCortexInds),3);
        strAllRegionsAllMusc = cat(1,regionCrossingStr{iAnimal,:,:});
        strAllRegionsAllMusc = mean(strAllRegionsAllMusc(:,:,nonINStriatunmInds),3);

        % adjust to baseline and get average/sem for shaded plot
        ctxAllRegionsAve = mean(ctxAllRegions - mean(ctxAllRegions(:,1:50),2))*1000;
        strAllRegionsAve = mean(strAllRegions - mean(strAllRegions(:,1:50),2))*1000;
        emgAllRegionsAve = mean(emgAllRegions - mean(emgAllRegions(:,1:50),2));

        ctxAllRegionsSem = std((ctxAllRegions - mean(ctxAllRegions(:,1:50),2))*1000)/sqrt(size(ctxAllRegions,1));
        strAllRegionsSem = std((strAllRegions - mean(strAllRegions(:,1:50),2))*1000)/sqrt(size(strAllRegions,1));
        emgAllRegionsSem = std((emgAllRegions - mean(emgAllRegions(:,1:50),2)))/sqrt(size(emgAllRegions,1));

        ctxAllRegionsAllMuscAve = mean(ctxAllRegionsAllMusc - mean(ctxAllRegionsAllMusc(:,1:50),2))*1000;
        strAllRegionsAllMuscAve = mean(strAllRegionsAllMusc - mean(strAllRegionsAllMusc(:,1:50),2))*1000;

        ctxAllRegionsAllMuscSem = std((ctxAllRegionsAllMuscAve - mean(ctxAllRegionsAllMuscAve(:,1:50),2))*1000)/sqrt(size(ctxAllRegionsAllMusc,1));
        strAllRegionsAllMuscSem = std((strAllRegionsAllMuscAve - mean(strAllRegionsAllMuscAve(:,1:50),2))*1000)/sqrt(size(strAllRegionsAllMusc,1));

        % for the controls, get the 95% confidence interval across all control shifts
        ctxAllRegionsControl = [];
        strAllRegionsControl = [];

        for iShift = 1:size(regionCrossingCtxControlNorm,4)

            % for both individual muscles (for example plots) and all
            % muscles comibned (for summary plots)
            ctxAllRegionsControl(iShift,:) = squeeze(mean(cat(1,regionCrossingCtxControlNorm{iAnimal,:,iMusc,iShift})));
            strAllRegionsControl(iShift,:) = squeeze(mean(cat(1,regionCrossingStrControlNorm{iAnimal,:,iMusc,iShift})));

            ctxAllRegionsAllMuscControl(iShift,:) = squeeze(mean(cat(1,regionCrossingCtxControlNorm{iAnimal,:,:,iShift})));
            strAllRegionsAllMuscControl(iShift,:) = squeeze(mean(cat(1,regionCrossingStrControlNorm{iAnimal,:,:,iShift})));

        end
        ctxConfLower = prctile(ctxAllRegionsControl,2.5)*1000;
        ctxConfUpper = prctile(ctxAllRegionsControl,97.5)*1000;
        strConfLower = prctile(strAllRegionsControl,2.5)*1000;
        strConfUpper = prctile(strAllRegionsControl,97.5)*1000;

        ctxAllMuscConfLower = prctile(ctxAllRegionsAllMuscControl,2.5)*1000;
        ctxAllMuscConfUpper = prctile(ctxAllRegionsAllMuscControl,97.5)*1000;
        strAllMuscConfLower = prctile(strAllRegionsAllMuscControl,2.5)*1000;
        strAllMuscConfUpper = prctile(strAllRegionsAllMuscControl,97.5)*1000;

        % get average activity in the post period for making summary plots
        % first for all regions
        avePopResponse{1,1}(iAnimal) = mean(mean(ctxAllRegionsAllMusc(:,periTransTimes(1)+1:end) - mean(ctxAllRegionsAllMusc(:,1:50),2)))*1000;
        avePopResponse{1,2}(iAnimal) = mean(mean(strAllRegionsAllMusc(:,periTransTimes(1)+1:end) - mean(strAllRegionsAllMusc(:,1:50),2)))*1000;

        % then get for each region individually
        for iRegion = 1:size(regionCrossingCtx,2)

            ctxAllMusc = cat(1,regionCrossingCtx{iAnimal,iRegion,:});
            strAllMusc = cat(1,regionCrossingStr{iAnimal,iRegion,:});

            avePopResponse{iRegion+1,1}(iAnimal) = mean(mean(mean(ctxAllMusc(:,periTransTimes(1)+1:end,:),3) - ...
                mean(mean(ctxAllMusc(:,1:50,:),3),2)))*1000;

            avePopResponse{iRegion+1,2}(iAnimal) = mean(mean(mean(strAllMusc(:,periTransTimes(1)+1:end,:),3) - ...
                mean(mean(strAllMusc(:,1:50,:),3),2)))*1000;

        end


        %next, get neural correlations across behaviors
        for iRegion1 = 1:size(regionCrossingCtx,2)
            for iRegion2 = 1:size(regionCrossingCtx,2)

                for iNeuron = 1:size(regionCrossingCtx{iAnimal,iRegion1,iMusc},3)
                    if iRegion1 == iRegion2
                        ctxCorrs{iRegion1,iRegion2,iAnimal,iMusc}(iNeuron) = 1;
                        strCorrs{iRegion1,iRegion2,iAnimal,iMusc}(iNeuron) = 1;
                    else
                        ctxCorrs{iRegion1,iRegion2,iAnimal,iMusc}(iNeuron) = corr(mean(regionCrossingCtx{iAnimal,iRegion1,iMusc}(:,:,iNeuron),1)',mean(regionCrossingCtx{iAnimal,iRegion2,iMusc}(:,:,iNeuron),1)');
                        strCorrs{iRegion1,iRegion2,iAnimal,iMusc}(iNeuron) = corr(mean(regionCrossingStr{iAnimal,iRegion1,iMusc}(:,:,iNeuron),1)',mean(regionCrossingStr{iAnimal,iRegion2,iMusc}(:,:,iNeuron),1)');
                    end
                end

                for iMusc2 = 1:size(regionCrossingEmg{iAnimal,iRegion1,iMusc},3)
                    if iRegion1 == iRegion2
                        emgCorrs{iRegion1,iRegion2,iAnimal}(iMusc2) = 1;
                    else
                        emgCorrs{iRegion1,iRegion2,iAnimal,iMusc}(iMusc2) = corr(mean(regionCrossingEmg{iAnimal,iRegion1,iMusc}(:,:,iMusc2),1)',mean(regionCrossingEmg{iAnimal,iRegion2,iMusc}(:,:,iMusc2),1)');
                    end
                end

            end

        end


        % now make individual animal plots
        figure;
        tiledlayout(1,2,"TileSpacing","tight","Padding","tight")
        nexttile
        shadedErrorBar(xTimes,emgAllRegionsAve,emgAllRegionsSem,'lineProps',{'color',[0.2 0.2 0.2],'linewidth',1.5})
        hold on
        shadedErrorBar(xTimes,ctxAllRegionsAve,ctxAllRegionsSem,'lineProps',{'color',plotColors(1,:),'linewidth',1.5})
        plot(xTimes,ctxConfLower,'color',plotColors(1,:),'LineStyle','--')
        plot(xTimes,ctxConfUpper,'color',plotColors(1,:),'LineStyle','--')
        line([0,0],get(gca,'YLim'),'linestyle','--','color',[0.2 0.2 0.2])

        ylabel('Population Rate (spks/s)')
        xlabel('Time (ms)')
        set(gca,'LineWidth',1.5)
        set(gca,'FontSize',14)
        set(gca,'TickDir','out')
        set(gca,'XColor',[0 0 0])
        set(gca,'YColor',plotColors(1,:))
        box off

        nexttile
        shadedErrorBar(xTimes,emgAllRegionsAve,emgAllRegionsSem,'lineProps',{'color',[0.2 0.2 0.2],'linewidth',1.5})
        hold on
        shadedErrorBar(xTimes,strAllRegionsAve,strAllRegionsSem,'lineProps',{'color',plotColors(2,:),'linewidth',1.5,'linestyle','-'})
        plot(xTimes,strConfLower,'color',plotColors(2,:),'LineStyle','--')
        plot(xTimes,strConfUpper,'color',plotColors(2,:),'LineStyle','--')
        line([0,0],get(gca,'YLim'),'linestyle','--','color',[0.2 0.2 0.2])

        ylabel('PL EMG (a.u.)')
        xlabel('Time (ms)')
        set(gca,'LineWidth',1.5)
        set(gca,'FontSize',14)
        set(gca,'TickDir','out')
        set(gca,'YAxisLocation','right')
        set(gca,'XColor',[0 0 0])
        set(gca,'YColor',[0 0 0])
        box off

        set(gcf,'color','w')

        title(['Animal ' num2str(iAnimal) ', Muscle ' num2str(iMusc)])

    end

end

% make plots for average activity comparing ctx and str
figure;
[barH] = barScatterPlot(avePopResponse,'none',ones(size(avePopResponse,1),2),[],[1 2]);
ylabel('Mean Population Rate (spks/s)')
set(gca,'XTickLabel',['All Behaviors' behvLabels])
legendH = legend(barH, 'Cortex','Striatum','box','off');
for iLabel = 1:length(legendH.String)
    legendH.String{iLabel} = ['\color[rgb]{' num2str(barH(iLabel).FaceColor) '} ' legendH.String{iLabel}];
end


% next, make example plot of EMG for different behaviors
exampleAnimal = 1;
figure
tiledlayout(1,size(regionCrossingEmg,2),"TileSpacing","none","Padding","tight")
for iRegion = 1:size(regionCrossingEmg,2)
    nexttile
    regionMeanEMG = mean(regionCrossingEmg{exampleAnimal,iRegion}(:,:,4));
    regionSemEMG = std(regionCrossingEmg{exampleAnimal,iRegion}(:,:,1))/sqrt(size(regionCrossingEmg{exampleAnimal,iRegion},1));
    shadedErrorBar(xTimes,regionMeanEMG,regionSemEMG,'lineProps',{'color',[0.2 0.2 0.2],'linewidth',2})
    ylim([0 155])

    title(behvLabels{iRegion})
    set(gca,'FontSize',14)
    if iRegion==1
        line([-99 1],[10 10],'linewidth',2,'color','k')
        ylabel('PL EMG (a.u.)')
        set(gca,'LineWidth',1.5)
        set(gca,'TickDir','out')
    else
        axis off
    end
end

set(gcf,'color','w')

%get raster plot and sort it by peak time of max activity
unNormNeurAve{iAnimal,iRegion} = squeeze(mean(regionCrossingCtx{iRegion},1))';
normNeurAve{iAnimal,iRegion} = unNormNeurAve{iAnimal,iRegion}(rastPerm,:)./allRegionMaxFRs(rastPerm);
unNormEmgAve{iAnimal,iRegion} = squeeze(mean(regionCrossingEmg{iRegion},1))';

%get correlations with muscle
for iNeuron = 1:size(regionCrossingCtx{iRegion},3)
    neurActivity = squeeze(regionCrossingCtx{iRegion}(:,:,iNeuron));
    neurEmgCorrs(iNeuron,iRegion) = corr(mean(neurActivity(:,:))',mean(regionCrossingEmg{iRegion}(:,:,threshChan))');
end

figure(popFigH)
nexttile
shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),regionAveEmg,regionStdEmg/sqrt(nSamplesEmg))
hold on;plotH(1) = shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),regionAveCtx,regionStdCtx/sqrt(nSamplesCtx),'lineProps',{'Color','r'});
hold on;plotH(2) = shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),regionAveStr,regionStdStr/sqrt(nSamplesStr),'lineProps',{'Color','b'});

line([0, 0], get(gca,'ylim'),'linestyle','--','color','r')
xlabel('Time (ms)')
title([behvRegionLabels{iRegion} ', ' num2str(size(regionCrossingStr{iRegion},1)) ' Trials']);

figure(rasterFigH)
nexttile
imagesc(normNeurAve{iAnimal,iRegion});
set(gca,'XTick',50:100:450)
set(gca,'XTickLabel',-100:100:300)
line([periTransTimes(1) periTransTimes(1)],get(gca,'YLim'),'color','r','linewidth',2)
title([behvRegionLabels{iRegion} ', ' num2str(size(regionCrossingStr{iRegion},1)) ' Trials']);

xlabel('Time (ms)')
ylabel('Cortical Neuron')



aveStr = mean(mean(periCrossingFR(:,:,1:length(striatumInds)),3)*1000,1) - mean(mean(mean(periCrossingFR(:,1:periTransTimes(1)-50,1:length(striatumInds)))))*1000;
aveCtx = mean(mean(periCrossingFR(:,:,length(striatumInds)+1:end),1)*1000,3) - mean(mean(mean(periCrossingFR(:,1:periTransTimes(1)-50,length(striatumInds)+1:end))))*1000;
aveEMG = mean(squeeze(periCrossingEMG(:,:,threshChan)))/100 - mean(mean(mean(periCrossingEMG(:,1:periTransTimes(1)-50,:))))/100;

% aveStrControl = mean(mean(periControlFR(:,:,1:length(striatumInds)),1)*1000,3) - mean(mean(mean(periControlFR(:,1:(sum(periTransTimes)+1),1:length(striatumInds)))))*1000;
% aveCtxControl = mean(mean(periControlFR(:,:,length(striatumInds)+1:end),1)*1000,3) - mean(mean(mean(periControlFR(:,1:(sum(periTransTimes)+1),length(striatumInds)+1:end))))*1000;
% aveEMGControl = mean(squeeze(periControlEMG(:,:,threshChan)))/100 - mean(mean(mean(periControlEMG(:,1:(sum(periTransTimes)+1),:))))/100;

% stdStr = std(reshape(permute(periCrossingFR(:,:,1:length(striatumInds))*1000,[2 1 3]),size(periCrossingFR,2),[]),[],2);
% stdCtx= std(reshape(permute(periCrossingFR(:,:,length(striatumInds+1))*1000,[2 1 3]),size(periCrossingFR,2),[]),[],2);
% stdEMG = std(squeeze(periCrossingEMG(:,:,threshChan))/100);

stdStr = std(mean(periCrossingFR(:,:,1:length(striatumInds))*1000,3));
stdCtx = std(mean(periCrossingFR(:,:,length(striatumInds)+1:end)*1000,3));
stdEMG = std(squeeze(periCrossingEMG(:,:,threshChan))/100);



% 
