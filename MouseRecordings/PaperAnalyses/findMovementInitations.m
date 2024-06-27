clear

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

allAnimalLabels = {'D020','D024','D026'};
behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};

% initation detection parameters
fs = 1000;
threshChan = 4;
thresh = [0.7 0.7 0.7 0.7;...
          0.5 0.5 0.5 0.7;...
          0.7 0.7 0.7 0.7];

baselineChanThresh = [0.5 0.5 0.5 0.5;...
                      0.5 0.5 0.5 0.5;...
                      0.5 0.5 0.5 0.5];
neurBinSize = 1; %in ms

preThreshDuration = 150; %in ms
preThreshbuffer = 100;

periTransTimes = [100 200]; %in ms

nControlResamples = 100; %number of times to get the controls

for iAnimal = 1:length(allDirs)

    baseDir = allDirs{iAnimal};
    load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'regionAssignmentsFiltered','reduction','origDownsampEMGInd','regionWatershedLabels')

    load(fullfile(baseDir,'ProcessedData','neuronDataStruct.mat'))
    load(fullfile(baseDir,'Neuropixels','artifactTimestamps.mat'))
    load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

    animalLabel = allAnimalLabels{iAnimal};
    behvAlignPerm = allBehvAlignPerms(iAnimal,:);
    regionAssignmentLabels = unique(regionAssignmentsFiltered);

    % load EMG data
    load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))

    % normalize by standard deviation
    downsampEMGNorm = downsampEMG./std(downsampEMG,[],2);

    % load in neural data
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'),'allFRs','cortexInds','striatumInds')

    for iThreshChan = 1:4

        % get all threshold crossings
        threshCrossings = find(downsampEMGNorm(iThreshChan,2:end) > thresh(iAnimal,iThreshChan) & downsampEMGNorm(iThreshChan,1:end-1) <= thresh(iAnimal,iThreshChan));

        % prune detected threshold crossings
        goodCrossings = zeros(1,length(threshCrossings));
        for iCross = 1:length(threshCrossings)

            % make sure it's not too close to the beginning or end of the session
            if (threshCrossings(iCross)-preThreshDuration < 0) || (threshCrossings(iCross)-1 > size(downsampEMG,2))
                continue
            end

            % get the baseline periods
            fullBaselineEMG = downsampEMGNorm(:,threshCrossings(iCross)-preThreshDuration : threshCrossings(iCross)-1);
            baselineEMG = downsampEMGNorm(:,threshCrossings(iCross)-preThreshDuration : threshCrossings(iCross)-preThreshbuffer);

            %There should be no threshhold crossings in the baseline period
            if any(threshCrossings > threshCrossings(iCross)-preThreshDuration & ...
                    threshCrossings < threshCrossings(iCross))
                continue
            end

            %Channel that we threshold to should have no activity in the baseline
            %period
            if any((baselineEMG(iThreshChan,:)) > baselineChanThresh(iAnimal,iThreshChan))
                continue
            end

            %if there are multiple zeros, then that means it was artifact, don't use
            if sum(baselineEMG(iThreshChan,:)==0) > 20
                continue
            end

            %baseline should have no flucations past a certain amount
            %     if any(max(baselineEMG,[],2)-min(baselineEMG,[],2) > baselineFlucThresh)
            %         continue
            %     end

            goodCrossings(iCross) = 1;

        end

        threshCrossings = threshCrossings(find(goodCrossings));

        % get the thresh crossing times in neural and umap indices
        [threshCrossings, reducInds, crossingNeurInds] = convertCrossInds(threshCrossings,origDownsampEMGInd,baseDir);


        % remove any emg crossings that is outside of neural data range
        outRangeCrossings = find(crossingNeurInds>size(allFRs,2)-periTransTimes(2));
        threshCrossings(outRangeCrossings) = [];
        crossingNeurInds(outRangeCrossings) = [];
        reducInds(outRangeCrossings) = [];

        % now save crossing times for this animal
        threshCrossOrigEMG{iAnimal,iThreshChan} = threshCrossings;
        threshCrossReduc{iAnimal,iThreshChan} = reducInds;
        threshCrossNeur{iAnimal,iThreshChan} = crossingNeurInds;

        % % also do shifts
        % for iShift = 1:100
        %
        %     shiftAmount(iShift) = randi(size(downsampEMG,2)-10000*2)+10000;
        %
        %     threshCrossingsShift = threshCrossings + shiftAmount(iShift);
        %     threshCrossingsShift(threshCrossingsShift > size(downsampEMG,2) - periTransTimes(2) - 1) = ...
        %         threshCrossingsShift(threshCrossingsShift > size(downsampEMG,2) - periTransTimes(2) - 1) - (size(downsampEMG,2) - periTransTimes(2) - 1);
        %
        %     threshCrossingsShift(threshCrossingsShift < periTransTimes(1)+1) = [];
        %     % find the corresponding neural index for each threshold crossing
        %     currentDir = pwd;
        %     cd(fullfile(baseDir,'ProcessedData'))
        %     crossingNeurIndsShift = round(NeurEMGSync(threshCrossingsShift*20,...
        %         frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
        %     cd(currentDir)
        %
        %     % remove any emg crossings that is outside of neural data range
        %     outRangeCrossings = find(crossingNeurIndsShift>size(allFRs,2)-periTransTimes(2));
        %     threshCrossingsShift(outRangeCrossings) = [];
        %     crossingNeurIndsShift(outRangeCrossings) = [];
        %
        %     periCrossingFRShift = zeros(length(threshCrossingsShift),sum(periTransTimes)+1,size(allFRs,1));
        %     periCrossingEMGShift = zeros(length(threshCrossingsShift),sum(periTransTimes)+1,size(downsampEMG,1));
        %
        %     for iCross = 1:length(threshCrossingsShift)
        %
        %         %get neural activity
        %         for iNeuron = 1:size(allFRs,1)
        %             periCrossingFRShift(iCross,:,iNeuron) = allFRs(iNeuron,...
        %                 crossingNeurIndsShift(iCross)-periTransTimes(1):crossingNeurIndsShift(iCross)+periTransTimes(2));
        %         end
        %
        %         %get EMG activity
        %         for iMuscle = 1:size(downsampEMG,1)
        %             periCrossingEMGShift(iCross,:,iMuscle) = downsampEMG(iMuscle,...
        %                 threshCrossingsShift(iCross)-periTransTimes(1):threshCrossingsShift(iCross)+periTransTimes(2));
        %         end
        %
        %     end
        %
        %     nanNeurCrossingsShift = find(any(isnan(periCrossingFRShift(:,:))'));
        %     threshCrossingsShift(nanNeurCrossingsShift) = [];
        %     crossingNeurIndsShift(nanNeurCrossingsShift) = [];
        %     periCrossingFRShift(nanNeurCrossingsShift,:,:) = [];
        %     periCrossingEMGShift(nanNeurCrossingsShift,:,:) = [];
        %
        %     aveStrShift(iShift,:) = mean(mean(periCrossingFRShift(:,:,1:length(striatumInds)),3)*1000,1) - mean(mean(mean(periCrossingFRShift(:,1:periTransTimes(1)-50,1:length(striatumInds)))))*1000;
        %     aveCtxShift(iShift,:) = mean(mean(periCrossingFRShift(:,:,length(striatumInds)+1:end),1)*1000,3) - mean(mean(mean(periCrossingFRShift(:,1:periTransTimes(1)-50,length(striatumInds)+1:end))))*1000;
        %     aveEMGShift(iShift,:) = mean(squeeze(periCrossingEMGShift(:,:,threshChan)))/100 - mean(mean(mean(periCrossingEMGShift(:,1:periTransTimes(1)-50,:))))/100;
        %
        % end

        % now get the actual data around the crossings
        [periCrossingFR, periCrossingEMG, crossRegion, nanNeurCrossings] = getDataFromCrossings(threshCrossings,crossingNeurInds,reducInds,...
            periTransTimes,allFRs,downsampEMG,regionAssignmentsFiltered,regionWatershedLabels);

        reducInds(nanNeurCrossings) = [];
        threshCrossings(nanNeurCrossings) = [];
        crossingNeurInds(nanNeurCrossings) = [];
        crossRegion(nanNeurCrossings) = [];
        periCrossingFR(nanNeurCrossings,:,:) = [];
        periCrossingEMG(nanNeurCrossings,:,:) = [];

        % now divide by behavioral region
        for iRegion = 1:length(behvAlignPerm)

            regionInd = behvAlignPerm(iRegion);

            % pull out the crossings in this region
            regionCrossInds = find(crossRegion==regionInd);

            regionCrossingStr{iAnimal,iRegion,iThreshChan}(:,:,:) = periCrossingFR(regionCrossInds,:,1:length(striatumInds));
            regionCrossingCtx{iAnimal,iRegion,iThreshChan}(:,:,:) = periCrossingFR(regionCrossInds,:,length(striatumInds)+1:end);
            regionCrossingEmg{iAnimal,iRegion,iThreshChan}(:,:,:) = periCrossingEMG(crossRegion==regionInd,:,:);

            % next do shift controls
            regionUMAPInds = find(regionAssignmentsFiltered == regionAssignmentLabels(regionInd));

            % for the controls only get segments that are all within the same
            % behavioral region
            regionSegsStarts = [regionUMAPInds(1) regionUMAPInds(find(diff(regionUMAPInds)>1)+1)];
            regionSegsEnds = [regionUMAPInds(find(diff(regionUMAPInds)>1)) regionUMAPInds(end)];
            segLengths = regionSegsEnds - regionSegsStarts;

            goodSegsCell = {};
            goodSegInd = 1;
            for iSeg = 1:length(regionSegsStarts)

                if segLengths(iSeg) < (sum(periTransTimes)+1)*2
                    continue
                end

                % also leave some padding at the beginning and end of each
                % segment so the whole extracted time period is within the region
                goodSegsCell{goodSegInd} = regionSegsStarts(iSeg)+periTransTimes(1) : regionSegsEnds(iSeg)-periTransTimes(2);
                goodSegInd = goodSegInd+1;
            end
            regionGoodInds = cat(2,goodSegsCell{:});

            % randomly sample controls
            for iControl = 1:nControlResamples

                gotGoodControls = false;
                while ~gotGoodControls

                    controlCrossings = origDownsampEMGInd(regionGoodInds(randperm(length(regionGoodInds),length(regionCrossInds))));

                    % get the thresh crossing times in neural and umap indices
                    [controlCrossings, reducIndsControl, crossingNeurIndsControl] = convertCrossInds(controlCrossings,origDownsampEMGInd,baseDir);

                    % remove any emg crossings that is outside of neural data range
                    outRangeCrossings = find(crossingNeurIndsControl>size(allFRs,2)-periTransTimes(2));
                    controlCrossings(outRangeCrossings) = [];
                    crossingNeurIndsControl(outRangeCrossings) = [];
                    reducIndsControl(outRangeCrossings) = [];

                    % sometimes if there are very few crossings, the random
                    % samples don't get good indices and there are no
                    % control crossings which causes errors, so keep
                    % resampling until there is same number of control
                    % crossings, or at least 5 valid control crossings
                    if length(controlCrossings) == length(regionCrossInds) || length(controlCrossings) > 5
                        gotGoodControls = true;
                    end

                end

                [controlCrossingFR, controlCrossingEMG, crossRegionControl, nanNeurCrossings] = getDataFromCrossings(...
                    controlCrossings,crossingNeurIndsControl,reducIndsControl, periTransTimes,allFRs,downsampEMG,regionAssignmentsFiltered,regionWatershedLabels);

                reducIndsControl(nanNeurCrossings) = [];
                controlCrossings(nanNeurCrossings) = [];
                crossingNeurIndsControl(nanNeurCrossings) = [];
                crossRegionControl(nanNeurCrossings) = [];
                controlCrossingFR(nanNeurCrossings,:,:) = [];
                controlCrossingEMG(nanNeurCrossings,:,:) = [];

                regionCrossIndsControls{iAnimal,iRegion,iThreshChan,iControl} = controlCrossings;

                % due to space constraints, save averages across neurons
                regionCrossingStrControlMean{iAnimal,iRegion,iThreshChan,iControl} = mean(controlCrossingFR(:,:,1:length(striatumInds)),3);
                regionCrossingCtxControlMean{iAnimal,iRegion,iThreshChan,iControl} = mean(controlCrossingFR(:,:,length(striatumInds)+1:end),3);
                regionCrossingEmgControl{iAnimal,iRegion,iThreshChan,iControl} = controlCrossingEMG(:,:,1:4);

                % but save 3 control samples for testing
                if iControl <= 3
                    regionCrossingStrControl{iAnimal,iRegion,iThreshChan,iControl} = controlCrossingFR(:,:,1:length(striatumInds));
                    regionCrossingCtxControl{iAnimal,iRegion,iThreshChan,iControl} = controlCrossingFR(:,:,length(striatumInds)+1:end);
                end

            end

        end

        % stdStrControl = std(reshape(permute(periControlFR(:,:,1:length(striatumInds))*1000,[2 1 3]),size(periControlFR,2),[]),[],2);
        % stdCtxControl = std(reshape(permute(periControlFR(:,:,length(striatumInds+1))*1000,[2 1 3]),size(periControlFR,2),[]),[],2);
        % stdEMGControl = std(squeeze(periCrossingEMG(:,:,threshChan))/100);

        % nSamplesStr = size(periControlFR,1);
        % nSamplesCtx = size(periControlFR,1);
        % nSamplesEmg = size(periCrossingEMG,1);

        % for iCross = 1:length(crossingNeurInds)
        %     periCrossingEMG(:,iCross) = downsampEMG(threshChan,...
        %         threshCrossings(iCross)-periTransTimes(1):threshCrossings(iCross)+periTransTimes(2));
        % end
        % periCrossingEMG = periCrossingEMG/50000;

        % controlFR = permute(controlFR,[2 1 3]);
        % periCrossingFR = permute(periCrossingFR(:,:,:),[2 1 3]);
        %
        % figure
        % shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),aveEMG,stdEMG/sqrt(nSamplesEmg))
        % hold on;plotH(1) = shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),aveCtx,stdCtx/sqrt(nSamplesCtx),'lineProps',{'Color','r'});
        % hold on;plotH(2) = shadedErrorBar((-1*periTransTimes(1)):periTransTimes(2),aveStr,stdStr/sqrt(nSamplesStr),'lineProps',{'Color','b'});
        % line([0 0],get(gca,'YLim'),'color','r','linewidth',2,'linestyle','--')
        % legend('EMG','Cortex','Striatum','box','off')
        % xlim([-1*periTransTimes(1) periTransTimes(2)])
        % xlabel('Time (ms)')
        % ylabel('Firing Rate (spks/s)')
        % set(gcf,'color','w')
        % title(animalLabel)

    end

end

% Now divide by region
for iAnimal = 1:3

    crossRegion = allCrossRegion{iAnimal};
    behvAlignPerm = allBehvAlignPerms(iAnimal,:);
    periCrossingFR = allPeriCrossingFR{iAnimal};
    periCrossingEMG = allPeriCrossingEMG{iAnimal};
    striatumInds = allStrInds{iAnimal};
    
    popFigH = figure('Color','w');
    popTileH = tiledlayout(2,4,'TileSpacing','tight');
    popTileH.Title.String = animalLabel;

    rasterFigH = figure('Color','w');
    rasterTileH = tiledlayout(2,4,'TileSpacing','tight');
    rasterTileH.Title.String = animalLabel;

    regionCrossingCtx = {};
    maxFRs = [];
    meanFRs = [];
    for iRegion = 1:length(unique(crossRegion))
        regionInd = behvAlignPerm(iRegion);
        regionCrossingCtx{iRegion}(:,:,:) = periCrossingFR(crossRegion==regionInd,:,length(striatumInds)+1:end);
        unNormNeurAve{iAnimal,iRegion} = squeeze(mean(regionCrossingCtx{iRegion},1))';
        maxFRs(:,iRegion) = max(unNormNeurAve{iAnimal,iRegion},[],2);
        meanFRs(:,iRegion) = mean(unNormNeurAve{iAnimal,iRegion}(:,1:300)');
    end
    allRegionMaxFRs = max(maxFRs,[],2);
    allRegionMeanFRs = max(meanFRs,[],2);
    goodNeurons = allRegionMeanFRs*1000 > 0.2;
    allRegionMaxFRs = allRegionMaxFRs(goodNeurons);

    orderRegion = 1;
    unNormNeurAve{iAnimal,iRegion} = squeeze(mean(regionCrossingCtx{orderRegion},1))';
    unNormNeurAve{iAnimal,iRegion} = unNormNeurAve{iAnimal,iRegion}(goodNeurons,:);
    [~, maxFRInds] = max(unNormNeurAve{iAnimal,iRegion},[],2);
    rastPerm = sortrows([maxFRInds (1:length(maxFRInds))']);
    rastPerm = rastPerm(:,2);

    regionCrossingStr = {};
    regionCrossingCtx = {};
    regionCrossingEmg = {};
    neurEmgCorrs = [];
    for iRegion = 1:length(unique(crossRegion))

        regionInd = behvAlignPerm(iRegion);

        %average activity across the population
        regionCrossingStr{iRegion}(:,:,:) = periCrossingFR(crossRegion==regionInd,:,1:length(striatumInds));
        regionCrossingCtx{iRegion}(:,:,:) = periCrossingFR(crossRegion==regionInd,:,length(striatumInds)+1:end);
        regionCrossingCtx{iRegion} = regionCrossingCtx{iRegion}(:,:,goodNeurons);
        regionCrossingEmg{iRegion}(:,:,:) = periCrossingEMG(crossRegion==regionInd,:,1:4);

        regionAveStr = mean(mean(regionCrossingStr{iRegion},1)*1000,3) - mean(mean(mean(regionCrossingStr{iRegion}(:,1:periTransTimes(1),:))))*1000;
        regionAveCtx = mean(mean(regionCrossingCtx{iRegion},1)*1000,3) - mean(mean(mean(regionCrossingCtx{iRegion}(:,1:periTransTimes(1),:))))*1000;
        regionAveEmg = mean(mean(regionCrossingEmg{iRegion},1)/25,3) - mean(mean(mean(regionCrossingEmg{iRegion}(:,1:periTransTimes(1),:))))/25;

        regionStdStr = std(reshape(permute(regionCrossingStr{iRegion}*1000,[2 1 3]),size(regionCrossingStr{iRegion},2),[]),[],2);
        regionStdCtx = std(reshape(permute(regionCrossingCtx{iRegion}*1000,[2 1 3]),size(regionCrossingCtx{iRegion},2),[]),[],2);
        regionStdEmg = std(squeeze(regionCrossingEmg{iRegion}(:,:,threshChan))/25);

        nSamplesStr = size(regionCrossingStr{iRegion},1)*size(regionCrossingStr{iRegion},3);
        nSamplesCtx = size(regionCrossingCtx{iRegion},1)*size(regionCrossingCtx{iRegion},3);
        nSamplesEmg = size(regionCrossingEmg{iRegion},1);

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
    end


end

for iAnimal = 1:3
    periCrossingFR = allPeriCrossingFR{iAnimal};
    periCrossingEMG = allPeriCrossingEMG{iAnimal};
    periCrossingFR = allPeriCrossingFR{iAnimal};
    behvAlignPerm = allBehvAlignPerms(iAnimal,:);
    goodNeurons = allGoodNeurons{iAnimal};
    crossRegion = allCrossRegion{iAnimal};
    allRegionMaxFRs = allMaxFRs{iAnimal};
    striatumInds = allStrInds{iAnimal};
    for iRegion = 1:7
        regionInd = behvAlignPerm(iRegion);
        %average activity across the population
        regionCrossingCtx{iAnimal,iRegion}(:,:,:) = periCrossingFR(crossRegion==regionInd,:,length(striatumInds)+1:end);
        regionCrossingCtx{iAnimal,iRegion} = regionCrossingCtx{iAnimal,iRegion}(:,:,goodNeurons);
        regionCrossingEmg{iAnimal,iRegion} = periCrossingEMG(crossRegion==regionInd,:,:);
        %get raster plot and sort it by peak time of max activity
        unNormNeurAve{iAnimal,iRegion} = squeeze(mean(regionCrossingCtx{iAnimal,iRegion},1))';
        unNormEmgAve{iAnimal,iRegion} = squeeze(mean(regionCrossingEmg{iAnimal,iRegion},1))';
        normNeurAve{iAnimal,iRegion} = unNormNeurAve{iAnimal,iRegion}./allRegionMaxFRs;
    end
end

% do pairwise comparison of correlation stability across regions
% don't use any neurons which have nans
goodCorrs = neurEmgCorrs(~any(isnan(neurEmgCorrs),2),:);
for iAnimal = 1:3

    for iRegion1 = 1:length(unique(crossRegion))
        for iRegion2 = 1:length(unique(crossRegion))
%             linReg = fitlm(goodCorrs(:,iRegion1),goodCorrs(:,iRegion2));
%             corrFit(iRegion1,iRegion2) = linReg.Rsquared.Ordinary;

            for iNeuron = 1:size(unNormNeurAve{iAnimal,1},1)
                neurCorrs{iRegion1,iRegion2,iAnimal}(iNeuron) = corr(unNormNeurAve{iAnimal,iRegion1}(iNeuron,:)',unNormNeurAve{iAnimal,iRegion2}(iNeuron,:)');
            end
            for iMusc = 1:size(unNormEmgAve{iAnimal,1},1)
                if iRegion1 == iRegion2
                    emgCorrs{iRegion1,iRegion2,iAnimal}(iMusc) = 1;
                else
                    emgCorrs{iRegion1,iRegion2,iAnimal}(iMusc) = corr(unNormEmgAve{iAnimal,iRegion1}(iMusc,:)',unNormEmgAve{iAnimal,iRegion2}(iMusc,:)');
                end
            end
        end

    end

end

neurCorrsMean = cellfun(@(x) nanmean(x), neurCorrs);
emgCorrsSingle = cellfun(@(x) x(threshChan), emgCorrs);

% plot hiearchy
figure;
tiledlayout(1,2)
nexttile
distMatrix = 1-neurCorrsMean(:,:,1);
% corrFit(corrFit==1) = nan;
plotDistMatrix = distMatrix;
plotDistMatrix(plotDistMatrix<0.0001) = nan;
imagesc(distMatrix);
set(gca,'XTick',1:7)
set(gca,'YTick',1:7)

set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
box off
set(gca,'FontSize',14)
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')

cH = colorbar;
cH.Label.String = '1 - Ave Correlation';
nexttile
[~, dendPerm] = customDendrogram(linkage(squareform(distMatrix)),[],gca,{'Linewidth',2,'Color','k'});
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
ylabel('1 - Ave Correlation')

set(gca,'xlim',[0 8])
box off
set(gca,'FontSize',14)
set(gca,'TickDir','out')
set(gca,'linewidth',2)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'Color','w')


figure;
tiledlayout(1,2)
nexttile
distMatrix = 1-emgCorrsSingle(:,:,1);
% corrFit(corrFit==1) = nan;
plotDistMatrix = distMatrix;
plotDistMatrix(plotDistMatrix<0.0001) = nan;
imagesc(distMatrix);
set(gca,'XTick',1:7)
set(gca,'YTick',1:7)

set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
box off
set(gca,'FontSize',14)
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')

cH = colorbar;
cH.Label.String = '1 - Ave Correlation';
set(gca,'clim',[0 0.9764])

nexttile
[~, dendPerm] = customDendrogram(linkage(squareform(distMatrix)),[],gca,{'Linewidth',2,'Color','k'});
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
ylabel('1 - Ave Correlation')

set(gca,'xlim',[0 8])
set(gca,'ylim',[0 1])
box off
set(gca,'FontSize',14)
set(gca,'TickDir','out')
set(gca,'linewidth',2)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'Color','w')


minRegionCrossings = min([sum(threshCrossingRegions==1), sum(threshCrossingRegions==2), sum(threshCrossingRegions==3)]);

figure('Units','pixels','OuterPosition',[100 100 1400 600],'Color','w');
tiledlayout(1,3);
for iRegion = 1:3
    nexttile
    regionInds = find(threshCrossingRegions == iRegion);
%     regionInds = regionInds(randperm(length(regionInds),minRegionCrossings));

    controlFRNeurCortex = squeeze(mean(controlFR(:,:,1:size(cortexFRs,1)),3));
    controlFRNeurStriatum = squeeze(mean(controlFR(:,:,size(cortexFRs,1)+1:end),3));

    xVals = -1*periTransTimes(1):periTransTimes(2);

    hold on
    cortexFR = squeeze(mean(periCrossingFR(:,regionInds,1:size(cortexFRs,1)),3));
    striatumFR = squeeze(mean(periCrossingFR(:,regionInds,size(cortexFRs,1)+1:end),3));
    shadedErrorBar(xVals,nanmean(cortexFR,2),nanstd(cortexFR,[],2)/sqrt(length(regionInds)),'lineProps',{'color','r'})
    shadedErrorBar(xVals,nanmean(striatumFR,2),nanstd(striatumFR,[],2)/sqrt(length(regionInds)),'lineProps',{'color','b'})
%     shadedErrorBar([],nanmean(controlFRNeurConcat,2),nanstd(controlFRNeurConcat,[],2)/sqrt(size(controlFRNeurConcat,2)))
    shadedErrorBar(xVals,nanmean(periCrossingEMG(:,regionInds),2),nanstd(periCrossingEMG(:,regionInds),[],2)/sqrt(length(regionInds)),'lineProps',{'color','k'})
%     line([periTransTimes(1) periTransTimes(1)],get(gca,'ylim'),'color','r')

    legend('Cortex','Striatum','PL','box','off')
    title(['Region ' num2str(iRegion)])
    xlabel('Time (ms)')
    set(gca,'fontsize',14)
    set(gca,'LineWidth',2)

end

figure;
tileH = tiledlayout(1,2);
tileH.Padding = 'compact';
tileH.TileSpacing = 'compact';

nexttile
yyaxis left
shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),aveCtx,semCtx,'lineProps',{'color',plotColors(1,:)})
hold on
plot(-1*periTransTimes(1):periTransTimes(2),mean(aveCtxShift)+2*std(aveCtxShift),'--','Color',plotColors(1,:),'LineWidth',2)
plot(-1*periTransTimes(1):periTransTimes(2),mean(aveCtxShift)-2*std(aveCtxShift),'--','Color',plotColors(1,:),'LineWidth',2)
% shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),mean(aveCtxShift),2*std(aveCtxShift),'lineProps',{'color',plotColors(1,:),'marker','none'})
xlim([-100 200])
ylim([-0.5 3.5])
ylabel('Population mean rate (spks/s)')
line([0 0],[-0.5 3.5],'color','k','linewidth',2,'linestyle','--')
yyaxis right
shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),aveEMG,semEMG,'lineProps',{'color',[0 0 0]})
ylim([-0.2 1.2])
ylabel('PL EMG (a.u.)')
set(gca,'Ycolor','k')
xlabel('Time (ms)')
set(gca,'LineWidth',2)
set(gca,'tickdir','out')
set(gca,'fontsize',13)
set(gca,'Xcolor','k')
set(gcf,'color','w')

nexttile

yyaxis left
shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),aveStr,semStr,'lineProps',{'color',plotColors(2,:)})
hold on
plot(-1*periTransTimes(1):periTransTimes(2),mean(aveStrShift)+2*std(aveStrShift),'--','Color',plotColors(2,:),'LineWidth',2)
plot(-1*periTransTimes(1):periTransTimes(2),mean(aveStrShift)-2*std(aveStrShift),'--','Color',plotColors(2,:),'LineWidth',2)
% shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),mean(aveStrShift),2*std(aveStrShift),'lineProps',{'color',plotColors(2,:),'marker','none'})
xlim([-100 200])
ylim([-0.5 3.5])
% ylim([-0.3 0.3])
% set(gca,'ytick',-0.3:0.15:0.3)
set(gca,'Ycolor',plotColors(2,:))
ylabel('Population mean rate (spks/s)')
line([0 0],[-0.5 3.5],'color','k','linewidth',2,'linestyle','--')
yyaxis right
shadedErrorBar(-1*periTransTimes(1):periTransTimes(2),aveEMG,semEMG,'lineProps',{'color',[0 0 0]})
ylim([-0.2 1.2])
ylabel('PL EMG (a.u.)')
set(gca,'Ycolor','k')
xlabel('Time (ms)')
set(gca,'LineWidth',2)
set(gca,'tickdir','out')
set(gca,'fontsize',13)
set(gca,'Xcolor','k')
set(gcf,'color','w')


% for iCross = 1:length(threshCrossings)-1
%     
%     figure;
%     crossEMG = emg(:,threshCrossings(iCross)-preThreshDuration : threshCrossings(iCross)+preThreshDuration); 
%     plot(crossEMG')
%     
%     %EMG data
%     emgTrialData(:,:,iCross) = emg(:,threshCrossings(iCross)-neurWindow(1)*20 : ...
%         threshCrossings(iCross)+neurWindow(2)*20);
% 
%     %get the video corresponding to the thresh crossing
%     for iVid = 1:length(frameEMGSamples{1})
%         
%         if frameEMGSamples{1}{iVid}(end) < threshCrossings(iCross)
%             continue
%         else
%             trialVid(iCross) = iVid;
%             trialFrame(iCross) = find(threshCrossings(iCross) <= frameEMGSamples{1}{iVid},1);
%             break;
%         end
%         
%     end
%     
%     %get the neural data surrounding the activation
%     neurSample = round((threshCrossings(iCross)-frameEMGSamples{1}{1}(1))/20*30) ...
%         + frameNeuropixelSamples{1}{1}(1);
%     
%     cortexTrialData(:,:,iCross) = cortexFRs(:,round((neurSample-neurWindow(1)*30)/(30*neurBinSize)) : ...
%         round((neurSample+neurWindow(2)*30)/(30*neurBinSize)))';
%     striatumTrialData(:,:,iCross) = striatumFRs(:,round((neurSample-neurWindow(1)*30)/(30*neurBinSize)) : ...
%         round((neurSample+neurWindow(2)*30)/(30*neurBinSize)))';
%     
%     if any(any(isnan(cortexTrialData(:,:,iCross)))) || any(any(isnan(striatumTrialData(:,:,iCross))))
%         trialHasNan(iCross) = 1;
%     else
%         trialHasNan(iCross) = 0;
%     end
%     
% end
% 
% plot(-1*neurWindow(1):neurBinSize:neurWindow(2), nanmean(nanmean(cortexTrialData,3)')*1000/neurBinSize)
% hold on
% yyaxis right
% plot(-1*neurWindow(1):neurBinSize:neurWindow(2),nanmean(nanmean(striatumTrialData,3)')*1000/neurBinSize)
% plot(-1*neurWindow(1):(1/20):neurWindow(2),nanmean(nanmean(emgTrialData,3))/80*1000/neurBinSize+1)
%
%


function [threshCrossings, reducInds, crossingNeurInds] = convertCrossInds(threshCrossings,origDownsampEMGInd,baseDir)
% convert the threshold crossings in EMG data to the corresponding indices
% in UMAP reduction data and neural data

% convert to UMAP time points, and make sure the initations exists in
% umap time points
badCross = [];
reducInds = [];
for iCross = 1:length(threshCrossings)
    thisReducInd = find(origDownsampEMGInd == threshCrossings(iCross));
    if isempty(thisReducInd)
        badCross(iCross) = 1;
        reducInds(iCross) = 0;
    else
        badCross(iCross) = 0;
        reducInds(iCross) = thisReducInd;
    end
end

threshCrossings(find(badCross)) = [];
reducInds(find(badCross)) = [];

% find the corresponding neural index for each threshold crossing
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
currentDir = pwd;
cd(fullfile(baseDir,'ProcessedData'))
crossingNeurInds = round(NeurEMGSync(threshCrossings*20,...
    frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
cd(currentDir)

% don't use nans
nanNeurInds = isnan(crossingNeurInds);
threshCrossings(nanNeurInds) = [];
reducInds(nanNeurInds) = [];
crossingNeurInds(nanNeurInds) = [];

end

% 


function [periCrossingFR, periCrossingEMG, crossRegion, nanNeurCrossings] = getDataFromCrossings(threshCrossings,crossingNeurInds,reducInds,...
    periTransTimes,allFRs,downsampEMG,regionAssignmentsFiltered,regionWatershedLabels)
% function to get neural, emg, and behavior class label data from each
% threshold crossing

periCrossingFR = zeros(length(threshCrossings),sum(periTransTimes)+1,size(allFRs,1));
periCrossingEMG = zeros(length(threshCrossings),sum(periTransTimes)+1,size(downsampEMG,1));

for iCross = 1:length(threshCrossings)

    %determine which UMAP region it's assigned to (based on plurality of
    %the surrounding time points)
    reducCrossInds = reducInds(iCross)-periTransTimes(1) : reducInds(iCross)+periTransTimes(2);
    regionDistribution = histcounts(regionAssignmentsFiltered(reducCrossInds),[regionWatershedLabels regionWatershedLabels(end)+1]);
    [~,highestRegion] = max(regionDistribution);

    %only use regions where at least half of the time points is in that
    %region
    if regionDistribution(highestRegion) <= round(length(reducCrossInds)/2)
        crossRegion(iCross) = 0;
    else
        crossRegion(iCross) = highestRegion;
    end

    %get neural activity
    for iNeuron = 1:size(allFRs,1)
        periCrossingFR(iCross,:,iNeuron) = allFRs(iNeuron,...
            crossingNeurInds(iCross)-periTransTimes(1):crossingNeurInds(iCross)+periTransTimes(2));
    end

    %get EMG activity
    for iMuscle = 1:size(downsampEMG,1)
        periCrossingEMG(iCross,:,iMuscle) = downsampEMG(iMuscle,...
            threshCrossings(iCross)-periTransTimes(1):threshCrossings(iCross)+periTransTimes(2));
    end

end

nanNeurCrossings = find(any(isnan(periCrossingFR(:,:))'));

end



% 
