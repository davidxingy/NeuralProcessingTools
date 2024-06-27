clear
close all

sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

iSess = 3;

load(fullfile(sessionDirs{iSess},'ProcessedData','BehaviorAnnotations','limbEventLabels.mat'))
load(fullfile(sessionDirs{iSess},'ProcessedData','videoSyncFrames.mat'))
load(fullfile(sessionDirs{iSess},'ProcessedData','emg1ms.mat'))
load(fullfile(sessionDirs{iSess},'ProcessedData','UMAP.mat'),'origDownsampEMGInd','regionAssignmentsFiltered','regionWatershedLabels')

channelNames = {'Right Biceps','Right Triceps','Right ECR','Right PL','Right TA','Right Gastr','Left Triceps','Left Biceps'};

downsampleFactor = 5;
nShifts = 100;

switch downsampleFactor
    case 5
        load(fullfile(sessionDirs{iSess},'ProcessedData','NeuralFiringRates5msBins30msGauss.mat'),'allFRs','striatumInds','cortexInds')
    case 1
        load(fullfile(sessionDirs{iSess},'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'),'allFRs','striatumInds','cortexInds')
end

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

regionWatershedLabelsPerm = regionWatershedLabels(allBehvAlignPerms(iSess,:));

usedBehvs = {'Walk - Right Paw Strike','Climb - Right Paw Off','Eat - Right Paw Move','Groom - Right Paw Lift','Rear - Right Paw Lift','Jump Down - Right Paw Strike'};
iUsedBehv = 1;

% downsample EMG if needed
nDownSampBins = floor(size(downsampEMG,2)/downsampleFactor);
downsampEMGReshape = reshape(downsampEMG(:,1:nDownSampBins*downsampleFactor),size(downsampEMG,1),downsampleFactor,nDownSampBins);
downSampMoreEMG = squeeze(nanmean(downsampEMGReshape,2));

% go through each tracked limb event
for iBehv = 1:length(limbEventFrames)

    iEvent = 1;
    periCrossingEMG{iBehv} = [];
    for iVid = 1:length(limbEventFrames{iBehv})
        for iFrame = 1:length(limbEventFrames{iBehv}{iVid})

            %get ind in emg index
            ind = round(frameEMGSamples{1}{iVid}(limbEventFrames{iBehv}{iVid}(iFrame))/20);
            limbEventInds{iBehv}(iEvent) = ind;

            %get emg, 200ms before and after
            periCrossingEMG{iBehv}(iEvent,:,:) = downsampEMG(:,ind-200:ind+200);
            iEvent = iEvent+1;

        end
    end

    periCrossingEMGMeans{iBehv} = squeeze(mean(periCrossingEMG{iBehv},1));

    % some limb events don't have any tracked times
    if isempty(periCrossingEMG{iBehv})
        continue
    end

%     %plot prealigned data
%     for iMusc = 1:8
%         if iMusc<=4
%             nexttile(iMusc)
%         else
%             nexttile(iMusc+4)
%         end
%         plot(-100:200,squeeze(periCrossingEMG{iBehv}(:,iMusc,101:end))','color',[0 0 1 0.2],'LineWidth',0.5)
%         hold on
%         plot(-100:200,periCrossingEMGMeans{iBehv}(iMusc,101:end),'k','LineWidth',2)
%     end

    %do some alignment by finding max correlation with the trial averaged
    for iEvent = 1:length(limbEventInds{iBehv})
        for iMusc = 1:8

            %calc correlation
            [xcf, lags] = crosscorr(periCrossingEMGMeans{iBehv}(iMusc,:),...
                periCrossingEMG{iBehv}(iEvent,iMusc,:),100);

            %use offset with max correlation
            [maxCorr(iMusc,iEvent), maxLagInd(iMusc,iEvent)] = max(xcf);
            preAlignedCorrs(iMusc,iEvent,iBehv) = xcf(101);

        end
    end

    [~, muscleToAlgin] = max(mean(maxCorr(1:4,:)'));

    % shift tracked index to the index with highest alignment
    for iEvent = 1:length(limbEventInds{iBehv})
        limbEventAlignedInds{iBehv}(iEvent) = limbEventInds{iBehv}(iEvent) + lags(maxLagInd(muscleToAlgin,iEvent));

        periCrossingEMGAligned{iBehv}(iEvent,:,:) = downsampEMG(:,...
            limbEventAlignedInds{iBehv}(iEvent)-100:limbEventAlignedInds{iBehv}(iEvent)+200);
    end

    % sometimes the alignment converged two events to the same time point,
    % do some post-filtering to remove events that are too close to each other
    badEvents = find(diff(limbEventAlignedInds{iBehv})<round(300/downsampleFactor));
    limbEventAlignedInds{iBehv}(badEvents) = [];
    periCrossingEMGAligned{iBehv}(badEvents,:,:) = [];

    periCrossingEMGAlignedMeans{iBehv} = squeeze(mean(periCrossingEMGAligned{iBehv},1));

    % get a sense of how correlated they are now
    for iEvent = 1:length(limbEventAlignedInds{iBehv})
        for iMusc = 1:8
            [xcf, lags] = crosscorr(periCrossingEMGAlignedMeans{iBehv}(iMusc,:),...
                periCrossingEMGAligned{iBehv}(iEvent,iMusc,:),100);

            postAlignedCorrs(iMusc,iEvent,iBehv) = xcf(101);
        end
    end

    % plot post-aligned trials
    figure('color','w')
    tlH = tiledlayout(2,4);
    tlH.Title.String = limbEventNames{iBehv};

    for iMusc = 1:8
        if iMusc<=4
            nexttile(iMusc)
        else
            nexttile(iMusc)
        end
        plot(-100:200,squeeze(periCrossingEMGAligned{iBehv}(:,iMusc,:))','color',[0 0 1 0.2],'LineWidth',0.5)
        hold on
        plot(-100:200,periCrossingEMGAlignedMeans{iBehv}(iMusc,:),'k','LineWidth',2)
        title(channelNames{iMusc})
    end

    % now for behaviors we want to use get neural and EMG data and save
    if any(contains(usedBehvs,limbEventNames{iBehv}))

        %align neur data to reduction
        currentDir = pwd;
        cd(fullfile(sessionDirs{iSess},'ProcessedData'))
        eventNeurInds = round(NeurEMGSync(limbEventAlignedInds{iBehv}*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
        maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/30)-1;
        neurIndsToUse = find(eventNeurInds < maxNeurSamples);
        cd(currentDir)

        emgInds = round(limbEventAlignedInds{iBehv}(neurIndsToUse)/downsampleFactor);
        neurInds = round(eventNeurInds(neurIndsToUse)/downsampleFactor);

        emgEvents{iUsedBehv} = [];
        striatumEvents{iUsedBehv} = [];
        cortexEvents{iUsedBehv} = [];

        for iEvent = 1:length(emgInds)

            % get 100ms before and 200ms after the event
            emgEvents{iUsedBehv}(iEvent,:,:) = downSampMoreEMG(1:4,emgInds(iEvent)-(100/downsampleFactor):emgInds(iEvent)+(200/downsampleFactor))';
            striatumEvents{iUsedBehv}(iEvent,:,:) = allFRs(1:length(striatumInds),neurInds(iEvent)-round(100/downsampleFactor):neurInds(iEvent)+round(200/downsampleFactor))';
            cortexEvents{iUsedBehv}(iEvent,:,:) = allFRs(length(striatumInds)+1:end,neurInds(iEvent)-round(100/downsampleFactor):neurInds(iEvent)+round(200/downsampleFactor))';

        end

        %now don't use any events where there are nans in the data
        nanInds = find(any(any(isnan(emgEvents{iUsedBehv}),2),3) | any(any(isnan(striatumEvents{iUsedBehv}),2),3) | any(any(isnan(cortexEvents{iUsedBehv}),2),3));
        emgEvents{iUsedBehv}(nanInds,:,:) = [];
        striatumEvents{iUsedBehv}(nanInds,:,:) = [];
        cortexEvents{iUsedBehv}(nanInds,:,:) = [];
        emgInds(nanInds) = [];
        neurInds(nanInds) = [];

        allEmgInds{iUsedBehv} = emgInds;
        allNeurInds{iUsedBehv} = neurInds;

        iUsedBehv = iUsedBehv + 1;

    end

%     % now that we have the template, go through and match with the rest of
%     % the time series
%     for iMusc = 1:4
%         for i=1:round(size(downsampEMG,2)/50)-400
% 
%             muscCorrs(iMusc,i) = corr(periCrossingEMGAlignedMeans{iBehv}(iMusc,1:201)',downsampEMG(iMusc,i:i+200)');
% 
%         end
%     end
% 
%     periCrossingEMGMatch{iBehv} = [];
%     periCrossingEMGMatchMeans{iBehv} = [];
%     matchedInd{iBehv} = [];
%     matchedUmapRegions{iBehv} = [];
% 
%     potentialInds = find(sum(muscCorrs>0.6)>3);
% 
%     if isempty(potentialInds)
%         continue
%     end
%     
%     potentialBlockStarts = [potentialInds(1) potentialInds(find(diff(potentialInds)>1)+1)];
%     potentialBlockEnds = [potentialInds(find(diff(potentialInds)>1)) potentialInds(end)];
% 
%     for iBlock = 1:length(potentialBlockStarts)
% 
%         [~, blockBestInd] = max(mean(muscCorrs(:,potentialBlockStarts(iBlock):potentialBlockEnds(iBlock))));
%         matchedInd{iBehv}(iBlock) = potentialBlockStarts(iBlock) + blockBestInd - 1;
% 
%     end
% 
%     matchedInd{iBehv}(diff(matchedInd{iBehv})<400) = [];
% 
%     for iEvent = 1:length(matchedInd{iBehv})
% 
%         % get EMG data
%         periCrossingEMGMatch{iBehv}(iEvent,:,:) = downsampEMG(:,...
%             matchedInd{iBehv}(iEvent):matchedInd{iBehv}(iEvent)+300);
% 
%         % get umap region
%         umapInd = find(origDownsampEMGInd==matchedInd{iBehv}(iEvent));
%         if isempty(umapInd)
%             matchedUmapRegions{iBehv}(iEvent) = 0;
%         else
%             matchedUmapRegions{iBehv}(iEvent) = find(regionWatershedLabelsPerm==regionAssignmentsFiltered(umapInd));
%         end
% 
%     end
% 
%     periCrossingEMGMatchMeans{iBehv} = squeeze(mean(periCrossingEMGMatch{iBehv},1));
% 
%     figure
%     tlH = tiledlayout(2,4);
%     tlH.Title.String = limbEventNames{iBehv};
% 
%     for iMusc = 1:8
%         if iMusc<=4
%             nexttile(iMusc)
%         else
%             nexttile(iMusc)
%         end
%         plot(-100:200,squeeze(periCrossingEMGMatch{iBehv}(:,iMusc,:))','color',[0 0 1 0.2],'LineWidth',0.5)
%         hold on
%         plot(-100:200,periCrossingEMGMatchMeans{iBehv}(iMusc,:),'k','LineWidth',2)
%     end

end

nTrials = cellfun(@(x) size(x,1),emgEvents);
behvLabelsDiv = {};
for iBehv = 1:length(nTrials)
    behvLabelsDiv{iBehv} = ones(1,nTrials(iBehv))*iBehv;
end
behvLabels = cat(2,behvLabelsDiv{:});

% divide up the trials into train, test, and validation sets
for iBehv=1:length(unique(behvLabels))

    behvInds{iBehv} = find(behvLabels==iBehv);
    behvInds{iBehv} = behvInds{iBehv}(randperm(length(behvInds{iBehv})));
    cvBlocks = divideBlocks(behvInds{iBehv},2);
    indBlocks(iBehv) = cvBlocks(1);

    validationInds{iBehv} = indBlocks{iBehv}(1:floor(length(indBlocks{iBehv})/2));
    testInds{iBehv} = indBlocks{iBehv}(floor(length(indBlocks{iBehv})/2)+1:floor(length(indBlocks{iBehv})/2)*2);

    trainInds{iBehv} = setdiff(behvInds{iBehv},indBlocks{iBehv});
    
end

validationTrials = cat(2,validationInds{:});
testTrials = cat(2,testInds{:});
trainTrials = cat(2,trainInds{:});

%also do shift controls
for iShift = 1:nShifts

    % don't shift EMG
    emgShiftEvents(:,iShift) = emgEvents;

    % we want to keep the same number of trials in the test and validation set so keep doing
    % shifts until we get one which doesn't have nans
    nansInShift = true;
    shiftAttemptcounter = 1;
    while nansInShift

        %min shift by 30s
        neurLength = size(allFRs,2);
        shiftAmount(iShift) = randi(neurLength - round(60000/downsampleFactor), 1) + round(30000/downsampleFactor);
        neurShiftInds = cat(2,allNeurInds{:}) + shiftAmount(iShift);

        %move points shifted past the end to the beginning (circularly
        %shifting)
        neurShiftInds(neurShiftInds >= neurLength) = neurShiftInds(neurShiftInds >= neurLength)-neurLength;

        %need space for pre and post data
        tooCloseToEdges = neurShiftInds <= round(100/downsampleFactor) | neurShiftInds >= neurLength - round(200/downsampleFactor);
        

        if any(tooCloseToEdges)
            shiftAttemptcounter = shiftAttemptcounter+1;
            continue
        end

        %now get the data around the shifted inds
        for iEvent = 1:length(neurShiftInds)

            allStriatumShiftEvents{iShift}(iEvent,:,:) = allFRs(1:length(striatumInds),neurShiftInds(iEvent)-round(100/downsampleFactor):neurShiftInds(iEvent)+round(200/downsampleFactor))';
            allCortexShiftEvents{iShift}(iEvent,:,:) = allFRs(length(striatumInds)+1:end,neurShiftInds(iEvent)-round(100/downsampleFactor):neurShiftInds(iEvent)+round(200/downsampleFactor))';

        end

        %now don't use any events where there are nans in the data 
        nanInds = any(any(isnan(allStriatumShiftEvents{iShift}),2),3) | any(any(isnan(allCortexShiftEvents{iShift}),2),3);

        if ~any(nanInds(testTrials))
            nansInShift = false;

            %split back into behaviors
            for iBehv = 1:length(unique(behvLabels))
            
                striatumShiftEvents{iBehv,iShift} = allStriatumShiftEvents{iShift}(behvLabels==iBehv,:,:);
                cortexShiftEvents{iBehv,iShift} = allCortexShiftEvents{iShift}(behvLabels==iBehv,:,:);
                [~, badTrainInds] = intersect(trainInds{iBehv}, find(nanInds));
                [~, badTestInds] = intersect(testInds{iBehv}, find(nanInds));
                [~, badValidationInds] = intersect(validationInds{iBehv}, find(nanInds));

                trainShiftInds{iBehv,iShift} = trainInds{iBehv};
                trainShiftInds{iBehv,iShift}(badTrainInds) = [];
                testShiftInds{iBehv,iShift} = testInds{iBehv};
                testShiftInds{iBehv,iShift}(badTestInds) = [];
                validationShiftInds{iBehv,iShift} = validationInds{iBehv};
                validationShiftInds{iBehv,iShift}(badValidationInds) = [];

            end

        else
            shiftAttemptcounter = shiftAttemptcounter+1;
        end

    end

end

% to keep things consistent, use the minimum value across all shifts
minNumTrain = min(cellfun(@length,trainShiftInds),[],2);
minNumTest = min(cellfun(@length,validationShiftInds),[],2);
minNumVal = min(cellfun(@length,validationShiftInds),[],2);
for iBehv = 1:length(unique(behvLabels))

    trainInds{iBehv} = trainInds{iBehv}(1:minNumTrain(iBehv));
    testInds{iBehv} = testInds{iBehv}(1:minNumTest(iBehv));
    validationInds{iBehv} = validationInds{iBehv}(1:minNumVal(iBehv));

    for iShift = 1:100
        trainShiftInds{iBehv,iShift} = trainShiftInds{iBehv,iShift}(1:minNumTrain(iBehv));
        testShiftInds{iBehv,iShift} = testShiftInds{iBehv,iShift}(1:minNumTest(iBehv));
        validationShiftInds{iBehv,iShift} = validationShiftInds{iBehv,iShift}(1:minNumVal(iBehv));
    end

end

for iShift = 1:100
    strShiftData{iShift} = cat(1,striatumShiftEvents{:,iShift});
    ctxShiftData{iShift} = cat(1,cortexShiftEvents{:,iShift});
    trainShiftTrials{iShift} = cat(2,trainShiftInds{:,iShift});
    testShiftTrials{iShift} = cat(2,testShiftInds{:,iShift});
    validationShiftTrials{iShift} = cat(2,validationShiftInds{:,iShift});
end

% concatenate data
trainTrials = cat(2,trainInds{:});
testTrials = cat(2,testInds{:});
validationTrials = cat(2,validationInds{:});

emgData = cat(1,emgEvents{:});
emgDataMean = squeeze(mean(emgData,1));
tmp = cellfun(@(x) mean(x,1),emgEvents,'UniformOutput',false);
emgDataRegionMean = cat(1,tmp{:});

strData = cat(1,striatumEvents{:});
strDataMean = squeeze(mean(strData,1));
tmp = cellfun(@(x) mean(x,1),striatumEvents,'UniformOutput',false);
strDataRegionMean = cat(1,tmp{:});

ctxData = cat(1,cortexEvents{:});
ctxDataMean = squeeze(mean(ctxData,1));
tmp = cellfun(@(x) mean(x,1),cortexEvents,'UniformOutput',false);
ctxDataRegionMean = cat(1,tmp{:});

%%


save('limbEventData5ms','activationCue','emgData','emgDataMean','emgDataRegionMean','strData','strDataMean','strDataRegionMean','behvLabels',...
    'usedBehvs','ctxData','ctxDataMean','ctxDataRegionMean','testInds','testTrials','trainInds','trainTrials','validationInds','validationTrials');

save('limbEventData5msShift','activationCue','emgData','emgDataMean','emgDataRegionMean','strShiftData','behvLabels',...
    'usedBehvs','ctxShiftData','testShiftInds','testShiftTrials','trainShiftInds','trainShiftTrials','validationShiftInds','validationShiftTrials','-v7.3');

% % regionWatershedLabelsPerm = regionWatershedLabels([1 2 3 4 5 6 7]);
% % for iRegion = 1:7
% %     regionInds{iRegion} = find(regionAssignmentsFiltered==regionWatershedLabelsPerm(iRegion));
% %     regionBoutStarts{iRegion} = [regionInds{iRegion}(1) regionInds{iRegion}(find(diff(regionInds{iRegion})>1)+1)];
% %     regionBoutEnds{iRegion} = [regionInds{iRegion}(find(diff(regionInds{iRegion})>1)) regionInds{iRegion}(end)];
% %     goodBouts = find((regionBoutEnds{iRegion} - regionBoutStarts{iRegion}) > 1000);
% %     regionBoutStarts{iRegion} = regionBoutStarts{iRegion}(goodBouts);
% %     regionBoutEnds{iRegion} = regionBoutEnds{iRegion}(goodBouts);
% %     for iBout = 1:length(regionBoutStarts{iRegion})
% %         startInd = origDownsampEMGInd(regionBoutStarts{iRegion}(iBout));
% %         endInd = origDownsampEMGInd(regionBoutEnds{iRegion}(iBout));
% %         if isempty(startInd) || isempty(endInd)
% %             regionBoutStartVids{iRegion}(iBout) = 0;
% %             regionBoutStartFrames{iRegion}(iBout) = 0;
% % 
% %             regionBoutEndVids{iRegion}(iBout) = 0;
% %             regionBoutEndFrames{iRegion}(iBout) = 0;
% %             continue
% %         end
% % 
% %         for iVid = 1:length(frameEMGSamples{1})
% %             if frameEMGSamples{1}{iVid}(end) < startInd*20
% %                 continue
% %             else
% %                 regionBoutStartVids{iRegion}(iBout) = iVid;
% %                 regionBoutStartFrames{iRegion}(iBout) = find(frameEMGSamples{1}{iVid}>=startInd*20,1,'first');
% %                 break;
% %             end
% %         end
% % 
% %         for iVid = 1:length(frameEMGSamples{1})
% %             if frameEMGSamples{1}{iVid}(end) < endInd*20
% %                 continue
% %             else
% %                 regionBoutEndVids{iRegion}(iBout) = iVid;
% %                 regionBoutEndFrames{iRegion}(iBout) = find(frameEMGSamples{1}{iVid}>=endInd*20,1,'first');
% %                 break;
% %             end
% %         end
% % 
% %     end
% % end


% 
