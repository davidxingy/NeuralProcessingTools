function optoTrialAnalysis(baseDirs)

if nargin < 1

    D036Sessions = {'X:\David\ArenaRecordings\D036-101623-ArenaRecording','X:\David\ArenaRecordings\D036-101723-ArenaRecording',...
        'X:\David\ArenaRecordings\D036-101823-ArenaRecording','X:\David\ArenaRecordings\D036-101923-ArenaRecording',...
        'X:\David\ArenaRecordings\D036-102023-ArenaRecording'};

    D040Sessions = {'X:\David\ArenaRecordings\D040-110223-ArenaRecording','X:\David\ArenaRecordings\D040-110323-ArenaRecording',...
        'X:\David\ArenaRecordings\D040-110423-ArenaRecording','X:\David\ArenaRecordings\D040-110623-ArenaRecording',...
        'X:\David\ArenaRecordings\D040-110723-ArenaRecording','X:\David\ArenaRecordings\D040-110823-ArenaRecording'};

    D041Sessions = {'X:\David\ArenaRecordings\D041-121123-ArenaRecording','X:\David\ArenaRecordings\D041-121223-ArenaRecording',...
        'X:\David\ArenaRecordings\D041-121323-ArenaRecording','X:\David\ArenaRecordings\D041-121423-ArenaRecording',...
        'X:\David\ArenaRecordings\D041-121523-ArenaRecording'};

    baseDirs = {D036Sessions, D040Sessions, D041Sessions};

end

emgTrigStims = false;

baselineDuration = 100;
postStimDuration = 200;

shiftControlOffsetMin = 600;
shiftControlOffsetRange = 200;

maxRandOffset = 200;

%initialize session data containers
nStims = [];
stimWindow = {};
stimWindowNoBaseine = {};
randControlWindow = {};
randControlWindowNoBaseine = {};
shiftControlWindow = {};
shiftControlWindowNoBaseine = {};

scrap(baseDirs,emgTrigStims)

for iSession = 1:length(baseDirs)

    if iSession == 1
%         load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'analyzedBehaviors','behvLabelsNoArt','regionBehvAssignments','regionWatershedLabels')
    end

%     load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'reduction','origDownsampEMGInd','regionAssignmentsNoBound')
    load(fullfile(baseDirs{iSession},'ProcessedData','EMG1ms'))

    sessTotalPulses = length(downsampleLaserOnsetInds);

    controlStimsShift = downsampleLaserOnsetInds(1:sessTotalPulses) + randi(shiftControlOffsetRange,1,sessTotalPulses)+shiftControlOffsetMin;
    % controlStimsRand = randperm(size(downsampEMG,2)-(baselineDuration+postStimDuration+2),nTotalPulses-1) + baselineDuration;
    
    usedPulseInds = downsampleLaserOnsetInds(1:sessTotalPulses);
    if iSession == 2
%         usedPulseInds = downsampleLaserOnsetInds(1:3900);
    end

    normalizedEMG = (downsampEMG-mean(downsampEMG,2))./std(downsampEMG,[],2);

    %Go through each UMAP behavior region
    for iRegion = 1:length(regionWatershedLabels)

        regionInds = find(regionAssignmentsNoBound == regionWatershedLabels(iRegion));

        %get names/labels for the region
        regionName{iRegion} = string();
        for iBehv = 1:length(regionBehvAssignments{iRegion})
            regionName{iRegion}{iBehv} = analyzedBehaviors{regionBehvAssignments{iRegion}(iBehv)};
        end
        regionName{iRegion} = join(regionName{iRegion},'\');

        %get continuous chunks of time within the umap region to sample random control stims from
        %don't use time points where there is stim window
        allStimWindowInds = {};
        for iPulse = 1:length(usedPulseInds)
            stimInd = usedPulseInds(iPulse);
            allStimWindowInds{iPulse} = stimInd-postStimDuration:stimInd+postStimDuration+baselineDuration;
        end
        regionIndsNoStim = setdiff(regionInds, cat(2,allStimWindowInds{:}));
        regionChunkStarts = [1 find(diff(regionIndsNoStim)>1)+1];
        regionChunkEnds = [find(diff(regionIndsNoStim)>1) length(regionIndsNoStim)];

        %get the start and stop times of each chunk for this region
        goodRegionChunks = find((regionChunkEnds - regionChunkStarts) > (baselineDuration + postStimDuration + maxRandOffset));
        regionChunkStarts = regionChunkStarts(goodRegionChunks);
        regionChunkEnds = regionChunkEnds(goodRegionChunks);

        %go through and sample random times in the chunks for control trials
        iChunk = 1;
        controlStimsRand = [];
        for iControlPulse = 1:length(usedPulseInds)

            if iControlPulse == 1
                controlStimsRand(iControlPulse) = regionIndsNoStim(regionChunkStarts(1))+baselineDuration+randperm(maxRandOffset,1);
            else
                if controlStimsRand(iControlPulse-1) + (baselineDuration + 2*postStimDuration + maxRandOffset) > ...
                        regionIndsNoStim(regionChunkEnds(iChunk))
                    iChunk = iChunk + 1;

                    if iChunk > length(goodRegionChunks)
                        break
                    end

                    controlStimsRand(iControlPulse) = regionIndsNoStim(regionChunkStarts(iChunk))+baselineDuration+randperm(maxRandOffset,1);
                else
                    controlStimsRand(iControlPulse) = controlStimsRand(iControlPulse-1) + (baselineDuration + postStimDuration + maxRandOffset);
                end
            end

        end

        nRandControlPulses = length(controlStimsRand);
        
        %now that we have the control pulses, go through each stim trial
        %and get EMG window

        %initialize data
        windowEMGCell = {};
        windowEMGNoBaselineCell = {};
        controlRandEMGCell = {};
        controlRandEMGNoBaselineCell = {};
        controlShiftEMGCell = {};
        controlShiftEMGNoBaselineCell = {};
        goodPulseInd = 1;

        for iPulse = 1:length(usedPulseInds)

            stimInd = usedPulseInds(iPulse);

            windowInds = stimInd-baselineDuration:stimInd+postStimDuration;
            controlRandInds = controlStimsRand(goodPulseInd)-baselineDuration:controlStimsRand(goodPulseInd)+postStimDuration;
            controlShiftInds = controlStimsShift(iPulse)-baselineDuration:controlStimsShift(iPulse)+postStimDuration;

            %if there is artifact in the stim EMG window, don't use
            if any(downsampleRemovedInds>=windowInds(1) & downsampleRemovedInds<=windowInds(end)) %|| ~any(stimInd==behvAllEMGInds)
                continue
            end

            %if there is no corresponding UMAP point to either the pulse
            %index, or the stop of the window, don't use
            if isempty(find(windowInds(1)==origDownsampEMGInd)) || isempty(find(windowInds(end)==origDownsampEMGInd)) || ...
                    isempty(find(windowInds(baselineDuration)==origDownsampEMGInd))
                continue
            end

            %only use stims where the stim timepoint and the stop of the
            %window are within the region
            if regionAssignmentsNoBound(find(windowInds(baselineDuration)==origDownsampEMGInd)) ~= regionWatershedLabels(iRegion) || ...
                    regionAssignmentsNoBound(find(windowInds(end)==origDownsampEMGInd)) ~= regionWatershedLabels(iRegion)
                continue
            end

%             if regionAssignmentsNoBound(controlShiftInds(1)) ~= regionWatershedLabels(iRegion) || regionAssignmentsNoBound(controlShiftInds(end)) ~= regionWatershedLabels(iRegion)
%                 continue
%             end

            %now get the actual data
            windowEMGCell{goodPulseInd} = normalizedEMG(:,windowInds);
            windowEMGNoBaselineCell{goodPulseInd} = windowEMGCell{goodPulseInd} - mean(windowEMGCell{goodPulseInd}(:,1:baselineDuration),2);

            controlRandEMGCell{goodPulseInd} = normalizedEMG(:,controlRandInds);
            controlRandEMGNoBaselineCell{goodPulseInd} = controlRandEMGCell{goodPulseInd} - mean(controlRandEMGCell{goodPulseInd}(:,1:baselineDuration),2);

            %just do a sanity check, if our random control index sampling
            %worked, this should never trigger
            if regionAssignmentsNoBound(controlRandInds(1)) ~= regionWatershedLabels(iRegion) || regionAssignmentsNoBound(controlRandInds(end)) ~= regionWatershedLabels(iRegion)
                warning('Rand Control inds not in region!')
            end

            controlShiftEMGCell{goodPulseInd} = normalizedEMG(:,controlShiftInds);
            controlShiftEMGNoBaselineCell{goodPulseInd} = controlShiftEMGCell{goodPulseInd} - mean(controlShiftEMGCell{goodPulseInd}(:,1:baselineDuration),2);

            goodPulseInd = goodPulseInd + 1;

            if goodPulseInd > nRandControlPulses
                warning('Not enough rand control trials!')
                break
            end

        end

        %now rearrange all the data
        windowEMG = cat(3,windowEMGCell{:});
        windowEMGNoBaseline = cat(3,windowEMGNoBaselineCell{:});
        controlRandEMG = cat(3,controlRandEMGCell{:});
        controlRandEMGNoBaseline = cat(3,controlRandEMGNoBaselineCell{:});
        controlShiftEMG = cat(3,controlShiftEMGCell{:});
        controlShiftEMGNoBaseline = cat(3,controlShiftEMGNoBaselineCell{:});

        %do outlier removal
        %arrange the baseline data to the right format
        stimBaseline = permute(windowEMG(:,1:baselineDuration,:), [3 2 1]);
        stimBaseline = stimBaseline(:,:);
        randControlBaseline = permute(controlRandEMG(:,1:baselineDuration,:), [3 2 1]);
        randControlBaseline = randControlBaseline(:,:);
        shiftControlBaseline = permute(controlShiftEMG(:,1:baselineDuration,:), [3 2 1]);
        shiftControlBaseline = randControlBaseline(:,:);

        %find the outliers (subfunction)
        randOutlierTrials = removeOptoOutliers(cat(1,stimBaseline, randControlBaseline));
        shiftOutlierTrials = removeOptoOutliers(cat(1,stimBaseline, shiftControlBaseline));

        %remove the outliers
        randOutlierStimTrials = randOutlierTrials(randOutlierTrials<=size(stimBaseline,1));
        randOutlierControlTrials = randOutlierTrials(randOutlierTrials>size(stimBaseline,1))-size(stimBaseline,1);
        shiftOutlierStimTrials = shiftOutlierTrials(shiftOutlierTrials<=size(stimBaseline,1));
        shiftOutlierControlTrials = shiftOutlierTrials(shiftOutlierTrials>size(stimBaseline,1))-size(stimBaseline,1);

        %for now use rand control trials
        outlierStimTrials = randOutlierStimTrials;
        windowEMG(:,:,outlierStimTrials) = [];
        controlRandEMG(:,:,randOutlierControlTrials) = [];
        controlShiftEMG(:,:,shiftOutlierControlTrials) = [];
        windowEMGNoBaseline(:,:,outlierStimTrials) = [];
        controlRandEMGNoBaseline(:,:,randOutlierControlTrials) = [];
        controlShiftEMGNoBaseline(:,:,shiftOutlierControlTrials) = [];

        %add data for the session
        nStims(iSession,iRegion) = size(windowEMG,3);
        stimWindow{iSession,iRegion} = windowEMG;
        stimWindowNoBaseine{iSession,iRegion} = windowEMGNoBaseline;
        randControlWindow{iSession,iRegion} = controlRandEMG;
        randControlWindowNoBaseine{iSession,iRegion} = controlRandEMGNoBaseline;
        shiftControlWindow{iSession,iRegion} = controlShiftEMG;
        shiftControlWindowNoBaseine{iSession,iRegion} = controlShiftEMGNoBaseline;

    end

end

%go through each region an consolidate all the trials across sessions
for iRegion = 1:length(regionWatershedLabels)
    allStimWindow{iRegion} = cat(3,stimWindow{:,iRegion});
    allStimWindowNoBaseine{iRegion} = cat(3,stimWindowNoBaseine{:,iRegion});
    allRandControlWindow{iRegion} = cat(3,randControlWindow{:,iRegion});
    allRandControlWindowNoBaseine{iRegion} = cat(3,randControlWindowNoBaseine{:,iRegion});
    allShiftControlWindow{iRegion} = cat(3,shiftControlWindow{:,iRegion});
    allShiftControlWindowNoBaseine{iRegion} = cat(3,shiftControlWindowNoBaseine{:,iRegion});

    %do some plotting
    figure('Units','pixels','OuterPosition',[50 50 1400 800])
    t = tiledlayout(2,4);
    title(t,[char(regionName{iRegion}) ', N Stims = ' num2str(sum(nStims(:,iRegion)))])
    for iChan = 1:8
        nexttile
        emgChan = iChan;
        chanData = squeeze(allStimWindowNoBaseine{iRegion}(emgChan,:,:));
        shiftData = squeeze(allShiftControlWindowNoBaseine{iRegion}(emgChan,:,:));
        randData = squeeze(allRandControlWindowNoBaseine{iRegion}(emgChan,:,:));
        shadedErrorBar(-baselineDuration:postStimDuration,mean(chanData,2),std(chanData,[],2)/sqrt(size(chanData,2)),'transparent',1,'lineProps',{'color','b'})
        hold on
        shadedErrorBar(-baselineDuration:postStimDuration,mean(shiftData,2),std(shiftData,[],2)/sqrt(size(shiftData,2)),'transparent',1,'lineProps',{'color','m'})
        shadedErrorBar(-baselineDuration:postStimDuration,mean(randData,2),std(randData,[],2)/sqrt(size(randData,2)),'transparent',1,'lineProps',{'color','k'})
        line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',2)
        title(channelNames{iChan})
        ylabel('Muscle Std.Dev')
        xlim([-baselineDuration postStimDuration])
        xlabel('Time (ms)')
    end
    set(gcf,'color','w')

end


nexttile
emgChan = 9;
test = squeeze(windowEMGNorm(emgChan,:,1:end));
testControlShift = squeeze(controlStimEMGNorm(emgChan,:,:));
testControlRand = squeeze(controlRandEMGNorm(emgChan,:,:));
shadedErrorBar(-baselineDuration:postStimDuration,mean(test,2)-mean(testControlShift,2),std(test,[],2)/sqrt(size(test,2)),'transparent',1,'lineProps',{'color','b'})
hold on
shadedErrorBar(-baselineDuration:postStimDuration,zeros(1,size(testControlShift,1)),std(testControlShift,[],2)/sqrt(size(testControlShift,2)),'transparent',1,'lineProps',{'color','r'})
line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',2)
title('PL')

nexttile
emgChan = 2;
test = squeeze(windowEMGNorm(emgChan,:,1:end));
testControlShift = squeeze(controlStimEMGNorm(emgChan,:,:));
testControlRand = squeeze(controlRandEMGNorm(emgChan,:,:));
shadedErrorBar(-baselineDuration:postStimDuration,mean(test,2)-mean(testControlShift,2),std(test,[],2)/sqrt(size(test,2)),'transparent',1,'lineProps',{'color','b'})
% hold on
shadedErrorBar(-baselineDuration:postStimDuration,zeros(1,size(testControlShift,1)),std(testControlShift,[],2)/sqrt(size(testControlShift,2)),'transparent',1,'lineProps',{'color','r'})
line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',2)
title('Biceps')
xlabel('Time (ms)')
ylabel('Muscle Std.Dev')

nexttile
emgChan = 6;
test = squeeze(windowEMGNorm(emgChan,:,1:end));
testControlShift = squeeze(controlStimEMGNorm(emgChan,:,:));
testControlRand = squeeze(controlRandEMGNorm(emgChan,:,:));
shadedErrorBar(-baselineDuration:postStimDuration,mean(test,2)-mean(testControlShift,2),std(test,[],2)/sqrt(size(test,2)),'transparent',1,'lineProps',{'color','b'})
hold on
shadedErrorBar(-baselineDuration:postStimDuration,zeros(1,size(testControlShift,1)),std(testControlShift,[],2)/sqrt(size(testControlShift,2)),'transparent',1,'lineProps',{'color','r'})
line([0 0],get(gca,'ylim'),'linestyle','--','linewidth',2)
title('ECR')
xlabel('Time (ms)')

title(t,['D034 ' analyzedBehaviors{iRegion} ', num stims = ' num2str(size(windowEMGNorm,3))])

set(gcf,'color','w')

end




function outlierTrials = removeOptoOutliers(baselineEMG)


nTrials = size(baselineEMG,1);
baselineDists = squareform(pdist(baselineEMG));
aveBaselineDists = mean(baselineDists);
cutoff = median(aveBaselineDists) + 2*std(aveBaselineDists);
outlierTrials = find(aveBaselineDists >= cutoff);

end


%

function scrap(baseDirs, emgTrigStims)

baselineDuration = 100;
postStimDuration = 200;

shiftControlOffsetMin = 600;
shiftControlOffsetRange = 200;

maxRandOffset = 200;

% animal 3 I had some bad behavior towards the end of each of the sessions
% since I put some honey on him to induce grooming, but he sort of freaked
% out when I did those, so don't use those time points
animal3Cutoffs = [4419351, 5707322, 4937315, 4479614, 4286971]; %in ms

%initialize session data containers
nStims = [];
stimWindow = {};
stimWindowNoBaseine = {};
randControlWindow = {};
randControlWindowNoBaseine = {};
shiftControlWindow = {};
shiftControlWindowNoBaseine = {};

% clear testEMG testEMGBaseline testInds
% nPulses = 4000;
% figNonNorm = figure;
% figNoBase = figure;
% 
% iRegion = 1;
% regionReducInds = find(regionAssignmentsNoBound == regionWatershedLabels(iRegion));
% regionInds = origDownsampEMGInd(regionReducInds);
% 
% goodTestInds = setdiff(1:length(normalizedEMG),cat(2,allStimWindowInds{:}));
% goodTestInds = intersect(goodTestInds,regionInds);
% goodTestInds((goodTestInds <= baselineDuration) | (goodTestInds >= length(normalizedEMG)-postStimDuration)) = [];
% 
% for iSamp = 1:10
%     testInds(iSamp,:) = goodTestInds(randperm(length(goodTestInds),nPulses));
%     for iPulse = 1:length(testInds(iSamp,:))
% 
%         testWindInds = testInds(iSamp,iPulse)-baselineDuration:testInds(iSamp,iPulse)+postStimDuration;
%         testEMG(:,:,iPulse,iSamp) = normalizedEMG(:,testWindInds);
%         testEMGBaseline(:,:,iPulse,iSamp) = normalizedEMG(:,testWindInds) - mean(normalizedEMG(:,testWindInds(1:baselineDuration)),2);
% 
%     end
% 
%     figure(figNonNorm)
%     hold on
%     emgChan = 1;
%     plotColors = turbo(10);
%     shadedErrorBar([],mean(testEMG(emgChan,:,iSamp,:),4),std(testEMG(emgChan,:,iSamp,:),[],4)/sqrt(length(testInds(iSamp,:))),'lineprops',{'Color',plotColors(iSamp,:)})
%
%     figure(figNoBase)
%     hold on
%     shadedErrorBar([],mean(testEMGBaseline(emgChan,:,iSamp,:),4),std(testEMGBaseline(emgChan,:,iSamp,:),[],4)/sqrt(length(testInds(iSamp,:))),'lineprops',{'Color',plotColors(iSamp,:)})
% end

for iAnimal = 1:length(baseDirs)

    animalDirs = baseDirs{iAnimal};

    for iSession = 1:length(animalDirs)

        if iSession == 1
            load(fullfile(animalDirs{iSession},'ProcessedData','UMAP'),'analyzedBehaviors','behvLabelsNoArt','regionBehvAssignments','regionWatershedLabels')
        end

        % load behavior parcellation
        load(fullfile(animalDirs{iSession},'ProcessedData','UMAP'),'reduction','origDownsampEMGInd','regionAssignmentsNoBound','regionAssignmentsFiltered')
        % load emg and laser stim
        load(fullfile(animalDirs{iSession},'ProcessedData','EMG1ms'))

        if iAnimal == 3
            badSectionStart = find(downsampleLaserOnsetInds > animal3Cutoffs(iSession),1);
            downsampleLaserOnsetInds(badSectionStart:end) = [];
            if emgTrigStims
                downsampleControlOnsetInds(badSectionStart:end) = [];
            end
        end

        sessTotalPulses = length(downsampleLaserOnsetInds);

        % get control stim pulses, defined as a time some random amount
        % after the actual pulses (defined by shiftControlOffsetMin and
        % shiftControlOffsetRange)
        controlStimsShift = downsampleLaserOnsetInds(1:sessTotalPulses) + randi(shiftControlOffsetRange,1,sessTotalPulses)+shiftControlOffsetMin;
        % controlStimsRand = randperm(size(downsampEMG,2)-(baselineDuration+postStimDuration+2),nTotalPulses-1) + baselineDuration;

        % don't use first or last few pulses since video and intan stopped
        % recording at different times
        usedPulseInds = downsampleLaserOnsetInds(2:sessTotalPulses-100);

        % for emg-triggered stims rather than random stims, we already have
        % experimentally defined control pulses
        if emgTrigStims
            controlPulseInds = downsampleControlOnsetInds(2:length(downsampleControlOnsetInds)-1);
        end

        if iSession == 1
            usedPulseInds = downsampleLaserOnsetInds(2:sessTotalPulses-10);
            %         controlPulseInds = downsampleControlOnsetInds(2:sessTotalPulses-25);
        end

        % z-score the EMG
        normalizedEMG = (downsampEMG)./std(downsampEMG,[],2);

        % now go through each pulse
        pulseRegion{iAnimal}{iSession} = [];
        for iPulse = 1:length(usedPulseInds)

            % get the window surrounding each pulse, get both the straight
            % up z-scored EMG levels, as well as a version where the
            % baseline (pre-stim) is subtracted off
            allWindInds = usedPulseInds(iPulse)-baselineDuration:usedPulseInds(iPulse)+postStimDuration;
            allStim{iAnimal}{iSession}(:,:,iPulse) = normalizedEMG(:,allWindInds);
            allStimBaseline{iAnimal}{iSession}(:,:,iPulse) = allStim{iAnimal}{iSession}(:,:,iPulse) - mean(allStim{iAnimal}{iSession}(:,1:baselineDuration,iPulse),2);

            % for controls, shift by a random amount
            shiftAmount = 600*(randi(2)*2-3);
            allStimShift{iAnimal}{iSession}(:,:,iPulse) = normalizedEMG(:,allWindInds+shiftAmount);
            allStimShiftBaseline{iAnimal}{iSession}(:,:,iPulse) = allStimShift{iAnimal}{iSession}(:,:,iPulse) - mean(allStimShift{iAnimal}{iSession}(:,1:baselineDuration,iPulse),2);

            %assigning regions to each pulse based on which region the window majority is in
            windowReducInds = find(origDownsampEMGInd>=usedPulseInds(iPulse)-baselineDuration & origDownsampEMGInd<= usedPulseInds(iPulse)+postStimDuration);
            regionDistribution = histcounts(regionAssignmentsFiltered(windowReducInds),[regionWatershedLabels regionWatershedLabels(end)+1])/301;
            
            [maxFraction,maxRegion] = max(regionDistribution);

            %for D041 (animal 3) I accidentally switched the order of climb
            %up and climb down, so fix that here
            if iAnimal == 3
                if maxRegion == 2
                    maxRegion = 1;
                elseif maxRegion == 1
                    maxRegion = 2;
                end
            end
            pulseRegion{iAnimal}{iSession}(iPulse) = maxRegion;
        end

        % get data for emg-triggered controls
        if emgTrigStims
            for iPulse = 1:length(controlPulseInds)
        
                allWindInds = controlPulseInds(iPulse)-baselineDuration:controlPulseInds(iPulse)+postStimDuration;
                allControl{iSession}(:,:,iPulse) = normalizedEMG(:,allWindInds);
                allControlBaseline{iSession}(:,:,iPulse) = allControl{iSession}(:,:,iPulse) - mean(allControl{iSession}(:,1:baselineDuration,iPulse),2);
        
            end
        end

        % now go through each region and get the data the corresponds to
        % that region
        for iRegion = 1:length(regionWatershedLabels)

            %get names/labels for the region
            regionName{iRegion} = string();
            for iBehv = 1:length(regionBehvAssignments{iRegion})
                regionName{iRegion}{iBehv} = analyzedBehaviors{regionBehvAssignments{iRegion}(iBehv)};
            end
            regionName{iRegion} = join(regionName{iRegion},'\');

            %use stim pulses that are within this region
            regionPulseInds = usedPulseInds(pulseRegion{iAnimal}{iSession}==iRegion);

            pulseInd = 1;
            for iPulse = 1:length(regionPulseInds)

                %get shift ind for control
                shiftAmount = 600*(randi(2)*2-3);
                shiftInd = regionPulseInds(iPulse)+shiftAmount;

                %make sure that the shift ind is also in the same behavior
                %region, if not, then don't use this control pulse (and the
                %corresponding real pulse)
                windowReducInds = find(origDownsampEMGInd>=shiftInd-baselineDuration & origDownsampEMGInd<= shiftInd+postStimDuration);
                regionDistribution = histcounts(regionAssignmentsFiltered(windowReducInds),[regionWatershedLabels regionWatershedLabels(end)+1])/301;
                [maxFraction,maxRegion] = max(regionDistribution);

                if iAnimal == 3
                    if maxRegion == 2
                        maxRegion = 1;
                    elseif maxRegion == 1
                        maxRegion = 2;
                    end
                end
                if maxRegion ~= iRegion
                    continue
                end

                % if both control and real pulse are in the desired region,
                % extract the data window
                regionWindInds = regionPulseInds(iPulse)-baselineDuration:regionPulseInds(iPulse)+postStimDuration;
                regionStim{iAnimal}{iSession,iRegion}(:,:,pulseInd) = normalizedEMG(:,regionWindInds);
                regionStimBaseline{iAnimal}{iSession,iRegion}(:,:,pulseInd) = regionStim{iAnimal}{iSession,iRegion}(:,:,pulseInd) - mean(regionStim{iAnimal}{iSession,iRegion}(:,1:baselineDuration,pulseInd),2);

                regionStimShift{iAnimal}{iSession,iRegion}(:,:,pulseInd) = normalizedEMG(:,regionWindInds+shiftAmount);
                regionStimShiftBaseline{iAnimal}{iSession,iRegion}(:,:,pulseInd) = regionStimShift{iAnimal}{iSession,iRegion}(:,:,pulseInd) - mean(regionStimShift{iAnimal}{iSession,iRegion}(:,1:baselineDuration,pulseInd),2);
                pulseInd = pulseInd+1;

            end

%             % remove outlier trials
%             for iMusc = 1:size(normalizedEMG,1)
%                 outlierTrials = removeOptoOutliers(squeeze(cat(3,regionStimBaseline{iAnimal}{iSession,iRegion}(iMusc,1:baselineDuration-1,:), ...
%                     regionStimShiftBaseline{iAnimal}{iSession,iRegion}(iMusc,1:baselineDuration-1,:)))');
%                 stimOutlierTrials = outlierTrials(outlierTrials <= size(regionStimBaseline{iAnimal}{iSession,iRegion},3));
%                 shiftOutlierTrials = outlierTrials(outlierTrials > size(regionStimBaseline{iAnimal}{iSession,iRegion},3)) - size(regionStimBaseline{iAnimal}{iSession,iRegion},3);
%                 regionStimBaseline{iAnimal}{iSession,iRegion}(:,:,unique([stimOutlierTrials shiftOutlierTrials])) = [];
%                 regionStimShiftBaseline{iAnimal}{iSession,iRegion}(:,:,unique([stimOutlierTrials shiftOutlierTrials])) = [];
%             end

        end

    end
end

regionLabels = {'Climb Up', 'Climb Down','Jump/Misc','Walk','Rear/Still/Misc','Groom','Eat'};

% combine across sessions and animals
sessionStimComb = {};
sessionShiftComb = {};
sessionStimOrigComb = {};
stimComb = {};
shiftComb = {};
stimOrigComb = {};
musclesToUse = 1:4;
for iRegion = 1:length(regionWatershedLabels)

    for iAnimal = 1:length(regionStimBaseline)

        tmpComb = cat(3,regionStimBaseline{iAnimal}{:,iRegion});
        tmpComb = permute(tmpComb(musclesToUse,:,:),[2 1 3]);
        sessionStimComb{iRegion,iAnimal} = tmpComb(:,:);

        tmpComb = cat(3,regionStimShiftBaseline{iAnimal}{:,iRegion});
        tmpComb = permute(tmpComb(musclesToUse,:,:),[2 1 3]);
        sessionShiftComb{iRegion,iAnimal} = tmpComb(:,:);

        tmpComb = cat(3,regionStim{iAnimal}{:,iRegion});
        tmpComb = permute(tmpComb(musclesToUse,:,:),[2 1 3]);
        sessionStimOrigComb{iRegion,iAnimal} = tmpComb(:,:);

    end

    stimComb{iRegion} = cat(2,sessionStimComb{iRegion,:});
    shiftComb{iRegion} = cat(2,sessionShiftComb{iRegion,:});
    stimOrigComb{iRegion} = cat(2,sessionStimOrigComb{iRegion,:});

    figure;
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(stimComb{iRegion},2),std(stimComb{iRegion},[],2)/sqrt(size(stimComb{iRegion},2)),'lineProps',{'r'})
    hold on
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(shiftComb{iRegion},2),std(shiftComb{iRegion},[],2)/sqrt(size(shiftComb{iRegion},2)),'lineProps',{'k'})
    line([0 0],get(gca,'ylim'),'color','k','linewidth',2,'linestyle','--')
    title([regionLabels{iRegion} ', ' num2str(size(stimComb{iRegion},2)) ' stims'])

    diffComb{iRegion} = shiftComb{iRegion} - stimComb{iRegion};
    figure;
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(diffComb{iRegion},2),std(diffComb{iRegion},[],2)/sqrt(size(diffComb{iRegion},2)),'lineProps',{'b'})
    line([0 0],get(gca,'ylim'),'color','k','linewidth',2,'linestyle','--')
    title([regionLabels{iRegion} ', ' num2str(size(diffComb{iRegion},2)) ' stims'])

    baseline{iRegion} = mean(diffComb{iRegion}(1:baselineDuration,:));
    first35{iRegion} = mean(diffComb{iRegion}(baselineDuration+1:baselineDuration+35,:));
    second35{iRegion} = mean(diffComb{iRegion}(baselineDuration+36:baselineDuration+70,:));
    third35{iRegion} = mean(diffComb{iRegion}(baselineDuration+71:baselineDuration+105,:));

    first100{iRegion} = mean(diffComb{iRegion}(baselineDuration+1:baselineDuration+100,:));
    fullDuration{iRegion} = mean(diffComb{iRegion}(baselineDuration:end,:));

end

% make figure from one animal 4 muscles, all behaviors
exampleFigH = figure;
exampleTileH = tiledlayout(1,4);
exampleTileH.Padding = 'compact'; exampleTileH.TileSpacing = 'Compact';
muscleLabels = {'Right Elbow Flexor','Right Elbow Extensor','Right Wrist Flexor','Right Wrist Extensor'};
muscPerm = [1 2 4 3];
for iMusc = 1:4
    nexttile
    hold on
    sampleStim = cat(3,regionStimBaseline{1}{:});
    sampleStim = squeeze(sampleStim(muscPerm(iMusc),:,:));
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(sampleStim,2),std(sampleStim,[],2)/sqrt(size(sampleStim,2)),'lineProps',{'color',lines(1)},'patchSaturation',0.5);
    sampleShift = cat(3,regionStimShiftBaseline{1}{:});
    sampleShift = squeeze(sampleShift(muscPerm(iMusc),:,:));
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(sampleShift,2),std(sampleShift,[],2)/sqrt(size(sampleShift,2)),'lineProps',{'color','k'},'patchSaturation',0.5);
    line([0 0],[-0.1 0.05],'color',lines(1),'linewidth',2,'linestyle','--')
    patchH = patch([0 0 50 50],[0.04 0.05 0.05 0.04],lines(1));
    patchH.FaceAlpha = 0.5;
    patchH.EdgeColor = 'none';
    set(gca,'FontSize',12)
    set(gca,'linewidth',1.5)
    set(gca,'TickDir','out')
    set(gca,'XColor','k')
    set(gca,'YColor','k')
    set(gca,'TickLength',[0.02 0.05])
    set(gcf,'Color','w')
    xlabel('Time (ms)')
    ylim([-0.10 0.05])
    title(muscleLabels{iMusc})
    if iMusc == 1
        ylabel({'Normalized EMG (std dev)'})
        legend('Laser','Control','box','off')
    end
end

% make bar graph showing overall difference in EMG
summaryFigH = figure('color','w');
plotColors = lines(2);
stimCat = cellfun(@(x) cat(3,x{:}), regionStimBaseline,'un',0);
allStimsCat = cat(3,stimCat{:});
allStimsCat = permute(allStimsCat(1:4,:,:),[2 1 3]);
allStimsCat = allStimsCat(:,:);
shiftCat = cellfun(@(x) cat(3,x{:}), regionStimShiftBaseline,'un',0);
allShiftsCat = cat(3,shiftCat{:});
allShiftsCat = permute(allShiftsCat(1:4,:,:),[2 1 3]);
allShiftsCat = allShiftsCat(:,:);

stimLaserMean = mean(mean(allStimsCat(101:end,:)));
stimLaserSem = std(mean(allStimsCat(101:end,:)))/sqrt(size(allStimsCat,2));
shiftLaserMean = mean(mean(allShiftsCat(101:end,:)));
shiftLaserSem = std(mean(allShiftsCat(101:end,:)))/sqrt(size(allShiftsCat,2));
bar(1,stimLaserMean*-1,0.6,'FaceColor',plotColors(1,:),'EdgeColor','none')
hold on
bar(2,shiftLaserMean*-1,0.6,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
errorbar([1,2],-1*[stimLaserMean shiftLaserMean],[stimLaserSem shiftLaserSem],'.','linewidth',1.5,'color',[0.2 0.2 0.2],'CapSize',16)
line([1 2],[0.045 0.045],'linewidth',1.5','color','k')
text(1.48,0.046,'*','FontSize',14)
ylim([-0.003 0.05])
xlim([0.5 2.5])
set(gca,'YTick',0:0.025:0.05)
set(gca,'YTickLabel',{'0','-0.025','-0.05'})
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'Laser','Control'})
set(gca,'TickLength',[0.02 0.05])
ylabel('Normalized EMG')
set(gca,'FontSize',12)
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')
box off

% do statistical test
[~,pAll] = ttest(mean(allStimsCat(101:end,:)),mean(allShiftsCat(101:end,:)))



combFigH = figure;
hold on;
sepFigH = figure;
tileH = tiledlayout(1,7);
tileH.Padding = 'compact'; tileH.TileSpacing = 'Compact';
plotColors = turbo(length(regionWatershedLabels)+1);
for iRegion = 1:length(regionWatershedLabels)
    figure(combFigH)
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(diffComb{iRegion},2),std(diffComb{iRegion},[],2)/sqrt(size(diffComb{iRegion},2)),'lineProps',{'color',plotColors(iRegion,:)},'patchSaturation',0.5);

    figure(sepFigH)
    nexttile
    shadedErrorBar(-1*baselineDuration:postStimDuration,mean(diffComb{iRegion},2),std(diffComb{iRegion},[],2)/sqrt(size(diffComb{iRegion},2)),'lineProps',{'color',lines(1)},'patchSaturation',0.5);

    xlim([-50,150])
%     if iRegion == 6
%         ylim([-0.07,0.32])
%         line([0 0],[-0.07,0.32],'color',lines(1),'linewidth',2,'linestyle','--')
%         patchH = patch([0 0 50 50],[0.28 0.32 0.32 0.28],lines(1));
%         patchH.FaceAlpha = 0.5;
%         patchH.EdgeColor = 'none';
%         line([-50 -50],[-0.07 0.03],'color','k','linewidth',2)
%         line([-50 0],[-0.07 -0.07],'color','k','linewidth',2)
%     else
        ylim([-0.03,0.18])
        line([0 0],[-0.03,0.18],'color',lines(1),'linewidth',2,'linestyle','--')
        patchH = patch([0 0 50 50],[0.165 0.18 0.18 0.165],lines(1));
        patchH.FaceAlpha = 0.5;
        patchH.EdgeColor = 'none';
        line([-50 -50],[-0.03 0.02],'color','k','linewidth',2)
        line([-50 0],[-0.03 -0.03],'color','k','linewidth',2)
%     end
    set(gca,'FontSize',12)
    set(gca,'linewidth',1.5)
    set(gca,'TickDir','out')
    set(gca,'XColor','k')
    set(gca,'YColor','k')
    set(gcf,'Color','w')
    xlabel('Time (ms)')
    ylabel({'Change over control','(std dev)'})
    title(regionLabels{iRegion})
%     axis off
end
figure(combFigH)
xlim([-50,150])
ylim([-0.06,0.3])
line([0 0],get(gca,'ylim'),'color','k','linewidth',2,'linestyle','--')
legend(regionLabels,'box','off','FontSize',12)
set(gca,'FontSize',14)
set(gca,'linewidth',1.5)
set(gca,'TickDir','out')
set(gcf,'Color','w')
xlabel('Time (ms)')
ylabel('Change over control (std dev)')

% plot bar graph of change in muscle activity relative to control for
% different time periods
figure
plotColors = lines(3);
hold on
bar(2:5:35,cellfun(@mean,[first35; first100]'),'EdgeColor','none')
errorbar(1.29:5:35,cellfun(@mean,first35),cellfun(@std,first35)./sqrt(cellfun(@length,first35)),'.','color','k','LineWidth',1.5)
errorbar(2.71:5:35,cellfun(@mean,first100),cellfun(@std,first100)./sqrt(cellfun(@length,first100)),'.','color','k','LineWidth',1.5)
set(gca,'XTick',2:5:35)
set(gca,'XTickLabel', regionLabels)
ylim([-0.02 0.2])
ylabel('Change over control (std dev)')
legend('35ms','100ms','box', 'off','fontsize',12)
set(gca,'FontSize',12)
set(gca,'linewidth',1)
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'Color','w')
set(gca,'XTickLabelRotation',0)

% do t-tests
for iBehv = 1:7
    [~,pBehvs35(iBehv)] = ttest(first35{iBehv});
    [~,pBehvs100(iBehv)] = ttest(first100{iBehv});
end



% plot change in muscle activity from baseline vs raw baseline muscle
% activity
rawBaseline = cellfun(@(x) mean(x(1:100,:),1) ,stimOrigComb, 'un',0);
fullDiff = cellfun(@(x) mean(x(101:end,:),1) ,diffComb, 'un',0);

figure
plotColors = turbo(7);
hold on;
for iRegion = 1:length(regionWatershedLabels)
    errH(iRegion) = errorbar(mean(rawBaseline{iRegion}),mean(fullDiff{iRegion}),...
        std(rawBaseline{iRegion})/sqrt(length(rawBaseline{iRegion})),std(fullDiff{iRegion})/sqrt(length(fullDiff{iRegion})),...
        'both','.','Color',plotColors(iRegion,:),'LineWidth',1.5);
end

legendH = legend(errH,regionLabels,'Box','off','FontSize',12);
for iLabel = 1:length(legendH.String)
    legendH.String{iLabel} = ['\color[rgb]{' num2str(errH(iLabel).Color) '} ' legendH.String{iLabel}];
end
xlabel('Baseline EMG Level (std. dev.)')
ylabel('Laser EMG change over control (std. dev.)')
set(gca,'FontSize',12)
set(gca,'linewidth',1.5)
set(gca,'TickDir','out')
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'Color','w')


% look at changes in stim response across time
stimMeans = cellfun(@(x) cellfun(@(y) mean(y,3),x,'un',0), allStimBaseline, 'un', 0);
stimSems = cellfun(@(x) cellfun(@(y) std(y,[],3)/sqrt(size(y,3)),x,'un',0), allStimBaseline, 'un', 0);
stimNums = cellfun(@(x) cellfun(@(y) size(y,3),x), allStimBaseline, 'un', 0);

% iRegion = 1;
% 
% stimMeans = cellfun(@(x) cellfun(@(y) mean(y,3),x(:,iRegion),'un',0), regionStimBaseline, 'un', 0);
% stimSems = cellfun(@(x) cellfun(@(y) std(y,[],3)/sqrt(size(y,3)),x(:,iRegion),'un',0), regionStimBaseline, 'un', 0);
% stimNums = cellfun(@(x) cellfun(@(y) size(y,3),x(:,iRegion)), regionStimBaseline, 'un', 0);

animalLabels = {'D036','D040','D041'};
muscleLabels = {'Elbow Flexor','Elbow Extensor','Wrist Flexor','Wrist Extensor'};
muscPlotOrder = [1 2 4 3];
plotColors = turbo(6);

tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for iAnimal = 1:3
    for iMusc = 1:4
        nexttile
        plotH = [];
        for iSess = 1:length(stimMeans{iAnimal})
            hold on
            shadedErrorBar(-100:200,stimMeans{iAnimal}{iSess}(muscPlotOrder(iMusc),:),stimSems{iAnimal}{iSess}(muscPlotOrder(iMusc),:),'lineprops',{'color',plotColors(iSess,:)});
%             ylim([-0.16 0.06])
            line([0 0],[-0.16 0.06],'color',[0.5 0.5 0.5 0.5],'linewidth',2,'linestyle','--')
            if iMusc == 1
                text(-90,-0.075-(iSess-1)*0.015,['Day ' num2str(iSess) ', ' num2str(stimNums{iAnimal}(iSess)) ' Stims'],'color',plotColors(iSess,:))
            end
        end
        title([animalLabels{iAnimal} ' ' muscleLabels{iMusc}])
        xlabel('Time (ms)')
        ylabel('Std Dev')
        set(gca,'LineWidth',1.5)
        set(gca,'FontSize',13)
        set(gca,'TickDir','out')
        box off
    end
end
set(gcf,'color','w');


% and do it with all animals consolidated
allStimBaseline5Days = allStimBaseline;
allStimBaseline5Days{2}(6) = [];
allStimBaselineCat = cat(1,allStimBaseline5Days{:});
for iSess = 1:size(allStimBaselineCat,2)
    allStimMeans{iSess} = mean(cat(3,allStimBaselineCat{:,iSess}),3);
    allStimNums(iSess) = size(cat(3,allStimBaselineCat{:,iSess}),3);
    allStimSems{iSess} = std(cat(3,allStimBaselineCat{:,iSess}),[],3)/sqrt(allStimNums(iSess));
end

figure;
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
for iMusc = 1:4
    nexttile
    plotH = [];
    for iSess = 1:length(stimMeans{iAnimal})
        hold on
        shadedErrorBar(-100:200,allStimMeans{iSess}(muscPlotOrder(iMusc),:),allStimSems{iSess}(muscPlotOrder(iMusc),:),'lineprops',{'color',plotColors(iSess,:)});
        ylim([-0.12 0.03])
        line([0 0],[-0.16 0.06],'color',[0.5 0.5 0.5 0.5],'linewidth',2,'linestyle','--')
        if iMusc == 1
            text(-45,-0.06-(iSess-1)*0.012,['Day ' num2str(iSess) ', ' num2str(allStimNums(iSess)) ' Stims'],'color',plotColors(iSess,:))
        end
    end
    xlim([-50 100])
    title(muscleLabels{iMusc})
    xlabel('Time (ms)')
    ylabel('Std Dev')
    set(gca,'LineWidth',1.5)
    set(gca,'FontSize',13)
    set(gca,'YColor','k')
    set(gca,'XColor','k')
    set(gca,'TickDir','out')
    box off
end
set(gcf,'color','w');


figure
tileH = tiledlayout(2,4,'TileSpacing','tight','Padding','tight');
plotColors = jet(length(baseDirs));
% plotColors = turbo(length(regionWatershedLabels));
% regionsToPlot = [1 4 5 7];    
for iChan = 1:8
    nexttile
    for iSession = 1:length(baseDirs)
%         nexttile
        shadedErrorBar((-1*baselineDuration):postStimDuration,mean(allStimBaseline{iSession}(iChan,:,:),3),std(allStimBaseline{iSession}(iChan,:,:),[],3)/sqrt(size(allStimBaseline{iSession},3)),'lineprops',{'Color',plotColors(iSession,:)})
        hold on
%         shadedErrorBar((-1*baselineDuration):postStimDuration,mean(allControlBaseline{iSession}(iChan,:,:),3),std(allControlBaseline{iSession}(iChan,:,:),[],3)/sqrt(size(allControlBaseline{iSession},3)),'lineprops',{'Color','k'})
% shadedErrorBar((-1*baselineDuration):postStimDuration,mean(allStimBaseline{iSession}(iChan,:,:),3) - mean(allControlBaseline{iSession}(iChan,:,:),3),...
%     std(allStimBaseline{iSession}(iChan,:,:),[],3)/sqrt(size(allStimBaseline{iSession},3)),'lineprops',{'Color',plotColors(iSession,:)})
%         shadedErrorBar((-1*baselineDuration):postStimDuration,mean(allStimShiftBaseline{iSession}(iChan,:,:),3),std(allStimShiftBaseline{iSession}(iChan,:,:),[],3)/sqrt(size(allStimShiftBaseline{iSession},3)),'lineprops',{'Color',plotColors(iSession,:)}) %plot shift controls
%         shadedErrorBar((-1*baselineDuration):postStimDuration,mean(regionStimBaselineComb{iSession}(iChan,:,:),3),std(regionStimBaselineComb{iSession}(iChan,:,:),[],3)/sqrt(size(regionStimBaselineComb{iSession},3)),'lineprops',{'Color',plotColors(iSession,:)}) %plot regions
%         shadedErrorBar((-1*baselineDuration):postStimDuration,mean(regionStimBaselineComb{regionsToPlot(iRegion)}(iChan,:,:),3),... %Plot individual regions
%             std(regionStimBaselineComb{regionsToPlot(iRegion)}(iChan,:,:),[],3)/sqrt(size(regionStimBaselineComb{regionsToPlot(iRegion)},3)),'lineprops',{'Color',plotColors(iSession,:)}) 
        if iChan == 1
            plotH(iSession) = plot(0,0,'.','color',plotColors(iSession,:),'MarkerSize',0.0001);
        end
%     end
    line([0 0],get(gca,'ylim'),'color','r','linestyle','--')
    title(channelNames{iChan})
    ylabel('Std Dev')
    xlabel('Time (ms)')
    ylim([-0.12 0.1])

    if iChan == 1
        legendNames = join([string(repmat({'Day '},1,length(baseDirs))); string(1:length(baseDirs))]','');
%         legendNames = {'Climb Up', 'Climb Down', 'Jump/Misc', 'Rear/Still/Misc', 'Walk', 'Groom', 'Eat'};
        legendH = legend(plotH,legendNames,'Box','off','FontSize',10);
        for iLabel = 1:length(legendH.String)
            legendH.String{iLabel} = ['\color[rgb]{' num2str(plotH(iLabel).Color) '} ' legendH.String{iLabel}];
        end
    end
    end
end

tileH.Title.String = 'D041';


for iRegion1 = 1:length(regionWatershedLabels)
    for iRegion2 = 1:length(regionWatershedLabels)

        for iChan = 1:4
            meanTrace1 = mean(regionStimBaselineComb{iRegion1}(iChan,:,:),3);
            meanTrace2 = mean(regionStimBaselineComb{iRegion2}(iChan,:,:),3);

            chanCorrs(iChan) = corr(meanTrace1',meanTrace2');

            maxEffect1 = max([meanTrace1 0]);
            minEffect1 = min([meanTrace1 0]);
            effectSize1(iChan) = max([abs(maxEffect1) abs(minEffect1)]);
            
            maxEffect2 = max([meanTrace2 0]);
            minEffect2 = min([meanTrace2 0]);
            effectSize2(iChan) = max([abs(maxEffect2) abs(minEffect2)]);
        end

        if iRegion1 == iRegion2
            correlations(iRegion1,iRegion2) = 0;
            effectSize(iRegion1,iRegion2) = 0;
        else
            correlations(iRegion1,iRegion2) = mean(chanCorrs);
            effectSize(iRegion1,iRegion2) = sqrt(sum(abs(effectSize1-effectSize2).^2));
        end

    end
end

figure; imagesc(effectSize)
box off
set(gcf,'color','w')

regionNames = {'Climb Up','Climb Down','Jumping/Misc','Rear/Still/Misc','Walk','Groom','Eat'};
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',1:7)
set(gca,'fontsize',14)
set(gca,'XTickLabel',regionNames)
set(gca,'YTickLabel',regionNames)

title('D040')

end
