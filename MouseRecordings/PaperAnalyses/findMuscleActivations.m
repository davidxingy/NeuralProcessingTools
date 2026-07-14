clear

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

sessionNames = {{'D020-062922-ArenaRecording'},{'D024-111022-ArenaRecording'},{'D026-032923-ArenaRecording'},...
    {'D043-020625-ArenaRecording','D043-020525-ArenaRecording','D043-020425-ArenaRecording','D043-020325-ArenaRecording'},...
    {'D047-090825-ArenaRecording','D047-090925-ArenaRecording','D047-091625-ArenaRecording','D047-091825-ArenaRecording'},...
    {'D050-120925-ArenaRecording','D050-121825-ArenaRecording','D050-120625-ArenaRecording','D050-120525-ArenaRecording'},...
    {'D054-011626-ArenaRecording','D054-012126-ArenaRecording','D054-011326-ArenaRecording','D054-012426-ArenaRecording'},...
    {'D056-012726-ArenaRecording','D056-020926-ArenaRecording','D056-012926-ArenaRecording','D056-013126-ArenaRecording'}};

nBrainRegions = [repmat({1},1,3) repmat({[2,2,2,2]},1,3) {[2 2 2]} {[2 2 2 2]}];
brainRegionNames = {'CFA', 'tjMC'};
allBehvAlignPerms = repmat(1:7,length(sessionNames),1);

usedThreshChannels = 1:6;
baselineDur = 50;

behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};

% initation detection parameters
fs = 1000;
exampleThreshChan = 4;
musclesForAnalysis = 1:4;

threshs = { {[0.8 0.8 0.8 0.8 0.8 0.8]}, {[2 0.6 0.8 0.8 0.8 0.8]}, {[0.8 0.8 0.8 0.8 1.5 1]}, ...
    {[2 1.5 1 1 2 1.5], [2 1.5 1 1 2 1.5], [2 1.5 1 0.8 2 1.5], [2 1 1 1 0.1 1.5]},...
    {[2 2 1 1.3 1.5 1], [2 2 1.2 1.2 1.2 1], [2 2 1 1 1.5 1], [2 2 1 1.3 1.5 1]},...
    {[2 1.5 1 1.2 1.6 1], [2 1.5 1.1 1.3 1.6 1], [2 1.5 1.2 1.2 1.6 1], [2 1.5 1 1.2 1.6 1]},...
    {[1.6 1.5 0.8 1.5 2 3], [1.6 0.08 1 1.5 3.5 3.5], [1.6 1.5 1 1.5 2 3], [1.6 1.5 0.8 1.5 2 3]},...
    {[1.8 0.5 3 1.2 1.5 1.2], [1.8 0.5 0.7 1.2 3.5 1.2], [1.8 0.25 3 1 3 1.2], [1.8 0.5 3 1 3 1.2]} };

% thresh = [2 2 1 1.5 2 2];
% 
% baselineChanThresh = [1 1 0.5 0.5 1 1];

% baselineChanThresh = [0.5 0.5 0.5 0.5;...
%                       0.5 0.5 0.5 0.5;...
%                       0.5 0.5 0.5 0.5];

baselineChanThresh = { {[0.6 0.6 0.6 0.6 0.6 0.6]}, {[1.5 0.5 0.6 0.6 0.6 0.6]}, {[0.6 0.6 0.6 0.6 1 0.6]}, ...
    {[1 1 0.6 0.6 1.1 1], [1 1 0.6 0.6 1.1 1], [1 0.8 0.6 0.6 1 1], [1 0.8 0.6 0.6 0.05 1]},...
    {[1.1 1 0.6 0.6 1 0.6], [1.1 1 0.7 0.7 0.8 0.5], [1.2 1.2 0.6 0.6 1 0.6], [1.1 1 0.6 0.6 1 0.6]},...
    {[1.2 0.7 0.7 0.6 1.1 0.6], [1.2 0.7 0.7 0.8 1.1 0.6], [1.2 0.7 0.8 0.6 1.1 0.6], [1.2 0.7 0.7 0.6 1.1 0.6]},...
    {[1.2 1 0.5 0.8 1.5 2], [1.2 0.04 1.5 0.8 3 3], [1.2 1 0.8 0.8 1.5 2], [1.2 1 0.5 0.8 1.5 2]},...
    {[1.2 0.3 2 0.6 1 0.7], [1.2 0.3 0.7 0.6 3 0.6], [1.2 0.15 2.2 0.6 2 0.7], [1.2 0.3 2.2 0.6 2 0.7]} };

neurBinSize = 1; %in ms

preThreshDuration = 500; %in ms %originally was 150 ms
preThreshbuffer = 50; %originally was 100 ms

periTransTimes = [500 250]; %in ms

nControlResamples = 100; %number of times to get the controls

for iAnimal = 1:length(sessionNames)

    animalSessions = sessionNames{iAnimal};

    for iSess = 1:length(animalSessions)

        % set the variables that we're going to save
        threshCrossOrigEMG = {};
        threshCrossReduc = {};
        threshCrossNeur = {};

        crossFRs = {};
        crossEMGs = {};
        crossUmapBehvs = {};
        crossClassifierBehvs = {};

        behvCrossingStr = [];
        behvCrossingCtx = [];
        behvCrossingEmg = [];

        % for shifts
        threshCrossOrigEMGShifts = {};
        threshCrossReducShifts = {};
        threshCrossNeurShifts = {};
        shiftAmounts = {};

        aveStrNoBaseShift = {};
        aveCtxNoBaseShift = {};
        aveEMGNoBaseShift = {};

        allFRShift = {};
        allEMGShift = {};
        crossUmapBehvsShift = {};
        crossClassifierBehvsShift = {};

        % and for behavior-specific
        aveStrNoBaseBehvsControl = [];
        aveCtxNoBaseBehvsControl = [];
        aveEMGNoBaseBehvsControl = [];

        allFRBehvControl = [];
        allEMGBehvControl = [];

        threshCrossOrigEMGBehvControls = [];

        sessionArtifacts = consolidateArtifactInds(animalSessions{iSess},'CFA');

        for iBrainRegion = 1:length(nBrainRegions{iAnimal}(iSess))

            filePaths = getMouseDataNames(animalSessions{iSess}(1:4),animalSessions{iSess},brainRegionNames{iBrainRegion});

            load(filePaths.UMAPFile,'regionAssignmentsFiltered','reduction','origDownsampEMGInd','regionWatershedLabels','classifierLabels')

            load(filePaths.neuronDataStruct)
            load(filePaths.npArtifactTimestamps)
            load(filePaths.VideoSyncFrames)

            animalLabel = animalSessions{iSess}(1:4);
            % behvAlignPerm = allBehvAlignPerms(iAnimal,:);
            behvAlignPerm = 1:7;
            regionAssignmentLabels = unique(regionAssignmentsFiltered);
            classifierLabelInds = unique(classifierLabels);

            % load EMG data
            load(filePaths.EMG1ms)

            % normalize by standard deviation
            downsampEMGNorm = downsampEMG./std(downsampEMG,[],2);

            % load in neural data
            load(filePaths.NeuralFiringRates1msBins10msGauss,'allFRs','cortexInds','striatumInds')

            for iThreshChan = 1:length(usedThreshChannels)

                thisThresh = threshs{iAnimal}{iSess}(iThreshChan);
                thisBaselineThresh = baselineChanThresh{iAnimal}{iSess}(iThreshChan);

                % get all threshold crossings
                threshCrossings = find(downsampEMGNorm(usedThreshChannels(iThreshChan),2:end) > thisThresh & ...
                    downsampEMGNorm(usedThreshChannels(iThreshChan),1:end-1) <= thisThresh);

                % prune detected threshold crossings
                goodCrossings = zeros(1,length(threshCrossings));
                for iCross = 1:length(threshCrossings)

                    % make sure it's not too close to the beginning or end of the session
                    if (threshCrossings(iCross)-periTransTimes(1) <= 0) || (threshCrossings(iCross)+periTransTimes(2) > size(downsampEMG,2))
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
                    if any((baselineEMG(usedThreshChannels(iThreshChan),:)) > thisBaselineThresh)
                        continue
                    end

                    %if there are multiple zeros, then that means it was artifact, don't use
                    if sum(baselineEMG(usedThreshChannels(iThreshChan),:)<1e-4) > 20
                        continue
                    end

                    %baseline should have no flucations past a certain amount
                    %     if any(max(baselineEMG,[],2)-min(baselineEMG,[],2) > baselineFlucThresh)
                    %         continue
                    %     end

                    % finally, don't use any crossings that have artifacts
                    if any(threshCrossings(iCross) == sessionArtifacts.allArtEmgInds1ms10msSmooth)
                        continue
                    end

                    goodCrossings(iCross) = 1;

                end

                threshCrossings = threshCrossings(find(goodCrossings));

                % get the thresh crossing times in neural and umap indices
                [threshCrossings, reducInds, crossingNeurInds] = convertCrossInds(threshCrossings,origDownsampEMGInd,filePaths);

                % remove any emg crossings that is outside of neural data range
                outRangeCrossings = find(crossingNeurInds>size(allFRs,2)-periTransTimes(2));
                threshCrossings(outRangeCrossings) = [];
                crossingNeurInds(outRangeCrossings) = [];
                reducInds(outRangeCrossings) = [];
                
                % now get the actual data around the crossings
                [periCrossingFR, periCrossingEMG, crossRegionUmap, crossRegionClassifier, nanNeurCrossings] = getDataFromCrossings(threshCrossings,crossingNeurInds,reducInds,...
                    periTransTimes,allFRs,downsampEMG,regionAssignmentsFiltered,regionWatershedLabels,classifierLabels,classifierLabelInds);

                reducInds(nanNeurCrossings) = [];
                threshCrossings(nanNeurCrossings) = [];
                crossingNeurInds(nanNeurCrossings) = [];
                crossRegionUmap(nanNeurCrossings) = [];
                crossRegionClassifier(nanNeurCrossings) = [];
                periCrossingFR(nanNeurCrossings,:,:) = [];
                periCrossingEMG(nanNeurCrossings,:,:) = [];

                % now save crossing times for this session
                threshCrossOrigEMG{iThreshChan} = threshCrossings;
                threshCrossReduc{iThreshChan} = reducInds;
                threshCrossNeur{iThreshChan} = crossingNeurInds;

                % save data
                crossFRs{iThreshChan} = periCrossingFR;
                crossEMGs{iThreshChan} = periCrossingEMG;
                crossUmapBehvs{iThreshChan} = crossRegionUmap;
                crossClassifierBehvs{iThreshChan} = crossRegionClassifier;

                % also do shifts
                for iShift = 1:100

                    shiftAmount(iShift) = randi(size(downsampEMG,2)-30000*2)+30000;

                    threshCrossingsShift = mod(threshCrossings + shiftAmount(iShift)-1,size(downsampEMG,2)-periTransTimes(2)-1)+1;
%                     threshCrossingsShift(threshCrossingsShift > size(downsampEMG,2) - periTransTimes(2) - 1) = ...
%                         threshCrossingsShift(threshCrossingsShift > size(downsampEMG,2) - periTransTimes(2) - 1) - ...
%                         (size(downsampEMG,2) - periTransTimes(2) - 1); %Old way of shifting before I learned about mod()

                    threshCrossingsShift(threshCrossingsShift < periTransTimes(1)+1) = [];
                    % find the corresponding neural index for each threshold crossing
                    [threshCrossingsShift, reducIndsShift, crossingNeurIndsShift] = convertCrossInds(threshCrossingsShift,origDownsampEMGInd,filePaths);

                    % remove any emg crossings that is outside of neural data range
                    outRangeCrossings = find(crossingNeurIndsShift>size(allFRs,2)-periTransTimes(2));
                    threshCrossingsShift(outRangeCrossings) = [];
                    crossingNeurIndsShift(outRangeCrossings) = [];
                    reducIndsShift(outRangeCrossings) = [];

% %                     periCrossingFRShift = zeros(length(threshCrossingsShift),sum(periTransTimes)+1,size(allFRs,1));
% %                     periCrossingEMGShift = zeros(length(threshCrossingsShift),sum(periTransTimes)+1,size(downsampEMG,1));
% % 
% %                     % get activity for the shifts
% %                     for iCross = 1:length(threshCrossingsShift)
% % 
% %                         %get neural activity
% %                         for iNeuron = 1:size(allFRs,1)
% %                             periCrossingFRShift(iCross,:,iNeuron) = allFRs(iNeuron,...
% %                                 crossingNeurIndsShift(iCross)-periTransTimes(1):crossingNeurIndsShift(iCross)+periTransTimes(2));
% %                         end
% % 
% %                         %get EMG activity
% %                         for iMuscle = 1:size(downsampEMG,1)
% %                             periCrossingEMGShift(iCross,:,iMuscle) = downsampEMG(iMuscle,...
% %                                 threshCrossingsShift(iCross)-periTransTimes(1):threshCrossingsShift(iCross)+periTransTimes(2));
% %                         end
% % 
% %                     end

                    % now get the actual data around the crossings
                    [periCrossingFRShift, periCrossingEMGShift, crossRegionUmapShift, crossRegionClassifierShift, nanNeurCrossingsShift] = ...
                        getDataFromCrossings(threshCrossingsShift,crossingNeurIndsShift,reducIndsShift,periTransTimes,allFRs,downsampEMG,...
                        regionAssignmentsFiltered,regionWatershedLabels,classifierLabels,classifierLabelInds);

                    nanNeurCrossingsShift = find(any(isnan(periCrossingFRShift(:,:))'));
                    threshCrossingsShift(nanNeurCrossingsShift) = [];
                    crossingNeurIndsShift(nanNeurCrossingsShift) = [];
                    periCrossingFRShift(nanNeurCrossingsShift,:,:) = [];
                    periCrossingEMGShift(nanNeurCrossingsShift,:,:) = [];
                    crossRegionUmapShift(nanNeurCrossingsShift) = [];
                    crossRegionClassifierShift(nanNeurCrossingsShift) = [];

                    % save the shifts indices
                    threshCrossOrigEMGShifts{iThreshChan,iShift} = threshCrossingsShift;
                    threshCrossReducShifts{iThreshChan,iShift} = reducIndsShift;
                    threshCrossNeurShifts{iThreshChan,iShift} = crossingNeurIndsShift;

                    % only save the average to save space
                    strShiftNoBaseline = periCrossingFRShift(:,:,1:length(striatumInds)) - nanmean(periCrossingFRShift(:,1:baselineDur,1:length(striatumInds)),2);
                    aveStrNoBaseShift{iThreshChan}(iShift,:) = nanmean(nanmean(strShiftNoBaseline,1),3);
                    ctxShiftNoBaseline = periCrossingFRShift(:,:,length(striatumInds)+1:end) - nanmean(periCrossingFRShift(:,1:baselineDur,length(striatumInds)+1:end),2);
                    aveCtxNoBaseShift{iThreshChan}(iShift,:) = nanmean(nanmean(ctxShiftNoBaseline,1),3);
                    emgShiftNoBaseline = periCrossingEMGShift(:,:,musclesForAnalysis) - nanmean(periCrossingEMGShift(:,1:baselineDur,musclesForAnalysis),2);
                    aveEMGNoBaseShift{iThreshChan}(iShift,:) = nanmean(nanmean(emgShiftNoBaseline,1),3);

                    % also save behavior labels
                    crossUmapBehvsShift{iThreshChan} = crossRegionUmapShift;
                    crossClassifierBehvsShift{iThreshChan} = crossRegionClassifierShift;

                    % but save 2 control samples for testing
                    if iShift <= 2
                        allFRShift{iThreshChan,iShift}(:,:,:) = periCrossingFRShift;
                        allEMGShift{iThreshChan,iShift}(:,:,:) = periCrossingEMGShift;
                    end

                end

                % save shift amount information
                shiftAmounts{iThreshChan} = shiftAmount;

                % now divide by behavioral region
                crossingBehvLabels = {crossRegionUmap, crossRegionClassifier};
                crossingBehvLabelsShift = {crossRegionUmapShift, crossRegionClassifierShift};
                allBehvLabels = {regionAssignmentsFiltered, classifierLabels};
                crossingBehvLabelInds = {regionWatershedLabels, classifierLabelInds};
                behvTypeNames = {'Umap','Classifier'};
                for iBehvType = 1:length(crossingBehvLabels)
                    for iBehv = 1:length(crossingBehvLabelInds{iBehvType})

                        behvLabelInd = crossingBehvLabelInds{iBehvType}(iBehv);

                        % pull out the crossings in this region
                        behvCrossInds = find(crossingBehvLabels{iBehvType}==iBehv);

                        if isempty(behvCrossInds)
                            continue
                        end

                        behvCrossingStr.(behvTypeNames{iBehvType}){iBehv,iThreshChan}(:,:,:) = periCrossingFR(behvCrossInds,:,1:length(striatumInds));
                        behvCrossingCtx.(behvTypeNames{iBehvType}){iBehv,iThreshChan}(:,:,:) = periCrossingFR(behvCrossInds,:,length(striatumInds)+1:end);
                        behvCrossingEmg.(behvTypeNames{iBehvType}){iBehv,iThreshChan}(:,:,:) = periCrossingEMG(behvCrossInds,:,:);

                        %don't need to save this, can easily do this in
                        %post, analysis code

                        % do shift controls
                        behvInds = find(allBehvLabels{iBehvType} == behvLabelInd);

                        % for the controls only get segments that are all within the same
                        % behavioral region
                        behvSegsStarts = [behvInds(1) behvInds(find(diff(behvInds)>1)+1)];
                        behvSegsEnds = [behvInds(find(diff(behvInds)>1)) behvInds(end)];
                        segLengths = behvSegsEnds - behvSegsStarts;

                        goodSegsCell = {};
                        goodSegInd = 1;
                        for iSeg = 1:length(behvSegsStarts)

                            if segLengths(iSeg) < (sum(periTransTimes)+1)*2
                                continue
                            end

                            % also leave some padding at the beginning and end of each
                            % segment so the whole extracted time period is within the region
                            goodSegsCell{goodSegInd} = behvSegsStarts(iSeg)+periTransTimes(1) : behvSegsEnds(iSeg)-periTransTimes(2);
                            goodSegInd = goodSegInd+1;
                        end
                        behvGoodInds = cat(2,goodSegsCell{:});

                        if isempty(behvGoodInds)
                            continue
                        end

                        % randomly sample controls
                        for iControl = 1:nControlResamples

                            gotGoodControls = false;
                            while ~gotGoodControls

                                controlCrossings = origDownsampEMGInd(behvGoodInds(randperm(length(behvGoodInds),length(behvCrossInds))));

                                % get the thresh crossing times in neural and umap indices
                                [controlCrossings, reducIndsControl, crossingNeurIndsControl] = convertCrossInds(controlCrossings,origDownsampEMGInd,filePaths);

                                % remove any emg crossings that is outside of neural data range
                                outRangeCrossings = find(crossingNeurIndsControl>size(allFRs,2)-periTransTimes(2));
                                controlCrossings(outRangeCrossings) = [];
                                crossingNeurIndsControl(outRangeCrossings) = [];
                                reducIndsControl(outRangeCrossings) = [];

                                % sometimes if there are very few crossings, the random
                                % samples don't get good indices and there are no
                                % control crossings which causes errors, so keep
                                % resampling until there is same number of control
                                % crossings, or at least 10 valid control crossings
                                if length(controlCrossings) == length(behvCrossInds) || length(controlCrossings) > 10
                                    gotGoodControls = true;
                                end

                            end

                            % get the data from these random time points
                            [controlCrossingFR, controlCrossingEMG, crossRegionUmapControl, crossRegionClassifierControl, nanNeurCrossings] = getDataFromCrossings(...
                                controlCrossings,crossingNeurIndsControl,reducIndsControl, periTransTimes,allFRs,downsampEMG,...
                                regionAssignmentsFiltered,regionWatershedLabels,classifierLabels,classifierLabelInds);

                            reducIndsControl(nanNeurCrossings) = [];
                            controlCrossings(nanNeurCrossings) = [];
                            crossingNeurIndsControl(nanNeurCrossings) = [];
                            crossRegionUmapControl(nanNeurCrossings) = [];
                            crossRegionClassifierControl(nanNeurCrossings) = [];
                            controlCrossingFR(nanNeurCrossings,:,:) = [];
                            controlCrossingEMG(nanNeurCrossings,:,:) = [];

                            % save crossing indices
                            threshCrossOrigEMGBehvControls.(behvTypeNames{iBehvType}){iThreshChan,iBehv,iControl} = controlCrossings;
                            
                            % only save the average to save space
                            strShiftNoBaseline = controlCrossingFR(:,:,1:length(striatumInds)) - nanmean(controlCrossingFR(:,1:baselineDur,1:length(striatumInds)),2);
                            aveStrNoBaseBehvsControl.(behvTypeNames{iBehvType}){iThreshChan,iBehv}(iControl,:) = nanmean(nanmean(strShiftNoBaseline,1),3);
                            ctxShiftNoBaseline = controlCrossingFR(:,:,length(striatumInds)+1:end) - nanmean(controlCrossingFR(:,1:baselineDur,length(striatumInds)+1:end),2);
                            aveCtxNoBaseBehvsControl.(behvTypeNames{iBehvType}){iThreshChan,iBehv}(iControl,:) = nanmean(nanmean(ctxShiftNoBaseline,1),3);
                            emgShiftNoBaseline = controlCrossingEMG(:,:,musclesForAnalysis) - nanmean(controlCrossingEMG(:,1:baselineDur,musclesForAnalysis),2);
                            aveEMGNoBaseBehvsControl.(behvTypeNames{iBehvType}){iThreshChan,iBehv}(iControl,:) = nanmean(nanmean(emgShiftNoBaseline,1),3);

                            % but save 2 control samples for testing
                            if iControl <= 2
                                allFRBehvControl.(behvTypeNames{iBehvType}){iThreshChan,iBehv,iControl}(:,:,:) = controlCrossingFR;
                                allEMGBehvControl.(behvTypeNames{iBehvType}){iThreshChan,iBehv,iControl}(:,:,:) = controlCrossingEMG;
                            end


                        end % of shift controls loop

                    end % of behaviors loop

                end % of behavior lable type loop

            end % of emg channel loop

            % save data for each session
            save(filePaths.EMGTrigInitiations, 'threshCrossOrigEMG', 'threshCrossReduc', 'threshCrossNeur',...
                'crossFRs', 'crossEMGs', 'crossUmapBehvs', 'crossClassifierBehvs', ...
                'threshCrossOrigEMGShifts', 'threshCrossReducShifts', 'threshCrossNeurShifts', 'shiftAmounts', ...
                'aveStrNoBaseShift', 'aveCtxNoBaseShift', 'aveEMGNoBaseShift', 'allFRShift', 'allEMGShift', ...
                'crossUmapBehvsShift','crossClassifierBehvsShift','aveStrNoBaseBehvsControl', 'aveCtxNoBaseBehvsControl',...
                'aveEMGNoBaseBehvsControl', 'allFRBehvControl', 'allEMGBehvControl','threshCrossOrigEMGBehvControls','-v7.3');

        end % of brain region loop

    end % of session loop

end % of animal loop

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


function [threshCrossings, reducInds, crossingNeurInds] = convertCrossInds(threshCrossings,origDownsampEMGInd,filePaths)
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
load(filePaths.VideoSyncFrames)
currentDir = pwd;
cd(filePaths.processedDataFolder)
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


function [periCrossingFR, periCrossingEMG, crossRegionUmap, crossRegionClassifier, nanNeurCrossings] = getDataFromCrossings(threshCrossings,crossingNeurInds,reducInds,...
    periTransTimes,allFRs,downsampEMG,regionAssignmentsFiltered,regionWatershedLabels,classifierAssignments,classifierLabelInds)
% function to get neural, emg, and behavior class label data from each
% threshold crossing

periCrossingFR = zeros(length(threshCrossings),sum(periTransTimes)+1,size(allFRs,1));
periCrossingEMG = zeros(length(threshCrossings),sum(periTransTimes)+1,size(downsampEMG,1));

for iCross = 1:length(threshCrossings)

    %determine which UMAP region it's assigned to (based on plurality of
    %the surrounding time points)
    reducCrossInds = reducInds(iCross)-periTransTimes(1) : reducInds(iCross)+periTransTimes(2);
    regionDistributionUmap = histcounts(regionAssignmentsFiltered(reducCrossInds),[regionWatershedLabels regionWatershedLabels(end)+1]);
    [~,highestRegionUmap] = max(regionDistributionUmap);

    %also do it for supervised classifier labels
    regionDistributionClassifier = histcounts(classifierAssignments(reducCrossInds),[classifierLabelInds classifierLabelInds(end)+1]);
    [~,highestRegionClassifier] = max(regionDistributionClassifier);

    %only use regions where at least half of the time points is in that
    %region
    if regionDistributionUmap(highestRegionUmap) <= round(length(reducCrossInds)/2)
        crossRegionUmap(iCross) = 0;
    else
        crossRegionUmap(iCross) = highestRegionUmap;
    end

    if regionDistributionClassifier(highestRegionClassifier) <= round(length(reducCrossInds)/2)
        crossRegionClassifier(iCross) = 0;
    else
        crossRegionClassifier(iCross) = highestRegionClassifier;
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
