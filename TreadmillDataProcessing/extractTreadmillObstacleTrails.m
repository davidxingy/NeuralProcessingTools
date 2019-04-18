function trials = extractTreadmillObstacleTrails(continuousData, eventsData, taskType, eventNames, nPreSteps, nPostSteps, standingStanceLength)
% trials = extractTreadmillObstacleTrails(continuousData, eventsData, taskType,...
%                                           eventNames, nPreSteps, nPostSteps)
% 
% Function to segement out data into trials/epochs for David's
% treadmill-obstacle recordings. Each type of motor task has specific ways
% of dividing the trials, normalizing, behavioral markers, ect.

switch taskType
    case 'Walk' %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %just basic walking trial divided into two segments, stance and
        %swing. Stance segment starts with hand/foot strike and ends with
        %hand/foot off. Swing is then hand/foot off to hand/foot strike
        strikes = eventsData.(eventNames.Strike);
        offs = eventsData.(eventNames.Off);
        
        %for walking, assume all cycles are basic walking trials unless
        %otherwise indicated
        nonTrialMarkers = [];
        for iMarker = 1:length(eventNames.NonTrialMarkers)
            nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames.NonTrialMarkers{iMarker})];
        end
        
        %first, the number strikes must be equal to or greater than the
        %number of offs
        if length(offs)>length(strikes)
            error('Can''t have more limb lift offs than limb strikes!')
        end
        
        %Divide into trials:
        iContinuousWalkingSeg = 1;
        mainEvents = [];
        continuousWalkingSeg = [];
        
        for iCycle = 1:length(strikes)-1
            
            if any(nonTrialMarkers>strikes(iCycle) & nonTrialMarkers<strikes(iCycle+1))
                %not a basic walking cycle
                iContinuousWalkingSeg = iContinuousWalkingSeg+1;
                continue
            end
            
            offInd = find(offs>strikes(iCycle) & offs<strikes(iCycle+1));
            if isempty(offInd)
                %not a gait cycle (no off defined)
                iContinuousWalkingSeg = iContinuousWalkingSeg+1;
                continue
            end
            
            mainEvents(end+1,:) = [strikes(iCycle) offs(offInd) strikes(iCycle+1)];
            continuousWalkingSeg(end+1) = iContinuousWalkingSeg;
            
        end
        
        %Event markers:
        %the data is divided by hand/foot strike, followed by hand/foot off
        trialDivisionNames = repmat({{'LimbStrike','LimbOff'}},1,size(mainEvents,1));
        
        %Normalization:
        %for walking, can normalize to gait cycle and swing. Gait cycle
        %normalization will assign [0-60) for stance, and [60-100) for swing.
        %Swing normalization will assign the whole swing phase as 0-100.
        gaitNormalizeEventInds = mainEvents(:,2)-mainEvents(:,1)+1;
        gaitNormalizeEventPercentages = repmat(60,size(mainEvents,1),1);
        usePreTrial = false;
        swingNormalizeEventInds = [mainEvents(:,2)-mainEvents(:,1)+1 mainEvents(:,3)-mainEvents(:,1)];
        swingNormalizeEventPercentages = repmat([1 100],size(mainEvents,1),1);
        
        %Peri-trial data:
        %get 500ms before and after
        preTrialSamples = 50;
        postTrialSamples = 50;
        
        
    case {'Obstacle', 'Reach'} %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %Standing obstacle step, or standing reaching(without eating). For 
        %both cases, there is only 1 segment, the "swing" phase. Starts with
        %hand/foot off, and ends with hand/foot strike.
        
        strikes = eventsData.(eventNames.Strike);
        offs = eventsData.(eventNames.Off);
        starts = eventsData.(eventNames.TrialStarts);
        stops = eventsData.(eventNames.TrialStops);
        
        %for these type of tasks, only get the trials that are marked as
        %obstacle or reach
        trialMarkers = eventsData.(eventNames.TrialMarkers);
        
        nonTrialMarkers = [];
        for iMarker = 1:length(eventNames.NonTrialMarkers)
            nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames.NonTrialMarkers{iMarker})];
        end
        
        %check that the are equal numbers of offs and strikes
        if length(offs) ~= length(strikes)
            error('Different number of foot strikes to foot offs!');
        end
        
        %Divide into trials:
        mainEvents = [];
        for iTrial = 1:length(offs)
            trialStart=starts(find(starts<offs(iTrial),1,'last'));
            trialStop=stops(find(stops>offs(iTrial),1));
            if any(nonTrialMarkers>trialStart & nonTrialMarkers<trialStop)
                continue
            end
            if any(trialMarkers>trialStart & trialMarkers<trialStop)
                mainEvents(end+1,:) = [offs(iTrial) strikes(iTrial)];
            end
        end
        
        %none of the trials should be continuous
        continuousWalkingSeg = 1:size(mainEvents,1);
        
        %Event markers:
        %the data shouldn't have any division, it just starts with
        %hand/foot off
        trialDivisionNames = repmat({{'LimbOff'}},1,size(mainEvents,1));
        
        %Normalization:
        %for Obstacle and Reaching, the stance will have to come from the
        %pre-trial data, and the amount to use will have to be user defined
        gaitNormalizeEventInds = repmat(standingStanceLength+1, size(mainEvents,1),1);
        gaitNormalizeEventPercentages = repmat(60, size(mainEvents,1),1);
        usePreTrial = true;
        swingNormalizeEventInds = zeros(size(mainEvents,1),0); %empty since will just use default 0-100
        swingNormalizeEventPercentages = zeros(size(mainEvents,1),0);
                    
        %Peri-trial data:
        %get 1000ms before and 500 ms after
        preTrialSamples = 100;
        postTrialSamples = 50;
        
        
    case 'ReachAndEat' %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %Standing reaching + eating the treat. Two segments: initial reach
        %starting with hand off to treat contact, followed by treat contact
        %to hand strike (with putting the treat in the mouth in between)
        
        strikes = eventsData.(eventNames.Strike);
        offs = eventsData.(eventNames.Off);
        treatContacts = eventsData.(eventNames.TreatContact);
        starts = eventsData.(eventNames.TrialStarts);
        stops = eventsData.(eventNames.TrialStops);
        
        %for these type of tasks, only get the trials that are marked as
        %both reach and eat
        if iscell(eventNames.TrialMarkers)
            for iMarker=1:length(eventNames.TrialMarkers)
                trialMarkers{iMarker} = eventsData.(eventNames.TrialMarkers{iMarker});
            end
        else
            trialMarkers{1} = eventsData.(eventNames.TrialMarkers);
        end
        
        %get bad trials
        nonTrialMarkers = [];
        for iMarker = 1:length(eventNames.NonTrialMarkers)
            nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames.NonTrialMarkers{iMarker})];
        end
        
        %check that the are equal numbers of offs and strikes
        if length(offs) ~= length(strikes)
            error('Different number of foot strikes to foot offs!');
        end
        
        %find trials
        mainEvents = [];
        for iTrial = 1:length(offs)
            
            %first determine which tracking trial this belongs to
            trackingTrialInd = find(offs(iTrial) > starts, 1, 'last');
            
            if any(nonTrialMarkers>starts(trackingTrialInd) & nonTrialMarkers<stops(trackingTrialInd))
                continue
            end
            
            %has to have 1 of each type of trial Markers
            markerFound = [];
            for iMarker = 1:length(trialMarkers)
                markerFound(iMarker) = any(trialMarkers{iMarker}>starts(trackingTrialInd) & ...
                    trialMarkers{iMarker}<stops(trackingTrialInd));
            end
            if all(markerFound)
                
                %get the treat contact 
                treatInd = find(treatContacts>offs(iTrial) & treatContacts<strikes(iTrial));
                
                %must have 1 treat contact event in this trial
                if length(treatInd) ~= 1
                    error('Did not find treat contact point, or more than 1 treat contact point!');
                end
                mainEvents(end+1,:) = [offs(iTrial) treatContacts(treatInd) strikes(iTrial)];
            end
        end
        
        %none of the trials should be continuous
        continuousWalkingSeg = 1:size(mainEvents,1);
        
        %Event markers:
        %the data is divided by hand off (to start), followed by
        %treat contact
        trialDivisionNames = repmat({{'LimbOff','TreatContact'}},1,size(mainEvents,1));
        
        %Normalization:
        %for reach and eat, can only normalize to the "swing" phase
        %since there is no stance phase analog, will count just the
        %reaching phase as part of the swing (eating is considered a
        %separate movement after treat contact).
        gaitNormalizeEventInds = repmat(standingStanceLength+1, size(mainEvents,1),1);
        gaitNormalizeEventPercentages = repmat(60, size(mainEvents,1),1);
        usePreTrial = true;
        swingNormalizeEventInds = mainEvents(:,2)-mainEvents(:,1); %treat contact will be end of phase (100%)
        swingNormalizeEventPercentages = repmat(100,size(mainEvents,1),1);
             
        %Peri-trial data:
        %get 1000ms before and 500 ms after
        preTrialSamples = 100;
        postTrialSamples = 50;
        
        
    case {'WalkingObstacle', 'WalkingReach', 'WalkingEat'} %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %Obstacle step or treat reach/eat (but not both) while walking.
        %Consists of the standard stance and swing phase as segments. Will
        %also get some number of gait cycles before and after the obstacle
        %or reach/each step.
        
        strikes = eventsData.(eventNames.Strike);
        offs = eventsData.(eventNames.Off);
        starts = eventsData.(eventNames.TrialStarts);
        stops = eventsData.(eventNames.TrialStops);
        
        %for these type of tasks, only get the trials that are marked as
        %obstacle or reach or eat
        trialMarkers = eventsData.(eventNames.TrialMarkers);
        
        nonTrialMarkers = [];
        for iMarker = 1:length(eventNames.NonTrialMarkers)
            nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames.NonTrialMarkers{iMarker})];
        end
        
        %find trials
        mainEvents = [];
        continuousWalkingSeg = [];
        mainEventOffset = [];
        
        iContinuousWalkingSeg = 1;
        for iTrial = 1:length(strikes)-1
            
            %first determine which tracking trial this belongs to
            trackingTrialInd = find(strikes(iTrial) > starts, 1, 'last');
            
            %don't get bad trials
            %TODO: be able to reject cycle specific markers (such as reach
            %or eat) in addition to tracking trial specific markers.
            if any(nonTrialMarkers>starts(trackingTrialInd) & nonTrialMarkers<stops(trackingTrialInd))
                continue
            end
            
            %see if this gait cycle is one of the wanted gait cycles
            if any(trialMarkers>strikes(iTrial) & trialMarkers<strikes(iTrial+1))
                
                %get the trial markers, along with the desired number of
                %previous and number of post step gait cycles
                for iAroundTrial = -1*nPreSteps:nPostSteps
                    
                    %get the off ind
                    offInd = find(...
                        offs>strikes(iTrial+iAroundTrial) & offs<strikes(iTrial+iAroundTrial+1));
                    
                    %must have 1 off event in this trial
                    if length(offInd) ~= 1
                        error('Did not find hand/foot off, or more than 1 hand/foot off!');
                    end
                    
                    mainEvents(end+1,:) = ...
                        [strikes(iTrial+iAroundTrial) offs(offInd) strikes(iTrial+iAroundTrial+1)];
                    
                    continuousWalkingSeg(end+1) = iContinuousWalkingSeg;
                    mainEventOffset(end+1) = iAroundTrial;
                end
                
                iContinuousWalkingSeg = iContinuousWalkingSeg+1;
                
            end
        end
        
        %Event Markers:
        %the data is just divided by hand/foot strike, followed by
        %hand/foot off
        trialDivisionNames = repmat({{'LimbStrike','LimbOff'}},1,size(mainEvents,1));

        %Normalization:
        %same as walking, assign [0-60) for stance, and [60-100) for swing.
        %Swing normalization will assign the whole swing phase as 0-100.
        gaitNormalizeEventInds = mainEvents(:,2)-mainEvents(:,1)+1;
        gaitNormalizeEventPercentages = repmat(60,size(mainEvents,1),1);
        usePreTrial = false;
        swingNormalizeEventInds = [mainEvents(:,2)-mainEvents(:,1)+1 mainEvents(:,3)-mainEvents(:,1)];
        swingNormalizeEventPercentages = repmat([1 100],size(mainEvents,1),1);
        
        %Peri-trial data:
        %get 500ms before and 500 ms after
        preTrialSamples = 50;
        postTrialSamples = 50;
        
    case {'WalkingReachAndEat' 'WalkingReachAndEatLeg'} %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %The special case of reaching+eating in the same gait cycle while
        %walking on the treadmill. Will have 3 segments, from hand off to
        %treat contact, followed by treat contact to hand strike (with
        %putting the treat in the mouth in the middle). Will also get the
        %surrounding gait cycles if desired.
        
        strikes = eventsData.(eventNames.Strike);
        offs = eventsData.(eventNames.Off);
        treatContacts = eventsData.(eventNames.TreatContact);
        starts = eventsData.(eventNames.TrialStarts);
        stops = eventsData.(eventNames.TrialStops);
        
        limbTouchesTreat=true;
        if strcmpi(taskType,'WalkingReachAndEatLeg')
            limbTouchesTreat=false;
        end
        
        %for these type of tasks, only get the trials that are marked as
        %both reach and eat
        if iscell(eventNames.TrialMarkers)
            for iMarker=1:length(eventNames.TrialMarkers)
                trialMarkers{iMarker} = eventsData.(eventNames.TrialMarkers{iMarker});
            end
        else
            trialMarkers{1} = eventsData.(eventNames.TrialMarkers);
        end
        
        nonTrialMarkers = [];
        for iMarker = 1:length(eventNames.NonTrialMarkers)
            nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames.NonTrialMarkers{iMarker})];
        end
        
        %find trials
        mainEvents = [];
        continuousWalkingSeg = [];
        mainEventOffset = [];
        trialDivisionNames = {};
        
        iContinuousWalkingSeg = 1;
        for iTrial = 1:length(strikes)-1
            
            %first determine which tracking trial this belongs to
            trackingTrialInd = find(strikes(iTrial) > starts, 1, 'last');
            
            %don't get bad trials
            if any(nonTrialMarkers>starts(trackingTrialInd) & nonTrialMarkers<stops(trackingTrialInd))
                continue
            end
            
            %see if this gait cycle is one of the wanted gait cycles
            %has to have 1 of each type of trial Markers
            markerFound = [];
            for iMarker = 1:length(trialMarkers)
                markerFound(iMarker) = any(...
                    trialMarkers{iMarker}>strikes(iTrial) & trialMarkers{iMarker}<strikes(iTrial+1));
            end
            
            if all(markerFound)
                
                %get the trial markers, along with the desired number of
                %previous and number of post step gait cycles
                for iAroundTrial = -1*nPreSteps:nPostSteps
                    
                    %get the off ind
                    offInd = find(...
                        offs>strikes(iTrial+iAroundTrial) & offs<strikes(iTrial+iAroundTrial+1));
                    
                    %must have 1 off event in this trial
                    if length(offInd) ~= 1
                        error('Did not find hand/foot off, or more than 1 hand/foot off!');
                    end
                    
                    %get the treat contact ind if it's not the previous or
                    %post gait cycles
                    if iAroundTrial==0 && limbTouchesTreat
                        treatInd = find(treatContacts>strikes(iTrial+iAroundTrial) & ...
                            treatContacts<strikes(iTrial+iAroundTrial+1));
                        
                        if length(treatInd) ~= 1
                            error('Did not find treat contact event, or more than 1 treat contact!');
                        end
                        
                        %get trial
                        mainEvents(end+1,:) = [strikes(iTrial+iAroundTrial) ...
                            offs(offInd) treatContacts(treatInd) strikes(iTrial+iAroundTrial+1)];
                        
                        %Event Markers:
                        %since this is treat picking trial, the events are
                        %hand strike, then hand off, then treat contact
                        trialDivisionNames{end+1}{1}='LimbStrike';
                        trialDivisionNames{end}{2}='LimbOff';
                        trialDivisionNames{end}{3}='TreatContact';                      
                        
                    else
                        %no treat contact for the pre and post gait cycles
                        mainEvents(end+1,:) = [strikes(iTrial+iAroundTrial) offs(offInd) ...
                            strikes(iTrial+iAroundTrial+1) strikes(iTrial+iAroundTrial+1)];
                        
                        %Event Markers:
                        %since this should be just normal walking, then
                        %events are just hand strike followed by hand off
                        trialDivisionNames{end+1}{1}='HandStrike';
                        trialDivisionNames{end}{2}='HandOff';
                        
                    end
                    
                    continuousWalkingSeg(end+1) = iContinuousWalkingSeg;
                    mainEventOffset(end+1) = iAroundTrial;
                end
                
                iContinuousWalkingSeg = iContinuousWalkingSeg+1;
                
            end
        end
        
        %Normalization:
        %same as walking, assign [0-60) for stance, and [60-100) for swing.
        %will not the eating as part of whole gait cycle.
        %Swing normalization will assign the whole swing phase as 0-100.
        %Will not include eating (so just up to treat contact).
        %For normal walking trials before/after just treat as walk trials
        gaitNormalizeEventInds = [mainEvents(:,2)-mainEvents(:,1)+1 mainEvents(:,3)-mainEvents(:,1)];
        gaitNormalizeEventPercentages = repmat([60 100],size(mainEvents,1),1);
        usePreTrial = false;
        swingNormalizeEventInds = [mainEvents(:,2)-mainEvents(:,1)+1 mainEvents(:,3)-mainEvents(:,1)];
        swingNormalizeEventPercentages = repmat([1 100],size(mainEvents,1),1);
        
        %Peri-trial data:
        %get 500ms before and 500 ms after
        preTrialSamples = 50;
        postTrialSamples = 50;
        
        
    case 'Baseline' %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'Cutaneous Touch' %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
end


% Now, do the actual trial extraction

dataSegs = segmentTrials(continuousData, mainEvents, preTrialSamples, postTrialSamples);

% save into the trials struct
for iTrial = 1:size(dataSegs,1)
    
    %parse pre and post data segments, and combine all the middle data
    %segments
    trials(iTrial).data = [dataSegs{iTrial,2:end-1}];
    trials(iTrial).preData = dataSegs{iTrial,1};
    trials(iTrial).postData = dataSegs{iTrial,end};
    
    %normalize data if desired
    if isempty(gaitNormalizeEventInds) | ~isnan(gaitNormalizeEventInds)
        if usePreTrial
            dataToNormalize = [trials(iTrial).preData(:,end-(standingStanceLength-1):end) trials(iTrial).data];
        else
            dataToNormalize = trials(iTrial).data;
        end
        trials(iTrial).dataGaitNormalized = timeNormalize(trials(iTrial).data,...
            gaitNormalizeEventInds(iTrial,:), gaitNormalizeEventPercentages(iTrial,:), 0:100);
    else
        trials(iTrial).dataGaitNormalized = [];
    end
    if isempty(swingNormalizeEventInds) |  ~isnan(swingNormalizeEventInds)
        trials(iTrial).dataSwingNormalized = timeNormalize(trials(iTrial).data,...
            swingNormalizeEventInds(iTrial,:), swingNormalizeEventPercentages(iTrial,:), 0:100);
    else
        trials(iTrial).dataSwingNormalized = [];
    end
    
    %save the division indices of the middle data sigments
    divisions=unique(cumsum(cellfun(@(x) size(x,2), dataSegs(iTrial,2:end-1)))+1);
    trials(iTrial).trialDivisions = [1 divisions(1:end-1)];
    
    %and save the name of what each division is
    trials(iTrial).trialDivisionNames=trialDivisionNames{iTrial};
    
    %save task type
    if exist('mainEventOffset') && mainEventOffset(iTrial) ~= 0
        trials(iTrial).task = [taskType ' ' num2str(mainEventOffset(iTrial))];
    else
        trials(iTrial).task = taskType;
    end
    
    %save continuous trials number
    trials(iTrial).continuousSegment = continuousWalkingSeg(iTrial);
    
end


% 
