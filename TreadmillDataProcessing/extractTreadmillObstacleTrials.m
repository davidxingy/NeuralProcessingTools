function [trialTimestamps, preTrialSamples, postTrialSamples] = extractTreadmillObstacleTrials(...
    Events, taskTypes, eventNames, nPreSteps, nPostSteps)
% trials = extractTreadmillObstacleTrials(continuousData, eventsData, taskType,...
%                                           eventNames, nPreSteps, nPostSteps)
%
% Function to segement out data into trials/epochs for David's
% treadmill-obstacle recordings. Each type of motor task has specific ways
% of dividing the trials, normalizing, behavioral markers, ect.

%get total number of blocks and make sure that this number is consistent
%across the inputs
assert(length(Events) == length(taskTypes) && length(Events) == length(eventNames), ...
    'The length of eventsData, taskTypes, and eventNames (i.e. the total number of blocks) must be the same!');
nBlocks = length(Events);

% make sure events has the simi events saved to it
assert(isfield(Events,'SimiEvents'), 'Events must have the simi tracked event timestamps!');

% go through each block and each type of task in that block
for iBlock = 1:nBlocks
    
    trialTimestamps{iBlock} = [];
    preTrialSamples{iBlock} = [];
    postTrialSamples{iBlock} = [];
    trialTypes{iBlock} = {};
    
    for iTask = 1:length(taskTypes{iBlock})
        
        currentTask = taskTypes{iBlock}{iTask};
        
        switch currentTask
            case 'Walk' %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %just basic walking trial divided into two segments, stance and
                %swing. Stance segment starts with hand/foot strike and ends with
                %hand/foot off. Swing is then hand/foot off to hand/foot strike
                strikes = Events(iBlock).SimiEvents.(eventNames{iBlock}.Strike);
                offs = Events(iBlock).SimiEvents.(eventNames{iBlock}.Off);
                
                %for walking, assume all cycles are basic walking trials unless
                %otherwise indicated
                nonTrialMarkers = [];
                for iMarker = 1:length(eventNames{iBlock}.NonTrialMarkers)
                    nonTrialMarkers = [nonTrialMarkers; ...
                        Events(iBlock).SimiEvents.(eventNames{iBlock}.NonTrialMarkers{iMarker})];
                end
                
                %first, the number strikes must be equal to or greater than the
                %number of offs
                if length(offs)>length(strikes)
                    error('Can''t have more limb lift offs than limb strikes!')
                end
                
                %Divide into trials:                
                for iCycle = 1:length(strikes)-1
                    
                    if any(nonTrialMarkers>strikes(iCycle) & nonTrialMarkers<strikes(iCycle+1))
                        %not a basic walking cycle
                        continue
                    end
                    
                    offInd = find(offs>strikes(iCycle) & offs<strikes(iCycle+1));
                    if isempty(offInd)
                        %not a gait cycle (no off defined)
                        continue
                    end
                    
                    trialTimestamps{iBlock}(end+1,:) = [strikes(iCycle) strikes(iCycle+1)];
                    trialTypes{iBlock}{end+1} = currentTask;
                    
                end
                
                %Peri-trial data:
                %get 500ms before and after for basic walking
                preTrialSamples{iBlock}(end+1) = 50;
                postTrialSamples{iBlock}(end+1) = 50;
                
                
            case {'Obstacle', 'Reach'} %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %Standing obstacle step, or standing reaching(without eating). For
                %both cases, there is only 1 segment, the "swing" phase. Starts with
                %hand/foot off, and ends with hand/foot strike.
                
                strikes = Events(iBlock).SimiEvents.(eventNames{iBlock}.Strike);
                offs = Events(iBlock).SimiEvents.(eventNames{iBlock}.Off);
                starts = Events(iBlock).SimiEvents.(eventNames{iBlock}.TrialStarts);
                stops = Events(iBlock).SimiEvents.(eventNames{iBlock}.TrialStops);
                
                %for these type of tasks, only get the trials that are marked as
                %obstacle or reach
                trialMarkers = Events(iBlock).SimiEvents.(eventNames{iBlock}.TrialMarkers);
                
                nonTrialMarkers = [];
                for iMarker = 1:length(eventNames{iBlock}.NonTrialMarkers)
                    nonTrialMarkers = [nonTrialMarkers; ...
                        Events(iBlock).SimiEvents.(eventNames{iBlock}.NonTrialMarkers{iMarker})];
                end
                
                %check that there are equal numbers of offs and strikes
                if length(offs) ~= length(strikes)
                    error('Different number of foot strikes to foot offs!');
                end
                
                %Divide into trials:
                for iTrial = 1:length(offs)
                    trialStart=starts(find(starts<offs(iTrial),1,'last'));
                    trialStop=stops(find(stops>offs(iTrial),1));
                    if any(nonTrialMarkers>trialStart & nonTrialMarkers<trialStop)
                        continue
                    end
                    if any(trialMarkers>trialStart & trialMarkers<trialStop)
                        trialTimestamps{iBlock}(end+1,:) = [offs(iTrial) strikes(iTrial)];
                        trialTypes{iBlock}{end+1} = currentTask;
                    end
                end
                
                %Peri-trial data:
                %get 1000ms before and 500 ms after
                preTrialSamples{iBlock}(end+1) = 100;
                postTrialSamples{iBlock}(end+1) = 50;
                
                
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
                if iscell(eventNames{iBlock}.TrialMarkers)
                    for iMarker=1:length(eventNames{iBlock}.TrialMarkers)
                        trialMarkers{iMarker} = eventsData.(eventNames{iBlock}.TrialMarkers{iMarker});
                    end
                else
                    trialMarkers{1} = eventsData.(eventNames{iBlock}.TrialMarkers);
                end
                
                %get bad trials
                nonTrialMarkers = [];
                for iMarker = 1:length(eventNames{iBlock}.NonTrialMarkers)
                    nonTrialMarkers = [nonTrialMarkers; eventsData.(eventNames{iBlock}.NonTrialMarkers{iMarker})];
                end
                
                %check that the are equal numbers of offs and strikes
                if length(offs) ~= length(strikes)
                    error('Different number of foot strikes to foot offs!');
                end

                %as well as starts and stops                
                if length(starts) ~= length(stops)
                    error('Different number of trial tracking starts to trial triacking stops!');
                end
                
                
                %find trials
                for iTrial = 1:length(offs)
                    
                    %first determine which tracking trial this belongs to
                    trackingTrialInd = find(offs(iTrial) > starts, 1, 'last');
                    
                    %if this trial has a non-trial marker, don't use it
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
                        trialTimestamps{iBlock}(end+1,:) = [offs(iTrial) strikes(iTrial)];
                        trialTypes{iBlock}{end+1} = currentTask;
                    end
                end
                
                
                %Peri-trial data:
                %get 1000ms before and 500 ms after
                preTrialSamples = 100;
                postTrialSamples = 50;
                
                
            case {'WalkingObstacle', 'WalkingReach', 'WalkingEat'} %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %Obstacle step or treat reach/eat (but not both) while walking.
                %Consists of the standard stance and swing phase as segments. Will
                %also get some number of gait cycles before and after the obstacle
                %or reach/each step.
                
                strikes = eventsData.(eventNames{iBlock}.Strike);
                offs = eventsData.(eventNames{iBlock}.Off);
                starts = eventsData.(eventNames{iBlock}.TrialStarts);
                stops = eventsData.(eventNames{iBlock}.TrialStops);
                
                %for these type of tasks, only get the trials that are marked as
                %obstacle or reach or eat
                trialMarkers = eventsData.(eventNames{iBlock}.TrialMarkers);
                
                nonTrialMarkers = [];
                for iMarker = 1:length(eventNames{iBlock}.NonTrialMarkers)
                    nonTrialMarkers = [nonTrialMarkers; eventsData{iBlock}.(eventNames.NonTrialMarkers{iMarker})];
                end
                
                %find trials
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
                        
                        %so we actually want to start not with the actual
                        %obstacle or reaching/eating step, but with some
                        %number of steps before it
                        
                        %get the off ind of the first step
                        for iAroundTrial = -1*nPreSteps{iBlock}{iTask}:nPostSteps{iBlock}{iTask}
                            
                            %get the off ind
                            offInd = find(...
                                offs>strikes(iTrial-nPreSteps{iBlock}{iTask}) & offs<strikes(iTrial+iAroundTrial+1));
                            
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
        
    end
    
    
    % Now, do the actual trial extraction
    
    dataSegs = segmentTrials(continuousData, mainEvents, preTrialSamples, postTrialSamples);
end

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
        trials(iTrial).dataGaitNormalized = timeNormalize(dataToNormalize,...
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
