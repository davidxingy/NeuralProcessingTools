function  [periEventFR, controlFR, transNeurInd, transControlInds, neuralRegionNames, removedEpoch] = findLabeledTransitions(baseDir, frFile, binSize, periTransTimes)
% [periEventFR, controlFR, transNeurInd, transControlInds, neuralRegionNames, removedEpoch] = ...
%   findLabeledTransitions(baseDir, frFile, binSize, periTransTimes)
% 
% Function to get transition times into labeled behaviors based on the
% human video annotation, and extract the neural firing rates around those
% transition times
% 
% 

load(fullfile(baseDir,'ProcessedData','BehaviorAnnotations','BehaviorLabels.mat'))
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
load(fullfile(baseDir,'ProcessedData',frFile),'cortexFRs','striatumFRs')
allFRs = cat(1,cortexFRs,striatumFRs);
regionFRs = {allFRs,cortexFRs,striatumFRs};
neuralRegionNames = {'Cortex_and_Striatum','Cortex','Striatum'};

frameRate = 1/mean(diff(frameNeuropixelSamples{1}{1})/30000);

% first get control trials
transControlInds = randperm(size(allFRs,2)-periTransTimes(2)-periTransTimes(1),1000)+periTransTimes(1);

transHasNan = zeros(1,length(transControlInds));
for iTrans = 1:length(transControlInds)
    for iRegion = 1:length(regionFRs)
        for iNeuron = 1:size(regionFRs{iRegion},1)

            %get control by randomly selecting points
            controlFR{iRegion}(iTrans,:,iNeuron) = regionFRs{iRegion}(iNeuron,...
                transControlInds(iTrans)-periTransTimes(1):transControlInds(iTrans)+periTransTimes(2));

        end

        %see if there's any nans
        if iRegion==1 && any(any(isnan(controlFR{iRegion}(iTrans,:,:))))
            transHasNan(iTrans) = true;
        end
    end
end

transHasNan = logical(transHasNan);
for iRegion = 1:length(regionFRs)
    controlFR{iRegion}(transHasNan,:,:) = [];
end
transControlInds(transHasNan) = [];

% now go through each behavior and get the transition trials
for iBehv = 1:length(behaviors)-1

        [behvEpochStartInds, behvEpochEndInds] = findChunksFromInds(behaviorFramesConcat{iBehv});
        behvEpochStarts = behaviorFramesConcat{iBehv}(behvEpochStartInds);
        behvEpochEnds = behaviorFramesConcat{iBehv}(behvEpochEndInds);

        %specifically for walking, don't use those epochs which
        %is the transition from grid walking into flat walking or vice versa
        badWalkTrans = zeros(1,length(behvEpochStarts));
        if strcmpi(behaviors{iBehv},'walkflat')
            %look 1 second before the start of the epoch
            for iEpoch = 1:length(behvEpochStarts)
                if any(behvEpochStarts(iEpoch)-round(frameRate) == behaviorFramesConcat{strcmpi(behaviors,'walkgrid')})
                    badWalkTrans(iEpoch) = true;
                end
            end

        elseif strcmpi(behaviors{iBehv},'walkgrid')
            for iEpoch = 1:length(behvEpochStarts)
                if any(behvEpochStarts(iEpoch)-round(2*frameRate) == behaviorFramesConcat{strcmpi(behaviors,'walkflat')})
                    badWalkTrans(iEpoch) = true;
                end
            end

        end

        %epoch starts are the transition times
        transHasNan = zeros(1,length(behvEpochStarts));
        for iTrans = 1:length(behvEpochStarts)

            %get the video the transition is from and the corresponding frame
            for iVid = 1:length(frameNeuropixelSamples{1})
                if sum(nFrames(1:iVid+1)) < behvEpochStarts(iTrans)
                    continue
                else
                    transVid = iVid;
                    transFrame = behvEpochStarts(iTrans) - sum(nFrames(1:iVid));
                    break;
                end
            end
            %get the neural recording index
            transNeurInd{iBehv}(iTrans) = round(frameNeuropixelSamples{1}{iVid}(transFrame)/(30*binSize));
            
            %sometimes datasets have sections of time where the data wasn't
            %saved, it will have han in the frame-index mapping, in this
            %case just use nans
            if isnan(transNeurInd{iBehv}(iTrans))
                for iRegion = 1:length(regionFRs)
                    periEventFR{iBehv}{iRegion}(iTrans,:,:) = nan;
                end
                transHasNan(iTrans) = true;
                continue
            end

            for iRegion = 1:length(regionFRs)
                %get the firing rates from all the different brain regions
                for iNeuron = 1:size(regionFRs{iRegion},1)
                    periEventFR{iBehv}{iRegion}(iTrans,:,iNeuron) = regionFRs{iRegion}(iNeuron,...
                        transNeurInd{iBehv}(iTrans)-periTransTimes(1):transNeurInd{iBehv}(iTrans)+periTransTimes(2));
                end

                %check if there are any nans and if there are, don't use
                %those transition trials
                if iRegion==1 && any(any(isnan(periEventFR{iBehv}{iRegion}(iTrans,:,:))))
                    transHasNan(iTrans) = true;
                end
            end

        end

        removedEpoch{iBehv} = logical(transHasNan) | logical(badWalkTrans);
        %remove any transition trials with nans 
        for iRegion = 1:length(regionFRs)
            periEventFR{iBehv}{iRegion}(removedEpoch{iBehv},:,:) = [];
        end
        transNeurInd{iBehv}(removedEpoch{iBehv}) = [];

end


% 
