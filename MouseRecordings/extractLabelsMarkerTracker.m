function extractLabelsMarkerTracker(sessionName,baseAnnotationFileName)

filePaths = getMouseDataNames(sessionName(1:4),sessionName,'CFA');

trackedVidFiles = string(ls(fullfile(filePaths.processedDataFolder,'BehaviorAnnotations')));
trackedVidFiles = trackedVidFiles(contains(trackedVidFiles,baseAnnotationFileName));

endOfExclusions = false; % just a variable to keep track of whether the last exclude has been marked

for iVid = 1:length(trackedVidFiles)
    
    vidNum = str2num(trackedVidFiles{iVid}(length(baseAnnotationFileName)+1:length(baseAnnotationFileName)+2));
    load(fullfile(filePaths.processedDataFolder,'BehaviorAnnotations', trackedVidFiles{iVid}));
        
    if iVid == 1
        %behaviors are labled as '*-start' and '*-stop'
        behaviorsAll = squeeze(split(string({MARKERTRACKERGUI_UserData.markersInfo(2:end).name}),'-'));
        epochMarkerNamesInds = contains(behaviorsAll(:,2),'start') | contains(behaviorsAll(:,2),'stop');
        behaviors = cellstr(unique(behaviorsAll(epochMarkerNamesInds,1)));
        
        %background is the catchall for all unlabeled frames
        behaviors{end+1} = 'background';
        
        for iBehv = 1:length(behaviors)-1
            behavMarkerInds(iBehv,1) = find(contains(string({MARKERTRACKERGUI_UserData.markersInfo.name}), [behaviors{iBehv} '-start']));
            behavMarkerInds(iBehv,2) = find(contains(string({MARKERTRACKERGUI_UserData.markersInfo.name}), [behaviors{iBehv} '-stop']));
        end

    end % of first video condition
    
    %get total frames in each video (for offseting)
    nFrames(1) = 0;
    nFrames(iVid+1) = size(MARKERTRACKERGUI_UserData.trackedData,2);
    
    for iBehv = 1:length(behaviors)-1
        
        if iVid == 1
            behaviorFramesConcat{iBehv} = [];
        end
        behaviorFrames{iBehv}{iVid} = [];
        
        %get all bouts
        boutStarts = find(~isnan(MARKERTRACKERGUI_UserData.trackedData(behavMarkerInds(iBehv,1),:,1)));
        boutStops = find(~isnan(MARKERTRACKERGUI_UserData.trackedData(behavMarkerInds(iBehv,2),:,1)));

        if endOfExclusions && strcmpi(behaviors{iBehv},'exclude') && ~isempty(boutStarts)
            % There is something wrong, since the last exclusion
            % has already been markedso there should be no more
            % exclusions
            error('More excludes after final exlucde marked!')
        end

        if strcmpi(behaviors{iBehv},'exclude') && length(boutStarts)-1 == length(boutStops)
            % for exclude, it's ok to have an extra bout start, because
            % that indicates the end of the session (and I move on to
            % treadmill walking)
            endOfExclusions = true;
            lastExclusionMarker = boutStarts(end);
            lastExclusionVid = iVid;
            boutStarts = boutStarts(1:end-1);
        end

        %need same number of bout starts and stops
        if length(boutStarts) ~= length(boutStops)
            error(['Number of starts does not equal number of stops for behavior ' behaviors{iBehv}])
        elseif any(boutStarts > boutStops)
            error(['Start marker does not precede stop marker for behavior ' behaviors{iBehv}])
        end

        behvFrames = {};
        for iBout = 1:length(boutStarts)
            behvFrames{iBout} = boutStarts(iBout):boutStops(iBout);
        end

        % if final exclude marked, all frames after that are considered
        % excluded
        if endOfExclusions && strcmpi(behaviors{iBehv},'exclude')
            if lastExclusionVid == iVid
                behvFrames{end+1} = lastExclusionMarker:nFrames(iVid+1);
            else
                behvFrames{end+1} = 1:nFrames(iVid+1);
            end
        end
        
        behaviorFramesConcat{iBehv} = [behaviorFramesConcat{iBehv} [behvFrames{:}]+sum(nFrames(1:iVid))];
        behaviorFrames{iBehv}{iVid} = [behvFrames{:}];
            

    end % of behaviors loop
    
end % of video files loop

% pull out exclude as it's own thing (not an actual behavior)
excludeBehvInd = find(strcmp(behaviors,'exclude'));
excludeFrames = behaviorFrames{excludeBehvInd};
excludeFramesConcat = behaviorFramesConcat{excludeBehvInd};

behaviorFrames(excludeBehvInd) = [];
behaviorFramesConcat(excludeBehvInd) = [];
behaviors(excludeBehvInd) = [];

% next, get exclusions time points in EMG indices and neural indices
load(filePaths.VideoSyncFrames)
excludeEMGInds = [];
excludeNPInds = [];
for iVid = 1:length(excludeFrames)

    frames = excludeFrames{iVid};
    if isempty(frames)
        disp(['No Exclusions found in video ' num2str(iVid)])
        continue
    end

    excludeStarts = [frames(1) frames(find(diff(frames)>1)+1)];
    excludeStops = [frames(find(diff(frames)>1)) frames(end)];

    for iBout = 1:length(excludeStarts)
        excludeEMGInds = [excludeEMGInds frameEMGSamples{1}{iVid}(excludeStarts(iBout)):frameEMGSamples{1}{iVid}(excludeStops(iBout))];
        excludeNPInds = [excludeNPInds frameNeuropixelSamples{1}{iVid}(excludeStarts(iBout)):frameNeuropixelSamples{1}{iVid}(excludeStops(iBout))];
    end

end

save(fullfile(filePaths.processedDataFolder,'BehaviorAnnotations','BehaviorLabels'),...
    'behaviorFrames','behaviorFramesConcat','excludeFrames','excludeFramesConcat','nFrames','behaviors','excludeEMGInds','excludeNPInds')



% 
