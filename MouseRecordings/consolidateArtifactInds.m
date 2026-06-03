function artInds = consolidateArtifactInds(sessionName, recordingRegion)
% function to get all the different sources of artifacts in recordings
% (neural, emg, bad experimental behaviors, ect) and output the indices in
% both neural, emg, umap, and video timeseries

filePaths = getMouseDataNames(sessionName(1:4),sessionName,recordingRegion);

% get EMG and NP timeseries limits based on video
load(filePaths.VideoSyncFrames)
load(filePaths.UMAPFile,'regionAssignmentsFiltered')

artInds.limitNeurInd = frameNeuropixelSamples{1}{end}(end);
artInds.limitEmgInd = frameEMGSamples{1}{end}(end);
artInds.limitUmap = length(regionAssignmentsFiltered);
artInds.limitFrame = sum(cellfun(@length,frameNeuropixelSamples{1}));

% first, do neural artifacts - Noise from Opto, random spikes in noise, and
% lick artifacts, would should be detected from my
% removeArtifactFromBin.m/removeOptoArtifacts.m scripts and outputted as
% artifactTimestamps.mat
load(filePaths.npArtifactTimestamps)
% also in some sessions, there were artifacts from baseline drifts, causing
% flatlines, these were detected by checkNPFlatlines.m and saved to a
% flatlineDetections file
if exist(filePaths.flatlineDetections)
    load(filePaths.flatlineDetections,'flatLineInds')
else
    flatLineInds = [];
end

% combine the two
npRecordingArtifacts = unique([artifactTS flatLineInds]);
artInds = convertAndAddArtInds(npRecordingArtifacts,artInds,filePaths,'neur');

% next get EMG artifacts (based on thresholding criterion from processEMG.m)
load(filePaths.emgMetaData,'removedInds')
artInds = convertAndAddArtInds(removedInds,artInds,filePaths,'emg');

% For UMAP, there were some sessions where the EMG features clustered into
% a "artifact" region, presumably do to some strangeness in the underlying EMG
% data which were not caught by my processEMG.m file.
umapArtInds = find(isnan(regionAssignmentsFiltered));
artInds = convertAndAddArtInds(umapArtInds,artInds,filePaths,'umap');

% Don't use time periods where I manually excluded based on
% looking at the video recording and there was out of place things
% happening, e.g. the animal getting tangled up with the cable
behvLabelsVars = load(filePaths.BehaviorAnnotationLabels);
if isfield(behvLabelsVars,'excludeFramesConcat')
    frameArtInds = behvLabelsVars.excludeFramesConcat;
else
    frameArtInds = [];
end
artInds = convertAndAddArtInds(frameArtInds,artInds,filePaths,'video');

% finally, consolidate all artifact sources together to get overall
% artifact indices

artTimeNames = {'NeurInds','NeurInds1ms0msSmooth','NeurInds1ms10msSmooth','NeurInds10ms0msSmooth','NeurInds10ms30msSmooth',...
    'NeurInds100ms0msSmooth','NeurInds100ms50msSmooth','EmgInds','EmgInds1ms0msSmooth','EmgInds1ms10msSmooth',...
    'EmgInds10ms0msSmooth','EmgInds10ms30msSmooth','EmgInds100ms0msSmooth','EmgInds100ms50msSmooth','UmapInds','FramesInds'};

for iArtTime = 1:length(artTimeNames)
    neurArts = artInds.(['neurArt' artTimeNames{iArtTime}]);
    emgArts = artInds.(['emgArt' artTimeNames{iArtTime}]);
    umapArts = artInds.(['umapArt' artTimeNames{iArtTime}]);
    videoArts = artInds.(['videoArt' artTimeNames{iArtTime}]);

    artInds.(['allArt' artTimeNames{iArtTime}]) = unique([neurArts(:); emgArts(:); umapArts(:); videoArts(:)]);
end


end % of main function


% Subfunctions------------------------

function artInds = convertAndAddArtInds(sigArtTimeStamps,artInds,filePaths,artSource)

load(filePaths.VideoSyncFrames)

switch artSource
    case 'neur'

        % add in neur time
        artInds = addArtInds(sigArtTimeStamps,artInds,artSource,'Neur');
        sigArtTimeStampsInNeur = sigArtTimeStamps;

        % do conversions for others
        convertNeur = false;
        convertEMG = true;
        convertUMAP = true;
        convertFrames = true;

    case 'emg'

        % add in emg time
        artInds = addArtInds(sigArtTimeStamps,artInds,artSource,'Emg');
        sigArtTimeStampsInEmg = sigArtTimeStamps;

        % do conversions for others
        convertNeur = true;
        convertEMG = false;
        convertUMAP = true;
        convertFrames = true;

    case 'umap'
        % add in umap time
        artInds = addArtInds(sigArtTimeStamps,artInds,artSource,'Umap');

        % then convert to EMG time and run through the rest of the
        % conversions
        if isempty(sigArtTimeStamps)
            sigArtTimeStampsInEmg = [];
        else
            load(filePaths.UMAPFile,'origDownsampEMGInd')
            origDownsampArtInds = origDownsampEMGInd(sigArtTimeStamps);

            % because this is in 1ms downsamples, have to convert it back
            % to the original time points (20 samples per ms)
            boutStarts = [origDownsampArtInds(1) origDownsampArtInds(find(diff(origDownsampArtInds)>1)+1)];
            boutEnds = [origDownsampArtInds(find(diff(origDownsampArtInds)>1)) origDownsampArtInds(end)];

            allBoutInds = {};
            for iBout = 1:length(boutStarts)
                allBoutInds{iBout} = boutStarts(iBout)*20:boutEnds(iBout)*20;
            end
            sigArtTimeStampsInEmg = [allBoutInds{:}];

        end
        artInds = addArtInds(sigArtTimeStampsInEmg,artInds,artSource,'Emg');

        % do conversions for others
        convertNeur = true;
        convertEMG = false;
        convertUMAP = false;
        convertFrames = true;

    case 'video'

        % add in video frames
        artInds = addArtInds(sigArtTimeStamps,artInds,artSource,'Frames');

        % then convert to EMG time and run through the rest of the
        % conversions
        if isempty(sigArtTimeStamps)
            sigArtTimeStampsInEmg = [];
        else
            allFramesEmgs = [frameEMGSamples{1}{:}];

            % because this is in video frames, should upsample back to
            % original EMG time points (20 samples per ms)
            boutStarts = [sigArtTimeStamps(1) sigArtTimeStamps(find(diff(sigArtTimeStamps)>1)+1)];
            boutEnds = [sigArtTimeStamps(find(diff(sigArtTimeStamps)>1)) sigArtTimeStamps(end)];

            allBoutInds = {};
            for iBout = 1:length(boutStarts)
                allBoutInds{iBout} = allFramesEmgs(boutStarts(iBout)):allFramesEmgs(boutEnds(iBout));
            end
            sigArtTimeStampsInEmg = [allBoutInds{:}];

        end
        artInds = addArtInds(sigArtTimeStampsInEmg,artInds,artSource,'Emg');

        % do conversions for others
        convertNeur = true;
        convertEMG = false;
        convertUMAP = true;
        convertFrames = false;

end

if convertNeur

    if isempty(sigArtTimeStampsInEmg)
        sigArtTimeStampsInNeur = [];
    else
        % convert using the sync function
        % because the emg points are in 20kHz sampling rate, have to fill in
        % the gaps of the higher sampling rate neural time series
        boutStarts = [sigArtTimeStampsInEmg(1) sigArtTimeStampsInEmg(find(diff(sigArtTimeStampsInEmg)>1)+1)];
        boutEnds = [sigArtTimeStampsInEmg(find(diff(sigArtTimeStampsInEmg)>1)) sigArtTimeStampsInEmg(end)];

        allBoutInds = {};
        currentDir = pwd;
        cd(filePaths.processedDataFolder)
        for iBout = 1:length(boutStarts)
            boutEdgesNeur = NeurEMGSync([boutStarts(iBout) boutEnds(iBout)],...
                frameEMGSamples, frameNeuropixelSamples, 'EMG');
            allBoutInds{iBout} = boutEdgesNeur(1):boutEdgesNeur(2);
        end
        cd(currentDir)

        sigArtTimeStampsInNeur = [allBoutInds{:}];
    end
    artInds = addArtInds(sigArtTimeStampsInNeur,artInds,artSource,'Neur');

end


if convertEMG
    % convert using the sync function
    currentDir = pwd;
    cd(filePaths.processedDataFolder)
    sigArtTimeStampsInEmg = unique(NeurEMGSync(sigArtTimeStampsInNeur,...
        frameEMGSamples, frameNeuropixelSamples, 'Neuropixel'));
    cd(currentDir)

    artInds = addArtInds(sigArtTimeStampsInEmg,artInds,artSource,'Emg');
end

if convertUMAP
    load(filePaths.UMAPFile,'origDownsampEMGInd')
    % find umap points that correspond to emg artifact inds
    [~,sigArtTimeStampsInUmap] = intersect(origDownsampEMGInd,unique(round(sigArtTimeStampsInEmg/20)));
    artInds = addArtInds(sigArtTimeStampsInUmap,artInds,artSource,'Umap');
end

if convertFrames
    allFramesEmgs = [frameEMGSamples{1}{:}];

    % I want to avoid doing a for loop going through each artifact ind an
    % comparing to the frame inds, so I will implement an algorithm using
    % the intersect function.

    % First will generate all EMG inds that are within the recording
    % session (all inds between first and last frame)
    frameEmgContSamps = allFramesEmgs(1):allFramesEmgs(end);

    % next get the inds from this set which are also artifact inds (using
    % intersect)
    [~,artContInds] = intersect(frameEmgContSamps,sigArtTimeStampsInEmg);

    % then intersect these inds with the inds corresponding to the frame
    % emg inds. However, there might be some brief artifact inds that are
    % in between the emg inds of two frames, so catch these, will shift all
    % of these inds up to the spacing between frames so that any artifact
    % inds will be caught by some frame. Note that this has an unintended
    % effect of possibling adding an adjacent frame to the artifact frames
    % which shouldn't be an artifact frame, but I'm willing to accept this
    % more conservative approach for the sake of speed.
    thisShiftFrameArts = {};
    for iShift = 1:max(diff(allFramesEmgs))
        [~,thisShiftFrameArts{iShift}] = intersect(allFramesEmgs,frameEmgContSamps(artContInds)-iShift+1);
    end
    sigArtTimeStampsInFrames = unique(cat(1,thisShiftFrameArts{:}));

    artInds = addArtInds(sigArtTimeStampsInFrames,artInds,artSource,'Frames');
end

end


function artInds = addArtInds(sigArtTimeStamps,artInds,artSource,artTime)
% add artifacts in different time bins

% bin size depends on the which signal type we're using
switch artTime
    case 'Emg'
        binSize1ms = 20;
        sigEndInd = artInds.limitEmgInd;
        doLargerBins = true;
        doSmooth = false;
    case 'Neur'
        binSize1ms = 30;
        sigEndInd = artInds.limitNeurInd;
        doLargerBins = true;
        doSmooth = true;
    case 'Umap'
        binSize1ms = 1;
        sigEndInd = artInds.limitUmap;
        doLargerBins = false;
        doSmooth = false;
    case 'Frames'
        binSize1ms = 1;
        sigEndInd = artInds.limitUmap;
        doLargerBins = false;
        doSmooth = false;
end

% add for raw sample indices
artInds.([artSource 'Art' artTime 'Inds']) = sigArtTimeStamps;

% then add for 1ms bin - 0ms smoothing, 1ms bin - 10ms smoothing,
% 10ms bin - 0ms smoothing, 10ms bin - 30ms smoothing,
% 100ms bins - 0ms smoothing, 100ms bins - 50ms smoothing
if doLargerBins
    countVals = histcounts(sigArtTimeStamps,1:binSize1ms:sigEndInd);
    countVals(countVals>=1) = 1;
    artInds.([artSource 'Art' artTime 'Inds1ms0msSmooth']) = find(countVals);
    if doSmooth
        convCountVals = find(convGauss(countVals, 1, 10, 0));
    else
        convCountVals = find(countVals);
    end
    artInds.([artSource 'Art' artTime 'Inds1ms10msSmooth']) = convCountVals;

    countVals = histcounts(sigArtTimeStamps,1:(binSize1ms*10):sigEndInd);
    countVals(countVals>=1) = 1;
    artInds.([artSource 'Art' artTime 'Inds10ms0msSmooth']) = find(countVals);
    if doSmooth
        convCountVals = find(convGauss(countVals, 10, 30, 0));
    else
        convCountVals = find(countVals);
    end
    artInds.([artSource 'Art' artTime 'Inds10ms30msSmooth']) = convCountVals;

    countVals = histcounts(sigArtTimeStamps,1:(binSize1ms*100):sigEndInd);
    countVals(countVals>=1) = 1;
    artInds.([artSource 'Art' artTime 'Inds100ms0msSmooth']) = find(countVals);
    if doSmooth
        convCountVals = find(convGauss(countVals, 100, 50, 0));
    else
        convCountVals = find(countVals);
    end
    artInds.([artSource 'Art' artTime 'Inds100ms50msSmooth']) = convCountVals;

end

end



%

