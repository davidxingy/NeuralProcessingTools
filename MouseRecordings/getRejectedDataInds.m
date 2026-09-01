function [rejectedInds rejectedBlockVids, rejectedBlockFrames] = getRejectedDataInds(frameNeuropixelSamples)
% function to get the indices in the neuropixel recordings which are
% removed (e.g. D024 which has some periods of dropped packets) and are
% identified by nans in the frame-neuropixel sync mappings
% 
% David updated 9/1/26 - also have it output the frames so outside the
% function I can use frameEMGSamples to get the timepoints in EMG time as
% using NeurEMGSync does not work

nFrames = cumsum(cellfun(@length, frameNeuropixelSamples{1}));
% concatenate across videos first
allFrameSyncInds = cat(2,frameNeuropixelSamples{1}{:});

nanFrames = find(isnan(allFrameSyncInds));

if isempty(nanFrames)
    rejectedInds = [];
    rejectedBlockVids = {};
    rejectedBlockFrames = {};
    return
end
    
nanBlocksStartFrames = nanFrames([1 find(diff(nanFrames)>1)+1]);
nanBlocksStopFrames = nanFrames([find(diff(nanFrames)>1) length(nanFrames)]);

for iNanBlock = 1:length(nanBlocksStartFrames)
    
    thisStartFrame = nanBlocksStartFrames(iNanBlock);
    thisStartVid = find(thisStartFrame < nFrames,1);
    thisStartVidFrames = thisStartFrame - nFrames(thisStartVid-1) - 1;
    
    thisStopFrame = nanBlocksStopFrames(iNanBlock);
    thisStopVid = find(thisStopFrame < nFrames,1);
    thisStopVidFrames = thisStopFrame - nFrames(thisStopVid-1) + 1;
    
    blockInds{iNanBlock} = frameNeuropixelSamples{1}{thisStartVid}(thisStartVidFrames) : ...
        frameNeuropixelSamples{1}{thisStopVid}(thisStopVidFrames);

    rejectedBlockVids{iNanBlock} = [thisStartVid thisStopVid];
    rejectedBlockFrames{iNanBlock} = [thisStartVidFrames thisStopVidFrames];
    
end

rejectedInds = cat(2,blockInds{:});


% 
