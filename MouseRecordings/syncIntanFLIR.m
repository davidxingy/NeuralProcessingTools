function frameSamples = syncIntanFLIR(emgSyncInds,risingEdgeFrames,vidFileDir, vidFileNames)
% synchronize EMG with video
% find the total number of sync trains

emgEdges = [emgSyncInds{:}];
trainEnds = [0 find(diff(emgEdges)>30000) length(emgEdges)];
for iSyncPulse = 1:length(trainEnds)-1
    
    trainPulseIndsEMG{iSyncPulse} = emgEdges(trainEnds(iSyncPulse)+1:trainEnds(iSyncPulse+1));
    
end

% now group according to each video file
syncTrainInd = 1;
for iVidChunk = 1:length(risingEdgeFrames)
    
    for iVidFile = 1:length(risingEdgeFrames{iVidChunk})
        
        if isempty(risingEdgeFrames{iVidChunk}{iVidFile})
            risingEdgeSamples{iVidChunk}{iVidFile}=[];
            continue
        end
        
        nSyncTrains = sum(diff(risingEdgeFrames{iVidChunk}{iVidFile})>50)+1;
        
        risingEdgeSamples{iVidChunk}{iVidFile} = [trainPulseIndsEMG{syncTrainInd:syncTrainInd+nSyncTrains-1}];
        syncTrainInd = syncTrainInd + nSyncTrains;
        
    end
end

% go through each video chunk and sync
for iVidChunk = 1:length(risingEdgeFrames)
    
    nFrames = [];
    edgeFrames={};
    nonEdgeFrames={};
    for iVidFile = 1:length(risingEdgeFrames{iVidChunk})
        
        %for each video of the chunk, offset by total number of frames of
        %all the files that came before it
        nOffsetFrames = sum(nFrames);
        v = VideoReader(fullfile(vidFileDir, vidFileNames{iVidChunk}{iVidFile}));
        nFrames(iVidFile) = round(v.FrameRate*v.Duration);
        
        edgeFrames{iVidFile} = risingEdgeFrames{iVidChunk}{iVidFile} + nOffsetFrames;
        nonEdgeFrames{iVidFile} = setdiff(1:nFrames(iVidFile),risingEdgeFrames{iVidChunk}{iVidFile}) + nOffsetFrames;
        
    end
    
    %concatenate across all vids
    allEdgeFrames = [edgeFrames{:}];
    allNonEdgeFrames = [nonEdgeFrames{:}];
    
    %do interpolation for all of the frames in the chunk
    allFrameSamples = round(sort([[risingEdgeSamples{iVidChunk}{:}],...
        interp1(allEdgeFrames,...
        [risingEdgeSamples{iVidChunk}{:}],allNonEdgeFrames,'linear','extrap')]));
    
    %divide chunk back up into separate files
    nFrames = [0 nFrames];
    for iVidFile = 1:length(risingEdgeFrames{iVidChunk})
        
        fileFrames = (1:nFrames(iVidFile+1)) + sum(nFrames(1:iVidFile));
        frameSamples{iVidChunk}{iVidFile} = allFrameSamples(fileFrames);
        
    end
            
end





% 

