function frameSamples = syncFLIR(emgSyncInds,risingEdgeFrames,pulseDuration,nPulsesPerTrain,fs,nVidFrames)
% synchronize EMG with video
% find the total number of sync trains

emgEdges = emgSyncInds;% emgEdges = [emgSyncInds{:}];
trainEnds = [0 find(diff(emgEdges)>pulseDuration*fs*1.5) length(emgEdges)];
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
        
        %determine how many different sync trains were applied in the video
        pulseDiffs = diff(risingEdgeFrames{iVidChunk}{iVidFile})>pulseDuration*40*1.5;
        
        %only keep sync trains that have the expected number of pulses
        %(more or less might be detected due to occusion of the LED)
        %get number of pulses per train
        [segStartInds, segLengths] = findConstants(double(pulseDiffs'));
        
        wholeTrain = ones(1,length(segStartInds));
        if pulseDiffs(1) == 1
            segStartInds = [1 segStartInds];
            segLengths = [1 segLengths];
            wholeTrain = [0 wholeTrain];
        end
        badTrains = find(segLengths+1 ~= nPulsesPerTrain);
        
        segStartInds = [segStartInds length(risingEdgeFrames{iVidChunk}{iVidFile})+1];
        
        iVidTrain = {};
        removeFrameInds = {};
        for iTrain = 1:length(segStartInds)-1
            
            if any(iTrain==badTrains)
                %don't use this pulse train, remove it from the frame inds
                removeFrameInds{iTrain} = segStartInds(iTrain):segStartInds(iTrain+1)-1;
                if wholeTrain(iTrain)
                    syncTrainInd = syncTrainInd + 1;
                end
                continue
            end
            
            %also make sure that the time between pulses is the same for each
            %train (somtimes the led is blocked briefly for the first pulse)
            pulseDurations = diff(risingEdgeFrames{iVidChunk}{iVidFile}(...
                segStartInds(iTrain):segStartInds(iTrain+1)-1));
            if length(unique(pulseDurations)) ~= 1
                if length(unique(pulseDurations)) > 2 || diff(unique(pulseDurations)) ~= 1
                    %don't use this pulse train, remove it from the frame inds
                    removeFrameInds{iTrain} = segStartInds(iTrain):segStartInds(iTrain+1)-1;
                    syncTrainInd = syncTrainInd + 1;
                    continue
                end
            end
            
        
            %add to the intan sample list of pulses
            iVidTrain{iTrain} = trainPulseIndsEMG{syncTrainInd};
            syncTrainInd = syncTrainInd + 1;
            
        end
        
        risingEdgeFrames{iVidChunk}{iVidFile}([removeFrameInds{:}]) = [];
        risingEdgeSamples{iVidChunk}{iVidFile} = [iVidTrain{:}];
        
    end
end

% go through each video chunk and sync
for iVidChunk = 1:length(risingEdgeFrames)
    
    edgeFrames={};
    nonEdgeFrames={};
    for iVidFile = 1:length(risingEdgeFrames{iVidChunk})
        
        %for each video of the chunk, offset by total number of frames of
        %all the files that came before it
        if iVidFile == 1    
            nOffsetFrames = 0;
        else
            nOffsetFrames = sum(nVidFrames{iVidChunk}(1:iVidFile-1));
        end
        
        edgeFrames{iVidFile} = risingEdgeFrames{iVidChunk}{iVidFile}' + nOffsetFrames;
        nonEdgeFrames{iVidFile} = setdiff(1:nVidFrames{iVidChunk}(iVidFile),...
            risingEdgeFrames{iVidChunk}{iVidFile}) + nOffsetFrames;
        
    end
    
    %concatenate across all vids
    allEdgeFrames = [edgeFrames{:}];
    allNonEdgeFrames = [nonEdgeFrames{:}];
    
    %do interpolation for all of the frames in the chunk
    allFrameSamples = round(sort([[risingEdgeSamples{iVidChunk}{:}],...
        interp1(allEdgeFrames,...
        [risingEdgeSamples{iVidChunk}{:}],allNonEdgeFrames,'linear','extrap')]));
    
    %plot for verification
    figure
    histogram(diff(allFrameSamples))
    
    %divide chunk back up into separate files
    for iVidFile = 1:length(risingEdgeFrames{iVidChunk})
        
        fileFrames = (1:nVidFrames{iVidChunk}(iVidFile)) + sum(nVidFrames{iVidChunk}(1:iVidFile-1));
        frameSamples{iVidChunk}{iVidFile} = allFrameSamples(fileFrames);
        
    end
            
end





% 

