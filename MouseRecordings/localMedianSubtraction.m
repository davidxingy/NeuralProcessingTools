function localMedianSubtraction(binDIR,metaFilename,binFilename,newBinFilename,nChans2Use,skipChannels,removedChans,chanReordering)
% localMedianSubtraction(binDIR,metaFilename,binFilename,newBinFilename,nChans2Use,skipChannels,removedChans,chanReordering)
% function to do some pre-processing on spikeGLX .bin data to do local
% median subtraction for all the channels (takes the median of the signals
% from the surrounding channels and subtracts it from the signal of that
% channel). This removes noise which may be present across multiple
% channels but with changing amplitudes across large number of channels.
% 
% Inputs:
% binDIR            - Path to the directory containin the spikeGLX data
%                     files. Output data file will also be saved here
% 
% metaFilename      - Name of the spikeGLX meta data file for the recording
% 
% binFilename       - Name of the spikeGLX .bin file to do the median
%                     subtraction on
% 
% newBinFilename    - Name of the output .bin file to save the preprocessed
%                     data as
% 
% nChans2Use        - Number of surrounding channels to use for the median
%                     calculation
% 
% skipChannels      - Any channels that you don't want to do the median
%                     subtraction for
% 
% removedChans      - This function also lets you "remove" channels by
%                     setting all it's values to zero. Specify which
%                     channels you want to do this for with this input.
%                     Just input empty vector of don't want any removed
% 
% chanReordering    - This function also lets you reorder channels (useful
%                     for when you change the spikeGLX imro mappping to
%                     shit the recording channels up or down and need to
%                     realign the channels by depth). Just input an empty
%                     vector of don't want any reordering.

% load in metadata
meta = ReadMeta(metaFilename, binDIR);

% total channel count
nChans = str2double(meta.nSavedChans);

% if chanReordering is empty, then there is no reordering (no shifting in
% the spikeGLX imro map)
if isempty(chanReordering)
    chanReordering = 1:nChans;
end


% total number of samples
totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));

% load in data, in segments of 1 000 000 samples  (to avoid running out of
% memory)
segLength = 3000000;
nSegs = ceil(totSamps/segLength);

%open file for reading
readFid = fopen(fullfile(binDIR, binFilename), 'rb');
writeFid = fopen(fullfile(binDIR, newBinFilename), 'wb');

for iSeg = 1:nSegs
    
    if iSeg == nSegs
        samples2Load = totSamps - (nSegs-1)*segLength;
    else
        samples2Load = segLength;
    end
    
    dataArray = fread(readFid, [nChans, samples2Load], 'int16=>double');
    
    processedArray = dataArray(chanReordering,:);
    
    tic
    
    %get local median
    for iChan = 1:size(dataArray,1)
        
        %don't replace channels that we want to skip
        if any(iChan==skipChannels)
            continue
        end
        
        localChanInds = max(1,iChan-round(nChans2Use/2)) : .../;'
            min(iChan+round(nChans2Use/2),size(dataArray,1));
        
        %don't use any channels that are being skipped in the calculation
        localChanInds = setdiff(localChanInds, skipChannels);
        
        %subtract the local median
        localMedian = median(dataArray(localChanInds,:),1);
        processedArray(iChan,:) = dataArray(iChan,:) - localMedian;
        
        %for channels we want to remove, just set to 0
        if any(iChan == removedChans)
            processedArray(iChan,:) = 0;
        end
        
        
    end
        
    meadSubTime = toc;
    
    tic
    
    %save currrent data block to bin file
    writtenCount = fwrite(writeFid,processedArray,'int16');
    if writtenCount ~= numel(processedArray)
        error('Not all datapoints written to bin file! @Seg %u', iSeg);
    end
    
    saveTime = toc;
    
end

% close file
fclose(readFid);
fclose(writeFid);


end


% 
