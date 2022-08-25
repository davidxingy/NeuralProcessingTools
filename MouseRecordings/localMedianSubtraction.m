function localMedianSubtraction(binDIR,metaFilename,binFilename,newBinFilename,nChans2Use,skipChannels,removedChans)

% load in metadata
meta = ReadMeta(metaFilename, binDIR);

% total channel count
nChans = str2double(meta.nSavedChans);

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
    
    processedArray = dataArray;
    
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
