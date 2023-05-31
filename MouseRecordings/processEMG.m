function [processedEMG,filteredEMG,syncInds] = processEMG(emgDir, baseFilename, varargin)

allfiles = string(ls(emgDir));
intanFiles = allfiles(cellfun(@(x) contains(x,'.rhd'),string(allfiles)));
channelInds = [1:16];
channelNames = {'Right Biceps','Right Triceps', 'Right ECR', 'Right PL', 'Right TA', 'Right Gastr', 'Left Triceps','Left Biceps'};

if nargin == 3
    fileBreakSize = varargin{1};
else
    fileBreakSize = 25000000;
end

iSeg = 1;
iSaveFile = 1;
for iFile = 1:length(intanFiles)
    
    outputData = read_Intan_RHD2000_file(emgDir,intanFiles{iFile});
    
    fileNumSamples(iFile) = size(outputData.amplifier_data,2);
    
    emgSig = outputData.amplifier_data(channelInds,:);
    outputData.board_adc_data = zeros(4,size(outputData.amplifier_data,2));
    %get sync pulse indices
    syncInds{iFile} = detectSyncPulse(outputData.board_adc_data(1,:),2.5)+sum(fileNumSamples)-fileNumSamples(end);
    
    %get lick activations (so to remove voltage artifacts from the touch
    %sensor)
    touchInds = detectTouchInds(outputData.board_adc_data(2:end,:),2.5);
    
    %remove data from lick activations and replace with line
    if ~isempty(touchInds)
        sigNoTouchArt = replaceLinear(emgSig, touchInds);
    else
        sigNoTouchArt = emgSig;
    end
        
    %even with the lick sensor, sometimes the animal appears to touch the
    %port without setting off the detector, so need some way to detect
    %artifacts based on the EMG signal
    artInds = detectLickArtifacts(sigNoTouchArt,[],2,10000);
    
    %remove data from lick activations and replace with line
    if ~isempty(artInds)
        sigNoEMGArt = replaceLinear(sigNoTouchArt, artInds);
    else
        sigNoEMGArt = sigNoTouchArt;
    end
    
    %save all the removed points (and offset by proper amount)
    removedInds{iFile} = [touchInds artInds]+sum(fileNumSamples)-fileNumSamples(1);
    
    %concatenate with last 5000 samples from last file (unless it's the
    %first file) to make sure the transition is smooth
    if iFile~=1
        emgSigCat = [prevDataTail sigNoEMGArt];
    else
        emgSigCat = sigNoEMGArt;
    end
    
    %process emg
    for iChan = 1:length(channelInds)
        [filtsig, env] = processEMGSig(emgSigCat(iChan,:));
        [b,a] = butter(3,[100/10000 3000/10000]);
        rawSig = filtfilt(b,a,emgSigCat(iChan,:));
        chanEnv{iSeg}(iChan,:) = env;
        chanFilt{iSeg}(iChan,:) = filtsig;
        chanFilt{iSeg}(iChan,:) = rawSig;
    end
    
    %if not first file, discard the concatentated samples, but keep 2000 of
    %them since the previous file discarded the last 2000
    if iFile~=1
        chanEnv{iSeg} = chanEnv{iSeg}(:,3001:end);
        chanFilt{iSeg} = chanFilt{iSeg}(:,3001:end);
    end
    
    %don't use last 2000 samples due to edge cases
    chanEnv{iSeg} = chanEnv{iSeg}(:,1:end-2000);
    chanFilt{iSeg} = chanFilt{iSeg}(:,1:end-2000);
    
    %save last 5000 samples for use with next file
    prevDataTail = emgSig(:,end-4999:end);
    
    %Save whenever we reach the sample limit
    totalSamples = sum(cellfun(@(x) size(x,2),chanEnv));
    
    if totalSamples >= fileBreakSize || iFile == length(intanFiles)
        
        processedEMG = [chanEnv{:}];    
        filteredEMG = [chanFilt{:}];
        
        chanEnv = {};
        chanFilt = {};
        
        if totalSamples > fileBreakSize
            
            chanEnv{1} = processedEMG(:,fileBreakSize+1:end);
            chanFilt{1} = filteredEMG(:,fileBreakSize+1:end);
            
            processedEMG = processedEMG(:,1:fileBreakSize);
            filteredEMG = filteredEMG(:,1:fileBreakSize);
            
            iSeg = 2;
        else
            iSeg = 1;
        end
        
        save([baseFilename '_ProcessedEMG_Block2_Part' num2str(iSaveFile)], 'processedEMG', 'filteredEMG')
        iSaveFile = iSaveFile+1;
        
        if iFile == length(intanFiles)
            removedInds = [removedInds{:}];
            syncInds = [syncInds{:}];
            save([baseFilename '_MetaData'],'removedInds','syncInds','fileNumSamples','channelNames');
        end
        
    else
        iSeg = iSeg + 1;
    end
    
end


for iChan = 1:length(channelInds)
    
    processedEMG(iChan,:) = processedEMG(iChan,:)/prctile(processedEMG(iChan,:),90);
    filteredEMG(iChan,:) = filteredEMG(iChan,:)/prctile(filteredEMG(iChan,:),90);
    
end

end



function risingEdges = detectSyncPulse(signal,thresh)

risingEdges = find(signal(2:end)>=thresh & signal(1:end-1)<thresh)+1;

end



function lickInds = detectTouchInds(signals, voltThresh)

% get the samples where the lick is detected
lickInds = sort(find(any(signals>voltThresh)));

if isempty(lickInds)
    return
end

% because it is noisy sometimes there are small flucuations. If any of the
% gaps between detection samples are smaller than some specified
% threshold, just fill in the gap (just call all those samples lick
% detected).
gapInds = find(diff(lickInds)~=1);

filledGapTimes = {};
for iGap = 1:length(gapInds)

    %add 1.5s to the end of the gap (due to the slow flucations)
    filledGapTimes{end+1} = lickInds(gapInds(iGap)):lickInds(gapInds(iGap))+30000;
    
    %also add 10ms before the start of the next gap due to the smalle delay
    %between the touch detection and artifact onset
    filledGapTimes{end+1} = lickInds(gapInds(iGap)+1)-2000:lickInds(gapInds(iGap)+1);
    
end

% also do it for the very first detection sample and last sample
filledGapTimes{end+1} = lickInds(1)-2000:lickInds(1);
filledGapTimes{end+1} = lickInds(end):lickInds(end)+30000;

lickInds = sort(unique([lickInds filledGapTimes{:}]));
% make sure we don't go under 1 or after the last index
lickInds(lickInds<1 | lickInds>size(signals,2)) = [];

end



function artInds = detectLickArtifacts(signals, thresh, minChans, minDuration)

% criteron will be that at least some number of channels maintains a 
% large voltage level constantly (i.e not oscilating) for a period of time
if isempty(thresh)
    %if no threshold given, then just use the 90th percentile
    threshMat = repmat(prctile(signals',90)', 1, size(signals,2));
else
    threshMat = repmat(thresh,size(signals,1),size(signals,2));
end
  
threshCrossings = abs(signals) > threshMat;

multiThreshCrossings = double(sum(threshCrossings,1) >= minChans);

[segStartInds, segLengths] = findConstants(multiThreshCrossings);

artSegs = segLengths >= minDuration & multiThreshCrossings(segStartInds+1)==1;
artStartInds = segStartInds(artSegs);
artLengths = segLengths(artSegs);

artSegInds = {};
for iSeg = 1:length(artStartInds)
    
    artSegInds{iSeg} = artStartInds(iSeg)-2000:artStartInds(iSeg)+artLengths(iSeg)+2000;
    
end

if isempty(artSegInds)
    artInds = [];
else
    artInds = [artSegInds{:}];
    artInds(artInds < 1 | artInds > size(signals,2)) = [];
end

end



function signalNoArt = replaceLinear(signal, lickInds)

% find the start and end of continuous blocks to remove
gaps = find(diff(lickInds)>1);
blockStarts = [lickInds(1) lickInds(gaps+1)];
blockEnds = [lickInds(gaps) lickInds(end)];

signalNoArt = signal;
for iBlock = 1:length(blockStarts)
    
    for iChan = 1:size(signal,1)
        lineReplacement = linspace(signal(iChan,blockStarts(iBlock)),signal(iChan,blockEnds(iBlock)), ...
            blockEnds(iBlock)-blockStarts(iBlock)+1);
        signalNoArt(iChan,blockStarts(iBlock):blockEnds(iBlock)) = lineReplacement;
    end
    
end


end



function [filtsig, envelope] = processEMGSig(signal)

% filter between 40 and 1000 hz
[b,a] = butter(3,[40/10000 1000/10000]);

filtsig = filtfilt(b,a,signal);

% rectify and do 40ms moving average filter
envelope = filtfilt(ones(200,1)/200,1,abs(filtsig));

end


% 
