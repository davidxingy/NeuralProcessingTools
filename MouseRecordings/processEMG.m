function [processedEMG,filteredEMG,syncInds] = processEMG(emgDir)

allfiles = string(ls(emgDir));
intanFiles = allfiles(cellfun(@(x) contains(x,'.rhd'),string(allfiles)));
channelInds = 1:4;
channelNames = {'PL','ECR', 'Bi', 'Tri'};

for iFile = 1:length(intanFiles)
    
    outputData = read_Intan_RHD2000_file(emgDir,intanFiles{iFile});
    
    fileNumSamples(iFile) = size(outputData.amplifier_data,2);
    
    amplifier_data = outputData.amplifier_data(channelInds,:);
    
    %get sync pulse indices
    syncInds{iFile} = detectSyncPulse(outputData.board_adc_data,2.5)+sum(fileNumSamples)-fileNumSamples(1);
    
    %concatenate with last 5000 samples from last file (unless it's the
    %first file) to make sure the transition is smooth
    if iFile~=1
        amplifier_data = [prevDataTail amplifier_data];
    end
    
    %process emg
    for iChan = 1:length(channelInds)
        [filtsig, env] = processEMGSig(amplifier_data(iChan,:));
        chanEnv{iFile}(iChan,:) = env;
        chanFilt{iFile}(iChan,:) = filtsig;
    end
    
    %if not first file, discard the concatentated samples, but keep 2000 of
    %them since the previous file discarded the last 2000
    if iFile~=1
        chanEnv{iFile} = chanEnv{iFile}(:,3001:end);
        chanFilt{iFile} = chanFilt{iFile}(:,3001:end);
    end
    
    %don't use last 2000 samples due to edge cases
    chanEnv{iFile} = chanEnv{iFile}(:,1:end-2000);
    chanFilt{iFile} = chanFilt{iFile}(:,1:end-2000);
    
    %save last 5000 samples for use with next file
    prevDataTail = amplifier_data(:,end-4999:end);
    
end

processedEMG = [chanEnv{:}];
filteredEMG = [chanFilt{:}];
for iChan = 1:length(channelInds)
    
    processedEMG(iChan,:) = processedEMG(iChan,:)/prctile(processedEMG(iChan,:),90);
    filteredEMG(iChan,:) = filteredEMG(iChan,:)/prctile(filteredEMG(iChan,:),90);
    
end

end


function risingEdges = detectSyncPulse(signal,thresh)

risingEdges = find(signal(2:end)>=thresh & signal(1:end-1)<thresh)+1;

end

function [filtsig, envelope] = processEMGSig(signal)

% filter between 40 and 1000 hz
[b,a] = butter(3,[40/10000 1000/10000]);

filtsig = filtfilt(b,a,signal);

% rectify and do 40ms moving average filter
envelope = filtfilt(ones(800,1)/1000,1,abs(filtsig));

end