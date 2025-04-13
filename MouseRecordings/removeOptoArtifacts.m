function removeOptoArtifacts(binDIR,metaFilename,binFilename,newBinFilename,laserInds)
% removeOptoArtifacts(binDIR,metaFilename,binFilename,newBinFilename)
% function to do some processing on spikeGLX .bin data to remove laser
% stimulation artifacts. It also detects time periods where large noise
% occurs across all channels (besides those time points around laser
% stimulation) and removes them. This code assumes that the laser onset is
% saved to the sync channel of the spikeGLX recording (i.e. channel 385 is
% the laser control signal).
% 
% Inputs:
% binDIR            - Path to the directory containin the spikeGLX data
%                     files. Output data file will also be saved here
% 
% metaFilename      - Name of the spikeGLX meta data file for the recording
% 
% binFilename       - Name of the spikeGLX .bin file to do the processing
%                     on
% 
% newBinFilename    - Name of the output .bin file to save the preprocessed
%                     data as
% 
% laserInds         - If the laser trigger pulses weren't sent to the sync
%                     port of the base station, input them here.

savePlots = true;

if nargin < 5 
    laserAsSyncChan = true;
else
    laserAsSyncChan = false;
end

% load in metadata
meta = ReadMeta(metaFilename, binDIR);

% total channel count
nChans = str2double(meta.nSavedChans);

% total number of samples
totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));
% totSamps = 156287391; nChans = 384;
% load in data, in segments of 1 000 000 samples  (to avoid running out of
% memory)
segLength = 3000000;
nSegs = ceil(totSamps/segLength);

% open file for reading
readFid = fopen(fullfile(binDIR, binFilename), 'rb');
writeFid = fopen(fullfile(binDIR, newBinFilename), 'wb');

% save plots
if savePlots
    plotHandles = [];
    currentH = 1;
    %make directory for figures
    if ~exist(fullfile(binDIR,'OptoArtifactRemovalPlots'),'dir')
        mkdir(fullfile(binDIR,'OptoArtifactRemovalPlots'))
    end
end

for iSeg = 1:nSegs

    tic

    if iSeg == nSegs
        samples2Load = totSamps - (nSegs-1)*segLength;
    else
        samples2Load = segLength;
    end
    
%     fseek(readFid,2310000000,-1)

    dataArray = fread(readFid, [nChans, samples2Load], 'int16=>double');
    
    processedArray = dataArray;
        
    preSamplesToRemoveOpto = 2;
    postSamplesToRemoveOpto = 20; %remove 2/3 ms after laser

    if laserAsSyncChan
        % get opto triggers (going to asssume it's the sync chan which is the
        % last channel)
        syncSig{iSeg} = dataArray(end,:);

        % first get opto artifacts
        optoTrigThresh = 50;

        optoTimes = find(syncSig{iSeg}(2:end)>optoTrigThresh & syncSig{iSeg}(1:end-1)<=optoTrigThresh);
        optoOffTimes = find(syncSig{iSeg}(2:end)<optoTrigThresh & syncSig{iSeg}(1:end-1)>=optoTrigThresh);
        allOptos = sort([optoTimes optoOffTimes]);

    else
        thisSegInds = (1:samples2Load) + (iSeg-1)*segLength;
        allOptos = intersect(thisSegInds,laserInds) - (iSeg-1)*segLength;
    end

    optoTimesSeg = zeros(1,size(dataArray,2));
    optoTimesSeg(allOptos) = 1;
    convWindowOpto= zeros(1,postSamplesToRemoveOpto*2+1);
    convWindowOpto(postSamplesToRemoveOpto-preSamplesToRemoveOpto:end) = 1;
    optoRemoveWindows = conv(optoTimesSeg,convWindowOpto,'same') > 0;

    %135.275, 2310000000

    % remove opto artifact data and replace with data right before each
    % opto artifact window
    for iOptoPulse = 1:length(allOptos)

        %due to imperfect syncing, the real laser artifact might be offset
        %from the index, find the real time, which should still be near the
        %index, by looking for sudden jump in all channels
        dataAroundIndex = processedArray(1:end-1,max(1,allOptos(iOptoPulse)-30):min(allOptos(iOptoPulse)+30,size(processedArray,2)));
        diffSums = sum(abs(diff(dataAroundIndex,[],2)));
        [~, artInd] = max(diffSums);

        numPreSamples = min(allOptos(iOptoPulse)+1,32);
        realArtInd = allOptos(iOptoPulse) + artInd - numPreSamples;

        processedArray(1:end-1,realArtInd-preSamplesToRemoveOpto:realArtInd+postSamplesToRemoveOpto) = ...
            processedArray(1:end-1,realArtInd-preSamplesToRemoveOpto*2-postSamplesToRemoveOpto-1:realArtInd-preSamplesToRemoveOpto-1);

    end

%     numOptoRemoveSamples = length(find(optoRemoveWindows));
%     optoRemoveData = repmat(baselineData,1,ceil(numOptoRemoveSamples/baselineSize));
%     processedArray(1:end-1,find(optoRemoveWindows)) = optoRemoveData(:,1:length(find(optoRemoveWindows)));

    % noise detection parameters:
    noiseParameters.noiseStdThresh = 15;
    noiseParameters.noiseStdLowThresh = 0.1;
    noiseParameters.noiseVoltThresh = 400;
    noiseParameters.maxStdThreshChans = 50;
    noiseParameters.maxStdLowThreshChans = 50;
    noiseParameters.maxVoltThreshChans = 10;
    noiseParameters.stdWindowLength = 600;
    noiseParameters.preSamplesToRemoveNoiseStd= 300; %10 ms
    noiseParameters.postSamplesToRemoveNoiseStd= 3000; %100 ms
    noiseParameters.preSamplesToRemoveNoiseValue = 3000; %100 ms
    noiseParameters.postSamplesToRemoveNoiseValue = 30000; %1 s

    %first filter 300-5000 hz
    [b,a] = butter(5,[300/15000 5000/15000],'bandpass');
    dataArrayFilt = filtfilt(b,a,dataArray(1:end-1,:)');
    dataArray(1:end-1,:) = dataArrayFilt';
    processedArrayFilt = filtfilt(b,a,processedArray(1:end-1,:)');
    processedArray(1:end-1,:) = processedArrayFilt';
    clear dataArrayFilt processedArrayFilt

    %get data to fill in removed points with
    if iSeg == 1

        % get times with large noise or electrical artifacts
        artifactWindows = findNoiseArtifacts(dataArray(1:end-1,:),noiseParameters);
        artifactWindowsInds = find(artifactWindows);

        %make sure no artifacts
        baselineSize = max([noiseParameters.preSamplesToRemoveNoiseStd,noiseParameters.preSamplesToRemoveNoiseValue]) + ...
            max([noiseParameters.postSamplesToRemoveNoiseStd,noiseParameters.postSamplesToRemoveNoiseValue]);
        if any(artifactWindowsInds < baselineSize)
            warning('Baseline has artifact, aborting')
            return
        end

        baselineData = dataArray(1:end-1,1:baselineSize);

    end
    
    % next get times with large noise or electrical artifacts and replace
    artifactWindows = findNoiseArtifacts(processedArray(1:end-1,:),noiseParameters);
    artifactWindowsInds = find(artifactWindows);
    artifactWindowsIndsCell{iSeg} = artifactWindowsInds + segLength * (iSeg-1);

    numNoiseRemoveSamples = length(artifactWindowsInds);
    noiseRemoveData = repmat(baselineData,1,ceil(numNoiseRemoveSamples/baselineSize));
    processedArray(1:end-1,artifactWindowsInds) = noiseRemoveData(:,1:length(artifactWindowsInds));
    
    
    % save currrent data block to bin file
    writtenCount = fwrite(writeFid,processedArray,'int16');
    if writtenCount ~= numel(processedArray)
        error('Not all datapoints written to bin file! @Seg %u', iSeg);
    end

    % save plots
    plotH = figure('Units','pixels','OuterPosition',[10 10 1000 800],'Color','w','Visible','off');
    ax(1) = subplot(2,1,1);
    plot(dataArray(10:29,:)')
    ax(2) = subplot(2,1,2);
    plot(artifactWindows,'k','linewidth',2)
    hold on
    plot(optoRemoveWindows*0.8,'b','linewidth',2)
    linkaxes(ax,'x')
    saveas(plotH,fullfile(binDIR,'OptoArtifactRemovalPlots',[binFilename '_Block' num2str(iSeg) '.png']),'png')
    close(plotH);

    disp(['Finished segment ' num2str(iSeg) ' of ' num2str(nSegs) ', time: ' num2str(toc) ' s'])
end

% close file
fclose(readFid);
fclose(writeFid);

% save sync signal
syncSignal = cat(2,syncSig{:});
save(fullfile(binDIR,'syncSignal.mat'),'syncSignal');

% save artifacts
artifactInds = cat(2,artifactWindowsIndsCell{:});
save(fullfile(binDIR,'artifactInds.mat'),'artifactInds');


end



function artifactWindows = findNoiseArtifacts(dataArray,noiseParameters)

% get std dev of voltages
dataStdDev = movstd(dataArray,noiseParameters.stdWindowLength,0,2);

% find periods that meet noise critereon
threshCrosses = sum(dataArray > noiseParameters.noiseVoltThresh,1) >= noiseParameters.maxVoltThreshChans;
noiseStdCrosses = sum(dataStdDev > noiseParameters.noiseStdThresh,1) >= noiseParameters.maxStdThreshChans;
noiseStdLows = sum(dataStdDev < noiseParameters.noiseStdLowThresh,1) >= noiseParameters.maxStdLowThreshChans;

% get windows of removal
convWindowValue = zeros(1,noiseParameters.postSamplesToRemoveNoiseValue*2+1);
convWindowValue(noiseParameters.postSamplesToRemoveNoiseValue-noiseParameters.preSamplesToRemoveNoiseValue:end) = 1;
convWindowStd = zeros(1,noiseParameters.postSamplesToRemoveNoiseStd*2+1);
convWindowStd(noiseParameters.postSamplesToRemoveNoiseStd-noiseParameters.preSamplesToRemoveNoiseStd:end) = 1;
artifactWindows = conv(threshCrosses,convWindowValue,'same') > 0 | conv(noiseStdCrosses | noiseStdLows,convWindowStd,'same');

end
