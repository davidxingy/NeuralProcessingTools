function neuronDataStruct = processNeurons(baseDIR,phyModified,binFile,startSample)

waveforms = readNPY(fullfile(baseDIR, 'templates.npy'));
electrodPos = readNPY(fullfile(baseDIR, 'channel_positions.npy'));

if phyModified
    
    % load in phy data
    phyClusters = readtable(fullfile(baseDIR,'cluster_info.tsv'), "FileType","text",'Delimiter', '\t');
    spikeClusterLabels = readNPY(fullfile(baseDIR,'spike_clusters.npy'));
    spikeTimes = readNPY(fullfile(baseDIR,'spike_times.npy'));
    spikeTimes = spikeTimes + startSample;

    % set up loading in raw bin file
    [binDIR, binFilename, binExt] = fileparts(binFile);
    meta = ReadMeta([binFilename binExt], binDIR);
    nChans = str2double(meta.nSavedChans);
    totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));
    nPreSamples = 40;
    nPostSamples =  60;
    sampleRate = 30000;

    readFid = fopen(fullfile(binDIR, [binFilename binExt]), 'rb');

    goodPhyClusters = find(contains(phyClusters.group,'good'));
    clusterIds = phyClusters.cluster_id(goodPhyClusters);

    for iUnit = 1:length(clusterIds)

        clusterSpikes = find(spikeClusterLabels == clusterIds(iUnit));
        neuronDataStruct(iUnit).timeStamps = sort(spikeTimes(clusterSpikes));
        %Remove any time stamps which are within 1ms of the previous time
        %stamp, which I sometimes do during Phy when kilosort does this
        %weird thing where it counts a spike multiple times with a small
        %offset
        neuronDataStruct(iUnit).timeStamps(find(diff(neuronDataStruct(iUnit).timeStamps) < 30)+1) = [];

        % now get the  waveforms (rather than using the templates,
        % use the averge waveforms from the data itself)

        % to save time, only do up to 2000 spikes (subsample if more than
        % 2000 spikes in a unit)
        maxSpikesToGet = 20;
        if length(neuronDataStruct(iUnit).timeStamps) > maxSpikesToGet
            spikeSampleInds = neuronDataStruct(iUnit).timeStamps(...
                randperm(length(neuronDataStruct(iUnit).timeStamps), maxSpikesToGet));
        else
            spikeSampleInds = neuronDataStruct(iUnit).timeStamps;
        end

        waveforms = zeros(nPreSamples+nPostSamples+1,nChans,length(spikeSampleInds));
        for iSpike = 1:length(spikeSampleInds)
            % got to read in data point for each spike individually

                startReadSample = spikeSampleInds(iSpike) - nPreSamples*2;
                fseek(readFid, startReadSample * 2 * nChans, 'bof');
                rawData = fread(readFid, [nChans, (nPreSamples+nPostSamples)*2], 'int16=>double');

                %sometimes spike is too close to beginning or end of file,
                %don't use in this case
                if size(rawData,2) ~= (nPreSamples+nPostSamples)*2
                    waveforms(:,:,iSpike) = nan;
                    continue
                end

                %bandpass filter to denoise and remove offsets
                [b,a] = butter(3,[300 3000]/(sampleRate/2),'bandpass');
                filtData = filtfilt(b, a, rawData');
                waveforms(:,:,iSpike) = filtData(nPreSamples:(nPostSamples+nPreSamples*2),:);

        end

        meanWaveforms = squeeze(nanmean(waveforms,3));
        % find the biggest channel
        [~, maxChan] = max(max(abs(meanWaveforms(:,1:end-1))));
        
        % get amplitude and peak-to-valley
        % first get the peak (whether positive or negative)
        [peakValue, peakInd] = max(abs(meanWaveforms(:,maxChan)));
        % next, get the next peak after the peak (in opposite polarity)
        if meanWaveforms(peakInd,maxChan) < 0 
            [valleyValue, valleyInd] = max(meanWaveforms(peakInd+1:end,maxChan));
            valleyInd = valleyInd + peakInd;
            peakValue = peakValue * -1;
        elseif meanWaveforms(peakInd,maxChan) > 0 
            [valleyValue, valleyInd] = max(-1*meanWaveforms(peakInd+1:end,maxChan));
            valleyInd = valleyInd + peakInd;
            valleyValue = -1*valleyValue;
        else
            warning('Neuron''s peak is 0, most likely error in getting waveform');
        end

        % save to data struct
        neuronDataStruct(iUnit).waveforms = meanWaveforms;
        neuronDataStruct(iUnit).biggestChan = maxChan;
        neuronDataStruct(iUnit).depth = electrodPos(maxChan,2);
        neuronDataStruct(iUnit).x = electrodPos(maxChan,1);
        neuronDataStruct(iUnit).amplitude = abs(peakValue - valleyValue);
        neuronDataStruct(iUnit).peakToValley = abs(valleyInd - peakInd);


    end

    % finally, sort by depth
    depths = [neuronDataStruct.depth];
    [~, sortInds] = sort(depths);
    neuronDataStruct = neuronDataStruct(sortInds);

else

    load(fullfile(baseDIR, 'rez.mat'));

    % no phy curation, just use ks output directly
    nUnits = size(waveforms,1);

    %get average waveform and also channel number of that waveform
    aveWaveform = mean(waveforms,3);
    [maxValues,maxInds] = max(waveforms,[],2);
    [minValues,minInds] = min(waveforms,[],2);
    [amplitudes, waveformChans] = max(squeeze(maxValues - minValues),[],2);


    for iUnit = 1:nUnits

        neuronDataStruct(iUnit).waveforms = squeeze(waveforms(iUnit,:,:));
        neuronDataStruct(iUnit).biggestChan = waveformChans(iUnit);
        neuronDataStruct(iUnit).depth = electrodPos(waveformChans(iUnit),2);
        neuronDataStruct(iUnit).amplitude = amplitudes(iUnit);
        spikeEventInds = rez.st3(:,2) == iUnit;
        neuronDataStruct(iUnit).timeStamps = sort(rez.st3(spikeEventInds,1));

    end

end


%
