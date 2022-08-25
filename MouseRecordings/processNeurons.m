function neuronDataStruct = processNeurons(baseDIR)

waveforms = readNPY(fullfile(baseDIR, 'templates.npy'));
electrodPos = readNPY(fullfile(baseDIR, 'channel_positions.npy'));
load(fullfile(baseDIR, 'rez.mat'));


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