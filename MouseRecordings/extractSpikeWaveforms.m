function extractSpikeWaveforms(baseDir,binNChans)

nWaveformChans = 12;
templateChansMethod = 'distance';

% filenames
binFile = 'temp_wh.dat';
spikeWaveformOutputFile = '_phy_spikes_subset.waveforms.npy';
spikeIndOutputFile = '_phy_spikes_subset.spikes.npy';
spikeChannelsOutputFile = '_phy_spikes_subset.channels.npy';
spikePCAFeaturesOutputFile = 'pc_features.npy';
spikePCAFeaturesChansOutputFile = 'pc_feature_ind.npy';

% load in phy data
load(fullfile(baseDir,'rez.mat'))
spikeTemplates = readNPY(fullfile(baseDir,'spike_templates.npy'));
spikeTimes = readNPY(fullfile(baseDir,'spike_times.npy'));
chanPhysCoords = readNPY(fullfile(baseDir,'channel_positions.npy'));
spikeTimes = readNPY(fullfile(baseDir,'spike_times.npy'));
templates = readNPY(fullfile(baseDir,'templates.npy'));

spikeTemplatesLabels = unique(spikeTemplates);

% see which channels correspond to which cluster templates
templateChanAmps = squeeze(max(templates,[],2)-min(templates,[],2));

% only save the <nWaveformChans> channels 
for iTemplate = 1:size(templateChanAmps,1)

    %can either get the channels associated with the template based on the
    %amplitude of the template waveforms, or based on the distance of the
    %channels from the main channel. Phy uses distance
    switch templateChansMethod
        case 'amplitude'
            [~, highestChanInds] = sort(templateChanAmps(iTemplate,:),'descend');
            templateUsedChans(iTemplate,:) = highestChanInds(1:nWaveformChans);
        case 'distance'
            % get the channel with the greatest amplitude.
            [~, maxAmpChan] = max(templateChanAmps(iTemplate,:));

            % get the closest channels
            dists = sqrt((rez.xcoords-rez.xcoords(maxAmpChan)).^2+(rez.ycoords-rez.ycoords(maxAmpChan)).^2);
            [~, closestChanInds] = sort(dists,'ascend');
            templateUsedChans(iTemplate,:) = closestChanInds(1:nWaveformChans);
    end
end

% _phy_spikes_subset.waveforms.npy should be a N x T x C matrix with N =
% the total number of spike instances, T is the 82 time points for the
% waveforms (centered at 41), and C are the 16 channels. Channels 13-16 are
% usually all just zeros. Because N is usually very large (can be > 4
% million spikes) will run into a memory problem if I try to hold all the
% data at once (>30 GB), so I will save it in chunks. Unfortunately, since
% _phy_spikes_subset.waveforms.npy saves the data along the first dimension
% first, then the second dimensions, then the third dimension, when I write
% the file I have to chunk along the third dim (the channels dim), which
% will be somewhat inefficient.

nChunks = 2; % how to split up the channels.
chunkNChans = ceil(nWaveformChans/nChunks);
nPreSamples = 41;
nPostSamples =  41;
totalSaveFileChans = 16;
nPCAChans = 6; %number of channels to do PCA on to get features, make sure this is less than chunkNChans
nPCAFeatures = 4;
rawDataFid = fopen(fullfile(baseDir, binFile), 'rb');

% initiate file to save waveforms to
header = constructNPYheader('double', [length(spikeTimes) nPreSamples+nPostSamples, totalSaveFileChans]);
waveformFid = fopen(spikeWaveformOutputFile, 'w');
fwrite(waveformFid, header, 'uint8');

% initate PC features matrix
pcFeatures = zeros(length(spikeTimes),nPCAFeatures,nPCAChans);

for iChunk = 1:nChunks

    % channels for this chunk
    chunkChanInds = chunkNChans*(iChunk-1)+1:min(chunkNChans*(iChunk),nWaveformChans);
    % initiate matrix for the waveforms of this chunk
    chunkWaveforms = zeros(length(spikeTimes),nPreSamples+nPostSamples,chunkNChans);

    tic
    for iSpike = 1:length(spikeTimes)
        
        % load in the data
        startReadSample = spikeTimes(iSpike) - nPreSamples;
        fseek(rawDataFid, startReadSample * 2 * binNChans, 'bof');
        rawData = fread(rawDataFid, [binNChans, (nPreSamples+nPostSamples)], 'int16=>double');

        % which channels to use
        thisSpikeTemplate = spikeTemplates(iSpike)+1; %python starts at 0
        thisSpikeChans = templateUsedChans(thisSpikeTemplate,chunkChanInds);
        
        % save data to chunk waveform matrix
        chunkWaveforms(iSpike,:,:) = rawData(thisSpikeChans,:)';

    end
    disp(['Loaded in data for chunk ' num2str(iChunk) ', time: ' num2str(toc) 's'])

    % save data
    tic
    fwrite(waveformFid, chunkWaveforms, 'double');
    disp(['Saved data for chunk ' num2str(iChunk) ', time: ' num2str(toc) 's'])

    % run PCA
    if iChunk == 1
        for iTemplate = 1:length(spikeTemplatesLabels)

            templateSpikeInds = spikeTemplates == spikeTemplatesLabels(iTemplate);
            for iChan = 1:nPCAChans
                [~, pcaPoints] = pca(squeeze(chunkWaveforms(templateSpikeInds,:,iChan)));
                pcFeatures(templateSpikeInds,1:nPCAFeatures,iChan) = pcaPoints(:,1:nPCAFeatures);
            end
        end
    end

end

% finally, save the zeros at the end
waveformZeros = zeros(length(spikeTimes),nPreSamples+nPostSamples,totalSaveFileChans-nWaveformChans);
fwrite(waveformFid, waveformZeros, 'double');

% also save spike index and channels
header = constructNPYheader('int64', [length(spikeTimes) 1]);
spikeIndsFid = fopen(spikeIndOutputFile, 'w');
fwrite(spikeIndsFid, header, 'uint8');
fwrite(spikeIndsFid, int64((1:length(spikeTimes))-1), 'int64'); %python starts at 0
fclose(spikeIndsFid);

% channels
header = constructNPYheader('int32', [length(spikeTimes) totalSaveFileChans]);
spikeChansFid = fopen(spikeChannelsOutputFile, 'w');
fwrite(spikeChansFid, header, 'uint8');
waveformChans = [templateUsedChans(spikeTemplates+1,:) zeros(length(spikeTimes),totalSaveFileChans-nWaveformChans)];
fwrite(spikeChansFid, int32(waveformChans-1), 'int32'); %python starts at 0
fclose(spikeChansFid);

% save PCA Features
header = constructNPYheader('double', [length(spikeTimes) nPCAFeatures nPCAChans]);
spikePCAFid = fopen(spikePCAFeaturesOutputFile, 'w');
fwrite(spikePCAFid, header, 'uint8');
fwrite(spikePCAFid, pcFeatures, 'double');
fclose(spikePCAFid);

header = constructNPYheader('int32', [size(templateUsedChans,1) nPCAChans]);
spikePCAChansFid = fopen(spikePCAFeaturesChansOutputFile, 'w');
fwrite(spikePCAChansFid, header, 'uint8');
fwrite(spikePCAChansFid, int32(templateUsedChans(:,1:nPCAChans)-1), 'int32'); %python starts at 0
fclose(spikePCAChansFid);


fclose(waveformFid);
fclose(rawDataFid);


% 
