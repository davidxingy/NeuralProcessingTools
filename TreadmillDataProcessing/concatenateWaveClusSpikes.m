function neuronInfo = concatenateWaveClusSpikes(folder,filenameStructure,varargin)
% neuronInfo = concatenateWaveClusSpikes(folder, filenameStructure, [SNRCutoff])
% 
% Function that takes the output of wave_clus-SWADE spike sorting (so
% times_....mat files that contain the timestamps of sorted neurons) and
% consolidates all the found units into a struct that contains all of the
% neurons (containing its waveform, timestamps, SNR, ect). Also allows you
% to optionally reject units that don't pass a certain SNR cutoff.
% 
% Inputs:
% folder -              String containing the path to folder containing the
%                       times_....mat files that were the output of
%                       WaveClus sorting. Make sure all your channels and
%                       blocks are contained in this folder.
% 
% filenameStructure -   String that defines the formatting of the
%                       times_...mat filenames. Put a '%c' where the
%                       channel number should be and a '%b' where the block
%                       number should be. Ex:
%                       'times_blahblah_Chan%c_Block%b.mat', will look for
%                       all files in the folder that follow that filename
%                       convention, such as times_blahblah_Chan96_Block1.mat
% 
% SNRCutoff -           Optional. If given, will only keep units that pass
%                       this SNR cutoff. SNR is calculated as the amplitude
%                       of the average waveform over 3x the standard
%                       deviation of the noise (noise is calculated from
%                       the first 5 samples of the waveform, across all
%                       waveforms). If not given, won't do any rejection of
%                       units.
% 
% Outputs:
% neuronInfo -          Nx1 struct containing the consolidated neurons. Has
%                       the following fields:
%                       SNR: signal to noise ration of the neuron, as
%                       defined in the description above in SNRCutoff
%                       Waveform: average waveform the neuron (across all
%                       the spikes of all the blocks)
%                       WaveformSTD: the standard deviation at each point
%                       of the waveform (across all the spikes of all the
%                       blocks)
%                       OriginChannel: The waveclus file channel that the
%                       neuron came from.
%                       OriginUnitNum: The unit label number that the
%                       neuron was assigned to in the waveclus file
%                       Timestamps: a cell array containing the timestamps
%                       (in whatever units that was output by WaveClus),
%                       with each cell corresponding to a block.
% 
% David Xing 4/10/2019

% check inputs
narginchk(2,3)

if nargin>2
    if isempty(varargin{1})
        SNRCutoff = 0;
    else
        SNRCutoff = varargin{1};
    end
else
    SNRCutoff = 0;
end

% first check that the folder with the spike data exists
if ~exist(folder,'dir')
    error('Folder containing the input spike data doesn''t exist!') 
end
% get the list of files in the folder to search fo the ones that follow the
% filename structure that contain the spikes
allFileNames = dir(folder);
allFileNames = string({allFileNames.name})';

% get the channel # location in the filename string
chanNumInd = strfind(filenameStructure,'%c');

% get the trail # location in the filename string
blockNumInd = strfind(filenameStructure,'%b');

% make sure at least channel # location was found
if isempty(chanNumInd)
    error('filenameStructure must containt ''%c'' to indicate where the channel num is!');
elseif isempty(blockNumInd)
    warning('''%b'' not found in filenameStructure, will assume there is only 1 block');
end

% change the filename structure to format that can be used by sscanf
filenameStructureMod = filenameStructure;
filenameStructureMod(chanNumInd+1) = 'u';
filenameStructureMod(blockNumInd+1) = 'u';

% determine whether the block # comes first or the channel #
if isempty(blockNumInd) || chanNumInd < blockNumInd
    chanIndFirst = true;
else
    chanIndFirst = false;
end

% now go through all the files and get the channel numbers and blocks of
% all the files that match the format
chanList = [];
blockList = [];
fileIndList = [];
for iFile = 1:length(allFileNames)
    
    %match current file to template
    foundInds = sscanf(allFileNames{iFile}, filenameStructureMod);
    
    %no match found, go to next filename
    if length(foundInds)~=1 && length(foundInds)~=2
        continue
    end
    
    %template found, assign channel and block number
    if isempty(blockNumInd)
        chanNum = foundInds;
        blockNum = 1;    
    elseif chanIndFirst
        chanNum = foundInds(1);
        blockNum = foundInds(2);
    else
        blockNum = foundInds(1);
        chanNum = foundInds(2);
    end
    
    chanList(end+1) = chanNum;
    blockList(end+1) = blockNum;
    fileIndList(end+1) = iFile;
    
end

% get list of all the unique channel numbers that are contained in the
% files
if isempty(blockNumInd)
    allBlocks = 1;
else
    allBlocks = unique(blockList);
end

% get list of all the unique channel numbers that are contained in the
% files
allChans = unique(chanList);

% and go through each channel to extract its neurons (for all blocks)
iNeuron = 1; %for keeping track of the number of neurons we are saving
for iChan = 1:length(allChans)
    
    dontUseChan = false;

    %go through each block and load in the data
    for iBlock = 1:length(allBlocks)
        
        filenameInd = find(chanList==allChans(iChan) & blockList==allBlocks(iBlock));
        
        %warn user if the block wasn't found
        if length(filenameInd)~=1
            warning(['No file found for Block %u of channel %u!' ...
                'Will not save neurons from this channel'], allBlocks(iBlock), allChans(iChan));
            dontUseChan = true;
            break;
        end
            
        
        %load in the data from the file
        try
            fileVars{allBlocks(iBlock)} = load(fullfile(folder, ...
                allFileNames{fileIndList(filenameInd)}));
        catch
            warning(['Unable to load channel %u file: %s. ' ...
                'Will not save neurons from this channel'], allChans(iChan), ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
        end
        
        %Make sure all the proper variables are in the file
        if ~isfield(fileVars{allBlocks(iBlock)},'cluster_class')
            warning(['File %s doesn''t have the ''cluster_class'' variable! ' ...
                'Will not save neurons from this channel'], ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
            
        elseif ~isfield(fileVars{allBlocks(iBlock)},'spikes')
            warning(['File %s doesn''t have the ''spikes'' variable! ', ...
                'Will not save neurons from this channel'], ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
        end

        %see how many neurons were found for this channel
        blockUnitLabels{allBlocks(iBlock)} = unique(fileVars{iBlock}.cluster_class(:,1));
        blockUnitLabels{allBlocks(iBlock)}(blockUnitLabels{allBlocks(iBlock)}==0) = [];
        
        %make sure the neuron labels are the same across blocks
        firstNonEmptyBlock = find(~cellfun(@isempty, fileVars),1);
        if any(blockUnitLabels{allBlocks(iBlock)}~=blockUnitLabels{firstNonEmptyBlock})
            warning(['Different neuron labels found for blocks %u and %u in channel %u! ' ...
                'Will not save neurons from this channel.'], allBlocks(iBlock),...
                firstNonEmptyBlock, allChans(iChan))
            dontUseChan = true;
            break
        end
    end
    
    if dontUseChan
        %Something wrong, don't use this channel
        continue
    end

    %if code reached this point, that means all the neuron labels were
    %consistent across blocks
    firstNonEmptyBlock = find(~cellfun(@isempty, fileVars),1);
    unitLabels = blockUnitLabels{firstNonEmptyBlock};

    %now go through each putative neuron
    for iUnit = 1:length(unitLabels)
        
        %calculate the SNR for each block
        for iBlock = 1:length(allBlocks)
            
            waveforms{iBlock} = fileVars{allBlocks(iBlock)}.spikes(...
                fileVars{allBlocks(iBlock)}.cluster_class(:,1)==unitLabels(iUnit),:);
            
            noise = waveforms{iBlock}(:,1:5);
            SNR(iBlock) = (max(nanmean(waveforms{iBlock}))-...
                min(nanmean(waveforms{iBlock})))/(3*nanstd(noise(:)));
            
        end
        
        %if the average SNR across blocks passes threshold, keep neuron
        if nanmean(SNR)<SNRCutoff
            continue
        else
            neuronInfo(iNeuron).SNR = nanmean(SNR);
            neuronInfo(iNeuron).Waveform = nanmean(cat(1,waveforms{:}));
            neuronInfo(iNeuron).WaveformStd = nanstd(cat(1,waveforms{:}));
            neuronInfo(iNeuron).OriginChannel = allChans(iChan);
            neuronInfo(iNeuron).OriginUnitNum = unitLabels(iUnit);
        end
        
        %now finally, get the timestamps
        for iBlock = 1:length(allBlocks)
            timestamps = fileVars{allBlocks(iBlock)}.cluster_class(...
                fileVars{allBlocks(iBlock)}.cluster_class(:,1)==unitLabels(iUnit),2);
            timestamps(isnan(timestamps)) = [];
            
            neuronInfo(iNeuron).Timestamps{allBlocks(iBlock)} = timestamps;
        end
        
        %increment neuron counter
        iNeuron = iNeuron+1;
    end
    
end


% 
