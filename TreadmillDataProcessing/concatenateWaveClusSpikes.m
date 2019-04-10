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
%                       trials are contained in this folder.
% 
% filenameStructure -   String that defines the formatting of the
%                       times_...mat filenames. Put a '%c' where the
%                       channel number should be and a '%t' where the trial
%                       number should be. Ex:
%                       'times_blahblah_Chan%c_Trial%t.mat', will look for
%                       all files in the folder that follow that filename
%                       convention, such as times_blahblah_Chan96_Trial1.mat
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
%                       the spikes of all the trials)
%                       WaveformSTD: the standard deviation at each point
%                       of the waveform (across all the spikes of all the
%                       trials)
%                       OriginChannel: The waveclus file channel that the
%                       neuron came from.
%                       OriginUnitNum: The unit label number that the
%                       neuron was assigned to in the waveclus file
%                       Timestamps: a cell array containing the timestamps
%                       (in whatever units that was output by WaveClus),
%                       with each cell corresponding to a trial.
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
trialNumInd = strfind(filenameStructure,'%t');

% make sure at least channel # location was found
if isempty(chanNumInd)
    error('filenameStructure must containt ''%c'' to indicate where the channel num is!');
elseif isempty(trialNumInd)
    warning('''%t'' not found in filenameStructure, will assume there is only 1 trial');
end

% change the filename structure to format that can be used by sscanf
filenameStructureMod = filenameStructure;
filenameStructureMod(chanNumInd+1) = 'u';
filenameStructureMod(trialNumInd+1) = 'u';

% determine whether the trial # comes first or the channel #
if isempty(trialNumInd) || chanNumInd < trialNumInd
    chanIndFirst = true;
else
    chanIndFirst = false;
end

% now go through all the files and get the channel numbers and trials of
% all the files that match the format
chanList = [];
trialList = [];
fileIndList = [];
for iFile = 1:length(allFileNames)
    
    %match current file to templte
    foundInds = sscanf(allFileNames{iFile}, filenameStructureMod);
    
    %no match found, go to next filename
    if length(foundInds)~=1 && length(foundInds)~=2
        continue
    end
    
    %template found, assign channel and trial number
    if isempty(trialNumInd)
        chanNum = foundInds;
        trialNum = 1;    
    elseif chanIndFirst
        chanNum = foundInds(1);
        trialNum = foundInds(2);
    else
        trialNum = foundInds(1);
        chanNum = foundInds(2);
    end
    
    chanList(end+1) = chanNum;
    trialList(end+1) = trialNum;
    fileIndList(end+1) = iFile;
    
end

% get list of all the unique channel numbers that are contained in the
% files
if isempty(trialNumInd)
    allTrials = 1;
else
    allTrials = unique(trialList);
end

% get list of all the unique channel numbers that are contained in the
% files
allChans = unique(chanList);

% and go through each channel to extract its neurons (for all trials)
iNeuron = 1; %for keeping track of the number of neurons we are saving
for iChan = 1:length(allChans)
    
    dontUseChan = false;

    %go through each trial and load in the data
    for iTrial = 1:length(allTrials)
        
        filenameInd = find(chanList==allChans(iChan) & trialList==allTrials(iTrial));
        
        %warn user if the trial wasn't found
        if length(filenameInd)~=1
            warning(['No file found for Trial %u of channel %u!' ...
                'Will not save neurons from this channel'], allTrials(iTrial), allChans(iChan));
            dontUseChan = true;
            break;
        end
            
        
        %load in the data from the file
        try
            fileVars{allTrials(iTrial)} = load(fullfile(folder, ...
                allFileNames{fileIndList(filenameInd)}));
        catch
            warning(['Unable to load channel %u file: %s. ' ...
                'Will not save neurons from this channel'], allChans(iChan), ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
        end
        
        %Make sure all the proper variables are in the file
        if ~isfield(fileVars{allTrials(iTrial)},'cluster_class')
            warning(['File %s doesn''t have the ''cluster_class'' variable! ' ...
                'Will not save neurons from this channel'], ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
            
        elseif ~isfield(fileVars{allTrials(iTrial)},'spikes')
            warning(['File %s doesn''t have the ''spikes'' variable! ', ...
                'Will not save neurons from this channel'], ...
                allFileNames{fileIndList(filenameInd)});
            dontUseChan = true;
            break;
        end

        %see how many neurons were found for this channel
        trialUnitLabels{allTrials(iTrial)} = unique(fileVars{iTrial}.cluster_class(:,1));
        trialUnitLabels{allTrials(iTrial)}(trialUnitLabels{allTrials(iTrial)}==0) = [];
        
        %make sure the neuron labels are the same across trials
        firstNonEmptyTrial = find(~cellfun(@isempty, fileVars),1);
        if any(trialUnitLabels{allTrials(iTrial)}~=trialUnitLabels{firstNonEmptyTrial})
            warning(['Different neuron labels found for trials %u and %u in channel %u! ' ...
                'Will not save neurons from this channel.'], allTrials(iTrial),...
                firstNonEmptyTrial, allChans(iChan))
            dontUseChan = true;
            break
        end
    end
    
    if dontUseChan
        %Something wrong, don't use this channel
        continue
    end

    %if code reached this point, that means all the neuron labels were
    %consistent across trials
    firstNonEmptyTrial = find(~cellfun(@isempty, fileVars),1);
    unitLabels = trialUnitLabels{firstNonEmptyTrial};

    %now go through each putative neuron
    for iUnit = 1:length(unitLabels)
        
        %calculate the SNR for each trial
        for iTrial = 1:length(allTrials)
            
            waveforms{iTrial} = fileVars{allTrials(iTrial)}.spikes(...
                fileVars{allTrials(iTrial)}.cluster_class(:,1)==unitLabels(iUnit),:);
            noise = waveforms{iTrial}(:,1:5);
            SNR(iTrial) = (max(nanmean(waveforms{iTrial}))-...
                min(nanmean(waveforms{iTrial})))/(3*nanstd(noise(:)));
            
        end
        
        %if the average SNR across trials passes threshold, keep neuron
        if mean(SNR)<SNRCutoff
            continue
        else
            neuronInfo(iNeuron).SNR = mean(SNR);
            neuronInfo(iNeuron).Waveform = nanmean(cat(1,waveforms{:}));
            neuronInfo(iNeuron).WaveformStd = nanstd(cat(1,waveforms{:}));
            neuronInfo(iNeuron).OriginChannel = allChans(iChan);
            neuronInfo(iNeuron).OriginUnitNum = unitLabels(iUnit);
        end
        
        %now finally, get the timestamps
        for iTrial = 1:length(allTrials)
            timestamps = fileVars{allTrials(iTrial)}.cluster_class(...
                fileVars{allTrials(iTrial)}.cluster_class(:,1)==unitLabels(iUnit),2);
            
            neuronInfo(iNeuron).Timestamps{allTrials(iTrial)} = timestamps;
        end
        
        %increment neuron counter
        iNeuron = iNeuron+1;
    end
    
end


% 
