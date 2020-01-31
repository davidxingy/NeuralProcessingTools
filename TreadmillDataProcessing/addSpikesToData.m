function [ContinuousData, ContinuousDataInfo, Events, EventsInfo] = addSpikesToData(...
    neuronInfos, outputLabels, simiFrameIndsLabels, ContinuousData, ContinuousDataInfo, Events, EventsInfo)
% [ContinuousData, ContinuousDataInfo, Events, EventsInfo] = addSpikesToData(...
%     neuronInfos, outputLabels, simiFrameIndsLabels, ContinuousData, ContinuousDataInfo, Events, EventsInfo)
% 
% Function to add spikes to the Events data structure and bin spikes and
% save it to the ContinuousData structure. Will be binned at 1ms windows as
% well as windows defined by the Simi camera frames (usually 10ms)
% 
% Inputs:
% neuronInfos -         1xN Cell array. Contains the consolidated
%                       data struct that's the output of the spike sorting
%                       process. Each of the N cells comes from a different
%                       MEA.
% 
% outputLabels -        1xN Cell array. The names of the different MEAs
%                       that we want to use in the data structs.
% 
% simiFrameIndsLabels - 1xN Cell array. I'm assuming each MEA goes to a
%                       separate NSP so each MEA will correspond to a
%                       different set of Simi trigger times. This input
%                       should contain the field names of the different
%                       Simi trigger times in the Events data struct.
% 
% ContinuousData -      1xP struct. Contains the continuous data to add the
%                       binned spike counts to. P is the total number of
%                       recording blocks
% 
% ContinuousDataInfo -  1x1 struct. Contains the metadata for the fields of
%                       ContinuousData
% 
% Events -              1xP struct. Contains the events timestamps to add
%                       the spike timestamps to.
% 
% EventsInfo -          1x1 struct. Contains the metadata for the fields of
%                       Events
%
% Outputs:
% ContinuousData -      Will return the inputted ContinuousData struct with
%                       the binned spike counts added to it.
%
% ContinuousDataInfo -  Will return the inputted ContinuousDataInfo struct
%                       with new metadata info about the new fields added to
%                       ContinuousData
% 
% Events -              Will return the inputted Events struct with the
%                       spike timestamps added to it.
%
% EventsInfo -          Will return the inputted EventsInfo struct with new
%                       metadata info about the new fied added to Events
%
% David Xing
% Last updated: 9/13/2019


% number of labels and number of arrays must match
assert(length(neuronInfos) == length(outputLabels), ...
    'Number of labels and number of neuronal data structs must match!');
assert(~isempty(neuronInfos), 'No data in neuronInfos!');
assert(~isempty(neuronInfos{1}), 'No data in neuronInfos!');
assert(~isempty(neuronInfos{1}(1).Timestamps), 'No data in neuronInfos!');

% number of blocks must match
assert(length(neuronInfos{1}(1).Timestamps) == length(Events), ...
    'Number of blocks in the sorted spikes and the data structs must match!');

% Simi frame times must be saved in the Events data struct
assert(isfield(Events, 'simiTriggers'), 'Events must have the simi camera triggers!')

% Number of arrays must match number of simi frame times (I'm assuming 1
% NSP per array)
assert(length(fieldnames(Events(1).simiTriggers)) == length(neuronInfos), ...
    'Each neuron data struct should correspond to a set set of Simi Frame times!');

nArrays = length(neuronInfos);
nBlocks = length(Events);

% go through each array
for iArray = 1:nArrays
    
    %go through each block
    for iBlock = 1:nBlocks

        %how long the recording is for this block
        recordingLength = Events(iBlock).simiTriggers.(simiFrameIndsLabels{iArray})(end);

        %save the timestamps to events
        Events(iBlock).SortedTimeStamps.(outputLabels{iArray}) = ...
            cellfun(@(x) x{iBlock}, {neuronInfos{iArray}.Timestamps}, 'un', 0);
        
        %add binned timestamps to continuous data
        nNeurons = length(neuronInfos{iArray});
        
        binEdges1ms = 0:30:recordingLength;
        for iNeuron = 1:nNeurons
            %bin timestamps by 1ms
            ContinuousData(iBlock).SpikeTrains1ms.(outputLabels{iArray})(iNeuron,:) = ...
                histcounts(neuronInfos{iArray}(iNeuron).Timestamps{iBlock}*30,binEdges1ms);
            
            if (iNeuron==1)
                ContinuousData(iBlock).SpikeTrains1msInds.(outputLabels{iArray}) = binEdges1ms(2:end);
            end
            
            %bin timestamps by simi frames (usually 10ms)
            ContinuousData(iBlock).SpikeCountsSimiTrig.(outputLabels{iArray})(iNeuron,:) = ...
                histcounts(neuronInfos{iArray}(iNeuron).Timestamps{iBlock}*30,...
                Events(iBlock).simiTriggers.(simiFrameIndsLabels{iArray}));
            
        end

    end
    
end

% finally, update the meta data for the ContinuousData and Events structs
ContinuousDataInfo.fields(end+1:end+3) = {'SpikeTrains1ms'; 'SpikeTrains1msInds'; ...
    'SpikeCountsSimiTrig'};
                         
ContinuousDataInfo.description(end+1:end+3) = {'Number of spikes in 1ms bins'; ...
    'Indices of the edges of the 1ms spike train bins'; ...
    'Number of spikes in 1ms bins'};

ContinuousDataInfo.units(end+1:end+3) = {'number of spikes'; ...
                            'NSP samples'; ...
                            'number of spikes'};

ContinuousDataInfo.processingStep(end+1:end+3) = repmat({'Spike Sorting'}, 3, 1);

EventsInfo.fields(end+1) = {'SortedTimeStamps'};
                 
EventsInfo.description(end+1) = {'Spike times for all the sorted neurons'};

EventsInfo.units(end+1) = {'miliseconds'};
EventsInfo.processingStep(end+1) = {'Spike Sorting'};
                    