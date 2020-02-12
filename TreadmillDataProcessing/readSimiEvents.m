function [Events, EventsInfo] = readSimiEvents(simiFileNames, markerNames, varargin)
% [Events, EventsInfo] = readSimiEvents(simiFileNames, markerNames, [Events], [EventsInfo])
%
% Function to read in Simi .p files which contain the marked events and
% extract the frames where the events are marked. Assumes that the .p files
% follow the structure of: header lines followed by a blank line, followed
% by the marker IDs, followed by the marker names, then the data as
% delimited values. If you want to append the extracted events to an
% existing events struct, you can give that events struct as an input.
%
% Inputs:
% simiFileNames -   1xN cell array of strings containing the path and
%                   filename of the simi .p files containing the marked
%                   events. Each file should correspond to each recording
%                   block in the Events struct. If a block doesn't have a
%                   file, just use an empty string, ''.
%
% markerNames -     The names in the labels in the .p file whos events 
%                   fields we want to extract. Has to be an exact match 
%                   (case sensitive)
%
% Events -          Optional. Nx1 struct, with fields corresponding to
%                   various events (e.g. simi triggers, tones, ect), where
%                   N is the total number of blocks of the data (recording
%                   sessions). If this isn't given, the function will
%                   generate a new blank Events struct and add the simi
%                   events to that
%
% EventsInfo -      Optional. 1x1 Struct, containing the metadata for the 
%                   Events struct. Like Events, if this isn't the function
%                   will generate a new blank EventsInfo struct.
%
% Outputs:
% Events -          Will return the inputted Events struct with the newly
%                   extracted Simi events appended as a new field
%
% EventsInfo -      Will return the inputted EventsInfo struct with new
%                   metadata info about the new field that was added to
%                   Events.
%
% David Xing
% Last updated: 9/13/2019

% check inputs and make defaults
narginchk(2,4);

if nargin>2
    Events=varargin{1};
else
    Events=[];
end

if nargin>3
    EventsInfo=varargin{2};
else
    EventsInfo.fields = {};
    EventsInfo.description = {};
    EventsInfo.units = {};
    EventsInfo.processingStep = {};
end

if isa(simiFileNames, 'char')
    simiFileNames = {simiFileNames};
end

% go through each recording block
for iBlock = 1:length(simiFileNames)
    
    %if a file's name is an empty string, then theres no data for this
    %block
    if isempty(simiFileNames{iBlock})
        continue
    end
    
    % first, find when the header stops (blank line)
    fid=fopen(simiFileNames{iBlock});
    
    if fid==-1
        error('Unable to open file!');
    end
    
    nHeaders=0;
    while true
        line=fgetl(fid);
        nHeaders=nHeaders+1;
        
        %keep reading until we get to the blank line (also errors return a
        %blank line so there also has to be no error), or until we reach the
        %end of the file
        if (isempty(line) && isempty(ferror(fid))) | line==-1
            break
        end
    end
    fclose(fid);
    
    nHeaders=nHeaders+2;
    
    % now read in the actual data from the file
    fileData=importdata(simiFileNames{iBlock}, '\t', nHeaders);
    
    % get the headers and the numeric data
    headers=fileData.textdata(end,:);
    data=fileData.data;
    
    % go through each marker name whos events we want to extract
    for iMarker=1:length(markerNames)
        markerInd=find(strcmp(headers,markerNames{iMarker}),1);
        
        if isempty(markerInd)
            warning('Unable to find %s!',markerNames{iMarker})
        end
        
        if markerInd<=size(data,2)
            %get the frames where this marker was marked and save to events
            Events(iBlock).SimiEvents.(markerNames{iMarker})=find(~isnan(data(:,markerInd)));
        else
            %there was no data for this marker, just give empty struct
            Events(iBlock).SimiEvents.(markerNames{iMarker})=[];
        end
    end
    
end

% finally, add metadata to EventsInfo
EventsInfo.fields(end+1) = {'SimiEvents'};
EventsInfo.description(end+1) = {'Behavioral Events obtained using the video recording'};
EventsInfo.units(end+1) = {'Video frame number'};
EventsInfo.processingStep(end+1) = {'Behavioral Event Extraction'};


%
