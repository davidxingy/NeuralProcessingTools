function Events = readSimiEvents(fullFileName, blockNum, markerNames, varargin)
% Events = readSimiEvents(fullFileName, blockNum, markerNames, [Events])
% 
% Function to read in Simi .p files which contain the marked events and
% extract the frames where the events are marked. Assumes that the .p files
% follow the structure of: header lines followed by a blank line, followed
% by the marker IDs, followed by the marker names, then the data as
% delimited values. If you want to append the extracted events to an
% existing events struct, you can give that events struct as an input.
% 
% Inputs:
% fullFileName -    String containing the path and filename of the simi .p
%                   file containing the marked events
% 
% blockNum -        Which block to add the new extracted events to (since
%                   the Events struct is arranged by blocks)
% 
% markerNames -     The names in the labels in the .p file whos events we
%                   want to extract. Has to be an exact match (case
%                   sensitive)
% 
% Events -          Optional. Nx1 struct, with fields corresponding to
%                   various events (e.g. simi triggers, tones, ect), where
%                   N is the total number of blocks of the data (recording
%                   sessions).
% 
% Outputs:
% Events -          If given as an input, will append the new extracted
%                   events into that input, otherwise, will just return a
%                   new Event struct which only contains the SimiEvents
%                   field.

% check inputs and make defaults
narginchk(3,4);

if nargin>3
    Events=varargin{1};
else
    Events=[];
end

% first, find when the header stops (blank line)
fid=fopen(fullFileName);

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
fileData=importdata(fullFileName, '\t', nHeaders);

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
        Events(blockNum).SimiEvents.(markerNames{iMarker})=find(~isnan(data(:,markerInd)));
    else
        %there was no data for this marker, just give empty struct
        Events(blockNum).SimiEvents.(markerNames{iMarker})=[];
    end
end



%
