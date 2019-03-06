function Events=readSimiEvents(fullFileName,trialNum,markerNames,varargin)
% 
% Function to read in Simi .p files which contain the marked events and
% extract the frames where the events are marked. Assumes that the .p files
% follow the structure of header lines followed by a blank line, followed
% by the marker IDs, followed by the marker names, then the data as
% delimited values. If you want to append the extracted events to an
% existing events struct, you can give that events struct as an input.
% 

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
        Events(trialNum).SimiEvents.(markerNames{iMarker})=find(~isnan(data(:,markerInd)));
    else
        %there was no data for this marker, just give empty struct
        Events(trialNum).SimiEvents.(markerNames{iMarker})=[];
    end
end



%
