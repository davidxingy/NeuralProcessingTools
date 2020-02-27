function simiEvents = extractTrackedSimiEvents(filenames, markerNames)
% simiEvents = extractTrackedSimiEvents(filenames, markerNames)
% 
% To get event times from camera recordings using Simi, I define some 
% markers to represent events (e.g. toe_off, treat_contact, ect). In the
% frames where an event happens, I click somewhere in the frame to add data
% for the marker that represents the event (the actual pixel locations
% don't matter, only that there is some data during that frame). This
% function reads in a bunch of Simi .p files and extract the frame times
% for all the markers that represent events by extracing the frames where
% data exists for those markers.
% 
% Inputs:
% filenames     - String, or cell array of strings that gives the full path
%                 to the .p files to extract the events from
% 
% markerNames   - Cell array of strings of the exact marker names in the .p
%                 files which represents events (as opposed to markers that
%                 represent continuous data such as kinematics).
% 
% Outputs:
% simiEvents -  a N-sized struct (where N is the number of files, same size
%               as filenames) containing the extracted events. The fields 
%               of the struct are the inputted markerNames, and the event
%               times are given as frame numbers.
% 
% David Xing
% Last updated: 9/10/2019

% check inputs
narginchk(2,3)

% if only one file was given (as a string), convert to cell array
if isa(filenames,'char')
    filenames = {filenames};
end

% go through each file
for iFile = 1:length(filenames)
    
    %open file
    fID = fopen(filenames{iFile});
    
    if fID==-1
        error('Unable to open file %s', filenames{iFile})
    end
    
    %read header
    header = string([]);
    while true
        header(end+1) = fgetl(fID);
        
        if isempty(header{end}) || strcmp(header{end},'-1')
            break
        end
    end
    
    %get markers name and ids (for some reason split with '\t'
    %doesn't work, have to use char(9))
    fileMarkerIDs = split(string(fgetl(fID)));
    fileMarkerIDs = fileMarkerIDs(1:2:end);
    fileMarkerNames = split(string(fgetl(fID)),char(9));
    
    %load the data
    fclose(fID);
    data = importdata(filenames{iFile},'\t',length(header)+2);
    data = data.data;
    
    %not all file markers have corresponding data
    fileMarkerNamesWithData = fileMarkerNames(1:size(data,2));
    
    %get the marker events that I want
    for iMarker = 1:length(markerNames)
        
        markerInd = find(strcmp(markerNames{iMarker}, fileMarkerNamesWithData),1);
        if isempty(markerInd)
            warning('%s wasn''t found in file %s!',markerNames{iMarker}, filenames{iFile})
            simiEvents(iFile).(markerNames{iMarker}) = [];
            continue
        end
        
        simiEvents(iFile).(markerNames{iMarker}) = find(~isnan(data(:,markerInd)));
    end
    
end

% save events to .mat file
% save('ExtractedEvents.mat', simiEvents)


% 
