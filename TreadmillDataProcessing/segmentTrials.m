function dataSegs = segmentTrials(continuousData,eventsTimestamps,varargin)
% dataSegs = segmentTrials(continuousData,eventsTimestamps,...
%   [preTrialsamples],[postTrialsamples],[segIncludesStartEvent])
% 
% Function to segment continuous data (such as neural spike counts) into
% trials based on inputted timestamps (in samples, though the ts can be a
% decimal even though it's in samples). Additionally, can get data
% segments that precedes the first event ts or after the last event ts by a
% certain amount. If segIncludesStartEvent is true, each segment will contain 
% the earliest datapoint that precedes the segment start timestamp and if the 
% ts is exactly on a sample, the segment will include that sample. The 
% segement will not include the segment end timestamp. If
% segIncludesStartEvent is false, then the starting ts sample isn't
% included in the segment, but the segment end ts is included.
% E.g.:
% segIncludesStartEvent=true
%       Seg start ts -?             ?-Seg start ts
% sample:   [1 2 3 4][5 6 7 8 9 10][11 12 13 ...
%             seg 1      seg 2          seg3
% 
% segIncludesStartEvent=false
%      Seg Start ts -?              ?-Seg start ts
% sample:   [1 2 3 4 5][6 7 8 9 10 11][12 13 ...
%             seg 1         seg 2        seg3
% 
% Inputs:
% continuousData -      NxM array containing the data that will be
%                       segmented .The segmentation will happen along the
%                       2nd dimension (M).
% 
% eventsTimestamps -    PxT array containing the timestamps (in samples)
%                       that demarks when the trial segments start and end.
%                       P is the total number of trials that will be
%                       segmented, and T-1 is the total number of segments
%                       per trial (The segments for each trial are
%                       continuous with each other).
% 
% preTrialsamples -     Optional. If you want to get the data that happens
%                       before the start of the trial, specify how many
%                       data points with this input. Default 0
% 
% postTrialsamples -    Optional. If you want to get the data that happens
%                       after the end of the trial, specify how many
%                       data points with this input. Default 0
% 
% segIncludesStartEvent - Optional. Determines how to segement the data. If
%                         true, then will get the data including the
%                         starting timestamp, up to, but not including the
%                         stop timestamp. If false, then will get the data
%                         not including the starting timestamp up to, and
%                         including, the stop timestamp. Default true.
% 
% Outputs:
% dataSegs -    Px(T+1) cell array. Each row is a trial, and each column is
%               a segment. The first and last columns are always the
%               segements preceding and following the actual trial segments
%               respectively.
% 
% David Xing 3/7/2019

% check inputs and set defaults
narginchk(2,5)

if nargin>2
    if isempty(varargin{1})
        preTrialsamples = [];
    else
        preTrialsamples = varargin{1};
    end
else
    preTrialsamples = [];
end

if nargin>3
    if isempty(varargin{2})
        postTrialsamples = [];
    else
        postTrialsamples = varargin{2};
    end
else
    postTrialsamples = [];
end

if nargin>4
    if isempty(varargin{3})
        segIncludesStartEvent = true;
    else
        segIncludesStartEvent = varargin{3};
    end
else
    segIncludesStartEvent = true;
end

% get how many trials and how many segments
nTrials = size(eventsTimestamps,1);
nSegs = size(eventsTimestamps,2)-1;

dataSegs = {};

% convert event timestamps to data indices
if segIncludesStartEvent
    eventsTimestamps = floor(eventsTimestamps);
else
    eventsTimestamps = floor(eventsTimestamps+1);
end

% go through and segment the data
for iTrial = 1:nTrials
    
    %first, get the pre-trial data
    if isempty(preTrialsamples) || preTrialsamples==0
        dataSegs{iTrial,1} = [];
    else 
        dataSegs{iTrial,1} = continuousData(:,...
            eventsTimestamps(iTrial,1)-preTrialsamples:eventsTimestamps(iTrial,1)-1);
    end
    
    %next go through all the events
    for iEvent = 1:nSegs
        dataSegs{iTrial,iEvent+1} = continuousData(:,...
            eventsTimestamps(iTrial,iEvent):eventsTimestamps(iTrial,iEvent+1)-1);
    end
    
    %finally, get the post-trial data
    if isempty(postTrialsamples) || postTrialsamples==0
        dataSegs{iTrial,nSegs+2} = [];
    else
        dataSegs{iTrial,nSegs+2} = continuousData(:,...
            eventsTimestamps(iTrial,end):eventsTimestamps(iTrial,end)+postTrialsamples-1);
    end
    
end


% 
