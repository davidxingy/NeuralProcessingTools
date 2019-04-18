function normalizedData = timeNormalize(data, eventInds, eventPercentages, desiredPercentages)
% normalizedData = timeNormalize(data, eventInds, eventPercentages, [desiredPercentages])
% function to normalize timeseries data of various lengths to 0-100
% percent. Will use linear normalization, and normalize in segments based
% on events markers (if any are given). Will always normalize the first
% data point to be 0% and the last data point to 100%, unless overwritten
% by event markers (e.g. you can have 1 in the eventInds input and 5% in
% the corresponding eventPercentages input and now the first index will be
% 5% rather than 0%).
% 
% Inputs:
% data -                MxN matrix, containing the time series data. Each
%                       row will be considered a time series.
% 
% eventInds -           A vector containing the indices (range of 1 to N)
%                       that indicate where events happen. Can be empty.
%                       Usually will normalize the first
%                       index to 0% and the last index to 100%, but you 
%                       can overwrite that here by putting 1 or N in 
%                       this input, and have the corresponding percentage
%                       in eventPercentages.
% 
% eventPercentages -    A vector containing the percentages corresponding
%                       to the indices in eventInds. Must be same length
%                       as eventInds. Usually will normalize the first 
%                       index to 0% and the last index to 100%, but you 
%                       can overwrite that here by putting 0 or 100 in 
%                       this input, and have the corresponding index in
%                       eventInds.
% 
% desiredPercentages -  Optional. A 1xK vector containing the percentages 
%                       at which you want the normalized data to be at. 
%                       If not given will use 0-100% in increments of 1%
%                       as default
% 
% Outputs:
% normalizedData -      MxK (default Mx100) matrix of the normalized time
%                       series data given by the data input.
% 
% David Xing 4/4/2019

narginchk(3,4)

if nargin==3
    %default output data in 1% increments
    desiredPercentages = 0:100;
end

% check that eventInds and eventPercentages are the same size
if length(eventInds)~=length(eventPercentages)
    error('eventInds and eventPercentages must have the same number of inputs!')
end

% if the first or last index isn't manually defined in eventInds, set them
% to be 0% and 100%
if ~any(eventInds==1) && ~any(eventPercentages==0)
    eventInds(end+1) = 1;
    eventPercentages(end+1) = 0;
end
if ~any(eventInds==size(data,2)) && ~any(eventPercentages==100)
    eventInds(end+1) = size(data,2);
    eventPercentages(end+1) = 100;
end

% get the "time points" (the indices) corresponding to the desired output
% percentages
desiredPercentageAsInds = interp1(eventPercentages, eventInds, desiredPercentages);

% now get the data values at those time points corresponding to the desired
% percentages
normalizedData = interp1(1:size(data,2), data', desiredPercentageAsInds)';


% 
