function [segStartInds, segLengths] = findConstants(data, varargin)
% [segStartInds, segLengths] = findConstants(data, [tolerance])
% 
% Find segments of data where the value is constant (within some tolerance)
% 
% Inputs:
% data -        A 1xN vector containing the data in which you want to find
%               the constant value segments. Should be a double
% 
% tolerance -   Optional. The tolerance for considering something constant.
%               Will be 0 as default (the subsequent values have to exactly
%               equal the previos value to be considered a constant
%               segment).
% 
% Outputs:
% segStartInds -    The indices of the start of each constant segment found
%                   in data.
% 
% segStartLengths - The lengths of each constant segment found in data.
% 
% David Xing 4/11/2019

% check inputs and set defaults
narginchk(1, 2)

if ~strcmpi(class(data), 'double')
    error('data input must be a double!');
end

if length(varargin)>1
    if isempty(varargin)
        tolerance = 0;
    else
        tolerance = varargin{1};
    end
else
    tolerance = 0;
end

% find differences in data
dataDiffs = [NaN diff(data)];

% snap to tolerance
dataDiffs(abs(dataDiffs)<tolerance) = 0;

% constants are where diffs are zeros
zeroInds = find(dataDiffs==0);

% constant segments are where the zero indices are incrementing by 1
% so the start of each segment is where it isn't incremented by 1
% also, always include the first zero ind
segStartZeroInds = [1 find(diff(zeroInds)~=1)+1];

% segment lengths are the differences in the start inds
% also, need to find the length of the last segment (from the last start
% ind to the length of the zero inds)
segLengths = [diff(segStartZeroInds) length(zeroInds)-segStartZeroInds(end)+1]+1;

% the actual start segments indices are the one before the segStartZeroInds
segStartInds = zeroInds(segStartZeroInds)-1;



% 
