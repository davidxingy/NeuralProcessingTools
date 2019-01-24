function histArray=addHistory(inputArray, nHist, varargin)
% function histArray=addHistory(inputArray, nHist, [dim], [shiftAmount], [keepOrder], [onlyHist])
% 
% This function shifts an input array along a specified dimension by some
% amount each time and then concatenated onto the array. The empty space 
% is padded by zeros. This is repeated nHist times. 
% 
% Inputs:
% inputArray -  a N x M 2D numerical array to add the history to
% 
% nHist -       the amount of shifting (history) to add. This doesn't
%               include the origianal, unshifted array, e.g. nHist=0
%               will return the original array, not an empty array
% 
% dim -         Optional. The dimension to shift the array along. The
%               shifted array will be concatenated along the other
%               dimension. Default 2.
% 
% shiftAmount - Optional. The number of elements to shift by for each
%               history. Positive is shift to the right. Default 1.
%
% keepOrder -   Optional. Setting this to true will add the history of each
%               row right underneath the original row/column, rather than
%               all the histories of the rows/columns as a block underneath
%               the block of original rows/columns. Default true.
% 
% onlyHist -    Optional. If this isn't true, the output will also contain
%               the original inputted data as a 0-history section. If you 
%               only want the time shifted values, set to true. Default
%               false.
% 
% Outputs:
% histArray -   a N*(nHist+1) x M (or N x M*(nHist+1) if dim==1) 2D
%               numerical array with the shifted history of inputArray.
% 
% Example:
% histArray=addHistory([1 2 3; 4 5 6], 2)
% will result in
% histArray = 
% [1    2   3
%  0    1   2
%  0    0   1
%  4    5   6
%  0    4   5
%  0    0   4]
% 
% whereas
% histArray=addHistory([1 2 3; 4 5 6], 2, [], -1, false)
% will result in 
% histArray = 
% [1    2   3
%  4    5   6
%  2    3   0
%  5    6   0
%  3    0   0
%  6    0   0]
% 
% 
% David Xing, last updated 8/9/2018

% parse inputs, set defaults
narginchk(2,6);
% dimension
if nargin>2
    if(isempty(varargin{1}))
        dim = 2;
    else
        dim = varargin{1};
    end
else
    dim = 2;
end
% shift amount
if nargin>3
    if(isempty(varargin{2}))
        shiftAmount = 1;
    else
        shiftAmount = varargin{2};
    end
else
    shiftAmount=1;
end
% keep row ordering
if nargin>4
    if(isempty(varargin{3}))
        keepOrder = true;
    else
        keepOrder = varargin{3};
    end
else
    keepOrder=true;
end
% keep unshifted (0-history) values
if nargin>5
    if(isempty(varargin{4}))
        onlyHist = false;
    else
        onlyHist = varargin{4};
    end
else
    onlyHist=false;
end


% calculate the final array size
if dim==1
    dim1Size=size(inputArray,1);
    dim2Size=(nHist+1)*size(inputArray,2);
    concatAmount=size(inputArray,2);
else
    dim1Size=(nHist+1)*size(inputArray,1);
    dim2Size=size(inputArray,2);
    concatAmount=size(inputArray,1);
end

%initate with zeroes
histArray=zeros(dim1Size,dim2Size);

% Add the first unshifted original array
if dim==1
    histArray(:,1:concatAmount)=inputArray;
else
    histArray(1:concatAmount,:)=inputArray;
end

% now go through and add each history
for iHist=1:nHist
    %shift to the right by number of timepoints (positive)
    trialShift=circshift(inputArray,iHist*shiftAmount,dim);
    
    %pad with 0's and append to complete history matrix
    if shiftAmount<0
        %add the zeros on the end
        if dim==1
            zeroInds=max(dim1Size-iHist*abs(shiftAmount)+1,1):dim1Size;
        else
            zeroInds=max(dim2Size-iHist*abs(shiftAmount)+1,1):dim2Size;
        end
    else
        if dim==1
            zeroInds=1:min(iHist*shiftAmount,dim1Size);
        else
            zeroInds=1:min(iHist*shiftAmount,dim2Size);
        end
    end
    
    if dim==1
        trialShift(zeroInds,:)=0;
        histArray(:,iHist*concatAmount+1:(iHist+1)*concatAmount)=trialShift;
    else
        trialShift(:,zeroInds)=0;
        histArray(iHist*concatAmount+1:(iHist+1)*concatAmount,:)=trialShift;
    end
    
end

% If we don't want the unshifted values, remove them
if(onlyHist)
    if dim==1
        histArray(:,1:concatAmount)=[];
    else
        histArray(1:concatAmount,:)=[];
    end
end

% for preserving order, get the proper indices where the history should
% be placed
histInd=reshape(1:concatAmount*(nHist+~onlyHist),concatAmount,nHist+~onlyHist);
histInd=histInd';
histInd=histInd(:);
    
% now finally reorder the rows to preserve the order if desired
if (keepOrder)
    if dim==1
        histArray=histArray(:,histInd);
    else
        histArray=histArray(histInd,:);
    end
end


% 
