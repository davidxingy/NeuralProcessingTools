function blocks=divideBlocks(vectorToBeSplit, nBlocks)
% function blocks=divideBlocks(vectorToBeSplit, nBlocks)
% 
% This function attempts to evenly divide a vector into blocks of
% smaller vectors. If the vector is evenly divisible by nBlocks, then
% all the blocks are the same size. Otherwise it will divide the blocks as 
% evenly as possible (no blocks will ever differ in size by more than 1).
% The order of the elements in the original vector will be preserved.
% 
% Inputs:
% vectorToBeSplit - the 1xN (or Nx1) vector that will be split into blocks
% 
% nBlocks -         the number of blocks to split inputTrials into.
% 
% Outputs:
% blocks -          a 1xnBlocks cell of smaller vectors whos union is inputTrials
% 
% Example:
% blocks=divideBlocks(1:20, 6)
% will result in 
% blocks = {[1 2 3], [4 5 6], [7 8 9], [10 11 12], [13 14 15 16], [17 18 19 20]}
% 
% 
% David Xing, last updated 7/9/2018


nValues=length(vectorToBeSplit);

% first check that there are at least as many elements as blocks
if nBlocks>nValues
    error('Can''t divide vector into more vectors than it has elements')
end

% See how many elements will be in each block
if floor(nValues/nBlocks)==nValues/nBlocks
    %we have perfect division of trials into the nfolds
    blockSizes{1}=nValues/nBlocks; %just one size needed
    numBlockWithSize{1}=nBlocks; %all blocks are that size
else
    %number of trials is not a multiple of number of
    %folds, have to split nBlocks into two, each with
    %different number of elements in the blocks
    blockSizes{1}=floor(nValues/nBlocks);
    blockSizes{2}=ceil(nValues/nBlocks);
    
    numBlockWithSize{2}=nValues-nBlocks*blockSizes{1}; %the remainder will be split amongst the remaining blocks
    numBlockWithSize{1}=nBlocks-numBlockWithSize{2};
end

% Now do the splitting into blocks
indCounter=1;
for iBlock=1:nBlocks
    
    %size of this block
    if iBlock<=numBlockWithSize{1}
        theBlockSize=blockSizes{1};
    else
        theBlockSize=blockSizes{2};
    end
    
    %move to output cell
    blocks{iBlock}=vectorToBeSplit(indCounter:indCounter+theBlockSize-1);
    indCounter=indCounter+theBlockSize;
end


% 
