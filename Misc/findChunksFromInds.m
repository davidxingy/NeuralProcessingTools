function [chunkStarts, chunkEnds] = findChunksFromInds(indices,varargin)
% [chunkStarts, chunkEnds] = findChunksFromInds(indices,[increaseAmount])
% 
% simple function to segment time series which is made of concatenated
% chunks of continuous data back into the individual chunks
% 
% Inputs:
% indices -         A 1xN vector containing the indices of the chunks
% 
% increaseAmount -  Optional. The amount the indices are increasing by for
%                   each time step. Default 1
% 
% Outputs:
% chunkStarts -     The indices of the start of each chunk
% 
% chunkEnds -       The indices of the end of each chunk
% 
% David Xing 1/16/2023

if length(varargin)>1
    if isempty(varargin)
        increaseAmount = 1;
    else
        increaseAmount = varargin{1};
    end
else
    increaseAmount = 1;
end

chunkBoundaries = find(diff(indices) > increaseAmount);
chunkStarts = [1 chunkBoundaries+1];
chunkEnds = [chunkBoundaries length(indices)];
