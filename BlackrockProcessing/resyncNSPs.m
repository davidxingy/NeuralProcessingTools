function [syncRemovalPoints, nResyncs]=resyncNSPs(triggerPoints, driftTolerance)
% [syncRemovalPoints, nResyncs]=resyncNSPs(triggerPoints, driftTolerance)
% 
% Function to re-synchronize two NSPs due to drifting of thier clocks. This
% method takes in the timestamps of sync points (points that are supposed
% to be simultaneously recorded on both NSPs, e.g. Simi triggers) from both
% NSPs and finds out which samples to remove to realign them when they
% drift over a specified amount.
% 
% Inputs:
% triggerPoints - 2xN array. Each row contains the time stamps (in samples)
%                 of each of the sync points for each NSPs.
% 
% driftTolerance - value which is the max amount of acceptable difference
%                  between the timestamps of a sync point (in samples).
% 
% Outputs:
% syncRemovalPoints - 2x1 cell array. Each cell contains an array of
%                     indicies which should be removed from the data of
%                     the NSP to realign the sync points. Each cell is for
%                     the corresponding row in the triggerPoints input.
% 
% nResyncs - the number of times the two NSPs desynced
%
% David Xing 1/31/2019

% check inputs
narginchk(2,2);
assert(size(triggerPoints,1)==2,'triggerPoints must be a 2xN array!');

syncRemovalPoints{1}=[];
syncRemovalPoints{2}=[];

%keep checking if there's any sync difference > the tolerance
nResyncs=0;
desyncFrame=find(abs(diff(triggerPoints))>driftTolerance,1);
while(~isempty(desyncFrame))
    
    %if there is see how much it's desynced by
    nDesyncSamples=abs(diff(triggerPoints(:,desyncFrame)));
    
    %and which NSP is faster
    if triggerPoints(1,desyncFrame)>triggerPoints(2,desyncFrame)
        removalNSP=1;
        unchangedTrigs=2;
    else
        removalNSP=2;
        unchangedTrigs=1;
    end
    
    %and flag the points in the faster NSP between the two triggers
    %(which there shouldn't be any if the triggers are exactly aligned)
    %to be removed
    syncRemovalPoints{removalNSP}=[syncRemovalPoints{removalNSP} ...
        (triggerPoints(unchangedTrigs,desyncFrame):triggerPoints(removalNSP,desyncFrame)-1)+...
        length(syncRemovalPoints{removalNSP})-...   %account for offset from
        length(syncRemovalPoints{unchangedTrigs})]; %previous resyncs
    
    %shift the trigger indices of the faster NSP so that it is
    %realigned with the other NSP at this location
    triggerPoints(removalNSP,desyncFrame:end)=...
        triggerPoints(removalNSP,desyncFrame:end)-nDesyncSamples;
    
    %add to counter, and then loop again to keep checking if there's
    %more desyncs
    nResyncs=nResyncs+1;
    desyncFrame=find(abs(diff(triggerPoints))>driftTolerance,1);
    
end


% 
