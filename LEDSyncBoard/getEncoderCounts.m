function count=getEncoderCounts(Avoltage, Bvoltage, thresh, varargin)
% count=getEncoderCounts(Avoltage, Bvoltage, thresh, [mode], [positiveIsALeading])
% 
% This function takes in qudrature encoder A and B voltage time series 
% signals and converts it into a count, where the count starts at 0 and is 
% changed by 1 upon edge detection. Can do X1, X2 or X4 modes. Whether the 
% count is incremented or decremented by 1 depends on whether B is leading
% or lagging A, and can be set using the positiveIsALeading input.
% 
% E.g (with positiveIsALeading set to true and x4 mode):
% A:        _______|-------|______|-------|______|------|_____|-------...
% B:        ___|-------|_____________|-------|______|------|______|---...
% Count:    0000-1-1-2-2-3-3-4-4-4-3-3-2-2-1-1000111222233344455556666...
% 
% Inputs:
% Avoltage -            a 1 x N numerical array which represents the voltage 
%                       readings from the A line of the quadrature encoder
% 
% Bvoltage -            a 1 x N numerical array which represents the voltage 
%                       readings from the B line of the quadrature encoder
% 
% thresh -              The threshold value for finding pulses in Avoltage
%                       and Bvoltage
% 
% mode -                Optional. An interger, 1, 2, or 4 to select between
%                       X1, X2, or X4 encoding modes. X1 will only 
%                       increment on rising edges of A when A leading and 
%                       decrement on falling edges of A when A lagging. X2 
%                       will increment/decrement on all edges of A, and X4 
%                       will increment/decrement on all edges of A and B.
%                       Default set to X4.
% 
% positiveIsALeading -  Optional. If this is set to true, then whenever A 
%                       is leading in time, the count will be incremented 
%                       by 1, if set to false, the count will be 
%                       decremented by 1. Default set to true.
% 
% Outputs:
% count -               a 1 x N array with starting value 0, representing
%                       the count of the encoder
% 
% 
% David Xing, last updated 8/13/2018


% parse inputs, set defaults
narginchk(3,5);
% mode (default 4)
if nargin>3
    if(isempty(varargin{1}))
        mode = 4;
    else
        mode = varargin{1};
        if ~any(mode==[1 2 4])
            error('Mode can only be 1, 2 or 4!')
        end
    end
else
    mode = 4;
end
% Increment for A leading rather than decrement (default true)
if nargin>4
    if(isempty(varargin{2}))
        positiveIsALeading = 1;
    else
        positiveIsALeading = varargin{2};
    end
else
    positiveIsALeading=1;
end


% find edge times
ARisingEdges=find(Avoltage(2:end)>thresh & Avoltage(1:end-1)<thresh)+1;
AFallingEdges=find(Avoltage(2:end)<thresh & Avoltage(1:end-1)>thresh)+1;
BRisingEdges=find(Bvoltage(2:end)>thresh & Bvoltage(1:end-1)<thresh)+1;
BFallingEdges=find(Bvoltage(2:end)<thresh & Bvoltage(1:end-1)>thresh)+1;


% For all edge times, determine whether A is leading or lagging.
% To do so, look at the next edge, and use the following rules:
% Current edge is Rising A:
%   Next edge is Rising B -     A is leading
%   Next edge is Falling B -    A is lagging
%   Next edge is Faling A -     A is lagging
% 
% Current edge is Falling A:
%   Next edge is Rising B -     A is lagging
%   Next edge is Falling B -    A is leading
%   Next edge is Rising A -     A is lagging
% 
% Current edge is Rising B:
%   Next edge is Rising A -     A is lagging
%   Next edge is Falling A -    A is leading
%   Next edge is Falling B -    A is leading
% 
% Current edge is Falling B:
%   Next edge is Rising A -     A is leading
%   Next edge is Falling A -    A is lagging
%   Next edge is Rising B -     A is leading


% put all of the edges into a cell so we can loop through all of them
edgesCell{1}=ARisingEdges;
edgesCellLabels{1}=ones(1,length(ARisingEdges)); %label for A rising is 1
edgesCell{2}=AFallingEdges;
edgesCellLabels{2}=2*ones(1,length(AFallingEdges)); %label for A falling is 2
edgesCell{3}=BRisingEdges;
edgesCellLabels{3}=3*ones(1,length(BRisingEdges)); %label for B rising is 3
edgesCell{4}=BFallingEdges;
edgesCellLabels{4}=4*ones(1,length(BFallingEdges)); %label for B falling is 4

% what the next edge must be to count A as leading, based on the table
% above
leadingConditions{1}=3;
leadingConditions{2}=4;
leadingConditions{3}=[2; 4];
leadingConditions{4}=[1; 3];

% loop through and assign leading vs lagging for all edge types
parfor iEdgeType=1:length(edgesCell)
    
    %get list of all edges that is not the current edge type
    allEdges=cat(2,edgesCell{setdiff(1:length(edgesCell),iEdgeType)})';
    allEdgesLabels=cat(2,edgesCellLabels{setdiff(1:length(edgesCell),iEdgeType)})';
    
%   Non-vectorized code just loop through all edges one by one~~~~~~~~~~~~~
% 
%     %go through for all edges of this edge type
%     for iEdge=1:length(edgesCell{iEdgeType})
%         %see if the next edge matches the condition for A to be leading
%         if any(allEdgesLabels(min(find(allEdges>edgesCell{iEdgeType}(iEdge))))==leadingConditions{iEdgeType})
%             leadingOrLagging{iEdgeType}(iEdge)=1; %use 1 for A leading
%         else
%             leadingOrLagging{iEdgeType}(iEdge)=0; %use 0 for A lagging
%         end
%     end
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    %divide into blocks for vectorization to make it faster (about 2x)
    blockSize=100; %increasing it by more doesn't seem to make it any faster
    
    %we want to get nextEdgeType, a vector listing the type of the next
    %edge for all the edges in the current edge type.
    nextEdgeType=zeros(1,length(edgesCell{iEdgeType}));
    for iBlock=1:ceil(length(edgesCell{iEdgeType})/blockSize)
        %get the block
        edgeBlock=edgesCell{iEdgeType}(1+(iBlock-1)*blockSize:min(iBlock*blockSize,length(edgesCell{iEdgeType})));
        %make matrices
        mat_AllEdges=repmat(int32(allEdges),1,size(edgeBlock,2));
        mat_blockEdges=repmat(int32(edgeBlock),size(allEdges,1),1);
        
        %find the difference between each of the edges in the block and all
        %of the edges of different types
        offset=mat_AllEdges-mat_blockEdges;
        
        %minimum of the (positive) difference is the closest next edge
        offset(offset<=0)=max(max(offset));
        [~,nextEdgeInd]=min(offset);
        
        %put the results of this block into nextEdgeType variable
        nextEdgeType(1+(iBlock-1)*blockSize:min(iBlock*blockSize,length(edgesCell{iEdgeType})))=allEdgesLabels(nextEdgeInd);
    end
    
    %now that we have the type of the next edge, use it to determine if the
    %current edges is A leading or A lagging
    leadingOrLagging{iEdgeType}=zeros(1,length(edgesCell{iEdgeType})); %use 0 for A lagging
    for iCondition=1:length(leadingConditions{iEdgeType}) %for any of the next edge types (if there are multiple)
        leadingOrLagging{iEdgeType}... %if it matches the condition, set to 1 (use 1 for A leading)
            (nextEdgeType==repmat(leadingConditions{iEdgeType}(iCondition),1,length(nextEdgeType)))=1;
    end
    
    leadingOrLagging{iEdgeType}=logical(leadingOrLagging{iEdgeType});
end

% Next, put 1 or -1 at edge times (which edges depends on the mode, and the
% sign of the 1 depends on positiveIsALeading)
dontFlipSign=2*positiveIsALeading-1;
countChanges=zeros(1,length(Avoltage));
switch mode
    case 1
        %only increase on rising edges of A when A leading, 
        %and decrease on falling edges of A when A lagging
        countChanges(ARisingEdges(leadingOrLagging{1}))=1*dontFlipSign;
        countChanges(AFallingEdges(~leadingOrLagging{2}))=-1*dontFlipSign;
    
    case 2
        %increase and decrease on all edges of A
        countChanges(ARisingEdges(leadingOrLagging{1}))=1*dontFlipSign;
        countChanges(ARisingEdges(~leadingOrLagging{1}))=-1*dontFlipSign;
        countChanges(AFallingEdges(leadingOrLagging{2}))=1*dontFlipSign;
        countChanges(AFallingEdges(~leadingOrLagging{2}))=-1*dontFlipSign;
        
    case 4
        %increase and decrease on all edges of A and B
        countChanges(ARisingEdges(leadingOrLagging{1}))=1*dontFlipSign;
        countChanges(ARisingEdges(~leadingOrLagging{1}))=-1*dontFlipSign;
        countChanges(AFallingEdges(leadingOrLagging{2}))=1*dontFlipSign;
        countChanges(AFallingEdges(~leadingOrLagging{2}))=-1*dontFlipSign;
        countChanges(BRisingEdges(leadingOrLagging{3}))=1*dontFlipSign;
        countChanges(BRisingEdges(~leadingOrLagging{3}))=-1*dontFlipSign;
        countChanges(BFallingEdges(leadingOrLagging{4}))=1*dontFlipSign;
        countChanges(BFallingEdges(~leadingOrLagging{4}))=-1*dontFlipSign;
end

% use cumsum to add up the counts
count=cumsum(countChanges);


% 
