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
% To do so, look at what B is, and use the following rules:
% Current edge is Rising A:
%   B is low    -   A is leading
%   B is high   -   A is lagging
% 
% Current edge is Falling A:
%   B is high   -   A is leading
%   B is low    -   A is lagging
% 
% Current edge is Rising B:
%   A is high   -   A is leading
%   A is low    -   A is lagging
% 
% Current edge is Falling B:
%   A is low    -   A is leading
%   A is high   -   A is lagging


% Assign each edge to either a leading or lagging edge based on the above
% rules:
leadingOrLagging{1}=Bvoltage(ARisingEdges)<thresh; %Rising A's leading if B is low
leadingOrLagging{2}=Bvoltage(AFallingEdges)>=thresh; %Falling A's are the opposite
leadingOrLagging{3}=Avoltage(BRisingEdges)>=thresh; %Continue following rules for rising B's
leadingOrLagging{4}=Avoltage(BFallingEdges)<thresh; %and falling B's

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
        countChanges(ARisingEdges(logical(leadingOrLagging{1})))=1*dontFlipSign;
        countChanges(ARisingEdges(logical(~leadingOrLagging{1})))=-1*dontFlipSign;
        countChanges(AFallingEdges(logical(leadingOrLagging{2})))=1*dontFlipSign;
        countChanges(AFallingEdges(logical(~leadingOrLagging{2})))=-1*dontFlipSign;
        countChanges(BRisingEdges(logical(leadingOrLagging{3})))=1*dontFlipSign;
        countChanges(BRisingEdges(logical(~leadingOrLagging{3})))=-1*dontFlipSign;
        countChanges(BFallingEdges(logical(leadingOrLagging{4})))=1*dontFlipSign;
        countChanges(BFallingEdges(logical(~leadingOrLagging{4})))=-1*dontFlipSign;
end

% use cumsum to add up the counts
count=cumsum(countChanges);


% 
