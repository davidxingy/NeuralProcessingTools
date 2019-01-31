function simiTriggers=extractSimiTriggers(trigSig, varargin)
% simiTriggers=extractSimiTriggers(trigSig, [maxSamplesBetweenPulses], [thresh], [plotFileName])
% 
% This function extracts the sample time indexs of the rising edges for a 
% Simi trigger recording. It'll only return the time point of edges between
% the start and stop of the video recording (i.e. there must be a pause
% before the first edge which indicates the start of the recording).
% 
% Inputs:
% trigSig - vector containing the recording of the Simi triggers.
% 
% maxSamplesBetweenPulses - Optional. Value which represents how many 
%                           samples you expect to be between two pulses 
%                           during continuous recording (beyond which, you 
%                           would consider to be a pause in the recording).
%                           Default 310 (appropriate for 100FPS video and 
%                           30kHz sampling rate).
% 
% thresh - Optional. Value to threshold the signal by to find rising edges.
%          Default 15000.
% 
% plotFileName - Optional. String containing the full path and file name to
%                save plots for verification that the triggers were 
%                correctly extracted. If not given, then plots won't be 
%                saved.
% 
% Outputs:
% simiTriggers - A 1xN array which contains the indices of the locations of
%                the detected trigger pulses. 
% 
% David Xing 1/31/2019

% parse inputs and set defaults
narginchk(1,4);

% threshold
if nargin>1
    if(isempty(varargin{1}))
        thresh = 1.5e4;
    else
        thresh = varargin{1};
    end
else
    thresh = 1.5e4;
end

% maxSamplesBetweenPulses
if nargin>2
    if(isempty(varargin{2}))
        maxSamplesBetweenPulses = 310;
    else
        maxSamplesBetweenPulses = varargin{2};
    end
else
    maxSamplesBetweenPulses = 310;
end

% plots file
if nargin>3
    if(isempty(varargin{3}))
        plotFileName = [];
    else
        plotFileName = varargin{3};
        
        %make directory if it doesn't already exist
        plotsDir =fileparts(plotFileName);
        if ~isempty(plotsDir) && ~exist(plotsDir,'dir')
            try
                dirSuccess = mkdir(plotsDir);
            catch
                error('Unable to access or make plots directory!');
            end
        end
    end
else
    plotFileName = [];
end



%find rising edges
threshCrossings=find(trigSig(2:end)>=thresh & trigSig(1:end-1)<thresh)+1;

%The camera may be getting triggers at the beginning or end of the
%recording (if I still had the live feed on), so I only want the
%triggers after and before the pause in the triggers, e.g.
% --|___|-----------------|___|---
% ?     ? actual frames  ?     ?
% |-Camera on but not recording-|
trigDiffs=diff([1 threshCrossings]);
trigPauses=find(trigDiffs(2:end)>maxSamplesBetweenPulses);

%can't have more than two pauses in the triggers, otherwise can't
%tell which block of triggers is the one that belongs to the video
%recording.
if length(trigPauses)>2
    error(['Multiple segments of threshold crossings found, ' ...
            'can''t distiguish which is the actual video recording!']);
end

if trigDiffs(1)<maxSamplesBetweenPulses
    %the recording started while the cameras were still on (though
    %recording video), find the pause where the cameras were turned
    %off; after which, the recording starts, e.g.:
    %---|____|-----------|___
    
    if isempty(trigPauses)
        %if there arn't any pauses, that means the video recording
        %already started when the NSP started recording, e.g.:
        %-------------|____ or
        %------------------
        error('Video recording started before NSP recording!');
    end
    
    %the video recording start is the first threshold crossing
    %after the first pause
    postStartTrig=trigPauses(1)+1;
else
    %the cameras were off when the NSP recording started, first
    %rising edge should be the first trigger of the video recording
    postStartTrig=1;
    
    %also, there can be at most 1 pause if the cameras were off at
    %the beginning of the NSP recording, e.g.:
    %___|----|____|----|__|-- is bad
    if length(trigPauses)>1
        error(['Multiple segments of threshold crossings found, ' ...
            'can''t distiguish which is the actual video recording!']);
    end
end

if postStartTrig~=1 && length(trigPauses)==2 ...
        || postStartTrig==1 && length(trigPauses)==1
    %there is an extra pause that happened after the video
    %recording stopped (and then the cameras were turned on again
    %before the NSP recording was stopped, e.g.:
    %____|---------|______|------ or
    %--|___|-----------|_____|---
    preEndTrig=trigPauses(end);
    
    %the camera recording cannot have stopped at the end of the NSP
    %recording, otherwise we have two blocks and it's unclear which
    %block is the video recording, e.g.:
    %___|---|____|---|____ or
    %-|_|---|____|---|____ is bad
    if length(trigSig)-threshCrossings(end)>maxSamplesBetweenPulses
        error('Multiple segments of threshold crossings found!');
    end
else
    %Cameras stayed off after the video recording was finished
    %until the NSP stopped recording, the last threshold crossing
    %is the last frame of the video, e.g.:
    %___|-------|____ or
    %--|__|-----|____
    preEndTrig=length(threshCrossings);
end

%now, only get the threshold crossings that correspond to the recorded
%video.
simiTriggers=threshCrossings(postStartTrig:preEndTrig);


%make plots for visual verification (histogram of trigger diffs)\
if ~isempty(plotFileName)
    nFrames=length(simiTriggers);
    plot_h=figure('Units','pixels','OuterPosition',...
        [100, 100, 1000, 700],'Visible','off');
    histogram(diff(simiTriggers))
    ylabel('# counts')
    xlabel('Time diff between threshold crossings (ms)')
    title(sprintf('%u total frames found',nFrames));
    set(gca,'FontSize',14);
    box off
    
    %save
    saveas(plot_h,plotFileName,'png');
    close(plot_h)
end


% 
