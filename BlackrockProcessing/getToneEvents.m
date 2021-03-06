function toneEvents=getToneEvents(micSignal,toneTypeEvents,varargin)
% toneEvents=getToneEvents(micSignal,toneTypeEvents,[thresholds], [plotsDirectory])
% This function extracts the time that the sounds are played in the
% treadmill based on a microphone recording. The sound type is determined
% based on the nearest NEV event received by the NSP from the control
% program.
% 
% Inputs:
% micSignal - vector containing the microphone recording from the treadmill
% 
% toneTypeEvents - Nx2 cell array. N is the total number of tones received
%                  and the first column is the name of the tone, while the
%                  second column is the timestamp of the tone in the
%                  micSignal vector
% 
% thresholds - Optional. Struct, with the tone type as field names. If you 
%              want to set a manual thesholds for finding the tones, use 
%              this input. If no input, then the threshold will be found 
%              automatically. The field name has to match the names of the
%              tones in toneTypeEvents
% 
% plotsDirectory - Optional. String containing the directory to save plots
%                  for verification that the events were correctly
%                  extracted. If not given, then plots won't be saved
% 
% Outputs:
% toneEvents - Struct. Each field is a name of the type of tone given in
%              toneTypeEvents/thresholds, and consists of a vector of
%              timestamps of when tones of that type was detected (in
%              # of samples) from micSignal
% 
% David Xing 1/12/2019

% parse inputs
narginchk(2,4);

% threshold
if nargin>2
    if(isempty(varargin{1}))
        thresholds = [];
    else
        thresholds = varargin{1};
    end
else
    thresholds = [];
end
% plots directory
if nargin>3
    if(isempty(varargin{2}))
        plotsDirectory = [];
    else
        plotsDirectory = varargin{2};
        if ~exist(plotsDirectory,'dir')
            try    
                dirSuccess = mkdir(plotsDirectory);
            catch 
                error('Unable to access or make plots directory!');
            end
        end
    end
else
    plotsDirectory = [];
end

% get the list of tone names
toneNames=unique(toneTypeEvents(:,1));
autoThreshTones={};

% calc thresholds if not manually given
if isempty(thresholds)
    %all the thresholds needs to be auto calculated
    autoThreshTones=toneNames;
else
    %check that the threshold names match the tone names in toneTypeEvents
    toneThreshNames=fieldnames(thresholds);
    for iName=1:length(toneNames)
        if ~any(strcmpi(toneThreshNames,toneNames{iName}))
            error('Field names in manual thresholds input must match with the names in toneTypeEvents!')
        end
        
        %if the threshold is empty for any of the tone names, calculated the
        %threshold automatically
        if isempty(thresholds.(toneNames{iName}))
            autoThreshTones(end+1)=toneNames(iName);
        end
    end
end

% go through each type of tone that needs auto-thresholding
for iTone=1:length(autoThreshTones)
    %get the max value of the 100 ms before the event for all the
    %events
    toneTimes=[toneTypeEvents{strcmp(toneTypeEvents, autoThreshTones{iTone}),2}];
    toneMaxValues=[];
    for iEvent=1:length(toneTimes)
        toneMaxValues(iEvent)=max(micSignal(max(1,toneTimes(iEvent)-3000):toneTimes(iEvent)));
    end
    
    %use 1.02x the 75th percentile of all those values as the threshold
    %since some of them may have high noise before the tone
    thresholds.(autoThreshTones{iTone})=prctile(toneMaxValues,75)*1.02;
end

% go through each type of tone
toneEvents=struct();
for iTone=1:length(toneNames)
    % find threshold crossings
    allToneStarts=find(micSignal(2:end)>thresholds.(toneNames{iTone}) & ...
        micSignal(1:end-1)<thresholds.(toneNames{iTone}))+1;
    
    % get the first threshold crossing after each event, and use that as
    % the time the tone is played
    toneTimes=[toneTypeEvents{strcmp(toneTypeEvents, toneNames{iTone}),2}];
    for iEvent=1:length(toneTimes)
        thisToneStart=min(allToneStarts(allToneStarts>toneTimes(iEvent)));
        
        %check that the time the tone is played is at most 500ms after the
        %event was sent (or that there actually was a tone played)
        if isempty(thisToneStart)
            warning(sprintf('No tone found for the %u event of the %s tone!',iEvent,toneNames{iTone}));
            continue
        elseif  thisToneStart>toneTimes(iEvent)+15000
            warning(sprintf('Event %u of the %s tone is %g.2 ms before the detected tone!',...
                iEvent,toneNames{iTone},thisToneStart-toneTimes(iEvent)));
        end
        
        %save the time
        toneEvents.(toneNames{iTone})(iEvent)=thisToneStart;
        
        %save plots if requested
        if ~isempty(plotsDirectory)
            %make plot
            plot_h=figure('Units','pixels','OuterPosition',[100, 100, 1000, 700],'Visible','off');
            plotSig=micSignal(max(1,toneTimes(iEvent)-10000):min(thisToneStart+30000,length(micSignal)));
            ax(1)=plot((1:length(plotSig))/30,plotSig);
            hold on
            ax(2)=line([10001 10001]/30,get(gca,'YLim'),'color','r');
            ax(3)=plot((length(plotSig)-30000)/30,plotSig(length(plotSig)-30000),'*');
            ylabel('Mic Signal')
            xlabel('Time (ms)')
            title(sprintf('%s Tone, %u Event',toneNames{iTone},iEvent))
            legend(ax,'Mic Signal','NEV Event Time','Detected Tone Time')
            set(gca,'FontSize',14);
            box off
            
            %save
            saveas(plot_h,fullfile(plotsDirectory,sprintf('%sTone_%u',toneNames{iTone},iEvent)),'png');
            close(plot_h)
            
        end
        
    end

end


% 
