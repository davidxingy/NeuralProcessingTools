function toneEvents=getToneEvents(micSignal,toneTypeEvents,varargin)
% toneEvents=getToneEvents(micSignal,toneTypeEvents,[threshold])
% This function extracts the time that the sounds are played in the
% treadmill based on a microphone recording. The sound type is determined
% based on the nearest NEV event received by the NSP from the control
% program.
% 
% Inputs:
% 
% micSignal - vector containing the voltage recording from the treadmill
%             mic. First second of the signal shouldn't contain any tones
%             as I use the first 30000 samples to estimate the noise level.
% 
% toneTypeEvents - Nx2 cell array. N is the total number of tones received
%                  and the first column is the name of the tone, while the
%                  second column is the timestamp of the tone in the
%                  micSignal vector
% 
% threshold - Optional. If you want to set a manual threshould for finding
%             the tones, use this input. If no input, then the threshould 
%             will be found automatically.
% 
% David Xing 1/12/2019

% parse inputs
threshold=[];
if nargin==3
    threshold=varargin{1};
end

% calc threshold if not manually given
% use first second of data (assuming 30ks/s) to get noise level
if isempty(threshold)
    maxNoiseValue=max(abs(micSignal(1:30000)));
    threshold=1.2*maxNoiseValue; %use 1.2x the max as threshold
end

% find threshold crossings
allToneStarts=micSignal(2:end)>threshold & micSignal(1:end-1)<threshold;

% start of tones should be at least 200ms apart
toneStarts=allToneStarts;
for iTone=2:length(allToneStarts)
    if (allToneStarts(iTone)-allToneStarts(iTone-1))<6000
        toneStarts(iTone)=NaN;
    end
end

toneStarts(iTone)=




% 
