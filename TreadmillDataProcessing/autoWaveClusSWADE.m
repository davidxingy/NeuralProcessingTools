function autoWaveClusSWADE(contDataFile, prevSpikesDataFile, savePlotFile)
% autoWaveClusSWADE(contDataFile, prevSpikesDataFile, savePlotFile)
% 
% Function to run spike detection and sorting automatically. The purpose 
% of this function is to run the template matching automatically, assuming
% that the templates have already been generated manually, and the new data
% that is being run on is from the same recording session as the manually
% templated data. (Ex, manual wave_clus/Carlos SWADE sorting on block 1 of
% a session, then running this automatic sorting code on the remaining
% blocks).
% Spike detection will be done with wave_clus's Get_spikes() function,
% while the sorting will be done using Carlos's SWADE algorithm. 
% 
% Inputs:
% contDataFile - full path to the data file containing the raw neural
%                recordings to run the spike sorting on
% 
% prevSpikesDataFile - full path to the output file from the manual
%                      wave_clus session (usually in the form of 
%                      'times_....mat' .The templates will be calculated 
%                      from the sorted spikes contained in this file
% 
% savePlotFile - full path to the file to save the verification plots to
% 
% The function will save a 'times_...' file in the same vein as the output
% files from wave_clus, containing the 'cluster_class', 'spikes', and 'par'
% variables as if these spikes were manually sorted in the wave_clus gui
% (however, note that there won't be any 'forced', 'gui_status', 'inspk',
% 'ipermut', or 'Temp' variables).
% Also note that the plots won't be the same as the wave_clus gui output
% plot (which is just a screenshot of the gui), but something similar,
% containing the waveforms of each neuron, and the ISI histogram (log)
% 
% David Xing 2/13/19

% first check that the input files exist
if ~exist(contDataFile, 'file')
    error('Input continuous data file doesn''t exist!')
end

if ~exist(prevSpikesDataFile,'file')
    error('Previous spike sorting file doesn''t exist!')
end

% next, load in continuous data file
contFileVars=load(contDataFile);
% has to contain the continuous data variable 'data'
if ~isfield(contFileVars, 'data')
    error('Continuous data file doesn''t have the ''data'' variable!');
end

% next, load in the manual wave_clus output file
waveclusVars=load(prevSpikesDataFile);
% has to contain the 'cluster_class', 'spikes', and 'par' field
if ~isfield(waveclusVars, 'cluster_class')
    error('Sorted spikes file doesn''t have the ''cluster_class'' variable!');
end
if ~isfield(waveclusVars, 'spikes')
    error('Sorted spikes file doesn''t have the ''spikes'' variable!');
end
if ~isfield(waveclusVars, 'par')
    error('Sorted spikes file doesn''t have the ''par'' variable!');
end

% make plots dir if it doesn't already exist
[plotDir, plotFile, plotExt]=fileparts(savePlotFile);
if ~exist(plotDir,'dir')
    if ~mkdir(plotDir)
        error('Unable to access or make output plot directory!');
    end
end

% Now, run the spike detection code
try
    Get_spikes(contDataFile,'par',waveclusVars.par);
catch ME
    fprintf('Wasn''t able to run spike detection!');
    rethrow(ME);
end

% the spike detection code should have saved (in the current working
% folder) another mat file with the same name with '_spikes' appended
[~, contFileName] = fileparts(contDataFile);
if ~exist([contFileName '_spikes.mat'], 'file')
    error('Couldn''t get output file from the spike detection!')
end

% next load in the spike data
detectedSpikesVars=load([contFileName '_spikes']);

% make sure that all the proper variables are there
if ~isfield(detectedSpikesVars,'spikes')
    error('Detected spikes file doesn''t have the ''spikes'' variable!')
elseif ~isfield(detectedSpikesVars,'threshold')
    error('Detected spikes file doesn''t have the ''threshold'' variable!')
elseif ~isfield(detectedSpikesVars,'index')
    error('Detected spikes file doesn''t have the ''index'' variable!')
elseif ~isfield(detectedSpikesVars,'psegment')
    error('Detected spikes file doesn''t have the ''psegment'' variable!')
end

%Now, run Carlos's template matching algorithm (SWADE)

%calculate templates from the previous wave_clus session sorted neurons
sortedClusterNums=sort(unique(waveclusVars.cluster_class(:,1)));
sortedClusterNums(sortedClusterNums==0)=[];
% go through each sorted neuron
for iTemplate=1:length(sortedClusterNums)
    templates(:,iTemplate)=nanmean(waveclusVars.spikes(waveclusVars.cluster_class(:,1)==...
        sortedClusterNums(iTemplate),:));
    templateNoiseValues{iTemplate}=waveclusVars.spikes(waveclusVars.cluster_class(:,1)==...
        sortedClusterNums(iTemplate),1:5);
    placeHolder{iTemplate}=[]; %just needs to have length == number of templates, not actually used
end

% calculated noise level for thresholding
if ~isnan(waveclusVars.par.SWADE_threshold_STDMultiplier)
    allNoiseValues=cat(1,templateNoiseValues{:});
    noiseSTD=nanstd(allNoiseValues(:));
    thresh=waveclusVars.par.SWADE_threshold_STDMultiplier*noiseSTD;
else
    thresh=par.SWADE_threshold;
end

% run the template matching
try
[sorted_timestamps, noise_waves, rec_waves, noise_amp, original_input_indices]= SWADE_XI(...
    detectedSpikesVars.spikes', detectedSpikesVars.index, [], placeHolder, [], [], templates,...
    [], thresh, waveclusVars.par.SWADE_average_threshold, waveclusVars.par.SWADE_overlap_threshold,...
    waveclusVars.par.SWADE_template_threshold, 0);
catch ME
    fprintf('Wasn''t able to run spike forcing!');
    rethrow(ME);
end
    
%sanity check, all original input indices have been assigned to
%a class (or noise)
allInds=sort(unique(cat(1,original_input_indices{:},noise_waves')));
if(length(allInds)~=length(detectedSpikesVars.index) || ...
        any(allInds'~=(1:length(detectedSpikesVars.index))))
    error('Not all input waves have been assigned to a neuron or noise! Something wrong with code')
end

%sanity check, the noise waveforms shouldn't appear in any of
%the classified spikes
if (any(cellfun(@(x) any(ismember(x,noise_waves)),original_input_indices)))
    error('Some noise waveforms have been assigned to a neuron somehow! Something wrong with code');
end

% now, we should make the verification plots
%make plot
plot_h=figure('Units','normalized','OuterPosition',[0, 0, 1, 1],'Visible','off');
nClasses=length(sorted_timestamps)+1;

colormap=hsv(nClasses-1);
colormap(end+1,:)=[0.2 0.2 0.2]; %rejected spikes are black

% top row should be the raw neural time series
subplot(3, nClasses, 1:nClasses)
% just plot 1 second (30k points)
voltageData=contFileVars.data(1:30000)-nanmean(contFileVars.data(1:30000));
plot((1:30000)/30,voltageData,'color',[0.1 0.1 0.1])
hold on

% make pretty
box off
set(gca,'FontSize',14)
set(gca,'TickDir','out')
xlabel('Time (ms)')
ylabel('Voltage (uV)')

% get y axis limits for waveform plots (the largest amplitude of non-noise
% waveforms)
ylimMax=max(max(detectedSpikesVars.spikes(...
    setdiff(1:size(detectedSpikesVars.spikes,1),noise_waves),:)))*1.1;
ylimMin=min(min(detectedSpikesVars.spikes(...
    setdiff(1:size(detectedSpikesVars.spikes,1),noise_waves),:)))*1.1;

% now go through all the sorted neurons, and plot spikes and histograms
for iNeuron=1:nClasses
    %first the waveforms
    subplot(3, nClasses, nClasses+iNeuron)
    %get a random subset of 2000 neurons
    if iNeuron~=nClasses
        sampleInds=randperm(size(rec_waves{iNeuron},2));
        sampleInds=sampleInds(1:min([2000, size(rec_waves{iNeuron},2)]));
        allWaveforms=rec_waves{iNeuron};
    else
        %for rejected waves, plot up to 5000 of them
        sampleInds=randperm(length(noise_waves));
        sampleInds=sampleInds(1:min([5000, length(noise_waves)]));
        allWaveforms=detectedSpikesVars.spikes(noise_waves,:)';
    end
    waveformsToPlot=allWaveforms(:,sampleInds);
    plot(waveformsToPlot,'Color',colormap(iNeuron,:));
    
    %plot the means and std
    meanWaveform=nanmean(allWaveforms');
    stdWaveform=nanstd(allWaveforms');
    hold on
    plot(meanWaveform,'k','LineWidth',2)
    plot(meanWaveform+stdWaveform,'k--')
    plot(meanWaveform-stdWaveform,'k--')
    
    %display the SNR 
    noiseLevel=nanstd(reshape(allWaveforms(1:5,:),1,numel(allWaveforms(1:5,:))))*3;
    signalLevel=max(meanWaveform)-min(meanWaveform);

    if iNeuron==1
        %display SNR definition (just once)
        title({sprintf('SNR = %.2f', signalLevel/noiseLevel),'(average waveform amplitude / 3 noise stds)'});
    else
        title(sprintf('SNR = %.2f', signalLevel/noiseLevel));
    end
    
    %make pretty
    box off
    set(gca,'FontSize',14)
    set(gca,'TickDir','out')
    xlabel('Sample')
    ylabel('Voltage (uV)')
    ylim([ylimMin ylimMax])
    xlim([0 60])
    
    %next, do histograms
    subplot(3, nClasses, 2*nClasses+iNeuron)
    %get ISIs
    if iNeuron~=nClasses
        neuronTimestamps=sorted_timestamps{iNeuron};
    else
        neuronTimestamps=detectedSpikesVars.index(noise_waves);
    end
    ISIs=diff(neuronTimestamps);
    
    %plot histogram
    histogram(ISIs,0:1:100)
    set(gca,'XScale','log')
    
    %show how many spikes there are and how many are < 2ms apart
    title({sprintf('%u total spikes',length(neuronTimestamps)),...
        sprintf('%u ISIs < 2ms',sum(ISIs<2))})
    
    %make pretty
    box off
    set(gca,'FontSize',14)
    set(gca,'TickDir','out')
    set(gca,'XTick',[1 2 3 4 5 10 50 100])
    xlabel('ISI (ms)')
    ylabel('Count')
    
    %finally add a marker to the time series plot
    subplot(3, nClasses, 1:nClasses)
    %since we're only plotting first 1s of data
    timestampInds=neuronTimestamps(neuronTimestamps<1000);
    plot(timestampInds,voltageData(round(timestampInds*30)),'*','color',colormap(iNeuron,:))
    
end

%title
[origFilePath, origFileName]=fileparts(contDataFile);
t=suptitle(sprintf('%s auto sorting',origFileName));
set(plot_h,'Visible','off');
set(t,'Interpreter','none');

%saveplot
saveas(plot_h,savePlotFile);
close(plot_h)

% finally, save the output file with the sorted spikes
% concatenate all spikes
spikes=[rec_waves{:} detectedSpikesVars.spikes(noise_waves,:)']';
 
% assign neuron numbers
classNums=cell(1,nClasses);
for iNeuron=1:nClasses-1
    classNums{iNeuron}=repmat(iNeuron,1,size(rec_waves{iNeuron},2));
end
% noise waveforms are assigned 0
classNums{end}=zeros(1,length(noise_waves));

cluster_class(:,1)=[classNums{:}]';

% finally save the timestamps
cluster_class(:,2)=[sorted_timestamps{:} detectedSpikesVars.index(noise_waves)]';

% parameters are just the parameters from the original wave_clus session
par=waveclusVars.par;

% save
save(fullfile(origFilePath,sprintf('times_%s',origFileName)),'spikes','cluster_class','par')


%
