% function runSALT(baseDir, binFile)

% load(fullfile(baseDir,'neuronDataStruct.mat'))
% load(fullfile(baseDir,'syncSignal.mat'))
% load(fullfile(baseDir,'artifactInds.mat'))

clear; baseDir = 'X:\David\ArenaRecordings\DT030-082425-OptoTest-Recording1\Neuropixels';
binFile = 'DT030-082425-OptoTest-Recording1_g0_t0_LocalMedianSubtr_OptoArtRemoved.imec0.ap.bin';
load(fullfile(baseDir,'syncSignal.mat'))

load(fullfile(baseDir,'..','ProcessedData','neuronDataStruct.mat'))
load(fullfile(baseDir,'artifactInds.mat'))
% load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

readFid = fopen(fullfile(baseDir,binFile));
nChans = 385;
useTotalSpikes = false;

optoStarts = find(syncSignal(2:end)>40 & syncSignal(1:end-1)<=40);
optoStarts = optoStarts(optoStarts > 1.11e8);
% load(fullfile(baseDir,'ProcessedData','D045-022225-ArenaRecording_ProcessedEMG_MetaData.mat'),'laserOnsetInds')
% optoStarts = NeurEMGSync(laserOnsetInds,frameEMGSamples,frameNeuropixelSamples,'emg');
% optoStarts = optoOns(601:1200);
nBaselineSamples = 200;
baselinePeriod = 400;
pulsePower = [0.1 0.5 3 5];
pulseDuration = [11];
pulseColor = {'Red','Blue'};
nReps = 200;

if length(optoStarts) ~= length(pulsePower)*length(pulseDuration)*length(pulseColor)*nReps
    error('Discrepency in the number of pulses found')
end

% assign pulses to each power level/duration
% will use protocol: first sweep through power levels, then repeat for
% different laser colors
colorCutoffs = length(optoStarts)/length(pulseColor) * (1:length(pulseColor));
powerCutoffs = (length(optoStarts)/length(pulseColor)) / length(pulsePower) * (1:length(pulsePower));

% get sample times for baseline control
baselineOffsets = randperm(baselinePeriod-20)+10;
baselineOffsets = baselineOffsets(1:nBaselineSamples);

for iNeuron = 1:length(neuronDataStruct)

    waveforms{iNeuron,1} = [];
    waveforms{iNeuron,2} = [];
    for iPulse = 1:length(optoStarts)

        allPulseColors(iPulse) = find(iPulse <= colorCutoffs,1);
        allPulsePowers(iPulse) = find((iPulse - (allPulseColors(iPulse)-1)*nReps*length(pulsePower)) <= powerCutoffs,1);

        % make sure pulse and baseline period to 10 ms after pulse doesn't have artifact
        if ~isempty(intersect(optoStarts(iPulse)-baselinePeriod*30:optoStarts(iPulse)+300,artifactInds))
            pulseHasArtifact(iPulse) = true;
        else
            pulseHasArtifact(iPulse) = false;
        end

        histEdges = optoStarts(iPulse)+30:30:optoStarts(iPulse)+300;
        spikeCounts{iNeuron}(iPulse,:) = histcounts(double(neuronDataStruct(iNeuron).timeStamps),histEdges);

        for iControl = 1:nBaselineSamples
            spikeCountsBaseline{iNeuron}(iPulse,:,iControl) = histcounts(double(neuronDataStruct(iNeuron).timeStamps),histEdges-baselineOffsets(iControl)*30);
        end

        %extract spike waveforms from raw data
        thisPulseTimes = intersect(double(neuronDataStruct(iNeuron).timeStamps), optoStarts(iPulse)+30:optoStarts(iPulse)+300);
        if ~isempty(thisPulseTimes)
            for iSpike = 1:length(thisPulseTimes)
                fseek(readFid,(thisPulseTimes(iSpike)-20)*nChans*2,'bof');
                voltDat = fread(readFid, [nChans, 50], 'int16=>double');
                waveforms{iNeuron,1}(end+1,:) = voltDat(neuronDataStruct(iNeuron).biggestChan,:);
            end
        end

        rasterTimes{iNeuron,iPulse} = intersect(double(neuronDataStruct(iNeuron).timeStamps), optoStarts(iPulse)-1500:optoStarts(iPulse)+3000) - optoStarts(iPulse);
        rasterTimes{iNeuron,iPulse}(rasterTimes{iNeuron,iPulse}>=0 & rasterTimes{iNeuron,iPulse} <= 30) = [];

        %         %load in 300 ms before and 100 ms after each opto pulse
        %         nPreSamples = 30*50;
        %         nPostSamples = 30*100;
        %         startReadSample = optoStarts(iPulse) - nPreSamples;
        %         fseek(readFid, startReadSample * 2 * nChans, 'bof');
        %         rawData = fread(readFid, [nChans, nPreSamples+nPostSamples], 'int16=>double');
        %
        %         pulseSigs(:,iNeuron,iPulse) = rawData(neuronDataStruct(iNeuron).biggestChan,:);

    end

    % get some waveforms from before the opto part of the session
    nOptoWaveforms = size(waveforms{iNeuron,1},1);
    timeCutoff = optoStarts(1)-30*30000; % only use time points at least 30 s before the first opto pulse
    nonOptoTimeStamps = neuronDataStruct(iNeuron).timeStamps(neuronDataStruct(iNeuron).timeStamps < timeCutoff);
    if length(nonOptoTimeStamps) <= nOptoWaveforms
        usedNonOptos = 1:length(nonOptoTimeStamps);
    else
        usedNonOptos = randperm(length(nonOptoTimeStamps),nOptoWaveforms);
    end
    for iSpike = 1:length(usedNonOptos)
        fseek(readFid,(nonOptoTimeStamps(usedNonOptos(iSpike))-20)*nChans*2,'bof');
        voltDat = fread(readFid, [nChans, 50], 'int16=>double');
        waveforms{iNeuron,2}(end+1,:) = voltDat(neuronDataStruct(iNeuron).biggestChan,:);
    end

    % for each power and duration compute the J-S divergence
    for iColor = 1:length(pulseColor)
        for iPower = 1:length(pulsePower)

            %get pulses corresponding to these parameters
            pulseInds = find(allPulsePowers == iPower & allPulseColors == iColor);

            pulseInds(pulseHasArtifact(pulseInds)) = [];

            spikeDistribution{iNeuron,iColor,iPower} = spikeCounts{iNeuron}(pulseInds,:);
            spikeControlDistribution{iNeuron,iColor,iPower} = spikeCountsBaseline{iNeuron}(pulseInds,:,:);

            % get J-S divergence between real distribution and control
            % distributions
            if useTotalSpikes
                for iControl = 1:nBaselineSamples
                    realDivs(iControl) = calcJSDivergence(spikeDistribution{iNeuron,iColor,iPower}, ...
                        spikeControlDistribution{iNeuron,iColor,iPower}(iControl,:));
                end

                controlPairs = nchoosek(1:nBaselineSamples,2);
                for iControlPair = 1:size(controlPairs,1)

                    controlDivs(iControlPair) = calcJSDivergence(spikeControlDistribution{iNeuron,iColor,iPower}(controlPairs(iControlPair,1),:), ...
                        spikeControlDistribution{iNeuron,iColor,iPower}(controlPairs(iControlPair,2),:));

                end

                % get p-value
                pVal(iNeuron,iColor,iPower) = (sum(nanmedian(realDivs) <= controlDivs)+1)/(sum(~isnan(controlDivs))+1);

                allControlDivs{iNeuron,iColor,iPower} = controlDivs;
                allRealDivs{iNeuron,iColor,iPower} = realDivs;

                controlDistributionPlot = spikeControlDistribution{iNeuron,iColor,iPower}';
                distributionPlot = spikeDistribution{iNeuron,iColor,iPower};

            else
                spt_baseline = spikeControlDistribution{iNeuron,iColor,iPower}(:,:);
                spt_test = spikeDistribution{iNeuron,iColor,iPower};
                [pVal(iNeuron,iColor,iPower), effSize(iNeuron,iColor,iPower), nhlsi, jsd] = salt(spt_baseline,spt_test,0.001,0.009);
                controlJSs = jsd(1:length(jsd)-1,1:length(jsd)-1);
                allControlDivs{iNeuron,iColor,iPower} = controlJSs(~isnan(controlJSs));
                allRealDivs{iNeuron,iColor,iPower} = jsd(1:end-1,end);
                controlDistributionPlot = nhlsi(:,1:200);
                distributionPlot = nhlsi(:,201);
            end

            % make plots
            if iPower == 1
                plotH = figure('Units','pixels','OuterPosition',[50 50 1800,1000],'color','w','Visible','on');
                tileH = tiledlayout(3,length(pulsePower),'Padding','compact','TileSpacing','compact');
            else
                figure(plotH)
            end
            nexttile(iPower + 0*length(pulsePower))
            rasterplot(cellfun(@(x) double(x)/30,rasterTimes(iNeuron,pulseInds),'un',0),'times','|')
            hold on;
            xlim([-10 10])
            title({['Neuron ' num2str(iNeuron) ', Pulse Color = ' pulseColor{iColor} ', Pulse Power=' num2str(pulsePower(iPower)) 'mW']; ...
                ['p=' num2str(round(pVal(iNeuron,iColor,iPower)*1e5)/1e5), 'effect size=' num2str(round(effSize(iNeuron,iColor,iPower)*1e5)/1e5)]})
            xlabel('Time since pulse on (ms)')
            ylabel('Trial')
            set(gca,'XTick',-10:1:10)
            nexttile(iPower + 1*length(pulsePower))
            plot(0:9,controlDistributionPlot)
            hold on
            plot(0:9,distributionPlot,'k','LineWidth',2)
            xlim([0 10])
            xlabel('Time bin (ms)')
            ylabel('# Spikes')
            histPlotH = nexttile(iPower + 2*length(pulsePower));
            histogram(allControlDivs{iNeuron,iColor,iPower},50)
            hold on
            line(repmat(nanmedian(allRealDivs{iNeuron,iColor,iPower}),2,1),get(gca,'ylim'),'color','k','linewidth',2)
            xlabel('J-S Divergence')
            ylabel('Count')

            if size(waveforms{iNeuron,1},1) >= 5 && size(waveforms{iNeuron,2},1) >= 5
                p = get(gca, 'Position');
                axes(plotH,Units="normalized",Position=[p(1)+0.005, p(2)+0.17, p(3)*0.2, p(4)*0.25])
                shadedErrorBar([],mean(waveforms{iNeuron,1}),std(waveforms{iNeuron,1}),'lineProps',{'Color','b'})
                ylims = get(gca,'YLim');
                axis off
                axes(plotH,Units="normalized",Position=[p(1)+0.05, p(2)+0.17, p(3)*0.2, p(4)*0.25])
                shadedErrorBar([],mean(waveforms{iNeuron,2}),std(waveforms{iNeuron,2}),'lineProps',{'Color','k'})
                set(gca,'YLim',ylims);
                axis off
            end
            axes(histPlotH)

            if iPower == length(pulsePower)
                saveas(plotH,fullfile(baseDir,'..','ProcessedData','SALT',['Neuron' num2str(iNeuron) '_Color' pulseColor{iColor} '.png']),'png')
            end


        end

        close(plotH);

    end

end

save('saltAnalysis','allPulseColors','allPulsePowers','pulseHasArtifact','rasterTimes','spikeDistribution','spikeControlDistribution',...
    'allControlDivs','allRealDivs','pVal','effSize','waveforms','-v7.3')

% end



function jsDiv = calcJSDivergence(distribution1, distribution2)

if sum(distribution1) == 0
    distr1Norm = ones(1,length(distribution1))/length(distribution1);
else
    distr1Norm = distribution1/sum(distribution1);
end

if sum(distribution2) == 0
    distr2Norm = ones(1,length(distribution2))/length(distribution2);
else
    distr2Norm = distribution2/sum(distribution2);
end

klDiv = @(x,y) nansum(x.*log(x./y));
jsDiv = 0.5*(klDiv(distr1Norm,0.5*(distr1Norm+distr2Norm)) + klDiv(distr2Norm,0.5*(distr1Norm+distr2Norm)));

end


%
