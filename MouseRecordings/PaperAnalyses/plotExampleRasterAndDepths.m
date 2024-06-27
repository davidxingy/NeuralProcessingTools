clear
close all

% all session locations
sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    7 6 3 5 4 2 1; ...
    1 2 4 5 3 6 7 ...
    ];

striatalDepthCutoffs = [2600 2400 2400];

plotColors = lines(7);
maxPlotTime = 30; %in seconds

for iSess = 1:length(sessionDirs)

    % load in data
    load(fullfile(sessionDirs{iSess},'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'cortexInds','striatumInds','allFRs');
    load(fullfile(sessionDirs{iSess},'ProcessedData','neuronDataStruct.mat'));

    % get depth and timestamps for each neuron
    depths = [neuronDataStruct([striatumInds cortexInds]).depth];
    spikeTimes = {neuronDataStruct([striatumInds cortexInds]).timeStamps};
    spikeTimes = cellfun(@(x) double(x(x<30000*maxPlotTime))/30000,spikeTimes,'UniformOutput',false);

    % get artifact times and remove spikes in those periods
    load(fullfile(sessionDirs{iSess},'Neuropixels','artifactTimestamps.mat'));
    artifactBoutStarts = [artifactTS(1) artifactTS(find(diff(artifactTS)>1)+1)]/30000;
    artifactBoutEnds = [artifactTS(find(diff(artifactTS)>1)) artifactTS(end)]/30000;
    for iBout = 1:length(artifactBoutEnds)
        for iNeuron = 1:length(spikeTimes)
            spikeTimes{iNeuron}(spikeTimes{iNeuron} > artifactBoutStarts(iBout) & spikeTimes{iNeuron} < artifactBoutEnds(iBout)) = [];
        end
    end

    % plot raster
    figure;
    tiledlayout(1,4,'TileSpacing','tight')
    axH(1) = nexttile(1,[1 3]);
    rasterplot(spikeTimes,'times','.',[],[],depths+randn(1,length(depths))*2)

    % add regions indicating artifact times
    hold on
    for iBout = 1:length(find(artifactBoutEnds<maxPlotTime))
        patchH = patch([repmat(artifactBoutStarts(iBout),1,2) repmat(artifactBoutEnds(iBout),1,2)],[0 4000 4000 0],'r');
        patchH.EdgeColor = 'none';
        patchH.FaceAlpha = 0.07;
    end

    line([0 5],[100 100],'color','k','linewidth',3)

    % make pretty
    set(axH(1),'Box','off')
    set(axH(1),'fontsize',14)
    set(axH(1),'tickdir','out')
    set(axH(1),'linewidth',1.5)
    set(axH(1),'XColor','k')
    set(axH(1),'YColor','k')
    set(axH(1),'ytick',0:1000:4000)
    set(axH(1),'yticklabel',4000:-1000:0)
    ylabel('Depth (mm)')
    xlabel('Time (s)')

    % plot depth histogram
    axH(2) = nexttile(4);
    histH = histogram(depths,0:100:4000,'Orientation','horizontal');
    hold on
    line(get(axH(2),'XLim'),repmat(striatalDepthCutoffs(iSess),1,2),'linewidth',2,'linestyle','--','color','k')
    histH.EdgeColor = 'none';
    set(axH(2),'Box','off')
    set(axH(2),'ytick',[])
    set(axH(2),'fontsize',14)
    set(axH(2),'tickdir','out')
    set(axH(2),'linewidth',1.5)
    set(axH(2),'XColor','k')
    set(axH(2),'YColor','k')
    xlabel('# of Neurons');
    linkaxes(axH,'y')

    set(gca,'YLim',[0 4000])
    set(gcf,'color','w')

    nStrCells(iSess) = length(striatumInds);
    nCtxCells(iSess) = length(cortexInds);

    %get average firing rate distributions
    strMeanFRs{iSess} = nanmean(allFRs(1:length(striatumInds),:),2)'*100;
    ctxMeanFRs{iSess} = nanmean(allFRs(length(striatumInds)+1:end,:),2)'*100;
    
    %plot distribution as cdf function
    cdfVals = 0:1:70;

    for iVal = 1:length(cdfVals)
        strCdfFreq{iSess}(iVal) = sum(strMeanFRs{iSess} <= cdfVals(iVal)) / length(strMeanFRs{iSess});
        ctxCdfFreq{iSess}(iVal) = sum(ctxMeanFRs{iSess} <= cdfVals(iVal)) / length(ctxMeanFRs{iSess});
    end

end

% plot number of cells
figure
tiledlayout(1,3,'Padding','tight','TileSpacing','compact')
scatterData{1} = nCtxCells;
scatterData{2} = nStrCells;

nexttile(1,[1 1])
barH = barScatterPlot(scatterData,'none',ones(1,2),repmat({[-0.05 0 0.05]},1,2),[]);
set(gca,'xtick',1:2)
set(gca,'xticklabels',{'Striatum','Cortex'})
xlim([0.5 2.5])

ylabel('# Neurons')

% plot firing rate distributions
nexttile(2,[1,2])
hold on
plotColors = lines(7);

% plot individual sessions
for iSess = 1:length(sessionDirs)

    firstMaxVal = find(strCdfFreq{iSess} == 1,1);
    plot(cdfVals(1:firstMaxVal),strCdfFreq{iSess}(1:firstMaxVal),'Color',[plotColors(2,:) 0.2],'LineWidth',1.5)
    firstMaxVal = find(ctxCdfFreq{iSess} == 1,1);
    plot(cdfVals(1:firstMaxVal),ctxCdfFreq{iSess}(1:firstMaxVal),'Color',[plotColors(1,:) 0.2],'LineWidth',1.5)

end

% plot combined across sessions
for iVal = 1:length(cdfVals)
    strAllCdfFreq(iVal) = sum([strMeanFRs{:}] <= cdfVals(iVal)) / length([strMeanFRs{:}]);
    ctxAllCdfFreq(iVal) = sum([ctxMeanFRs{:}] <= cdfVals(iVal)) / length([ctxMeanFRs{:}]);
end

firstMaxVal = find(strAllCdfFreq == 1,1);
plot(cdfVals(1:firstMaxVal),strAllCdfFreq(1:firstMaxVal),'Color',plotColors(2,:),'LineWidth',3)
firstMaxVal = find(ctxAllCdfFreq == 1,1);
plot(cdfVals(1:firstMaxVal),ctxAllCdfFreq(1:firstMaxVal),'Color',plotColors(1,:),'LineWidth',3)

line(repmat(mean([strMeanFRs{:}]),1,2),[0 1],'color',plotColors(2,:),'linewidth',1,'linestyle','--')
line(repmat(mean([ctxMeanFRs{:}]),1,2),[0 1],'color',plotColors(1,:),'linewidth',1,'linestyle','--')

xlim([0 40])
xlabel('Firing Rate (Hz)')
ylabel('Frequency')

box off
set(gca,'fontsize',14)
set(gca,'tickdir','out')
set(gca,'linewidth',1.5)
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gcf,'color','w')

% 
