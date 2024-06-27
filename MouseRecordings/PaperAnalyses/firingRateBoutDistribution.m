clear
close all

sessionNames = {'pcaAlignment_D020.mat','pcaAlignment_D024.mat','pcaAlignment_D026.mat'};
sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

behvRegionLabels = {'Climb Up','Climb Down','Jump Down/Across','Walk Flat/Grid','Rearing/Still','Grooming','Eating'};
minBoutLength = [1000 500 300 1000 500 1000 1000]; %in ms
minBoutSpacing = [1000 500 1000 1000 500 1000 500]; %in ms

exampleNeurons = [129 65 174 144 142 52 18;...
                  1 1 1 1 1 1 1;...
                  119 112 87 104 50 110 39];

for iSess = 1:length(sessionNames)

    behvAlignPerm = allBehvAlignPerms(iSess,:);
    baseDir = sessionDirs{iSess};

    %load in spiking data
    % use 1ms bins
    load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'),'allFRs','striatumInds','cortexInds');

     % load UMAP projection
    load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','regionAssignmentsFiltered','regionWatershedLabels')

    %align neur data to reduction
    currentDir = pwd;
    cd(fullfile(baseDir,'ProcessedData'))
    reducNeurInds = round(NeurEMGSync(origDownsampEMGInd*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
    maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/30)-1;
    reducIndsToUse = find(reducNeurInds < maxNeurSamples);
    cd(currentDir)

    % get neural data corresponding to the UMAP time points
    reducFRs = allFRs(:,reducNeurInds(reducIndsToUse));
    clear allFRs

    % don't use any points where there's nans in the FR data
    nanReducPoints = any(isnan(reducFRs),1);
    reducIndsToUse(nanReducPoints) = [];

    reducFRs = reducFRs(:,~nanReducPoints);

    %get the behavior labels
    behvLabels = regionAssignmentsFiltered(reducIndsToUse);
    behvClasses = unique(behvLabels);

    %now get the neurons that are behavior specific
    load(fullfile(baseDir,'ProcessedData','UMAPFRs','NeuronRegionProps.mat'))
    regionAveSigsShuff = regionAveSigsShuff(:,behvAlignPerm,:);
    regionAveSigs = regionAveSigs(:,behvAlignPerm);

    %determine neurons that are modulated for behavior
    shuffMean = mean(regionAveSigsShuff,3);
    shuffStd = std(regionAveSigsShuff,[],3);
    modulation = regionAveSigs >= shuffMean + 3*shuffStd;
    behvSpecNeurons = find(sum(modulation,2) == 1);
    singleBehvSpec = zeros(size(modulation,1),1);
    for iNeuron = 1:length(behvSpecNeurons)
        singleBehvSpec(behvSpecNeurons(iNeuron)) = find(modulation(behvSpecNeurons(iNeuron),:));
    end
    
    %also just get multi-behavioral specifity for finding interesting
    %neurons to plot
    multiBehvSpec = sqrt(sum((regionAveSigs./sum(regionAveSigs,2)).^2,2));

    %go through each behavior region
    for iBehv = 1:length(behvClasses)

        thisBehvInds = find(behvLabels==behvClasses(iBehv));
        boutStarts = [thisBehvInds(1) thisBehvInds(find(diff(thisBehvInds)>1)+1)];
        boutEnds = [thisBehvInds(diff(thisBehvInds)>1) thisBehvInds(end)];

        boutDurations = boutEnds - boutStarts;
        
        %only use bouts of a certain duration, and a minimum time away from
        %the previous bout
        goodBouts = find((boutDurations > minBoutLength(iBehv)) & [1 (boutStarts(2:end)-boutEnds(1:end-1) > minBoutSpacing(iBehv))]);

        if length(goodBouts) < 10
            warning(['Less than 10 bouts found in session ' num2str(iSess) ', behavior ' behvRegionLabels{iBehv}])
        end

        nBouts(iSess,iBehv) = length(goodBouts);

        boutStarts = boutStarts(goodBouts);
        boutEnds = boutEnds(goodBouts);
        boutDurations = boutDurations(goodBouts);

        boutFRs = {};
        boutInds = {};
        for iBout = 1:length(boutStarts)
            boutFRs{iBout} = reducFRs(:,boutStarts(iBout):boutEnds(iBout));
            boutInds{iBout} = boutStarts(iBout):boutEnds(iBout);
        end

        %make plots
        plotExtAmount = round(mean(boutDurations));

        figure
        plotBouts = 1:min(length(boutStarts),10);
        tiledlayout(2,length(plotBouts),'TileSpacing','compact','Padding','compact')
        exampleBehvNeuron = exampleNeurons(iSess,iBehv);

        for iBout = 1:length(plotBouts)-1
            axH(1) = nexttile(iBout);
            plot(reducFRs(exampleBehvNeuron,boutStarts(iBout+1)-plotExtAmount:boutEnds(iBout+1)+plotExtAmount),'LineWidth',2)
            axH(2) = nexttile(iBout + length(plotBouts));
            plotLabel = zeros(1,length(boutStarts(iBout+1)-plotExtAmount:boutEnds(iBout+1)+plotExtAmount));
            plotLabel(plotExtAmount:end-plotExtAmount) = 1;
            plot(plotLabel,'r','LineWidth',2);
        end

        title(behvRegionLabels{iBehv})

        allBoutInds = cat(2,boutInds{:});
        allBoutFRs = cat(2,boutFRs{:});
        maxFRs = prctile(allBoutFRs,95,2);

        %go through each bout and time warp to percentage of bout
        normalizedBoutFR = [];
        for iBout = 1:length(boutStarts)
            
            desiredPercentageAsInds = interp1([1 100], [1 size(boutFRs{iBout},2)], 1:100, 'linear', 'extrap');
            normalizedBoutFR(:,:,iBout) = interp1(1:size(boutFRs{iBout},2), (boutFRs{iBout}./maxFRs)', desiredPercentageAsInds, 'cubic')';

        end

        singleSpecNeurons = find(singleBehvSpec==iBehv);
        modulatedNeurons = find(modulation(:,iBehv));
        singleSpecBoutDist{iSess,iBehv} = squeeze(mean(normalizedBoutFR(singleSpecNeurons,:,:),3));
        singleSpecBoutDist{iSess,iBehv}(any(isnan(singleSpecBoutDist{iSess,iBehv}),2),:) = [];
        modulatedBoutDist{iSess,iBehv} = squeeze(mean(normalizedBoutFR(modulatedNeurons,:,:),3));
        modulatedBoutDist{iSess,iBehv}(any(isnan(modulatedBoutDist{iSess,iBehv}),2),:) = [];

    end

end





% 
