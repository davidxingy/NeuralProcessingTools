clear
close all

% all session locations
sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording\ProcessedData',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording\ProcessedData'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    7 6 3 5 4 2 1; ...
    1 2 4 5 3 6 7 ...
    ];

striatalDepthCutoffs = [2600 2400 2400];

plotColors = lines(7);

for iSess = 1:length(sessionDirs)

    % go through each session and get neural activity distribution across
    % behaviors
    selectivity = [];
    behvNeurs = {};
    behvAlignPerm = allBehvAlignPerms(iSess,:);

    load(fullfile(sessionDirs{iSess},'UMAPFRs','NeuronRegionProps.mat'));
    load(fullfile(sessionDirs{iSess},'NeuralFiringRates1msBins10msGauss.mat'),'cortexInds','striatumInds');
    load(fullfile(sessionDirs{iSess},'neuronDataStruct.mat'));

    normAveSigs = regionAveSigs(:,behvAlignPerm)./sum(regionAveSigs,2);
    depths = [neuronDataStruct([striatumInds cortexInds]).depth];

    % calculate selectivity, and also get neurons that are only modulated
    % by eating
    regionAveSigsShuffMean = mean(regionAveSigsShuff(:,behvAlignPerm,:),3);
    regionAveSigsShuffStd = std(regionAveSigsShuff(:,behvAlignPerm,:),[],3);
    regionAveSigsZscore = (regionAveSigs(:,behvAlignPerm,:) - regionAveSigsShuffMean)./regionAveSigsShuffStd;

    for iBehv = 1:size(regionAveSigs,2)

        % old method, used selectivity as the way to find behavior specific
        % neurons
        selectivity(:,iBehv) = normAveSigs(:,iBehv)./max(normAveSigs(:,setdiff(1:7,iBehv)),[],2);
        behvNeurs{iBehv} = find(selectivity(:,iBehv)>4);
        behvDepths{iBehv,iSess} = depths(behvNeurs{iBehv});

        % only use neurons with a minimum firing rate
        [~,badNeurs] = intersect(behvNeurs{iBehv}, find([all(regionAveSigs < 0.5,2)]));
        behvNeursGood{iBehv} = behvNeurs{iBehv};
        behvDepthsGood{iBehv,iSess} = behvDepths{iBehv,iSess};

        behvDepthsGood{iBehv,iSess}(badNeurs) = [];
        behvNeursGood{iBehv}(badNeurs) = [];

        % better method, use statistical method (3 std above mean of
        % shuffles) to find behavior specific neurons
        behvShuffNeurs{iBehv} = find(regionAveSigsZscore(:,iBehv) > 3 & all(regionAveSigsZscore(:,setdiff(1:7,iBehv)) <= 3,2));
        behvShuffDepths{iBehv,iSess} = depths(behvShuffNeurs{iBehv});

        % again, only use neurons which are above a minimim firing rate
        [~,badNeurs] = intersect(behvShuffNeurs{iBehv}, find([all(regionAveSigs < 0.5,2)]));
        behvShuffNeursGood{iBehv} = behvShuffNeurs{iBehv};
        behvShuffDepthsGood{iBehv,iSess} = behvShuffDepths{iBehv,iSess};

        behvShuffDepthsGood{iBehv,iSess}(badNeurs) = [];
        behvShuffNeursGood{iBehv}(badNeurs) = [];
    end

    % now get the histogram counts across striatal depth
    depthBins = 0:100:2600;
    counts{iSess} = hist(striatalDepthCutoffs(iSess) - behvShuffDepthsGood{7,iSess},depthBins); 
    counts{iSess} = counts{iSess}/sum(counts{iSess}(2:end));

    % find the cutoff of the eating region using 80% threshold
    % (depth below which contains 80% of the eating striatum neurons)
    strEatNeurons = sort(behvShuffDepthsGood{7,iSess}(behvShuffDepthsGood{7,iSess} <= striatalDepthCutoffs(iSess)));
    eatBoundaries(iSess) = strEatNeurons(ceil(length(strEatNeurons)*0.8));

end

aveCounts = mean(cat(1,counts{:}));

hold on
for iSess = 1:3
    plot(counts{iSess}(2:end),depthBins(2:end),'color',[plotColors(iSess,:) 0.2],'LineWidth',1.5)
    line([0 0.2],repmat(striatalDepthCutoffs(iSess) - eatBoundaries(iSess),1,2),...
        'linestyle','--','color',[plotColors(iSess,:) 0.7],'LineWidth',1.5)
end

plot(aveCounts(2:end),depthBins(2:end),'color','k','LineWidth',2.5)
box off
set(gcf,'Color','w')
set(gca,'TickDir','out')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',14)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
ylabel('Striatum Depth (um)')
xlabel('Fraction of eating-specific cells')
set(gca,'YDir','reverse')



% 

