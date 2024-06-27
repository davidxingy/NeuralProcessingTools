clear
close all

sessionNames = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

sessionNames = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording',...
    'X:\David\ArenaRecordings\D036-101623-ArenaRecording',...
    'X:\David\ArenaRecordings\D040-110223-ArenaRecording'};


% do for human annotated behaviors as well as for UMAP behavior regions

annotatedbehaviorFieldnames = {...
'climbup',...
'climbdown',...
'jumpdown',...
'jumping',...
'walkflat',...
'walkgrid',...
'rearing',...
'still',...
'eating',...
'grooming'};

annotatedbehaviorLabels = {...
'Climb Up',...
'Climb Down',...
'Jump Down',...
'Jump Across',...
'Walk Flat',...
'Walk Grid',...
'Rear',...
'Still',...
'Groom',...
'Eat'};

umapBehaviorLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rearing/Still','Groom','Eat'};
allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7; ...
    1 2 3 5 4 6 7; ...
    1 2 3 5 4 6 7 ...
    ];

for iSess = 1:length(sessionNames)

    if iSess <= 3
        % load video annotationed data
        load(fullfile(sessionNames{iSess},'ProcessedData','EpochedData10ms.mat'))
        structFieldNames = fieldnames(behavioralData);

        % total time points annotated
        totalDurationAnnotated(iSess) = sum(cellfun(@(x) size(x.allBoutFRs,2),struct2cell(behavioralData)))/100;
        totalNumBoutsAnnotated(iSess) = sum(cellfun(@(x) length(x.boutFRs),struct2cell(behavioralData)));

        % get behavior specifc durations and n bouts
        for iBehv = 1:length(annotatedbehaviorFieldnames)

            behvInd = find(strcmpi(structFieldNames,annotatedbehaviorFieldnames{iBehv}));
            behvDurationsAnnotated(iSess,iBehv) = size(behavioralData.(annotatedbehaviorFieldnames{iBehv}).allBoutFRs,2)/100;
            behvDurationsFractionAnnotated(iSess,iBehv) = behvDurationsAnnotated(iSess,iBehv)/totalDurationAnnotated(iSess);

            behvNumBoutsAnnotated(iSess,iBehv) = length(behavioralData.(annotatedbehaviorFieldnames{iBehv}).boutFRs);
            behvNumBoutsFractionAnnotated(iSess,iBehv) = behvNumBoutsAnnotated(iSess,iBehv)/totalNumBoutsAnnotated(iSess);

        end
    end

    % now do for UMAP annotated
    load(fullfile(sessionNames{iSess},'ProcessedData','UMAP.mat'),'regionAssignmentsFiltered','regionWatershedLabels')
    regionWatershedLabels = regionWatershedLabels(allBehvAlignPerms(iSess,:));

    % total time points and bouts
    boutStartInds = [1 find(diff(regionAssignmentsFiltered)~=0)+1];
    boutEndInds = [find(diff(regionAssignmentsFiltered)~=0) length(regionAssignmentsFiltered)];
    boutDurations = boutEndInds - boutStartInds;

    % don't count bouts that last less than 50ms (too short)
    goodBouts = boutDurations>50;
    boutStartInds = boutStartInds(goodBouts);
    boutEndInds = boutEndInds(goodBouts);
    boutDurations = boutDurations(goodBouts);
    boutRegions = regionAssignmentsFiltered(boutStartInds);

    totalDurationUMAP(iSess) = sum(boutDurations)/1000;
    totalNumBoutsUMAP(iSess) = length(boutStartInds);

    for iBehv = 1:length(umapBehaviorLabels)

        regionBoutInds = find(boutRegions==regionWatershedLabels(iBehv));
        behvDurationsUMAP(iSess,iBehv) = sum(boutDurations(regionBoutInds))/1000;
        behvDurationsFractionUMAP(iSess,iBehv) = behvDurationsUMAP(iSess,iBehv) / totalDurationUMAP(iSess);

        behvNumBoutsUMAP(iSess,iBehv) = length(regionBoutInds);
        behvNumBoutsFractionUMAP(iSess,iBehv) = behvNumBoutsUMAP(iSess,iBehv) / totalNumBoutsUMAP(iSess);

    end

end

figure;
tiledlayout(2,1,'padding','compact','tileSpacing','compact')

scatterData = mat2cell(behvDurationsFractionUMAP,length(sessionNames),ones(1,length(umapBehaviorLabels)))';
plotJitter = repmat({[-0.10:0.05:0.1]},length(scatterData),1);

nexttile
barH = barScatterPlot(scatterData,'none',1,plotJitter,[]);
ylabel('Fraction of time')
set(gca,'XTickLabel',umapBehaviorLabels)
barH.BarWidth = 0.6;

scatterData = mat2cell(behvNumBoutsFractionUMAP,length(sessionNames),ones(1,length(umapBehaviorLabels)))';
nexttile
barH = barScatterPlot(scatterData,'none',1,plotJitter,[]);
ylabel('Fraction of Total # of Bouts')
set(gca,'XTickLabel',umapBehaviorLabels)
barH.BarWidth = 0.6;



figure;
tiledlayout(2,1,'padding','compact','tileSpacing','compact')

scatterData = mat2cell(behvDurationsFractionAnnotated,size(behvDurationsAnnotated,1),ones(1,length(annotatedbehaviorLabels)))';
plotJitter = repmat({[-0.1:0.1:0.1]},length(scatterData),1);

nexttile
barH = barScatterPlot(scatterData,'none',1,plotJitter,[]);
ylabel('Fraction of time')
set(gca,'XTickLabel',annotatedbehaviorLabels)
barH.BarWidth = 0.6;

scatterData = mat2cell(behvNumBoutsFractionAnnotated,size(behvDurationsAnnotated,1),ones(1,length(annotatedbehaviorLabels)))';
nexttile
barH = barScatterPlot(scatterData,'none',1,plotJitter,[]);
ylabel('Fraction of Total # of Bouts')
set(gca,'XTickLabel',annotatedbehaviorLabels)
barH.BarWidth = 0.6;


% 
