clear
close all

sessionNames = {'D020-062922-ArenaRecording',...
    'D024-111022-ArenaRecording',...
    'D026-032923-ArenaRecording',...
    'D036-101623-ArenaRecording',...
    'D040-110223-ArenaRecording',...
    'D041-121123-ArenaRecording',...
    'D050-121825-ArenaRecording',...
    'D054-012426-ArenaRecording',... %maybe better umap?
    'D056-020926-ArenaRecording'};

% do for human annotated behaviors as well as for UMAP behavior regions
annotatedbehaviorLabels = {...
'Climb up',...
'Climb down',...
'Eat',...
'Groom',...
'Jump down',...
'Jump',...
'Rear',...
'Still',...
'Walk flat',...
'Walk Grid'};

umapBehaviorLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};

allUmapAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    2 1 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7 ...
    ];
allAnnotAlignPerms = [
    [2 1 3:10]; ...
    [2 1 3:10]; ...
    [2 1 3:6 8:11]; ...
    [2 1 3:10]; ...
    [2 1 3:10]; ...
    [2 1 3:6 8:11]; ...
    [2 1 3:10]; ...
    [2 1 3:10]; ...
    [2 1 3:10] ...
    ];

behvDurations = struct();
behvDurationsFraction = struct();
behvNumBouts = struct();
behvNumBoutsFraction = struct();
behvGoodInds = struct();
labelTypeNames = {'UMAPRegions','HumanAnnotation'};

for iSess = 1:length(sessionNames)

    dataNames = getMouseDataNames(sessionNames{iSess}(1:4),sessionNames{iSess},'CFA');
    
%     if iSess <= 3
%         % load video annotationed data
%         load(fullfile(dataNames.processedDataFolder,'EpochedData10ms.mat'))
%         structFieldNames = fieldnames(behavioralData);
% 
%         % total time points annotated
%         totalDurationAnnotated(iSess) = sum(cellfun(@(x) size(x.allBoutFRs,2),struct2cell(behavioralData)))/100;
%         totalNumBoutsAnnotated(iSess) = sum(cellfun(@(x) length(x.boutFRs),struct2cell(behavioralData)));
% 
%         % get behavior specifc durations and n bouts
%         for iBehv = 1:length(annotatedbehaviorFieldnames)
% 
%             behvInd = find(strcmpi(structFieldNames,annotatedbehaviorFieldnames{iBehv}));
%             behvDurationsAnnotated(iSess,iBehv) = size(behavioralData.(annotatedbehaviorFieldnames{iBehv}).allBoutFRs,2)/100;
%             behvDurationsFractionAnnotated(iSess,iBehv) = behvDurationsAnnotated(iSess,iBehv)/totalDurationAnnotated(iSess);
% 
%             behvNumBoutsAnnotated(iSess,iBehv) = length(behavioralData.(annotatedbehaviorFieldnames{iBehv}).boutFRs);
%             behvNumBoutsFractionAnnotated(iSess,iBehv) = behvNumBoutsAnnotated(iSess,iBehv)/totalNumBoutsAnnotated(iSess);
% 
%         end
%     end

    % now do for UMAP annotated
    load(dataNames.UMAPFile,'regionAssignmentsFiltered','regionWatershedLabels','behvLabelsNoArt')
    regionWatershedLabels = regionWatershedLabels(allUmapAlignPerms(iSess,:));
    

    labels = {regionAssignmentsFiltered, behvLabelsNoArt};
    labelNames = {umapBehaviorLabels, annotatedbehaviorLabels};
    labelIds = {regionWatershedLabels, allAnnotAlignPerms(iSess,:)};
    
    for iLabel = 1:length(labels)
        
        usedLabels = labels{iLabel};
        
        % total time points and bouts
        boutStartInds = [1 find(diff(usedLabels)~=0)+1];
        boutEndInds = [find(diff(usedLabels)~=0) length(usedLabels)];
        boutDurations = boutEndInds - boutStartInds + 1;
        boutRegions = usedLabels(boutStartInds);

        % don't count bouts that last less than 50ms (too short)
        goodBouts = boutDurations>50 & .... %and also only use bouts that correspond to a used behavior
            any(repmat(boutRegions,length(labelIds{iLabel}),1) == repmat(labelIds{iLabel}',1,length(boutRegions)));
        boutStartInds = boutStartInds(goodBouts);
        boutEndInds = boutEndInds(goodBouts);
        boutDurations = boutDurations(goodBouts);
        boutRegions = boutRegions(goodBouts);
        
        totalDurationUMAP.(labelTypeNames{iLabel})(iSess) = sum(boutDurations)/1000;
        totalNumBoutsUMAP.(labelTypeNames{iLabel})(iSess) = length(boutStartInds);
        
        behvGoodInds.(labelTypeNames{iLabel}) = zeros(length(labelNames{iLabel}), length(usedLabels));
        for iBehv = 1:length(labelNames{iLabel})
            
            behvBoutInds = find(boutRegions==labelIds{iLabel}(iBehv));
            
            %get indices corresponding to all bouts of this behavior
            for iBout = 1:length(behvBoutInds)
                behvGoodInds.(labelTypeNames{iLabel})(iBehv,boutStartInds(behvBoutInds(iBout)):boutEndInds(behvBoutInds(iBout))) = 1;
            end
            
            behvDurations.(labelTypeNames{iLabel})(iSess,iBehv) = sum(boutDurations(behvBoutInds))/1000;
            behvDurationsFraction.(labelTypeNames{iLabel})(iSess,iBehv) = ...
                behvDurations.(labelTypeNames{iLabel})(iSess,iBehv) / totalDurationUMAP.(labelTypeNames{iLabel})(iSess);
            
            behvNumBouts.(labelTypeNames{iLabel})(iSess,iBehv) = length(behvBoutInds);
            behvNumBoutsFraction.(labelTypeNames{iLabel})(iSess,iBehv) = ...
                behvNumBouts.(labelTypeNames{iLabel})(iSess,iBehv) / totalNumBoutsUMAP.(labelTypeNames{iLabel})(iSess);
            
        end
        
    end
    
    % now, calculate the correspondence between the umap regions and the
    % human annotations
    for iUmap = 1:length(labelIds{1})
        for iAnnot = 1:length(labelIds{2})
            thisUmapInds = logical(behvGoodInds.(labelTypeNames{1})(iUmap,:));
            thisAnnotInds = logical(behvGoodInds.(labelTypeNames{2})(iAnnot,:));
            nAnnotPerUmap(iUmap,iAnnot,iSess) = sum(thisUmapInds(thisAnnotInds));
            nUmapPerAnnot(iAnnot,iUmap,iSess) = sum(thisAnnotInds(thisUmapInds));
        end
    end

end

% make behavior region frequency plots
figure;
tiledlayout(2,1,'padding','compact','tileSpacing','compact')

scatterData = mat2cell(behvDurationsFraction.UMAPRegions,length(sessionNames),ones(1,length(umapBehaviorLabels)))';
plotJitter = repmat({linspace(-0.1,0.1,length(scatterData{1}))},length(scatterData),1);

nexttile
barH = barScatterPlot(scatterData,'none',ones(size(scatterData)),plotJitter,[]);
ylabel('Fraction of time')
set(gca,'XTickLabel',umapBehaviorLabels)
barH.BarWidth = 0.6;

scatterData = mat2cell(behvNumBouts.UMAPRegions,length(sessionNames),ones(1,length(umapBehaviorLabels)))';
nexttile
barH = barScatterPlot(scatterData,'none',ones(size(scatterData)),plotJitter,[]);
ylabel('Fraction of Total # of Bouts')
set(gca,'XTickLabel',umapBehaviorLabels)
barH.BarWidth = 0.6;


% make umap region - human annotation correspondence plots
figure
fracAnnotPerUmap = nAnnotPerUmap./sum(nAnnotPerUmap,2);
tiledlayout(size(fracAnnotPerUmap,1),1,'padding','compact','tileSpacing','compact')

for iRegion = 1:size(fracAnnotPerUmap,1)
    
    nexttile
    scatterData = mat2cell(squeeze(fracAnnotPerUmap(iRegion,:,:))',length(sessionNames),ones(1,size(fracAnnotPerUmap,2)))';
    plotJitter = repmat({linspace(-0.1,0.1,length(scatterData{1}))},length(scatterData),1);
    barH = barScatterPlot(scatterData,'none',ones(size(scatterData)),plotJitter,[]);
    ylim([0 1])
    set(gca,'ytick',[0 1])
    barH.BarWidth = 0.6;
    if iRegion == size(fracAnnotPerUmap,1)
        set(gca,'XTickLabelRotation',45)
        set(gca,'XTickLabel',annotatedbehaviorLabels)
    else
        set(gca,'xtick',[])
    end
    ylabel(umapBehaviorLabels{iRegion})
    
end


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
