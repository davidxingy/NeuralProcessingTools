clear

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

allAnnotatedBehvAlignPerms = [
    2 1 5 9 10 7 4 3; ...
    2 1 5 9 10 7 4 3; ...
    2 1 5 10 11 8 4 3; ...
    ];

allAnimalLabels = {'D020','D024','D026'};
behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};
annotatedBehvLabels = {'Climb Up','Climb Down','Jump Down','Walk Flat','Walk Grid','Rear','Groom','Eat'};

onlyHumanAnnotated = false;

for iAnimal = 1:length(allDirs)

    baseDir = allDirs{iAnimal};
    load(fullfile(baseDir,'ProcessedData','UMAP'),'reduction','freqData','behvLabelsNoArt','regionAssignmentsFiltered','regionWatershedLabels','origDownsampEMGInd','analyzedBehaviors')
%     load(fullfile(baseDir,'ProcessedData','TSNE'),'reduction','behvLabelsNoArt','origDownsampEMGInd','analyzedBehaviors')
    behvLabelsNoArt = behvLabelsNoArt(behvLabelsNoArt~=0);
    behvLabelsNoArt = behvLabelsNoArt(1:50:end);
    load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))
    timeEmgData = downsampEMG(:,origDownsampEMGInd)';

    % get time points where there's only human labels (that corresponds to
    % the behaviors that we want)
    annotatedBehvClassLabels = allAnnotatedBehvAlignPerms(iAnimal,:);
    annotatedInds = find(any(repmat(behvLabelsNoArt,size(allAnnotatedBehvAlignPerms,2),1) == repmat(annotatedBehvClassLabels',1,size(behvLabelsNoArt,2))));

%     for iBehv = 1:length(annotatedBehvClassLabels)
%         %calculate silhouette scores
%         clustScores{iAnimal,iBehv} = simplifiedSilhouette(freqData(annotatedInds,:),behvLabelsNoArt(annotatedInds));
%     end

    if onlyHumanAnnotated
        freqData = freqData(annotatedInds,:);
        timeEmgData = timeEmgData(annotatedInds,:);
        reduction = reduction(annotatedInds,:);
        behvLabelsNoArt = behvLabelsNoArt(annotatedInds);
        nBehvRegions = length(annotatedBehvClassLabels);
    else
        nBehvRegions = size(allBehvAlignPerms,2);
    end

%     [projs,trajs,vaf] = pca(freqData);

    behvAlignPerm = allBehvAlignPerms(iAnimal,:);
    regionWatershedLabels = regionWatershedLabels(behvAlignPerm);

    % go through each pair of behavior regions and calculate distance
    for iRegion1 = 1:nBehvRegions
        for iRegion2 = 1:nBehvRegions

            % calculate distance in UMAP space

            if onlyHumanAnnotated
                region1Inds = behvLabelsNoArt == annotatedBehvClassLabels(iRegion1);
                region2Inds = behvLabelsNoArt == annotatedBehvClassLabels(iRegion2);
            else
                region1Inds = regionWatershedLabels(iRegion1) == regionAssignmentsFiltered;
                region2Inds = regionWatershedLabels(iRegion2) == regionAssignmentsFiltered;
            end

            umapDists(iAnimal,iRegion1,iRegion2) = sqrt(sum((median(reduction(region1Inds,:)) - median(reduction(region2Inds,:))).^2));

            tsneDists(iAnimal,iRegion1,iRegion2) = sqrt(sum((median(reduction(region1Inds,:)) - median(reduction(region2Inds,:))).^2));

            %as well as in full high dimensional frequency space
            freqDists(iAnimal,iRegion1,iRegion2) = sqrt(sum((median(freqData(region1Inds,:)) - median(freqData(region2Inds,:))).^2));

            %also look at it in time domain EMG activity space
            timeEmgDists(iAnimal,iRegion1,iRegion2) = sqrt(sum((median(timeEmgData(region1Inds,:)) - median(timeEmgData(region2Inds,:))).^2));

            %also look at it in PCA projected freq space
            pcaDists(iAnimal,iRegion1,iRegion2) = sqrt(sum((median(trajs(region1Inds,1:3)) - median(trajs(region2Inds,1:3))).^2));

        end
    end
    

end

% 
