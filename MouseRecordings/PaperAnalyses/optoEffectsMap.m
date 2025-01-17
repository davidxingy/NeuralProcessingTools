% function allChansFiltData = optoUMAPOverlay(baseDir)
clearvars -except D040Sessions
% baseDir = 'X:\David\ArenaRecordings\D040-110223-ArenaRecording';

D040Sessions = {'X:\David\ArenaRecordings\D040-110223-ArenaRecording','X:\David\ArenaRecordings\D040-110323-ArenaRecording'};

baseDirs = D040Sessions;

for iSession = 1:length(baseDirs)

    % load UMAP projection
    if iSession == 1
        load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'analyzedBehaviors','behvLabelsNoArt','regionBehvAssignments','regionWatershedLabels',...
            'gridXInds','gridYInds','watershedRegions','annotatedBehvLabels');
        % load in sync data
        load(fullfile(baseDirs{iSession},'ProcessedData','VideoSyncFrames.mat'))
    end

    load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'reduction','origDownsampEMGInd','regionAssignmentsNoBound')

    % load EMG data
    load(fullfile(baseDirs{iSession},'ProcessedData','EMG1ms'))

    %get meta data as well
    sessionFiles = string(ls(fullfile(baseDirs{iSession},'ProcessedData')));
    metaDataFile = sessionFiles(contains(sessionFiles,'ProcessedEMG_MetaData'));
    load(fullfile(baseDirs{iSession},'ProcessedData', metaDataFile))

    sessTotalPulses = length(downsampleLaserOnsetInds);

    if iSession == 2
        %         usedPulseInds = downsampleLaserOnsetInds(1:3900);
    end

    % get control pulse times by taking random time between 200-400ms
    % before the actual laser pulses
    controlOffsets = randi(200,1,length(downsampleLaserOnsetInds))+200;
    downsampleControlOnsetInds = downsampleLaserOnsetInds - controlOffsets;
    
    % convert to UMAP reduction time inds
    usePulseCounter = 1;
    for iPulse = 1:length(downsampleLaserOnsetInds)
        
        %get time windows for the laser and control pulses that we need to
        %use (-10ms to 30ms) and make sure those points are all in the UMAP
        %reduction
        laserWindow = downsampleLaserOnsetInds(iPulse)-10:downsampleLaserOnsetInds(iPulse)+30;
        controlWindow = downsampleControlOnsetInds(iPulse)-10:downsampleControlOnsetInds(iPulse)+30;
        
        if ~isempty(setdiff([laserWindow controlWindow],origDownsampEMGInd))
            usePulse(iPulse) = false;
            continue
        end
        
        laserReducInds(usePulseCounter) = find(downsampleLaserOnsetInds(iPulse) == origDownsampEMGInd);
        controlReducInds(usePulseCounter) = find(downsampleControlOnsetInds(iPulse) == origDownsampEMGInd);
        
        usePulseCounter = usePulseCounter+1;
        
    end
    
    % now go through each used pulse and get the EMG window around it and
    % the reduction
    usedPulses = find(usePulse);
    for iPulse = 1:length(laserReducInds)
        
        thisLaserInd = downsampleLaserOnsetInds(usedPulses(iPulse));
        sessionLaserEmgs{iSession}(iPulse,:,:) = downsampEMG(:,thisLaserInd-10:thisLaserInd+30);
        thisControlInd = downsampleControlOnsetInds(usedPulses(iPulse));
        sessionControlEmgs{iSession}(iPulse,:,:) = downsampEMG(:,thisControlInd-10:thisControlInd+30);
        
        sessionLaserReducs{iSession}(iPulse,:) = reduction(laserReducInds(iPulse)-10:laserReducInds(iPulse)+10);
        sessionControlReducs{iSession}(iPulse,:) = reduction(controlReducInds(iPulse)-10:controlReducInds(iPulse)+10);
        
    end

end
