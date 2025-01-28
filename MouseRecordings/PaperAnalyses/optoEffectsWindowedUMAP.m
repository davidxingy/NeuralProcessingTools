% function allChansFiltData = optoUMAPOverlay(baseDir)
% clearvars -except D040Sessions
% baseDir = 'X:\David\ArenaRecordings\D040-110223-ArenaRecording';

D036Sessions = {'X:\David\ArenaRecordings\D036-101623-ArenaRecording','X:\David\ArenaRecordings\D036-101723-ArenaRecording',...
    'X:\David\ArenaRecordings\D036-101823-ArenaRecording','X:\David\ArenaRecordings\D036-101923-ArenaRecording',...
    'X:\David\ArenaRecordings\D036-102023-ArenaRecording'};

D040Sessions = {'X:\David\ArenaRecordings\D040-110223-ArenaRecording','X:\David\ArenaRecordings\D040-110323-ArenaRecording',...
    'X:\David\ArenaRecordings\D040-110423-ArenaRecording','X:\David\ArenaRecordings\D040-110623-ArenaRecording',...
    'X:\David\ArenaRecordings\D040-110723-ArenaRecording','X:\David\ArenaRecordings\D040-110823-ArenaRecording'};

D041Sessions = {'X:\David\ArenaRecordings\D041-121123-ArenaRecording','X:\David\ArenaRecordings\D041-121223-ArenaRecording',...
    'X:\David\ArenaRecordings\D041-121323-ArenaRecording','X:\David\ArenaRecordings\D041-121423-ArenaRecording',...
    'X:\David\ArenaRecordings\D041-121523-ArenaRecording'};

baseDirs = D040Sessions;
nShuffs = 100;
runUMAP = true;
if runUMAP

    for iSession = 1:length(baseDirs)

        tic

        % load behavior labels
        if iSession == 1 || iSession == length(baseDirs)
            load(fullfile(baseDirs{iSession},'ProcessedData','UMAP'),'behvLabels');

            % load in sync data
            load(fullfile(baseDirs{iSession},'ProcessedData','VideoSyncFrames.mat'))

            %load in EMG channel names
            emgMetaFile = dir(fullfile(baseDirs{1},'ProcessedData','*_MetaData.mat'));
            if ~isempty (emgMetaFile)
                load(fullfile(baseDirs{1},'ProcessedData',emgMetaFile.name),'channelNames')
            end
        end

        % load EMG data
        load(fullfile(baseDirs{iSession},'ProcessedData','EMG1ms'))

        %get meta data as well
        sessionFiles = string(ls(fullfile(baseDirs{iSession},'ProcessedData')));
        metaDataFile = sessionFiles(contains(sessionFiles,'ProcessedEMG_MetaData'));
        load(fullfile(baseDirs{iSession},'ProcessedData', metaDataFile))

        artifactInds = unique(downsampleRemovedInds);
        sessTotalPulses = length(downsampleLaserOnsetInds);

        if iSession == 2
            %         usedPulseInds = downsampleLaserOnsetInds(1:3900);
        end

        % get control pulse times by taking random time between 200-400ms
        % before the actual laser pulses
        control1Offsets = randi(200,1,length(downsampleLaserOnsetInds))+100;
        downsampleControl1OnsetInds = downsampleLaserOnsetInds - control1Offsets;

        control2Offsets = randi(200,1,length(downsampleLaserOnsetInds))+300;
        downsampleControl2OnsetInds = downsampleLaserOnsetInds - control2Offsets;

        % See which pulses can be used
        usePulseCounter = 1;
        usePulse = true(1,length(downsampleLaserOnsetInds));
        laserReducInds = [];
        control1ReducInds = [];
        control2ReducInds = [];
        for iPulse = 1:length(downsampleLaserOnsetInds)

            %get time windows for the laser and control pulses that we need to
            %use (-10ms to 30ms) and make sure those points don't have
            %artifacts

            laserWindowInds = downsampleLaserOnsetInds(iPulse)-80:downsampleLaserOnsetInds(iPulse)+90;
            control1WindowInds = downsampleControl1OnsetInds(iPulse)-80:downsampleControl1OnsetInds(iPulse)+90;
            control2WindowInds = downsampleControl2OnsetInds(iPulse)-80:downsampleControl2OnsetInds(iPulse)+90;


            if ~isempty(intersect(unique([laserWindowInds control1WindowInds control2WindowInds]), artifactInds))
                usePulse(iPulse) = false;
                continue
            end
            %
            %         laserWindowInds{iSession}()
            %         laserReducInds(usePulseCounter) = find(downsampleLaserOnsetInds(iPulse) == origDownsampEMGInd);
            %         control1ReducInds(usePulseCounter) = find(downsampleControlOnsetInds(iPulse) == origDownsampEMGInd);
            %
            %         for iShuff = 1:nShuffs
            %             shuffReducInds(iShuff,usePulseCounter) = find(downsampleShuffOnsetInds(iShuff,iPulse) == origDownsampEMGInd);
            %         end
            %
            usePulseCounter = usePulseCounter+1;

        end

        % now go through each used pulse and get the EMG window around it and
        % the behavior label
        usedPulses = find(usePulse);
        for iPulse = 1:length(usedPulses)

            thisLaserInd = downsampleLaserOnsetInds(usedPulses(iPulse));
            sessionLaserEmgs{iSession}(iPulse,:,:) = downsampEMG(:,thisLaserInd-80:thisLaserInd+90);
            thisControl1Ind = downsampleControl1OnsetInds(usedPulses(iPulse));
            sessionControl1Emgs{iSession}(iPulse,:,:) = downsampEMG(:,thisControl1Ind-80:thisControl1Ind+90);
            thisControl2Ind = downsampleControl2OnsetInds(usedPulses(iPulse));
            sessionControl2Emgs{iSession}(iPulse,:,:) = downsampEMG(:,thisControl2Ind-80:thisControl2Ind+90);

            if iSession == 1 %|| iSession == length(baseDirs)
                sessionLaserBehvs{iSession}(iPulse,:) = behvLabels(thisLaserInd-80:thisLaserInd+90);
                sessionControl1Behvs{iSession}(iPulse,:) = behvLabels(thisControl1Ind-80:thisControl1Ind+90);
                sessionControl2Behvs{iSession}(iPulse,:) = behvLabels(thisControl2Ind-80:thisControl2Ind+90);
            else
                sessionLaserBehvs{iSession}(iPulse,:) = zeros(1,171);
                sessionControl1Behvs{iSession}(iPulse,:) = zeros(1,171);
                sessionControl2Behvs{iSession}(iPulse,:) = zeros(1,171);
            end

            % also convert to UMAP features
            [features, ~, behvs] = extractFeatures(squeeze(sessionLaserEmgs{iSession}(iPulse,:,:)),sessionLaserBehvs{iSession}(iPulse,:));
            sessionLaserFeatures{iSession}(iPulse,:,:) = features;
            sessionLaserFeatureBehvs{iSession}(iPulse,:,:) = behvs;
            [features, ~, behvs] = extractFeatures(squeeze(sessionControl1Emgs{iSession}(iPulse,:,:)),sessionControl1Behvs{iSession}(iPulse,:));
            sessionControl1Features{iSession}(iPulse,:,:) = features;
            sessionControl1FeatureBehvs{iSession}(iPulse,:,:) = behvs;
            [features, ~, behvs] = extractFeatures(squeeze(sessionControl2Emgs{iSession}(iPulse,:,:)),sessionControl2Behvs{iSession}(iPulse,:));
            sessionControl2Features{iSession}(iPulse,:,:) = features;
            sessionControl2FeatureBehvs{iSession}(iPulse,:,:) = behvs;

        end

        disp(['Loaded session ' num2str(iSession) ', ' num2str(toc) 's'])

    end

    % combine sessions
    allLaserEMGs = cat(1,sessionLaserEmgs{:});
    allControl1EMGs = cat(1,sessionControl1Emgs{:});
    allControl2EMGs = cat(1,sessionControl2Emgs{:});

    allLaserFeatures = cat(1,sessionLaserFeatures{:});
    allControl1Features = cat(1,sessionControl1Features{:});
    allControl2Features = cat(1,sessionControl2Features{:});

    allLaserBehvs = cat(1,sessionLaserBehvs{:});
    allControl1Behvs = cat(1,sessionControl1Behvs{:});
    allControl2Behvs = cat(1,sessionControl2Behvs{:});

    allLaserFeatureBehvs = cat(1,sessionLaserFeatureBehvs{:});
    allControl1FeatureBehvs = cat(1,sessionControl1FeatureBehvs{:});
    allControl2FeatureBehvs = cat(1,sessionControl2FeatureBehvs{:});

    % calc UMAP
    allFeatures = cat(3,permute(allLaserFeatures,[2 3 1]),permute(allControl1Features,[2 3 1]),permute(allControl2Features,[2 3 1]));
    allFeatures = allFeatures(:,:)';
    distKL = @(x,y) log(y./sum(y,2))*(x/sum(x))'*-1+repmat((x/sum(x))*log((x/sum(x))'),size(y,1),1);
    % inputFeatures = featuresData(intersect(indsToUse,[annotatedBehvLabels]),:);

    projDownSamp = 1;
    nUMAPDims = 2;
    nUMapNeighbors = 30;

    [reduction,umap,clusterIdentifiers,extras] = run_umap(allFeatures(1:projDownSamp:end,:),...
        'n_components',nUMAPDims,'n_neighbors',nUMapNeighbors,'save_template_file', 'projUMAPTemplate.mat','Distance','Euclidean');

    % save umap output and data
    save('allSessionsUMAP','allLaserEMGs','allControl1EMGs','allControl2EMGs','allLaserFeatures','allControl1Features','allControl2Features',...
        'allLaserBehvs','allControl1Behvs','allControl2Behvs','projDownSamp','nUMAPDims','nUMapNeighbors','reduction','umap',...
        'allLaserFeatureBehvs','allControl1FeatureBehvs','allControl2FeatureBehvs','-v7.3')

else

    load('allSessionsUMAP')

end


% divide up reduction back up
groupSize = size(allLaserFeatures,1)*size(allLaserFeatures,3);
tempArray = reduction(1:groupSize,:);
laserReductions = permute(reshape(tempArray,size(allLaserFeatures,3),size(allLaserFeatures,1),2), [3,1,2]);
tempArray = reduction(groupSize+1:groupSize*2,:);
control1Reductions = permute(reshape(tempArray,size(allLaserFeatures,3),size(allLaserFeatures,1),2), [3,1,2]);
tempArray = reduction(groupSize*2+1:groupSize*3,:);
control2Reductions = permute(reshape(tempArray,size(allLaserFeatures,3),size(allLaserFeatures,1),2), [3,1,2]);

% plot UMAP for visualization
plotColors = turbo(10);
figure;
hold on;
for iBehv = 1:10

    for iWindow = 1:size(laserReductions,1)

        tmpH = plot(squeeze(laserReductions(1,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),...
            squeeze(laserReductions(2,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),'.','Color',plotColors(iBehv,:),'MarkerSize',0.1);

        plot(squeeze(control1Reductions(1,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),...
            squeeze(control1Reductions(2,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),'.','Color',plotColors(iBehv,:),'MarkerSize',0.1)

        plot(squeeze(control2Reductions(1,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),...
            squeeze(control2Reductions(2,iWindow,allLaserFeatureBehvs(:,iWindow)==iBehv)),'.','Color',plotColors(iBehv,:),'MarkerSize',0.1)


    end
    
    plotH(iBehv) = tmpH;
end

legH = legend(plotH,analyzedBehaviors);
legH.Box = "off";
legH.FontSize = 12;
for iLabel = 1:length(legH.String)
legH.String{iLabel} = ['\color[rgb]{' num2str(plotH(iLabel).Color) '} ' legH.String{iLabel}];
end



% z-score EMGs for calculating effect sizes
flattenedLaserEmgs = permute(allLaserEMGs,[2 3 1]);
flattenedLaserEmgs = flattenedLaserEmgs(:,:);
emgMeans = mean(flattenedLaserEmgs,2);
emgStds = std(flattenedLaserEmgs,[],2);

allLaserEMGs = (allLaserEMGs - repmat(emgMeans',size(allLaserEMGs,1),1,size(allLaserEMGs,3)))./...
    repmat(emgStds',size(allLaserEMGs,1),1,size(allLaserEMGs,3));
allControl1EMGs = (allControl1EMGs - repmat(emgMeans',size(allControl1EMGs,1),1,size(allControl1EMGs,3)))./...
    repmat(emgStds',size(allControl1EMGs,1),1,size(allControl1EMGs,3));
allControl2EMGs = (allControl2EMGs - repmat(emgMeans',size(allControl2EMGs,1),1,size(allControl2EMGs,3)))./...
    repmat(emgStds',size(allControl2EMGs,1),1,size(allControl2EMGs,3));

% make grid points
nXGridPoints = 100;
nYGridPoints = 100;

UMAPextents = [max(reduction(:,1)) min(reduction(:,1)) max(reduction(:,2)) min(reduction(:,2))];

xGridPoints = linspace(UMAPextents(2)*1.05,UMAPextents(1)*1.05,nXGridPoints);
yGridPoints = linspace(UMAPextents(4)*1.05,UMAPextents(3)*1.05,nXGridPoints);

% make control shuffs for null
combControlReductions = cat(3,control1Reductions,control2Reductions);
combControlEMGs = cat(1,allControl1EMGs,allControl2EMGs);
for iShuff = 1:nShuffs
    
    null1Inds(iShuff,:) = randperm(size(combControlReductions,3),size(combControlReductions,3)/2);
    null2Inds(iShuff,:) = setdiff(1:size(combControlReductions,3),null1Inds(iShuff,:));
    
end


% quest version, get all the weights first, then export and let quest do
% the weighted averaging which takes the most time. Then import back in and
% calculate effect sizes (fast)
if useQuest
    save('gridCalcVars','laserReductions','control1Reductions','combControlReductions','xGridPoints','yGridPoints',...
        'nShuffs','null1Inds','null2Inds','allLaserEMGs','allControl1EMGs','combControlEMGs','UMAPextents');
else
    
    for iGridX = 1:nXGridPoints
        for iGridY = 1:nYGridPoints
            
            tic
            
            %calculate the distance between the grid point and each of the
            %laser and control points
            
            %first just use UMAP location at 10ms
            laserDists = squeeze(sqrt((laserReductions(1,5,:) - xGridPoints(iGridX)).^2 + ...
                (laserReductions(2,5,:) - yGridPoints(iGridY)).^2));
            controlDists = squeeze(sqrt((control1Reductions(1,5,:) - xGridPoints(iGridX)).^2 + ...
                (control1Reductions(2,5,:) - yGridPoints(iGridY)).^2));
            
            for iShuff = 1:nShuffs
                null1Dists(iShuff,:) = sqrt((combControlReductions(1,5,null1Inds(iShuff,:)) - xGridPoints(iGridX)).^2 + ...
                    (combControlReductions(2,5,null1Inds(iShuff,:)) - yGridPoints(iGridY)).^2);
                null2Dists(iShuff,:) = sqrt((combControlReductions(1,5,null2Inds(iShuff,:)) - xGridPoints(iGridX)).^2 + ...
                    (combControlReductions(2,5,null2Inds(iShuff,:)) - yGridPoints(iGridY)).^2);
            end
            
            %change distance to a weight, using gaussian
            gaussKernal = mean([UMAPextents(1)-UMAPextents(2),UMAPextents(3)-UMAPextents(4)])/10;
            laserWeights = exp(laserDists.^2/(-2*gaussKernal^2));
            controlWeights = exp(controlDists.^2/(-2*gaussKernal^2));
            
            null1Weights = exp(null1Dists.^2/(-2*gaussKernal^2));
            null2Weights = exp(null2Dists.^2/(-2*gaussKernal^2));
            
            %save sum of weights
            laserWeightsSum(iGridX,iGridY) = sum(laserWeights);
            controlWeightsSum(iGridX,iGridY) = sum(controlWeights);
            null1WeightsSum(iGridX,iGridY,:) = sum(null1Weights,2);
            null2WeightsSum(iGridX,iGridY,:) = sum(null2Weights,2);
            
            % now get weighted average EMG time series
            laserEMG(iGridX,iGridY,:,:) = squeeze(pagemtimes(laserWeights'/sum(laserWeights),allLaserEMGs));
            controlEMG(iGridX,iGridY,:,:) = squeeze(pagemtimes(controlWeights'/sum(controlWeights),allControl1EMGs));
            
            for iShuff = 1:nShuffs
                
                thisNull1EMG = squeeze(pagemtimes(null1Weights(iShuff,:)/sum(null1Weights(iShuff,:)), ...
                    squeeze(combControlEMGs(null1Inds(iShuff,:),:,:))));
                
                thisNull2EMG = squeeze(pagemtimes(null2Weights(iShuff,:)/sum(null2Weights(iShuff,:)), ...
                    squeeze(combControlEMGs(null2Inds(iShuff,:),:,:))));
                
                %just save the first shuffle since don't want to take up too
                %much memory
                if iShuff==1
                    null1EMG(iGridX,iGridY,:,:) = thisNull1EMG;
                    null2EMG(iGridX,iGridY,:,:) = thisNull2EMG;
                end
                
            end
            
            disp(['X Grid point ' num2str(iGridX) ', Y Grid Point ' num2str(iGridY) ', ' num2str(toc) 's'])
            
        end
    end
    
end

% calculate slopes and effect sizes
if useQuest
    load('questOutputVars')
end

effectSizes = laserSlope - controlSlope;
effectSizesLong = laserSlopeLong - controlSlopeLong;

for iShuff = 1:nShuffs
    effectSizesShuff(iShuff,:,:,:) = squeeze(null2Slope(iShuff,:,:,:)) - squeeze(null1Slope(iShuff,:,:,:));
    effectSizesLongShuff(iShuff,:,:,:) = squeeze(null2SlopeLong(iShuff,:,:,:)) - squeeze(null1SlopeLong(iShuff,:,:,:));
end

% generate p-values
for iGridX = 1:size(effectSizes,1)
    for iGridY = 1:size(effectSizes,2)
        for iMusc = 1:size(effectSizes,3)
            
            pVal(iGridX,iGridY,iMusc) = sum(abs(effectSizesShuff(:,iGridX,iGridY,iMusc)) > abs(effectSizes(iGridX,iGridY,iMusc)))/size(effectSizesShuff,1);
            pValLong(iGridX,iGridY,iMusc) = sum(abs(effectSizesLongShuff(:,iGridX,iGridY,iMusc)) > abs(effectSizesLong(iGridX,iGridY,iMusc)))/size(effectSizesLongShuff,1);

        end
    end
end

% save
save('umapEffectSizes','effectSizes','laserSlope','controlSlope','effectSizesShuff','null1Slope','null2Slope',...
    'laserEMG','controlEMG','null1EMG','null2EMG','laserWeightsSum','controlWeightsSum','null1WeightsSum',...
    'null2WeightsSum','gaussKernal','null1Inds','null2Inds');

% % get watershedding region splitting
% [regionBoundaryIndsY, regionBoundaryIndsX] = find(watershedRegions==0);

% make mask for grid points that have sufficient closeby points
alphaMask = logisticTransparency(laserWeightsSum,1000,50)';

% make some plots
animalTitle = 'D040';
figure
tH = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
title(tH,[animalTitle ', Fast Latency (10-30ms)'])

for iMusc = 1:size(effectSizes,3)

    nexttile
    imagesc(xGridPoints,yGridPoints,effectSizes(:,:,iMusc)','AlphaData',alphaMask);
%     hold on
%     plot(UMAPGridXInds(regionBoundaryIndsX),UMAPGridYInds(regionBoundaryIndsY),'.k')
    set(gca,'YDir','normal')
    set(gcf,'color','w')
    title(channelNames{iMusc})
    axis off
    colorbar

end

for iMusc = 1:size(effectSizes,3)

    nexttile
    imagesc(xGridPoints,yGridPoints,pVal(:,:,iMusc)','AlphaData',alphaMask);
%     hold on
%     plot(UMAPGridXInds(regionBoundaryIndsX),UMAPGridYInds(regionBoundaryIndsY),'.k')
    set(gca,'YDir','normal')
    set(gca,'CLim',[0 0.05])
    set(gcf,'color','w')
    title(channelNames{iMusc})
    axis off
    colorbar

end

% 
figure
tH = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
title(tH,[animalTitle ', Long Latency (70-90ms)'])

for iMusc = 1:size(effectSizes,3)

    nexttile
    imagesc(xGridPoints,yGridPoints,effectSizesLong(:,:,iMusc)','AlphaData',alphaMask);
%     hold on
%     plot(UMAPGridXInds(regionBoundaryIndsX),UMAPGridYInds(regionBoundaryIndsY),'.k')
    set(gca,'YDir','normal')
    set(gcf,'color','w')
    title(channelNames{iMusc})
    axis off
    colorbar

end

for iMusc = 1:size(effectSizes,3)

    nexttile
    imagesc(xGridPoints,yGridPoints,pValLong(:,:,iMusc)','AlphaData',alphaMask);
%     hold on
%     plot(UMAPGridXInds(regionBoundaryIndsX),UMAPGridYInds(regionBoundaryIndsY),'.k')
    set(gca,'YDir','normal')
    set(gca,'CLim',[0 0.05])
    set(gcf,'color','w')
    title(channelNames{iMusc})
    axis off
    colorbar

end


% p-value histogram:
figure('color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
histogram(pVal(:))
box off
xlabel('p-values')
ylabel('# grid points')
title([animalTitle ', short latency (10-30ms)'])

nexttile
histogram(pValLong(:))
box off
xlabel('p-values')
ylabel('# grid points')
title([animalTitle ', long latency (70-90ms)'])


figure
tH = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
title(tH,'Null Control')

shuffToPlot = 5;
for iMusc = 1:size(effectSizesInst,3)

    nexttile
    imagesc(xGridPoints,yGridPoints,squeeze(effectSizesShuffInst(shuffToPlot,:,:,iMusc))','AlphaData',alphaMask)
    hold on
    plot(UMAPGridXInds(regionBoundaryIndsX),UMAPGridYInds(regionBoundaryIndsY),'.k')
    set(gca,'YDir','normal')
    set(gcf,'color','w')
    title(channelNames{iMusc})
    colorbar

end



function [features, outputInds, featureBehvs] = extractFeatures(signal,behaviorLabels)

% overall, divide into 50ms windows, sliding every 10ms
% however, each of the windows will average 5 ms bins

windowSize = 50;
slideAmount = 10;
binSize = 5;


% only use the points that fit into 10 ms slidings
nSlides = floor((size(signal,2)-1)/slideAmount);
nWindows = nSlides - windowSize/slideAmount + 1;

signalToUse = signal(:,1:slideAmount*nSlides);

% also get the derrivative (here difference) of the input sig
signalDiffs = diff(signal,1,2);
signalToUseDiffs = signalDiffs(:,1:slideAmount*nSlides);

% now get average in 5 ms bins
meanSig = zeros(size(signalToUse,1),size(signalToUse,2)/binSize);
meanDiffSig = zeros(size(signalToUse,1),size(signalToUse,2)/binSize);
for iTime = 1:binSize

    meanSig = meanSig+signalToUse(:,iTime:binSize:end);
    meanDiffSig = meanDiffSig+signalToUseDiffs(:,iTime:binSize:end);

end

meanSig = meanSig/binSize;
meanDiffSig = meanDiffSig/binSize;

% each window contains 50 ms total, so shift 10 times
meanSigFeatures = addHistory(meanSig,windowSize/binSize-1,2,-1);
meanDiffSigFeatures = addHistory(meanDiffSig,windowSize/binSize-1,2,-1);

% final output features are the combined signal and signal derrivative
features = [meanSigFeatures(:,1:2:end-windowSize/binSize+1) ; meanDiffSigFeatures(:,1:2:end-windowSize/binSize+1)];

outputInds = (0:size(features,2)-1)*10+1;

% calculate which behavior is labeled most in the 50 ms window of the
% features
for iInd = 1:windowSize
    allIndsLabels(:,iInd) = behaviorLabels(iInd:slideAmount:end-windowSize+iInd);
end
featureBehvs = mode(allIndsLabels,2);


end



% 

