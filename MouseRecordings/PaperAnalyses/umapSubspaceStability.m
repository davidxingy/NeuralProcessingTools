clear

% do analysis for each of the datasets
recordingSessions = {
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording', ...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording' ...
};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];


inputData = 'umapregions';
use10msBins = true;
nShifts = 1;

for iSess = 1:length(recordingSessions)

    baseDir = recordingSessions{iSess};
    behvAlignPerm = allBehvAlignPerms(iSess,:);

    switch lower(inputData)

        case 'umapregions'
            % -----Use data divded by the 7 umap regions-----

            behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rearing/Still','Groom','Eat'};

            if use10msBins
                load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'), 'cortexInds', 'striatumInds','allFRs')
            else
                load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds','allFRs')
            end

            load(fullfile(baseDir,'ProcessedData','UMAP.mat'), 'regionAssignmentsFiltered',...
                'origDownsampEMGInd','subRegionAssignments','analyzedBehaviors','regionBehvAssignments','regionWatershedLabels')
            load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))
            load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

            % annotated behavior stuff
            load(fullfile(baseDir,'ProcessedData','EpochedData10ms.mat'))
            behvNames = fieldnames(behavioralData);

            % Sync EMG and neural data
            % outputInds = NeurEMGSync(inputInds, frameEMGSamples, frameNeuropixelSamples, inputType);

            % normalize data
            normalizedFRs = (allFRs-nanmean(allFRs,2))./nanstd(allFRs,[],2);
            % normalizedFRs(isinf(normalizedFRs)) = nan;

            % if we want to downsample emg beyond 1ms bins, take the average across bins
            if use10msBins
                emg1ms = downsampEMG;
                downsampleFactor = 10;
                clear downsampEMG
                n1msBins = size(emg1ms,2);
                ndownSampBins = floor(n1msBins/downsampleFactor);

                downsampEMG = reshape(emg1ms(:,1:ndownSampBins*downsampleFactor),size(emg1ms,1),downsampleFactor,ndownSampBins);
                downsampEMG = squeeze(mean(downsampEMG,2));
                if mod(n1msBins,downsampleFactor) ~= 0
                    downsampEMG(:,end+1) = mean(emg1ms(:,ndownSampBins*downsampleFactor+1:end),2);
                end
            end

            normalizedEMG = (downsampEMG-nanmean(downsampEMG,2))./nanstd(downsampEMG,[],2);

            %align behavior regions across animals
            regionWatershedLabels = regionWatershedLabels(behvAlignPerm);
            nRegions = length(regionWatershedLabels);
            %     regionAssignmentsFiltered(:) = 1;
            %go through each region and divide data
            for iRegion = 1:nRegions

                [regionEMGInds{iRegion}, regionNeurInds{iRegion}, regionEMGs{iRegion}, regionFRs{iRegion}, regionFRsOrig{iRegion}] = ...
                    getBehvRegionData(regionAssignmentsFiltered, regionWatershedLabels(iRegion), frameEMGSamples, frameNeuropixelSamples, use10msBins, ...
                    origDownsampEMGInd, baseDir, normalizedEMG, normalizedFRs, allFRs);

                % do shift controls to break behavior relationship with neur/emg
                for iShift = 1:nShifts
                    %shift by at least 30 seconds
                    if use10msBins
                        shiftAmount(iShift) = randi(length(regionAssignmentsFiltered)-3000*2)+3000;
                    else
                        shiftAmount(iShift) = randi(length(regionAssignmentsFiltered)-30000*2)+30000;
                    end
                    regionAssignmentsShifted = circshift(regionAssignmentsFiltered,shiftAmount(iShift));

                    [regionEMGIndsShift{iRegion,iShift}, regionNeurIndsShift{iRegion,iShift}, regionEMGsShift{iRegion,iShift},...
                        regionFRsShift{iRegion,iShift}, ~] = ...
                        getBehvRegionData(regionAssignmentsShifted, regionWatershedLabels(iRegion), frameEMGSamples, frameNeuropixelSamples, use10msBins, ...
                        origDownsampEMGInd, baseDir, normalizedEMG, normalizedFRs, allFRs);
                end

            end

        case 'humanannotated'
            % -----Use data divded by human annotations-----

            behvRegionFieldNames = {'climbup','climbdown','jumpdown','walkflat','walkgrid','rearing','grooming','eating'};
            behvRegionLabels = {'Climb Up','Climb Down','Jump Down','Walk Flat','Walk Grid','Rear','Groom','Eat'};

            load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds')
            if use10msBins
                load(fullfile(baseDir,'ProcessedData','EpochedData10ms.mat'))
            else
                load(fullfile(baseDir,'ProcessedData','EpochedData1ms.mat'))
            end

            sessStructNames = fieldnames(behavioralData);
            nRegions = length(behvRegionFieldNames);

            %get data for each behavior
            for iBehv = 1:nRegions

                behvStructInd(iBehv) = find(strcmpi(sessStructNames,behvRegionFieldNames{iBehv}));
                regionFRsOrig{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).allBoutFRs;
                regionEMGsOrig{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).allBoutEMGs;
                
                %remove nans
                nanInds = unique([find(any(isnan(regionFRsOrig{iBehv}),1)) find(any(isnan(regionEMGsOrig{iBehv}),1))]);
                regionFRsOrig{iBehv}(:,nanInds) = [];
                regionEMGsOrig{iBehv}(:,nanInds) = [];

                %also get information for doing the shift controls
                behvBoutFRs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutFRs;
                behvBoutEMGs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutEMGs;

                %get the time of each bout, to put in order for getting the shift controls
                behvBoutStartInds{iBehv} = cellfun(@(x) x(1),allNeurInds(behvStructInd(iBehv),~cellfun(@isempty,allNeurInds(behvStructInd(iBehv),:))));
                behvBoutLabels{iBehv} = cellfun(@(x) repmat(iBehv,1,size(x,2)), behvBoutFRs{iBehv}, 'un', 0);

            end

            %normalize based on activity across all behaviors
            allFRs = cat(2,regionFRsOrig{:});
            allEMGs = cat(2,regionEMGsOrig{:});

            allFRMeans = mean(allFRs,2);
            allFRStds = std(allFRs,[],2);
            allEMGMeans = mean(allEMGs,2);
            allEMGStds = std(allEMGs,[],2);

            regionFRs = cellfun(@(x) (x-allFRMeans)./allFRStds,regionFRsOrig,'un',0);
            regionEMGs = cellfun(@(x) (x-allEMGMeans)./allEMGStds,regionEMGsOrig,'un',0);

            % do shift controls to break behavior relationship with neur/emg
            % first combine all behaviors (sorting bouts by time)
            [~, sortPerm] = sort([behvBoutStartInds{:}]);
            allBoutFRs = [behvBoutFRs{:}];
            allBoutFRs = allBoutFRs(sortPerm);
            allBoutEMGs = [behvBoutEMGs{:}];
            allBoutEMGs = allBoutEMGs(sortPerm);
            allBoutLabels = [behvBoutLabels{:}];
            allBoutLabels = allBoutLabels(sortPerm);

            allBoutFRsCat = cat(2,allBoutFRs{:});
            allBoutEMGsCat = cat(2,allBoutEMGs{:});
            allBoutLabelsCat = cat(2,allBoutLabels{:});

            % don't use points with nans
            catFRNans = find(any(isnan(allBoutFRsCat)));
            catEMGNans = find(any(isnan(allBoutEMGsCat)));

            allBoutFRsCat(:,unique([catFRNans catEMGNans])) = [];
            allBoutEMGsCat(:,unique([catFRNans catEMGNans])) = [];
            allBoutLabelsCat(:,unique([catFRNans catEMGNans])) = [];

            % do shift then divide into separate behaviors again
            for iShift = 1:nShifts

                %shift by at least 30 seconds
                if use10msBins
                    shiftAmount(iShift) = randi(length(allBoutLabelsCat)-3000*2)+3000;
                else
                    shiftAmount(iShift) = randi(length(allBoutLabelsCat)-30000*2)+30000;
                end
                labelsShift = circshift(allBoutLabelsCat,shiftAmount(iShift));

                % get FRs and EMGs based on new shifted labels
                for iBehv = 1:length(behvRegionFieldNames)
                    shiftBehvInds = find(labelsShift == iBehv);
                    regionFRsShift{iBehv,iShift} = allBoutFRsCat(:,shiftBehvInds);
                    regionEMGsShift{iBehv,iShift} = allBoutEMGsCat(:,shiftBehvInds);
                end

            end

        case 'muscletrig'
            % -----Use data divded by the 7 umap regions, but only during the muscle activation periods-----
            

    end


    %get the least number of points to keep the same data size across
    %behavior regions
    leastPoints = min(cellfun(@length,regionFRs));

    %calculate average firing rates for cells within each region
    regionMeanFRsCell = cellfun(@(x) nanmean(x,2),regionFRsOrig,'UniformOutput',false);
    regionMeanFRs = cat(2,regionMeanFRsCell{:});

    %change to spks/sec rather than ms
    if use10msBins
        regionMeanFRs = regionMeanFRs * 100;
    else
        regionMeanFRs = regionMeanFRs * 1000;
    end

    %only use neurons above a certain firing rate within at least one behavior
    firingRateCutoff = 0.2;                                                                                                          
    goodNeuronsAll = find(any(regionMeanFRs > firingRateCutoff,2));

    %more strict criterion we can use - use only neurons above a certain
    %firing rate across all behaviors
    reallyGoodNeuronsAll = find(all(regionMeanFRs > firingRateCutoff,2));

    %also divide into striatum and cortex
    goodNeuronsStr = intersect(goodNeuronsAll,1:length(striatumInds));
    goodNeuronsCtx = intersect(goodNeuronsAll,length(striatumInds)+1:size(allFRs,1));

%     goodNeuronsStr = goodNeuronsStr(randperm(length(goodNeuronsStr),length(goodNeuronsCtx)));

    reallyGoodNeuronsStr = intersect(reallyGoodNeuronsAll,1:length(striatumInds));
    reallyGoodNeuronsCtx = intersect(reallyGoodNeuronsAll,length(striatumInds)+1:size(allFRs,1));

    % do PCA on all time points
    allRegionsFR = cat(2,regionFRs{:});
    allRegionsEMG = cat(2,regionEMGs{:});
    allRegionsInds = [0 cumsum(cellfun(@(x) size(x,2), regionFRs))];

    [pcaProjAllRegionsStr{iSess}, pcaTrajAllRegionsStr{iSess}, varExpAllRegionsStr{iSess}] = pca(allRegionsFR(goodNeuronsStr,:)');
    [pcaProjAllRegionsCtx{iSess}, pcaTrajAllRegionsCtx{iSess}, varExpAllRegionsCtx{iSess}] = pca(allRegionsFR(goodNeuronsCtx,:)');
    [pcaProjAllRegionsEmg{iSess}, pcaTrajAllRegionsEmg{iSess}, varExpAllRegionsEmg{iSess}] = pca(allRegionsEMG');

    %calculate dimensionality
    [lbmleDimAllRegions(1), paDimAllRegions(iSess,1), cutoffDimAllRegions(iSess,1)] = dimEst(allRegionsFR(goodNeuronsStr,:)', 5);
    [lbmleDimAllRegions(2), paDimAllRegions(iSess,2), cutoffDimAllRegions(iSess,2)] = dimEst(allRegionsFR(goodNeuronsCtx,:)', 5);
    [lbmleDimAllRegions(3), paDimAllRegions(iSess,3), cutoffDimAllRegions(iSess,3)] = dimEst(allRegionsEMG', 5);

%     % do CCA on all time points
%     [ccaNeurProjCombRegStr{iSess},ccaEMGProjCombRegStr{iSess},cannonCorrsCombRegStr(iSess,:),ccaNeurTrajCombRegStr{iSess},ccaEMGTrajCombRegStr{iSess}] = ...
%         canoncorr(allRegionsFR(goodNeuronsStr,:)',allRegionsEMG(1:4,:)');
%     [ccaNeurProjCombRegCtx{iSess},ccaEMGProjCombRegCtx{iSess},cannonCorrsCombRegCtx(iSess,:),ccaNeurTrajCombRegCtx{iSess},ccaEMGTrajCombRegCtx{iSess}] = ...
%         canoncorr(allRegionsFR(goodNeuronsCtx,:)',allRegionsEMG(1:4,:)');
% 
%     %do CCA using only top PC's
%     [ccaNeurProjCombRegLowDimStr{iSess},ccaEMGProjCombRegLowDimStr{iSess},cannonCorrsCombRegLowDimStr(iSess,:),...
%         ccaNeurTrajCombRegLowDimStr{iSess},ccaEMGTrajCombRegLowDimStr{iSess}] = ...
%         canoncorr(pcaTrajAllRegionsStr{iSess}(:,1:max(paDimAllRegions(iSess,1),4)),allRegionsEMG(1:4,:)');
%     [ccaNeurProjCombRegLowDimCtx{iSess},ccaEMGProjCombRegLowDimCtx{iSess},cannonCorrsCombRegLowDimCtx(iSess,:),...
%         ccaNeurTrajCombRegLowDimCtx{iSess},ccaEMGTrajCombRegLowDimCtx{iSess}] = ...
%         canoncorr(pcaTrajAllRegionsCtx{iSess}(:,1:max(paDimAllRegions(iSess,2),4)),allRegionsEMG(1:4,:)');
% 
%     %shift EMG and neural data relative to each other for controls for CCA
%     for iShift = 1:100
% 
%         %shift by at least 5 seconds
%         if use10msBins
%             emgShiftAmount(iShift) = randi(size(allRegionsEMG,2)-500*2)+500;
%         else
%             emgShiftAmount(iShift) = randi(size(allRegionsEMG,2)-5000*2)+5000;
%         end
% 
%         shiftEMG = circshift(allRegionsEMG(1:4,:),emgShiftAmount(iShift),2)';
%         [~,~,cannonCorrsCombRegStrShift(iSess,iShift,:)] = canoncorr(allRegionsFR(goodNeuronsStr,:)',shiftEMG);
%         [~,~,cannonCorrsCombRegCtxShift(iSess,iShift,:)] = canoncorr(allRegionsFR(goodNeuronsCtx,:)',shiftEMG);
% 
%         [~,~,cannonCorrsCombRegLowDimStrShift(iSess,iShift,:)] = canoncorr(pcaTrajAllRegionsStr{iSess}(:,1:max(paDimAllRegions(iSess,1),4)),shiftEMG);
%         [~,~,cannonCorrsCombRegLowDimCtxShift(iSess,iShift,:)] = canoncorr(pcaTrajAllRegionsCtx{iSess}(:,1:max(paDimAllRegions(iSess,2),4)),shiftEMG);
% 
%     end

    %Free up some memory
    clear allFRs
    clear normalizedFRs
    clear regionFRsOrig
    clear allRegionsFR allRegionsEMG

    %make pseudo-data, take random subset, of the top 50 PC's,
    %rotate it then project back up to
%     for iShuff = 1:nShifts
%     
%         for iRegion1 = 1:length(regionWatershedLabels)
%             for iRegion2 = 1:length(regionWatershedLabels)
% 
%                 %do random rotation
%     %             randRotMat = randn(length(goodNeuronsStr),alignmentDim);
%                 nTotNeurons = size(pcaTrajCombRegionsAll{iSess},2);
%                 randRotMat1 = randn(nTotNeurons,nTotNeurons);
%                 randRotMat1 = orth(randRotMat1);
%                 randRotMat2 = randn(nTotNeurons,nTotNeurons);
%                 randRotMat2 = orth(randRotMat2);
%     
%                 psuedoData1Inds = randperm(size(pcaTrajCombRegionsAll{iSess},1), leastPoints);
%                 psuedoData2Inds = randperm(size(pcaTrajCombRegionsAll{iSess},1), leastPoints);
%     
%                 psuedoData1 = pcaTrajCombRegionsAll{iSess}(psuedoData1Inds,1:nTotNeurons)*randRotMat1*pcaProjCombRegionsAll{iSess}(:,1:nTotNeurons)';
%                 psuedoData2 = pcaTrajCombRegionsAll{iSess}(psuedoData2Inds,1:nTotNeurons)*randRotMat2*pcaProjCombRegionsAll{iSess}(:,1:nTotNeurons)';
%     
%                 [psuedoProj1, ~, psuedoVarExp1] = pca(psuedoData1);
%                 [psuedoProj2, ~, psuedoVarExp2] = pca(psuedoData2);
%     
%                 pcaNeurAligment10DimAll(iRegion1,iRegion2,iShuff) = calcPCAAlignment(psuedoData1,psuedoData2,psuedoProj1,psuedoProj2,10);
%     
%             end
%         end
%     
%     end

    % now go through each region and run PCA and CCA
    
    for iRegion = 1:nRegions

        tic

        %downsample so same number of points in all regions
        timeIndsToUse{iSess,iRegion} = randperm(size(regionFRs{iRegion},2),leastPoints);

        % do PCA
        [pcaProjStr{iSess,iRegion}, pcaTrajStr{iRegion}, varExpStr{iSess,iRegion}] = pca(regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})');
        [pcaProjCtx{iSess,iRegion}, pcaTrajCtx{iRegion}, varExpCtx{iSess,iRegion}] = pca(regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})');
        [pcaProjEmg{iSess,iRegion}, pcaTrajEmg{iRegion}, varExpEmg{iSess,iRegion}] = pca(regionEMGs{iRegion}(:,timeIndsToUse{iSess,iRegion})');

        %calculate dimensionality
        [lbmleDim(iRegion,1), paDim(iSess,iRegion,1), cutoffDim(iSess,iRegion,1)] = dimEst(regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})', 5);
        [lbmleDim(iRegion,2), paDim(iSess,iRegion,2), cutoffDim(iSess,iRegion,2)] = dimEst(regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})', 5);
        [lbmleDim(iRegion,3), paDim(iSess,iRegion,3), cutoffDim(iSess,iRegion,3)] = dimEst(regionEMGs{iRegion}(:,timeIndsToUse{iSess,iRegion})', 5);

%         %do CCA
%         [ccaNeurProjStr{iRegion},ccaEMGProjStr{iRegion},cannonCorrsStr(iSess,iRegion,:),ccaNeurTrajStr{iRegion},ccaEMGTrajStr{iRegion}] = ...
%             canoncorr(regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})',regionEMGs{iRegion}(1:4,timeIndsToUse{iSess,iRegion})');
%         [ccaNeurProjCtx{iRegion},ccaEMGProjCtx{iRegion},cannonCorrsCtx(iSess,iRegion,:),ccaNeurTrajCtx{iRegion},ccaEMGTrajCtx{iRegion}] = ...
%             canoncorr(regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})',regionEMGs{iRegion}(1:4,timeIndsToUse{iSess,iRegion})');
% 
%         %do CCA using only top PC's
%         [ccaNeurProjLowDimStr{iRegion},ccaEMGProjLowDimStr{iRegion},cannonCorrsLowDimStr(iSess,iRegion,:),ccaNeurTrajLowDimStr{iRegion},ccaEMGTrajLowDimStr{iRegion}] = ...
%             canoncorr(pcaTrajStr{iRegion}(:,1:max(paDim(iSess,iRegion,1),4)),regionEMGs{iRegion}(1:4,timeIndsToUse{iSess,iRegion})');
%         [ccaNeurProjLowDimCtx{iRegion},ccaEMGProjLowDimCtx{iRegion},cannonCorrsLowDimCtx(iSess,iRegion,:),ccaNeurTrajLowDimCtx{iRegion},ccaEMGTrajLowDimCtx{iRegion}] = ...
%             canoncorr(pcaTrajCtx{iRegion}(:,1:max(paDim(iSess,iRegion,2),4)),regionEMGs{iRegion}(1:4,timeIndsToUse{iSess,iRegion})');

        % do PCA with the behavior shift controls
        for iShift = 1:nShifts
            
            shiftLeastPoints = min(cellfun(@(x) size(x,2),regionFRsShift(:,iShift)));
            shiftTimeIndsToUse{iSess,iRegion,iShift} = randperm(size(regionFRsShift{iRegion,iShift},2),shiftLeastPoints);
            [pcaProjStrShift{iSess,iRegion,iShift}, ~, ~] = ...
                pca(regionFRsShift{iRegion,iShift}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iRegion,iShift})');
            [pcaProjCtxShift{iSess,iRegion,iShift}, ~, ~] = ...
                pca(regionFRsShift{iRegion,iShift}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iRegion,iShift})');
            [pcaProjEmgShift{iSess,iRegion,iShift}, ~, ~] = ...
                pca(regionEMGsShift{iRegion,iShift}(:,shiftTimeIndsToUse{iSess,iRegion,iShift})');
            
            [lbmleDimShift(iRegion,iShift,1), paDimShift(iSess,iRegion,iShift,1), cutoffDimShift(iSess,iRegion,iShift,1)] = ...
                dimEst(regionFRsShift{iRegion,iShift}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iRegion,iShift})', 5);
            [lbmleDimShift(iRegion,iShift,2), paDimShift(iSess,iRegion,iShift,2), cutoffDimShift(iSess,iRegion,iShift,2)] = ...
                dimEst(regionFRsShift{iRegion,iShift}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iRegion,iShift})', 5);
            [lbmleDimShift(iRegion,iShift,3), paDimShift(iSess,iRegion,iShift,3), cutoffDimShift(iSess,iRegion,iShift,3)] = ...
                dimEst(regionEMGsShift{iRegion,iShift}(:,shiftTimeIndsToUse{iSess,iRegion,iShift})', 5);
            
        end
        
%         %shift EMG and neural data relative to each other for controls for CCA
%         for iShift = 1:100
% 
%             %shift by at least 5 seconds
%             if use10msBins
%                 emgShiftAmount(iShift) = randi(size(regionEMGs{iRegion},2)-500*2)+500;
%             else
%                 emgShiftAmount(iShift) = randi(size(regionEMGs{iRegion},2)-5000*2)+5000;
%             end
% 
%             shiftEMG = circshift(regionEMGs{iRegion}(1:4,:),emgShiftAmount(iShift),2)';
%             [~,~,cannonCorrsStrShift(iSess,iRegion,iShift,:)] = canoncorr(regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})',shiftEMG(timeIndsToUse{iSess,iRegion},:));
%             [~,~,cannonCorrsCtxShift(iSess,iRegion,iShift,:)] = canoncorr(regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})',shiftEMG(timeIndsToUse{iSess,iRegion},:));
% 
%             [~,~,cannonCorrsLowDimStrShift(iSess,iRegion,iShift,:)] = canoncorr(pcaTrajStr{iRegion}(:,1:max(paDim(iSess,iRegion,1),4)),shiftEMG(timeIndsToUse{iSess,iRegion},:));
%             [~,~,cannonCorrsLowDimCtxShift(iSess,iRegion,iShift,:)] = canoncorr(pcaTrajCtx{iRegion}(:,1:max(paDim(iSess,iRegion,2),4)),shiftEMG(timeIndsToUse{iSess,iRegion},:));
%             
%         end

%         %do just correlations between neurons and EMG
%         corrEMGAll{iRegion} = corr(regionFRs{iRegion}(reallyGoodNeuronsAll,timeIndsToUse{iRegion})',regionEMGs{iRegion}(:,timeIndsToUse{iRegion})');
%         corrEMGStr{iRegion} = corr(regionFRs{iRegion}(reallyGoodNeuronsStr,timeIndsToUse{iRegion})',regionEMGs{iRegion}(:,timeIndsToUse{iRegion})');
%         corrEMGCtx{iRegion} = corr(regionFRs{iRegion}(reallyGoodNeuronsCtx,timeIndsToUse{iRegion})',regionEMGs{iRegion}(:,timeIndsToUse{iRegion})');

        %do pairwise neuron correlations in addition to PCA
%         corrMatrixStr{iRegion} = corr(regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})',regionFRs{iRegion}(goodNeuronsStr,timeIndsToUse{iSess,iRegion})');
%         corrMatrixStr{iRegion}(logical(diag(ones(length(corrMatrixStr{iRegion}),1),0))) = 0;
%         corrMatrixCtx{iRegion} = corr(regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})',regionFRs{iRegion}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion})');
%         corrMatrixCtx{iRegion}(logical(diag(ones(length(corrMatrixCtx{iRegion}),1),0))) = 0;
%         corrMatrixEmg{iRegion} = corr(regionEMGs{iRegion}(:,timeIndsToUse{iSess,iRegion})',regionEMGs{iRegion}(:,timeIndsToUse{iSess,iRegion})');
%         corrMatrixEmg{iRegion}(logical(diag(ones(length(corrMatrixEmg{iRegion}),1),0))) = 0;

        disp(['Get region activity: ' num2str(toc)])
    end

    alignmentDim = 10;
    % now for each pair of subregions, do PCA and CCA and get alignment
    for iRegion1 = 1:nRegions
        for iRegion2 = 1:nRegions

            if iRegion1 < iRegion2 %since the metrics are symmetric, save time by just calculating upper triangle

                tic

%                 [ccaPrinAngleCtx(iRegion1,iRegion2), ccaAlignmentCtx(iRegion1,iRegion2)] = calcCCAAlignmentMetrics(...
%                     regionFRs{iRegion1}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion1}),regionFRs{iRegion2}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion2}),...
%                     regionEMGs{iRegion1}(:,timeIndsToUse{iSess,iRegion1}),regionEMGs{iRegion2}(:,timeIndsToUse{iSess,iRegion2}),...
%                     ccaNeurProjCtx{iRegion1},ccaNeurProjCtx{iRegion2},...
%                     ccaEMGProjCtx{iRegion1},ccaEMGProjCtx{iRegion2},cannonCorrsCtx{iRegion1},cannonCorrsCtx{iRegion2});
% 
%                 [ccaPrinAngleStr(iRegion1,iRegion2), ccaAlignmentStr(iRegion1,iRegion2)] = calcCCAAlignmentMetrics(...
%                     regionFRs{iRegion1}(goodNeuronsStr,timeIndsToUse{iSess,iRegion1}),regionFRs{iRegion2}(goodNeuronsStr,timeIndsToUse{iSess,iRegion2}),...
%                     regionEMGs{iRegion1}(:,timeIndsToUse{iSess,iRegion1}),regionEMGs{iRegion2}(:,timeIndsToUse{iSess,iRegion2}),...
%                     ccaNeurProjStr{iRegion1},ccaNeurProjStr{iRegion2},...
%                     ccaEMGProjStr{iRegion1},ccaEMGProjStr{iRegion2},cannonCorrsStr{iRegion1},cannonCorrsStr{iRegion2});

                %               ccaNeurAligment(iRegion1,iRegion2) = calcPCAAlignment(regionFRs{iRegion1}(neurIndsToUse,:)',...
                %                   regionFRs{iRegion2}(neurIndsToUse,:)',ccaNeurProj{iRegion1},ccaNeurProj{iRegion2},4);
                %               pcaNeurAligment10DimAll(iRegion1,iRegion2) = calcPCAAlignment(regionFRs{iRegion1}(goodNeuronsAll,timeIndsToUse{iRegion1})',...
                %                   regionFRs{iRegion2}(goodNeuronsAll,timeIndsToUse{iRegion2})',pcaProjAll{iRegion1},pcaProjAll{iRegion2},'auto');

                [pcaDiverStr(iSess,iRegion1,iRegion2), pcaAlignStr(iSess,iRegion1,iRegion2), pcaAngleStr(iSess,iRegion1,iRegion2)] = ...
                    calcPCAAlignment(regionFRs{iRegion1}(goodNeuronsStr,timeIndsToUse{iSess,iRegion1})',...
                    regionFRs{iRegion2}(goodNeuronsStr,timeIndsToUse{iSess,iRegion2})',pcaProjStr{iSess,iRegion1},pcaProjStr{iSess,iRegion2},...
                    paDim(iSess,iRegion1,1),paDim(iSess,iRegion2,1));

                [pcaDiverCtx(iSess,iRegion1,iRegion2), pcaAlignCtx(iSess,iRegion1,iRegion2), pcaAngleCtx(iSess,iRegion1,iRegion2)] = ...
                    calcPCAAlignment(regionFRs{iRegion1}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion1})',...
                    regionFRs{iRegion2}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion2})',pcaProjCtx{iSess,iRegion1},pcaProjCtx{iSess,iRegion2},...
                    paDim(iSess,iRegion1,2),paDim(iSess,iRegion2,2));

                [pcaDiverEmg(iSess,iRegion1,iRegion2), pcaAlignEmg(iSess,iRegion1,iRegion2), pcaAngleEmg(iSess,iRegion1,iRegion2)] = ...
                    calcPCAAlignment(regionEMGs{iRegion1}(:,timeIndsToUse{iSess,iRegion1})',...
                    regionEMGs{iRegion2}(:,timeIndsToUse{iSess,iRegion2})',pcaProjEmg{iSess,iRegion1},pcaProjEmg{iSess,iRegion2},...
                    paDim(iSess,iRegion1,3),paDim(iSess,iRegion2,3));


                %do shuff                
                allRegDataStr = [regionFRs{iRegion1}(goodNeuronsStr,timeIndsToUse{iSess,iRegion1}) regionFRs{iRegion2}(goodNeuronsStr,timeIndsToUse{iSess,iRegion2})];
                allRegDataCtx = [regionFRs{iRegion1}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion1}) regionFRs{iRegion2}(goodNeuronsCtx,timeIndsToUse{iSess,iRegion2})];
                allRegDataEmg = [regionEMGs{iRegion1}(:,timeIndsToUse{iSess,iRegion1}) regionEMGs{iRegion2}(:,timeIndsToUse{iSess,iRegion2})];
                for iShuff = 1:nShifts

                    %for orth shuff, just do random rotation
                    individualProjsRot = 1;

                    if individualProjsRot
                        rotTrajStr1 = pcaTrajStr{iRegion1};
                        rotTrajStr2 = pcaTrajStr{iRegion2};
                        rotTrajCtx1 = pcaTrajCtx{iRegion1};
                        rotTrajCtx2 = pcaTrajCtx{iRegion2};
                        rotTrajEmg1 = pcaTrajEmg{iRegion1};
                        rotTrajEmg2 = pcaTrajEmg{iRegion2};

                        rotProjStr1 = pcaProjStr{iRegion1};
                        rotProjStr2 = pcaProjStr{iRegion2};
                        rotProjCtx1 = pcaProjCtx{iRegion1};
                        rotProjCtx2 = pcaProjCtx{iRegion2};
                        rotProjEmg1 = pcaProjEmg{iRegion1};
                        rotProjEmg2 = pcaProjEmg{iRegion2};
                    else
                        %                       rotTrajStr1 = pcaTrajAllRegionsStr(allRegionsInds(iRegion1)+1 : allRegionsInds(iRegion1+1),:);
                        %                       rotTrajStr2 = pcaTrajAllRegionsStr(allRegionsInds(iRegion2)+1 : allRegionsInds(iRegion2+1),:);
                        %                       rotTrajCtx1 = pcaTrajAllRegionsCtx(allRegionsInds(iRegion1)+1 : allRegionsInds(iRegion1+1),:);
                        %                       rotTrajCtx2 = pcaTrajAllRegionsCtx(allRegionsInds(iRegion2)+1 : allRegionsInds(iRegion2+1),:);
                        %                       rotTrajEmg1 = pcaTrajAllRegionsEmg(allRegionsInds(iRegion1)+1 : allRegionsInds(iRegion1+1),:);
                        %                       rotTrajEmg2 = pcaTrajAllRegionsEmg(allRegionsInds(iRegion2)+1 : allRegionsInds(iRegion2+1),:);

                        randInds1 = randperm(size(pcaTrajAllRegionsStr{iSess},1),leastPoints);
                        randInds2 = randperm(size(pcaTrajAllRegionsStr{iSess},1),leastPoints);

                        rotTrajStr1 = pcaTrajAllRegionsStr{iSess}(randInds1,:);
                        rotTrajStr2 = pcaTrajAllRegionsStr{iSess}(randInds2,:);
                        rotTrajCtx1 = pcaTrajAllRegionsCtx{iSess}(randInds1,:);
                        rotTrajCtx2 = pcaTrajAllRegionsCtx{iSess}(randInds2,:);
                        rotTrajEmg1 = pcaTrajAllRegionsEmg{iSess}(randInds1,:);
                        rotTrajEmg2 = pcaTrajAllRegionsEmg{iSess}(randInds2,:);

                        rotProjStr1 = pcaProjAllRegionsStr{iSess};
                        rotProjStr2 = pcaProjAllRegionsStr{iSess};
                        rotProjCtx1 = pcaProjAllRegionsCtx{iSess};
                        rotProjCtx2 = pcaProjAllRegionsCtx{iSess};
                        rotProjEmg1 = pcaProjAllRegionsEmg{iSess};
                        rotProjEmg2 = pcaProjAllRegionsEmg{iSess};

                    end

                    if iShuff == 1
                        %since calculating the pa dimensionality takes so
                        %long, just calculate it once for the first
                        %shuffle, since it doesn't change much across
                        %shuffles anyways
                        [pcaDiverRotStr(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotStr(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotStr(iSess,iRegion1,iRegion2,iShuff), rotPaDim1Str, rotPaDim2Str] = ...
                            randomRotationAlignment(rotTrajStr1, rotTrajStr2, rotProjStr1, rotProjStr2,[],[]);

                        [pcaDiverRotCtx(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotCtx(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotCtx(iSess,iRegion1,iRegion2,iShuff), rotPaDim1Ctx, rotPaDim2Ctx] = ...
                            randomRotationAlignment(rotTrajCtx1, rotTrajCtx2, rotProjCtx1, rotProjCtx2,[],[]);

                        [pcaDiverRotEmg(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotEmg(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotEmg(iSess,iRegion1,iRegion2,iShuff), rotPaDim1Emg, rotPaDim2Emg] = ...
                            randomRotationAlignment(rotTrajEmg1, rotTrajEmg2, rotProjEmg1, rotProjEmg2,[],[]);
                    else
                        [pcaDiverRotStr(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotStr(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotStr(iSess,iRegion1,iRegion2,iShuff)] = ...
                            randomRotationAlignment(rotTrajStr1, rotTrajStr2, rotProjStr1, rotProjStr2,rotPaDim1Str,rotPaDim2Str);

                        [pcaDiverRotCtx(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotCtx(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotCtx(iSess,iRegion1,iRegion2,iShuff)] = ...
                            randomRotationAlignment(rotTrajCtx1, rotTrajCtx2, rotProjCtx1, rotProjCtx2,rotPaDim1Ctx,rotPaDim2Ctx);

                        [pcaDiverRotEmg(iSess,iRegion1,iRegion2,iShuff), pcaAlignRotEmg(iSess,iRegion1,iRegion2,iShuff), pcaAngleRotEmg(iSess,iRegion1,iRegion2,iShuff)] = ...
                            randomRotationAlignment(rotTrajEmg1, rotTrajEmg2, rotProjEmg1, rotProjEmg2,rotPaDim1Emg,rotPaDim2Emg);
                    end

                    % now for behavioral label shifts
                    [pcaDiverShiftStr(iSess,iRegion1,iRegion2,iShuff), pcaAlignShiftStr(iSess,iRegion1,iRegion2,iShuff), pcaAngleShiftStr(iSess,iRegion1,iRegion2,iShuff)] = ...
                        calcPCAAlignment(regionFRsShift{iRegion1,iShuff}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iRegion1,iShuff})',...
                        regionFRsShift{iRegion2,iShuff}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iRegion2,iShuff})',pcaProjStrShift{iSess,iRegion1,iShuff},...
                        pcaProjStrShift{iSess,iRegion2,iShuff},paDimShift(iSess,iRegion1,iShuff,1),paDimShift(iSess,iRegion2,iShuff,1));

                    [pcaDiverShiftCtx(iSess,iRegion1,iRegion2,iShuff), pcaAlignShiftCtx(iSess,iRegion1,iRegion2,iShuff), pcaAngleShiftCtx(iSess,iRegion1,iRegion2,iShuff)] = ...
                        calcPCAAlignment(regionFRsShift{iRegion1,iShuff}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iRegion1,iShuff})',...
                        regionFRsShift{iRegion2,iShuff}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iRegion2,iShuff})',pcaProjCtxShift{iSess,iRegion1,iShuff},...
                        pcaProjCtxShift{iSess,iRegion2,iShuff},paDimShift(iSess,iRegion1,iShuff,2),paDimShift(iSess,iRegion2,iShuff,2));

                    [pcaDiverShiftEmg(iSess,iRegion1,iRegion2,iShuff), pcaAlignShiftEmg(iSess,iRegion1,iRegion2,iShuff), pcaAngleShiftEmg(iSess,iRegion1,iRegion2,iShuff)] = ...
                        calcPCAAlignment(regionEMGsShift{iRegion1,iShuff}(:,shiftTimeIndsToUse{iSess,iRegion1,iShuff})',...
                        regionEMGsShift{iRegion2,iShuff}(:,shiftTimeIndsToUse{iSess,iRegion2,iShuff})',pcaProjEmgShift{iSess,iRegion1,iShuff},...
                        pcaProjEmgShift{iSess,iRegion2,iShuff},paDimShift(iSess,iRegion1,iShuff,3),paDimShift(iSess,iRegion2,iShuff,3));


%                     %permutation shuff
%                     shuffPerm = randperm(leastPoints*2);
%                     [pcaProjShuff1Str, pcaTrajShuff1Str, varExpShuff1Str] = pca(allRegDataStr(:,shuffPerm(1:leastPoints))');
%                     [pcaProjShuff2Str, pcaTrajShuff2Str, varExpShuff2Str] = pca(allRegDataStr(:,shuffPerm(leastPoints+1:end))');
% 
%                     pcaNeurAligment10DimShuffStr(iRegion1,iRegion2,iShuff) = calcPCAAlignment(allRegDataStr(:,shuffPerm(1:leastPoints))',...
%                         allRegDataStr(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Str,pcaProjShuff2Str,'auto');
% 
%                     [pcaProjShuff1Ctx, pcaTrajShuff1Ctx, varExpShuff1Ctx] = pca(allRegDataCtx(:,shuffPerm(1:leastPoints))');
%                     [pcaProjShuff2Ctx, pcaTrajShuff2Ctx, varExpShuff2Ctx] = pca(allRegDataCtx(:,shuffPerm(leastPoints+1:end))');
% 
%                     pcaNeurAligment10DimShuffCtx(iRegion1,iRegion2,iShuff) = calcPCAAlignment(allRegDataCtx(:,shuffPerm(1:leastPoints))',...
%                         allRegDataCtx(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Ctx,pcaProjShuff2Ctx,'auto');
% 
%                     [pcaProjShuff1Emg, pcaTrajShuff1Emg, varExpShuff1Emg] = pca(allRegDataEmg(:,shuffPerm(1:leastPoints))');
%                     [pcaProjShuff2Emg, pcaTrajShuff2Emg, varExpShuff2Emg] = pca(allRegDataEmg(:,shuffPerm(leastPoints+1:end))');
% 
%                     pcaNeurAligment10DimShuffEmg(iRegion1,iRegion2,iShuff) = calcPCAAlignment(allRegDataEmg(:,shuffPerm(1:leastPoints))',...
%                         allRegDataEmg(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Emg,pcaProjShuff2Emg,'auto');
% 
% 
%                     [ccaNeurProjShuff1Ctx,ccaEMGProjShuff1Ctx,cannonCorrsShuff1Ctx,ccaNeurTrajShuff1Ctx,ccaEMGTrajShuff1Ctx] = ...
%                         canoncorr(allRegDataCtx(:,shuffPerm(1:leastPoints))',allRegDataEmg(:,shuffPerm(1:leastPoints))');
%                     [ccaNeurProjShuff2Ctx,ccaEMGProjShuff2Ctx,cannonCorrsShuff2Ctx,ccaNeurTrajShuff2Ctx,ccaEMGTrajShuff2Ctx] = ...
%                         canoncorr(allRegDataCtx(:,shuffPerm(leastPoints+1:end))',allRegDataEmg(:,shuffPerm(leastPoints+1:end))');
% 
%                     [ccaPrinAngleShuffCtx(iRegion1,iRegion2,iShuff), ccaAlignmentShuffCtx(iRegion1,iRegion2,iShuff)] = calcCCAAlignmentMetrics(...
%                         allRegDataCtx(:,shuffPerm(1:leastPoints)),allRegDataCtx(:,shuffPerm(leastPoints+1:end)),...
%                         allRegDataEmg(:,shuffPerm(1:leastPoints)),allRegDataEmg(:,shuffPerm(leastPoints+1:end)),...
%                         ccaNeurProjShuff1Ctx,ccaNeurProjShuff2Ctx,...
%                         ccaEMGProjShuff1Ctx,ccaEMGProjShuff2Ctx,cannonCorrsShuff1Ctx,cannonCorrsShuff2Ctx);
% 
%                     [ccaNeurProjShuff1Str,ccaEMGProjShuff1Str,cannonCorrsShuff1Str,ccaNeurTrajShuff1Str,ccaEMGTrajShuff1Str] = ...
%                         canoncorr(allRegDataStr(:,shuffPerm(1:leastPoints))',allRegDataEmg(:,shuffPerm(1:leastPoints))');
%                     [ccaNeurProjShuff2Str,ccaEMGProjShuff2Str,cannonCorrsShuff2Str,ccaNeurTrajShuff2Str,ccaEMGTrajShuff2Str] = ...
%                         canoncorr(allRegDataStr(:,shuffPerm(leastPoints+1:end))',allRegDataEmg(:,shuffPerm(leastPoints+1:end))');
% 
%                     [ccaPrinAngleShuffStr(iRegion1,iRegion2,iShuff), ccaAlignmentShuffStr(iRegion1,iRegion2,iShuff)] = calcCCAAlignmentMetrics(...
%                         allRegDataStr(:,shuffPerm(1:leastPoints)),allRegDataStr(:,shuffPerm(leastPoints+1:end)),...
%                         allRegDataEmg(:,shuffPerm(1:leastPoints)),allRegDataEmg(:,shuffPerm(leastPoints+1:end)),...
%                         ccaNeurProjShuff1Str,ccaNeurProjShuff2Str,...
%                         ccaEMGProjShuff1Str,ccaEMGProjShuff2Str,cannonCorrsShuff1Str,cannonCorrsShuff2Str);
% 
% 
%                     reg1Proj = pcaProjCtx{iRegion1}(:);
%                     reg1Proj = reg1Proj(randperm(length(reg1Proj),length(reg1Proj)));
%                     reg1Proj = reshape(reg1Proj,size(pcaProjCtx{iRegion1},1),size(pcaProjCtx{iRegion1},2));
% 
%                     reg2Proj = pcaProjCtx{iRegion2}(:);
%                     reg2Proj = reg2Proj(randperm(length(reg2Proj),length(reg2Proj)));
%                     reg2Proj = reshape(reg2Proj,size(pcaProjCtx{iRegion2},1),size(pcaProjCtx{iRegion2},2));
%                     reg1Proj = pcaProjCtx{iRegion1}(:,randperm(size(pcaProjCtx{iRegion1},1),size(pcaProjCtx{iRegion1},1)));
%                     reg2Proj = pcaProjCtx{iRegion2}(:,randperm(size(pcaProjCtx{iRegion1},1),size(pcaProjCtx{iRegion1},1)));
% 
%                     pcaNeurAligment10DimCtxShuff(iRegion1,iRegion2,iShuff) = calcPCAAlignment(regionFRs{iRegion1}(goodNeuronsCtx,timeIndsToUse{iRegion1})',...
%                         regionFRs{iRegion2}(goodNeuronsCtx,timeIndsToUse{iRegion2})',reg1Proj,reg2Proj,alignmentDim);
                end

%                 %do single neuron correlations
%                 corrAlignmentStr(iRegion1,iRegion2) = corr(squareform(corrMatrixStr{iRegion1})',squareform(corrMatrixStr{iRegion2})');
%                 corrAlignmentCtx(iRegion1,iRegion2) = corr(squareform(corrMatrixCtx{iRegion1})',squareform(corrMatrixCtx{iRegion2})');
%                 corrAlignmentEmg(iRegion1,iRegion2) = corr(squareform(corrMatrixEmg{iRegion1})',squareform(corrMatrixEmg{iRegion2})');

            elseif iRegion1 == iRegion2

                % same behaviors set to 0 divergence/1 alignemnt
                pcaDiverStr(iSess,iRegion1,iRegion2) = 0; pcaDiverCtx(iSess,iRegion1,iRegion2) = 0; pcaDiverEmg(iSess,iRegion1,iRegion2) = 0;
                pcaAlignStr(iSess,iRegion1,iRegion2) = 1; pcaAlignCtx(iSess,iRegion1,iRegion2) = 1; pcaAlignEmg(iSess,iRegion1,iRegion2) = 1;
                pcaAngleStr(iSess,iRegion1,iRegion2) = 0; pcaAngleCtx(iSess,iRegion1,iRegion2) = 0; pcaAngleEmg(iSess,iRegion1,iRegion2) = 0;

                pcaDiverRotStr(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaDiverRotCtx(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaDiverRotEmg(iSess,iRegion1,iRegion2,1:nShifts) = 0;
                pcaAlignRotStr(iSess,iRegion1,iRegion2,1:nShifts) = 1; pcaAlignRotCtx(iSess,iRegion1,iRegion2,1:nShifts) = 1; pcaAlignRotEmg(iSess,iRegion1,iRegion2,1:nShifts) = 1;
                pcaAngleRotStr(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaAngleRotCtx(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaAngleRotEmg(iSess,iRegion1,iRegion2,1:nShifts) = 0;

                pcaDiverShiftStr(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaDiverShiftCtx(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaDiverShiftEmg(iSess,iRegion1,iRegion2,1:nShifts) = 0;
                pcaAlignShiftStr(iSess,iRegion1,iRegion2,1:nShifts) = 1; pcaAlignShiftCtx(iSess,iRegion1,iRegion2,1:nShifts) = 1; pcaAlignShiftEmg(iSess,iRegion1,iRegion2,1:nShifts) = 1;
                pcaAngleShiftStr(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaAngleShiftCtx(iSess,iRegion1,iRegion2,1:nShifts) = 0; pcaAngleShiftEmg(iSess,iRegion1,iRegion2,1:nShifts) = 0;

            end

            disp(['Calc Alignment: ' num2str(toc)])

        end
    end

    %fill in lower triangle of the matrix
    pcaDiverStr(iSess,:,:) = triu(squeeze(pcaDiverStr(iSess,:,:)),1)' + triu(squeeze(pcaDiverStr(iSess,:,:)));
    pcaDiverCtx(iSess,:,:) = triu(squeeze(pcaDiverCtx(iSess,:,:)),1)' + triu(squeeze(pcaDiverCtx(iSess,:,:)));
    pcaDiverEmg(iSess,:,:) = triu(squeeze(pcaDiverEmg(iSess,:,:)),1)' + triu(squeeze(pcaDiverEmg(iSess,:,:)));

    pcaAlignStr(iSess,:,:) = triu(squeeze(pcaAlignStr(iSess,:,:)),1)' + triu(squeeze(pcaAlignStr(iSess,:,:)));
    pcaAlignCtx(iSess,:,:) = triu(squeeze(pcaAlignCtx(iSess,:,:)),1)' + triu(squeeze(pcaAlignCtx(iSess,:,:)));
    pcaAlignEmg(iSess,:,:) = triu(squeeze(pcaAlignEmg(iSess,:,:)),1)' + triu(squeeze(pcaAlignEmg(iSess,:,:)));

    pcaAngleStr(iSess,:,:) = triu(squeeze(pcaAngleStr(iSess,:,:)),1)' + triu(squeeze(pcaAngleStr(iSess,:,:)));
    pcaAngleCtx(iSess,:,:) = triu(squeeze(pcaAngleCtx(iSess,:,:)),1)' + triu(squeeze(pcaAngleCtx(iSess,:,:)));
    pcaAngleEmg(iSess,:,:) = triu(squeeze(pcaAngleEmg(iSess,:,:)),1)' + triu(squeeze(pcaAngleEmg(iSess,:,:)));

    for iShuff = 1:nShifts

        pcaDiverRotStr(iSess,:,:,iShuff) = triu(squeeze(pcaDiverRotStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotStr(iSess,:,:,iShuff)));
        pcaDiverRotCtx(iSess,:,:,iShuff) = triu(squeeze(pcaDiverRotCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotCtx(iSess,:,:,iShuff)));
        pcaDiverRotEmg(iSess,:,:,iShuff) = triu(squeeze(pcaDiverRotEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotEmg(iSess,:,:,iShuff)));

        pcaAlignRotStr(iSess,:,:,iShuff) = triu(squeeze(pcaAlignRotStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotStr(iSess,:,:,iShuff)));
        pcaAlignRotCtx(iSess,:,:,iShuff) = triu(squeeze(pcaAlignRotCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotCtx(iSess,:,:,iShuff)));
        pcaAlignRotEmg(iSess,:,:,iShuff) = triu(squeeze(pcaAlignRotEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotEmg(iSess,:,:,iShuff)));

        pcaAngleRotStr(iSess,:,:,iShuff) = triu(squeeze(pcaAngleRotStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotStr(iSess,:,:,iShuff)));
        pcaAngleRotCtx(iSess,:,:,iShuff) = triu(squeeze(pcaAngleRotCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotCtx(iSess,:,:,iShuff)));
        pcaAngleRotEmg(iSess,:,:,iShuff) = triu(squeeze(pcaAngleRotEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotEmg(iSess,:,:,iShuff)));

        pcaDiverShiftStr(iSess,:,:,iShuff) = triu(squeeze(pcaDiverShiftStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftStr(iSess,:,:,iShuff)));
        pcaDiverShiftCtx(iSess,:,:,iShuff) = triu(squeeze(pcaDiverShiftCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftCtx(iSess,:,:,iShuff)));
        pcaDiverShiftEmg(iSess,:,:,iShuff) = triu(squeeze(pcaDiverShiftEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftEmg(iSess,:,:,iShuff)));

        pcaAlignShiftStr(iSess,:,:,iShuff) = triu(squeeze(pcaAlignShiftStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftStr(iSess,:,:,iShuff)));
        pcaAlignShiftCtx(iSess,:,:,iShuff) = triu(squeeze(pcaAlignShiftCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftCtx(iSess,:,:,iShuff)));
        pcaAlignShiftEmg(iSess,:,:,iShuff) = triu(squeeze(pcaAlignShiftEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftEmg(iSess,:,:,iShuff)));

        pcaAngleShiftStr(iSess,:,:,iShuff) = triu(squeeze(pcaAngleShiftStr(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftStr(iSess,:,:,iShuff)));
        pcaAngleShiftCtx(iSess,:,:,iShuff) = triu(squeeze(pcaAngleShiftCtx(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftCtx(iSess,:,:,iShuff)));
        pcaAngleShiftEmg(iSess,:,:,iShuff) = triu(squeeze(pcaAngleShiftEmg(iSess,:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftEmg(iSess,:,:,iShuff)));
        
    end

end

commonSaveVarsPCA = {'use10msBins','shiftAmount','pcaProjAllRegionsStr','pcaTrajAllRegionsStr','varExpAllRegionsStr',...
    'pcaProjAllRegionsCtx','pcaTrajAllRegionsCtx','varExpAllRegionsCtx',...
    'pcaProjAllRegionsEmg','pcaTrajAllRegionsEmg','varExpAllRegionsEmg','timeIndsToUse','shiftTimeIndsToUse',...
    'pcaProjStr','pcaTrajStr','varExpStr','pcaProjCtx','pcaTrajCtx','varExpCtx',...
    'pcaProjEmg','pcaTrajEmg','varExpEmg','paDim','cutoffDim','pcaProjStrShift',...
    'pcaProjCtxShift','pcaProjEmgShift','paDimShift','cutoffDimShift',...
    'pcaDiverStr','pcaAlignStr','pcaAngleStr','pcaDiverCtx','pcaAlignCtx','pcaAngleCtx',...
    'pcaDiverEmg','pcaAlignEmg','pcaAngleEmg','rotPaDim1Str','rotPaDim2Str',...
    'pcaDiverRotStr','pcaAlignRotStr','pcaAngleRotStr','pcaDiverRotCtx','pcaAlignRotCtx','pcaAngleRotCtx',...
    'pcaDiverRotEmg','pcaAlignRotEmg','pcaAngleRotEmg','pcaDiverShiftStr','pcaAlignShiftStr','pcaAngleShiftStr',...
    'pcaDiverShiftCtx','pcaAlignShiftCtx','pcaAngleShiftCtx','pcaDiverShiftEmg','pcaAlignShiftEmg','pcaAngleShiftEmg'};

commonSaveVarsCCA = {'ccaNeurProjCombRegStr','ccaEMGProjCombRegStr','cannonCorrsCombRegStr','ccaNeurTrajCombRegStr',...
    'ccaEMGTrajCombRegStr','ccaNeurProjCombRegCtx','ccaEMGProjCombRegCtx','cannonCorrsCombRegCtx','ccaNeurTrajCombRegCtx',...
    'ccaEMGTrajCombRegCtx','ccaNeurProjCombRegLowDimStr','ccaEMGProjCombRegLowDimStr','cannonCorrsCombRegLowDimStr',...
    'ccaNeurTrajCombRegLowDimStr','ccaEMGTrajCombRegLowDimStr','ccaNeurProjCombRegLowDimCtx','ccaEMGProjCombRegLowDimCtx',...
    'cannonCorrsCombRegLowDimCtx','ccaNeurTrajCombRegLowDimCtx','ccaEMGTrajCombRegLowDimCtx','cannonCorrsCombRegStrShift',...
    'cannonCorrsCombRegCtxShift','cannonCorrsCombRegLowDimStrShift','cannonCorrsCombRegLowDimCtxShift','ccaNeurProjStr',...
    'ccaEMGProjStr','cannonCorrsStr','ccaNeurTrajStr','ccaEMGTrajStr','ccaNeurProjCtx','ccaEMGProjCtx','cannonCorrsCtx',...
    'ccaNeurTrajCtx','ccaEMGTrajCtx','ccaNeurProjLowDimStr','ccaEMGProjLowDimStr','cannonCorrsLowDimStr','ccaNeurTrajLowDimStr',...
    'ccaEMGTrajLowDimStr','ccaNeurProjLowDimCtx','ccaEMGProjLowDimCtx','cannonCorrsLowDimCtx','ccaNeurTrajLowDimCtx',...
    'ccaEMGTrajLowDimCtx','timeIndsToUse','paDim','cutoffDim'};

switch lower(inputData)
    case 'humanannotated'
        save('PCASubspacesAnnotated','behvRegionLabels',commonSaveVarsPCA{:},'-v7.3');
        save('CCAAnalysisAnnotated','behvRegionLabels',commonSaveVarsCCA{:},'-v7.3')
    case 'umapregions'
        save('PCASubspaces','behvRegionLabels',commonSaveVarsPCA{:},'-v7.3');
        save('CCAAnalysis','behvRegionLabels',commonSaveVarsCCA{:},'-v7.3')
end

% ccaAlignmentStrShift(:,:,iShift) = ccaAlignmentStr;
% ccaAlignmentCtxShift(:,:,iShift) = ccaAlignmentCtx;
% ccaPrinAngleStrShift(:,:,iShift) = ccaPrinAngleStr;
% ccaPrinAngleCtxShift(:,:,iShift) = ccaPrinAngleCtx;
% 
% pcaAlignmentStrShift(:,:,iShift) = pcaNeurAligment10DimStr;
% pcaAlignmentCtxShift(:,:,iShift) = pcaNeurAligment10DimCtx;
% pcaAlignmentEmgShift(:,:,iShift) = pcaNeurAligment10DimEmg;
% 
% end

plotColors = lines(3);

% plot CCA
plotColors = lines(2);

% treat each region as a separate plot
tmp = permute(cannonCorrsCtx,[3 1 2]);
allRegionCtx = tmp(:,:);

tmp = permute(cannonCorrsStr,[3 1 2]);
allRegionStr = tmp(:,:);

for iRegion = 1:7
    figure
    hold on
    plot(squeeze(cannonCorrsCtx(iRegion,:,:))','Color',[plotColors(1,:), 0.5]);
    plot(squeeze(cannonCorrsStr(iRegion,:,:))','Color',[plotColors(2,:) 0.5])
end

% make example alignment
exBehv1 = 7;
exBehv2 = 6;
exBehv3 = 1;

dat1 = regionFRs{exBehv1}(goodNeuronsStr,timeIndsToUse{iSess,exBehv1})';
dat2 = regionFRs{exBehv2}(goodNeuronsStr,timeIndsToUse{iSess,exBehv2})';
dat3 = regionFRs{exBehv3}(goodNeuronsStr,timeIndsToUse{iSess,exBehv3})';
proj1 = pcaProjStr{exBehv1};
proj2 = pcaProjStr{exBehv2};
proj3 = pcaProjStr{exBehv3};

var1Self = var(dat1 * proj1);
var1Cross2 = var(dat1 * proj2);
var1Cross3 = var(dat1 * proj3);

div2 = calcPCAAlignment(dat1,dat2,proj1,proj2,'auto');
div3 = calcPCAAlignment(dat1,dat3,proj1,proj3,'auto');

figure
tiledlayout(2,1,'Padding','compact','TileSpacing','compact')
nexttile    
plot(var1Self,'-o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(1,:),'Color',plotColors(1,:),'LineWidth',1.5)
hold on
plot(var1Cross2,'-o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(2,:),'Color',plotColors(2,:),'LineWidth',1.5)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickDir','out')
box off
xlim([1,length(var1Self)])
xlabel('Principal Component')
ylabel('Variance')
text(110,6.5,'Climb Up projected on climb up','FontSize',13,'Color',plotColors(1,:))
text(110,5.8,'Climb Up projected on climb down','FontSize',13,'Color',plotColors(2,:))
text(5,6.5,['Divergence = ' num2str(div2)],'fontsize',13)

nexttile
plot(var1Self,'-o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(1,:),'Color',plotColors(1,:),'LineWidth',1.5)
hold on
plot(var1Cross3,'-o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor',plotColors(2,:),'Color',plotColors(2,:),'LineWidth',1.5)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickDir','out')
box off
xlim([1,length(var1Self)])
xlabel('Principal Component')
ylabel('Variance')
text(110,6.5,'Climb Up projected on climb up','FontSize',13,'Color',plotColors(1,:))
text(110,5.8,'Climb Up projected on eating','FontSize',13,'Color',plotColors(2,:))
text(5,6.5,['Divergence = ' num2str(div3)],'fontsize',13)

% make matrices and dendrograms across all behaviors
figure;
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
plotTitle = 'D020 Striatum';
imagesc(pcaNeurAligment10DimStr(behvAlignPerm,behvAlignPerm,1))
set(gca,'XTick',1:length(regionBehvAssignments))
set(gca,'YTick',1:length(regionBehvAssignments))
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gca,'XTickLabelRotation',30)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
cH = colorbar;
% caxis([0 0.55])
cH.Label.String = 'Subspace Divergence';
title(plotTitle)

nexttile
plotTitle = 'D020 Cortex';
imagesc(pcaNeurAligment10DimCtx(behvAlignPerm,behvAlignPerm,1))
set(gca,'XTick',1:length(regionBehvAssignments))
set(gca,'YTick',1:length(regionBehvAssignments))
set(gca,'XTickLabelRotation',30)
set(gca,'XTickLabel',behvRegionLabels)
set(gca,'YTickLabel',behvRegionLabels)
set(gcf,'Color','w')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',13)
set(gca,'TickLength',[0 0])
cH = colorbar;
% caxis([0 0.55])
cH.Label.String = 'Subspace Divergence';
title(plotTitle)

% plot dendrograms
nexttile
plotMatrix = pcaNeurAligment10DimStr(behvAlignPerm,behvAlignPerm,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});

plotMatrix = pcaNeurAligment10DimShuffStr(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(1,:) 0.3]})

plotMatrix = pcaNeurAligment10DimOrthStr(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)-0.15,gca,{'Linewidth',2,'color',[plotColors(2,:) 0.3]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
ylabel('Subspace Divergence')
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')
ylim([0 0.35])
axH = gca;

nexttile
plotMatrix = pcaNeurAligment10DimCtx(behvAlignPerm,behvAlignPerm,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
[~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});

plotMatrix = pcaNeurAligment10DimShuffCtx(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(1,:) 0.3]})

plotMatrix = pcaNeurAligment10DimOrthCtx(:,:,1);
plotMatrix(logical(diag(ones(length(plotMatrix),1),0))) = 0;
customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)-0.15,gca,{'Linewidth',2,'color',[plotColors(2,:) 0.3]})

set(gca,'TickDir','out')
set(gca,'LineWidth',2)
set(gca,'XColor',[0 0 0])
set(gca,'YColor',[0 0 0])
set(gca,'FontSize',14)
ylabel('Subspace Divergence')
set(gca,'XTickLabel',behvRegionLabels(dendPerm))
set(gca,'XTickLabelRotation',30)
set(gca,'TickDir','out')
ylim([0 0.35])
axH = gca;

% plotTitle = 'D020 Cortex';
% imagesc(pcaNeurAligment10DimCtx)
% set(gca,'XTick',1:length(behvNames))
% set(gca,'YTick',1:length(behvNames))
% set(gca,'XTickLabel',behvNames)
% set(gca,'YTickLabel',behvNames)
% set(gcf,'Color','w')
% caxis([0 1])
% colorbar
% title(plotTitle)
% 
% figure;
% dendrogram(linkage(1-squareform(pcaNeurAligment10DimCtx)))
% dendOrder = str2num(get(gca,'XTickLabel'));
% ylim([0 1])
% set(gca,'XTickLabel',behvNames(dendOrder))
% set(gcf,'Color','w')
% title(plotTitle)


% cluster
cgObj = clustergram(alignment);
clustOrdering = cellfun(@(x) str2num(x),cgObj.RowLabels);


% 


function [regionEMGInds, regionNeurInds, regionEMGs, regionFRs, regionFRsOrig] = getBehvRegionData(regionAssignments, regionLabel,...
    frameEMGSamples, frameNeuropixelSamples, use10msBins, origDownsampEMGInd, baseDir, normalizedEMG, normalizedFRs, allFRs)
%get the indices in the UMAP reduction for a region

regionTimeInds = find(regionAssignments == regionLabel);
regionTimeInds(regionTimeInds > floor(frameEMGSamples{1}{end}(end)/20)) = [];

%get emg data for this region
if use10msBins
    regionEMGInds = unique(round(origDownsampEMGInd(regionTimeInds)/10));
    regionEMGInds(regionEMGInds==0) = [];
else
    regionEMGInds = origDownsampEMGInd(regionTimeInds);
end

regionEMGs = normalizedEMG(:,regionEMGInds);


%get neural data for this region
if use10msBins

    %sync umap(EMG) indices to neural indices
    currentDir = pwd;
    cd(fullfile(baseDir,'ProcessedData'))
    regionNeurInds = round(NeurEMGSync(regionEMGInds*200,...
        frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
    cd(currentDir)

    %don't use any time points past the video data
    outOfBoundInds = find(regionNeurInds > floor(frameNeuropixelSamples{1}{end}(end)/300)-1);

    %also don't use points with nans or if the ind is 0
    nanInds = find(isnan(regionNeurInds));
    zeroInds = find(regionNeurInds==0);

    %remove those bad time points
    regionNeurInds([outOfBoundInds nanInds zeroInds]) = [];
    regionEMGs(:,[outOfBoundInds nanInds zeroInds]) = [];

else
    %same as above but without downsampling
    currentDir = pwd;
    cd(fullfile(baseDir,'ProcessedData'))
    regionNeurInds = round(NeurEMGSync(regionEMGInds(regionTimeInds)*20,...
        frameEMGSamples, frameNeuropixelSamples, 'EMG')/30);
    cd(currentDir)

    outOfBoundInds = regionNeurInds > floor(frameNeuropixelSamples{1}{end}(end)/30);
    nanInds = find(isnan(regionNeurInds));
    zeroInds = find(regionNeurInds==0);

    regionNeurInds([outOfBoundInds nanInds zeroInds]) = [];
    regionEMGs(:,[outOfBoundInds nanInds zeroInds]) = [];
end

regionFRs = normalizedFRs(:,regionNeurInds);
%save the un-normalized data as well for calculating firing rates
regionFRsOrig = allFRs(:,regionNeurInds);

%remove Nans
nonNanNeurons = find(~all(isnan(normalizedFRs),2));

nanEMGInds = find(any(isnan(regionEMGs)));
nanFRInds = find(any(isnan(regionFRs(nonNanNeurons,:))));
regionEMGs(:,unique([nanEMGInds nanFRInds])) = [];
regionFRs(:,unique([nanEMGInds nanFRInds])) = [];
regionFRsOrig(:,unique([nanEMGInds nanFRInds])) = [];

end


% 


function alignment = calcPCAAlignmentOld(dat1,dat2,proj1,proj1Cross,proj2,proj2Cross,dims)

if strcmpi(dims,'auto')

    traj1 = dat1*proj1;
    vaf1 = cumsum(var(traj1))/sum(var(traj1));
    [~, dim1] = find(vaf1>0.80,1);

    traj2 = dat2*proj2;
    vaf2 = cumsum(var(traj2))/sum(var(traj2));
    [~, dim2] = find(vaf2>0.80,1);

else

    dim1 = dims;
    dim2 = dims;

end

alignment1 = sum(trace(proj1Cross(:,1:dim1)'*cov(dat1)*proj1Cross(:,1:dim1))) / sum(trace(proj1(:,1:dim1)'*cov(dat1)*proj1(:,1:dim1)));
alignment2 = sum(trace(proj2Cross(:,1:dim2)'*cov(dat2)*proj2Cross(:,1:dim2))) / sum(trace(proj2(:,1:dim2)'*cov(dat2)*proj2(:,1:dim2)));

% alignment1d = sum(var((dat1-mean(dat1))*proj2(:,1:dims)))/sum(var((dat1-mean(dat1))*proj1(:,1:dims)));
% alignment2d = sum(var((dat2-mean(dat2))*proj1(:,1:dims)))/sum(var((dat2-mean(dat2))*proj2(:,1:dims)));

alignment = mean([alignment1 alignment2]);

end


% 


function [divMetric, elsayedAlign, prinAngle] = calcPCAAlignment(dat1,dat2,proj1,proj2,dim1,dim2)

traj1 = (dat1-mean(dat1))*proj1;
vaf1 = cumsum(var(traj1))/sum(var(traj1));

traj2 = (dat2-mean(dat2))*proj2;
vaf2 = cumsum(var(traj2))/sum(var(traj2));

divDim1 = length(proj1);
divDim2 = length(proj2);

% if size(dat1,2) <= 8
%     dim1 = 4;
%     dim2 = 4;
% else
%     dim1 = 10;
%     dim2 = 10;
% end


divergence1 = var(traj1(:,1:divDim1)).^1/sum(var(traj1(:,1:divDim1)).^1) * (abs(diag(proj2(:,1:divDim1)'*cov(dat1)*proj2(:,1:divDim1)) - diag(proj1(:,1:divDim1)'*cov(dat1)*proj1(:,1:divDim1))) ./ ...
    (diag(proj2(:,1:divDim1)'*cov(dat1)*proj2(:,1:divDim1)) + diag(proj1(:,1:divDim1)'*cov(dat1)*proj1(:,1:divDim1))));

divergence2 = var(traj2(:,1:divDim2)).^1/sum(var(traj2(:,1:divDim2)).^1) * (abs(diag(proj1(:,1:divDim2)'*cov(dat2)*proj1(:,1:divDim2)) - diag(proj2(:,1:divDim2)'*cov(dat2)*proj2(:,1:divDim2))) ./ ...
    (diag(proj1(:,1:divDim2)'*cov(dat2)*proj1(:,1:divDim2)) + diag(proj2(:,1:divDim2)'*cov(dat2)*proj2(:,1:divDim2))));

% dim1 = 10; 
% dim2 = 10;
alignment1 = sum(trace(proj2(:,1:dim1)'*cov(dat1)*proj2(:,1:dim1))) / sum(trace(proj1(:,1:dim1)'*cov(dat1)*proj1(:,1:dim1)));
alignment2 = sum(trace(proj1(:,1:dim2)'*cov(dat2)*proj1(:,1:dim2))) / sum(trace(proj2(:,1:dim2)'*cov(dat2)*proj2(:,1:dim2)));

divMetric = mean([divergence1 divergence2]);
elsayedAlign = mean([alignment1 alignment2]);

[~,S,~] = svd(proj1(:,1:dim1)'*proj2(:,1:dim2));
prinAngle = acosd(S(1));

% prinAngle = subspace(proj1(:,1:dim1),proj2(:,1:dim2));

end


% 


function [divMetricRot, elsayedAlignRot, prinAngleRot, paDim1, paDim2] = randomRotationAlignment(pcaTraj1, pcaTraj2, pcaProj1, pcaProj2,dim1,dim2)

% region1Dim = find(cumsum(var(pcaTraj1))/sum(var(pcaTraj1))>0.8,1);
% region2Dim = find(cumsum(var(pcaTraj2))/sum(var(pcaTraj2))>0.8,1);
% 
% region1Dim = length(pcaProj1);
% region2Dim = length(pcaProj2);
% 
% randRotMat1 = randn(region1Dim,region1Dim);
% randRotMat1 = orth(randRotMat1);
% psuedoData1 = pcaTraj1(:,1:region1Dim)*randRotMat1 * pcaProj1';
% 
% randRotMat2 = randn(region2Dim,region2Dim);
% randRotMat2 = orth(randRotMat2);
% psuedoData2 = pcaTraj2(:,1:region2Dim)*randRotMat2  * pcaProj2';
% 
% [psuedoProj1, ~, psuedoVarExp1] = pca(psuedoData1);
% [psuedoProj2, ~, psuedoVarExp2] = pca(psuedoData2);
% 
% randRotAlignment = calcPCAAlignment(psuedoData1,psuedoData2,psuedoProj1,psuedoProj2,10);

region1Dim = find(cumsum(var(pcaTraj1))/sum(var(pcaTraj1))>0.9,1);
region2Dim = find(cumsum(var(pcaTraj2))/sum(var(pcaTraj2))>0.9,1);

region1Dim = length(pcaProj1);
region2Dim = length(pcaProj2);

randRotMat1 = randn(region1Dim,region1Dim);
randRotMat1 = orth(randRotMat1);
psuedoData1 = pcaTraj1(:,1:region1Dim)*randRotMat1 * pcaProj1';

randRotMat2 = randn(region2Dim,region2Dim);
randRotMat2 = orth(randRotMat2);
psuedoData2 = pcaTraj2(:,1:region2Dim)*randRotMat2  * pcaProj2';

[psuedoProj1, ~, psuedoVarExp1] = pca(psuedoData1);
[psuedoProj2, ~, psuedoVarExp2] = pca(psuedoData2);

if isempty(dim1)
    [~, paDim1, ~] = dimEst(psuedoData1, 5);
    [~, paDim2, ~] = dimEst(psuedoData2, 5);
else
    paDim1 = dim1;
    paDim2 = dim2;
end

[divMetricRot, elsayedAlignRot, prinAngleRot] = calcPCAAlignment(psuedoData1,psuedoData2,psuedoProj1,psuedoProj2,paDim1,paDim2);

end


% 


function [prinAngle ccAlignment] = calcCCAAlignmentMetrics(neur1,neur2,emg1,emg2,neurProj1,neurProj2,emgProj1,emgProj2,cc1,cc2)

%calc principle angle
%orthonormalize the projection matrices
neurProjOrtho1 = gson(neurProj1);
neurProjOrtho2 = gson(neurProj2);
[~, S, ~] = svd(neurProjOrtho1'*neurProjOrtho2);
Svalues = diag(S);
allPrinAngles = acosd(Svalues);

prinAngle = allPrinAngles(1);


%calc CC Alignment
centeredNeur1 = neur1 - mean(neur1,2);
centeredEMG1 = emg1 - mean(emg1,2);

trajCrossNeur1 = centeredNeur1'*neurProj2;
trajCrossEMG1 = centeredEMG1'*emgProj2;
crossCC1 = corr(trajCrossNeur1,trajCrossEMG1);
crossCC1 = abs(diag(crossCC1));
% cc1Ratio = sum(crossCC1)/sum(cc1);
cc1Ratio = (cc1/sum(cc1))*(abs(crossCC1'-cc1)./(crossCC1'+cc1))';

centeredNeur2 = neur2 - mean(neur2,2);
centeredEMG2 = emg2 - mean(emg2,2);

trajCrossNeur2 = centeredNeur2'*neurProj1;
trajCrossEMG2 = centeredEMG2'*emgProj1;
crossCC2 = corr(trajCrossNeur2,trajCrossEMG2);
crossCC2 = abs(diag(crossCC2));
% cc2Ratio = sum(crossCC2)/sum(cc2);
cc2Ratio = (cc2/sum(cc2))*(abs(crossCC2'-cc2)./(crossCC2'+cc2))';

ccAlignment = mean([cc1Ratio cc2Ratio]);

end


% 
