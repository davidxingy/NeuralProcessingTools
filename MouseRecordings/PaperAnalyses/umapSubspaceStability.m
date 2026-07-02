clear

% do analysis for each of the datasets
sessionNames = {'D020-062922-ArenaRecording', 'D024-111022-ArenaRecording', ...
    'D026-032923-ArenaRecording', 'D050-121825-ArenaRecording', 'D054-012126-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7 ...
    ];

inputData = 'umapregions';
binSize = 100;
useSmoothed = false;
runCCA = false;
runPCA = false;
runDyn = true;
runPLDS = false;
runSubregions = true;
subBehvMinPoints = 30000/binSize;
nShifts = 100;
runDimCalc = false;
minLabelPoints = 30000/binSize; %Don't limit number of points for calculating things within behaviors to less than this value

dynNoId = false;
dynNoiseStd = 0.01;

exampleControl = 2;

rng(2025)

for iSess = 1:length(sessionNames)

    filePaths = getMouseDataNames(sessionNames{iSess}(1:4),sessionNames{iSess},'CFA');
    sessionArtifacts = consolidateArtifactInds(sessionNames{iSess},'CFA');
    behvAlignPerm = allBehvAlignPerms(iSess,:);

    % load in labels to split the groups (behavior groups, muscle
    % groups, ect)
    load(filePaths.UMAPFile, 'regionAssignmentsFiltered','behvLabelsNoArt','classifierLabels','musclePercentileGroupLabels',...
        'origDownsampEMGInd','subRegionAssignments','analyzedBehaviors','regionBehvAssignments','regionWatershedLabels')

    % load in neural activity
    switch binSize
        case 100
            if useSmoothed
                load(filePaths.NeuralFiringRates100msBins50msGauss, 'cortexInds', 'striatumInds','allFRs')
            else
                load(filePaths.NeuralFiringRates100msBins0msGauss, 'cortexInds', 'striatumInds','allFRs')
            end
        case 10
            if useSmoothed
                load(filePaths.NeuralFiringRates10msBins30msGauss, 'cortexInds', 'striatumInds','allFRs')
            else
                load(filePaths.NeuralFiringRates10msBins0msGauss, 'cortexInds', 'striatumInds','allFRs')
            end
        case 1
            if useSmoothed
                load(filePaths.NeuralFiringRates1msBins10msGauss, 'cortexInds', 'striatumInds','allFRs')
            else
                load(filePaths.NeuralFiringRates1msBins0msGauss, 'cortexInds', 'striatumInds','allFRs')
            end
    end
    % normalize data
    normalizedFRs = (allFRs-nanmean(allFRs,2))./nanstd(allFRs,[],2);
    % normalizedFRs(isinf(normalizedFRs)) = nan;

    % load in EMG activity
    load(filePaths.EMG1ms)
    load(filePaths.VideoSyncFrames)

    % if we want to downsample emg beyond 1ms bins, take the average across bins
    if binSize > 1
        emg1ms = downsampEMG;
        downsampleFactor = binSize;
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

    % do it for different label types
    labelTypeNames = {'UmapRegions','ClassifierBehvs','MusclePercentile','HumanLabels'};

    allUsedLabels = {regionAssignmentsFiltered, classifierLabels, musclePercentileGroupLabels, behvLabelsNoArt};
    if iscell(regionWatershedLabels)
        regionWatershedLabelsCat = cat(2,regionWatershedLabels{:});
    else
        regionWatershedLabelsCat = regionWatershedLabels;
    end
    allUsedLabelIds = {regionWatershedLabelsCat, 1:10, 1:8, 1:10};

    % loop through each behavior label type
    behvAveSigs = [];
    behvPrctSigs = [];
    behvAveSigsResamp = [];
    behvPrctSigsResamp = [];
    behvAveSigsShuff = [];
    behvPrctSigsShuff = [];

    for iLabelType = 1:length(allUsedLabels)

        usedLabels = allUsedLabels{iLabelType};
        usedLabelIds = allUsedLabelIds{iLabelType};
        nBehvs = length(usedLabelIds);

        % specific processing for each label type:
        switch lower(labelTypeNames{iLabelType})

            case 'umapregions'
                % -----Use data divded by the 7 umap regions-----

                behvLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rearing/Still','Groom','Eat'};

                exampleSess = 3;
                exampleBaseBehv = 1;
                exampleLowDivBehv = 2;
                exampleHighDivBehv = 7;
                crossRegionAllTimes = false;
                exampleCCATimes = 6000:7000;

                if runSubregions
                    load(filePaths.UMAPFile, 'regionAssignmentsFineNoBound','fineRegionCorrespondence')
                end
                labelTypeRunSubregions = runSubregions;

                %align behavior regions across animals
                usedLabelIds = usedLabelIds(behvAlignPerm);

            case 'classifierbehvs'

                behvLabels = {'Climb Down','Climb Up','Eat','Groom','Jump Down','Jump','Rear','Still','Walk Flat','Walk Grid'};

                exampleSess = 3;
                exampleBaseBehv = 2;
                exampleLowDivBehv = 1;
                exampleHighDivBehv = 3;
                crossRegionAllTimes = false;
                exampleCCATimes = 6000:7000;

                labelTypeRunSubregions = false;

            case 'musclepercentile'

                behvLabels = {'EMG 12.5th Percentile','EMG 25th Percentile','EMG 37.5th Percentile','EMG 50th Percentile',...
                    'EMG 62.5th Percentile','EMG 75th Percentile','EMG 87.5th Percentile','EMG 100th Percentile'};

                exampleSess = 3;
                exampleBaseBehv = 1;
                exampleLowDivBehv = 2;
                exampleHighDivBehv = 8;
                crossRegionAllTimes = false;
                exampleCCATimes = 6000:7000;

                labelTypeRunSubregions = false;

            case 'humanlabels'

                behvLabels = {'Climb Down','Climb Up','Eat','Groom','Jump Down','Jump','Rear','Still','Walk Flat','Walk Grid'};

                exampleSess = 3;
                exampleBaseBehv = 2;
                exampleLowDivBehv = 1;
                exampleHighDivBehv = 3;
                crossRegionAllTimes = false;
                exampleCCATimes = 6000:7000;

                labelTypeRunSubregions = false;



                % % case 'humanannotated'
                % %     % -----Use data divded by human annotations-----
                % %
                % %     behvRegionFieldNames = {'climbup','climbdown','jumpdown','walkflat','walkgrid','rearing','grooming','eating'};
                % %     behvRegionLabels = {'Climb Up','Climb Down','Jump Down','Walk Flat','Walk Grid','Rear','Groom','Eat'};
                % %
                % %     exampleSess = 3;
                % %     exampleLowDivRegion = 2;
                % %     exampleHighDivRegion = 8;
                % %
                % %     load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'), 'cortexInds', 'striatumInds')
                % %     if use10msBins
                % %         load(fullfile(baseDir,'ProcessedData','EpochedData10ms.mat'))
                % %     else
                % %         load(fullfile(baseDir,'ProcessedData','EpochedData1ms.mat'))
                % %     end
                % %
                % %     sessStructNames = fieldnames(behavioralData);
                % %     nRegions = length(behvRegionFieldNames);
                % %
                % %     %get data for each behavior
                % %     for iBehv = 1:nRegions
                % %
                % %         behvStructInd(iBehv) = find(strcmpi(sessStructNames,behvRegionFieldNames{iBehv}));
                % %         regionFRsOrig{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).allBoutFRs;
                % %         regionEMGsOrig{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).allBoutEMGs;
                % %
                % %         %remove nans
                % %         nanInds = unique([find(any(isnan(regionFRsOrig{iBehv}),1)) find(any(isnan(regionEMGsOrig{iBehv}),1))]);
                % %         regionFRsOrig{iBehv}(:,nanInds) = [];
                % %         regionEMGsOrig{iBehv}(:,nanInds) = [];
                % %
                % %         %also get information for doing the shift controls
                % %         behvBoutFRs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutFRs;
                % %         behvBoutEMGs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutEMGs;
                % %
                % %         %get the time of each bout, to put in order for getting the shift controls
                % %         behvBoutStartInds{iBehv} = cellfun(@(x) x(1),allNeurInds(behvStructInd(iBehv),~cellfun(@isempty,allNeurInds(behvStructInd(iBehv),:))));
                % %         behvBoutLabels{iBehv} = cellfun(@(x) repmat(iBehv,1,size(x,2)), behvBoutFRs{iBehv}, 'un', 0);
                % %
                % %     end
                % %
                % %     %normalize based on activity across all behaviors
                % %     allFRs = cat(2,regionFRsOrig{:});
                % %     allEMGs = cat(2,regionEMGsOrig{:});
                % %
                % %     allFRMeans = mean(allFRs,2);
                % %     allFRStds = std(allFRs,[],2);
                % %     allEMGMeans = mean(allEMGs,2);
                % %     allEMGStds = std(allEMGs,[],2);
                % %
                % %     regionFRs = cellfun(@(x) (x-allFRMeans)./allFRStds,regionFRsOrig,'un',0);
                % %     regionEMGs = cellfun(@(x) (x-allEMGMeans)./allEMGStds,regionEMGsOrig,'un',0);
                % %
                % %     % do shift controls to break behavior relationship with neur/emg
                % %     % first combine all behaviors (sorting bouts by time)
                % %     [~, sortPerm] = sort([behvBoutStartInds{:}]);
                % %     allBoutFRs = [behvBoutFRs{:}];
                % %     allBoutFRs = allBoutFRs(sortPerm);
                % %     allBoutEMGs = [behvBoutEMGs{:}];
                % %     allBoutEMGs = allBoutEMGs(sortPerm);
                % %     allBoutLabels = [behvBoutLabels{:}];
                % %     allBoutLabels = allBoutLabels(sortPerm);
                % %
                % %     allBoutFRsCat = cat(2,allBoutFRs{:});
                % %     allBoutEMGsCat = cat(2,allBoutEMGs{:});
                % %     allBoutLabelsCat = cat(2,allBoutLabels{:});
                % %
                % %     % don't use points with nans
                % %     catFRNans = find(any(isnan(allBoutFRsCat)));
                % %     catEMGNans = find(any(isnan(allBoutEMGsCat)));
                % %
                % %     allBoutFRsCat(:,unique([catFRNans catEMGNans])) = [];
                % %     allBoutEMGsCat(:,unique([catFRNans catEMGNans])) = [];
                % %     allBoutLabelsCat(:,unique([catFRNans catEMGNans])) = [];
                % %
                % %     % do shift then divide into separate behaviors again
                % %     for iShift = 1:nShifts
                % %
                % %         %shift by at least 30 seconds
                % %         if use10msBins
                % %             shiftAmount(iShift) = randi(length(allBoutLabelsCat)-3000*2)+3000;
                % %         else
                % %             shiftAmount(iShift) = randi(length(allBoutLabelsCat)-30000*2)+30000;
                % %         end
                % %         labelsShift = circshift(allBoutLabelsCat,shiftAmount(iShift));
                % %
                % %         % get FRs and EMGs based on new shifted labels
                % %         for iBehv = 1:length(behvRegionFieldNames)
                % %             shiftBehvInds = find(labelsShift == iBehv);
                % %             regionFRsShift{iBehv,iShift} = allBoutFRsCat(:,shiftBehvInds);
                % %             regionEMGsShift{iBehv,iShift} = allBoutEMGsCat(:,shiftBehvInds);
                % %         end
                % %
                % %     end
                % %



        end

        %go through each region and divide data
        behvEMGs = {};
        behvFRs = {};
        behvFRsOrig = {};
        usedDynInds = {};
        behvNeurIndsShift = {};
        behvEMGsShift = {};
        usedDynIndsShift = {};
        subRegionFRs = {};
        subRegionFRsOrig = {};
        subRegionEMGs = {};
        usedDynSubRegionInds = {};
        usedSubregions = {};
        subRegionLeastPoints = [];
        subRegionDynLeastPoints = [];

        for iBehv = 1:nBehvs

            [~, behvNeurInds, behvEMGs{iBehv}, behvFRs{iBehv}, behvFRsOrig{iBehv}] = ...
                getBehvData(usedLabels, usedLabelIds(iBehv), binSize, useSmoothed, ...
                origDownsampEMGInd, sessionNames{iSess}, sessionArtifacts, normalizedEMG, normalizedFRs, allFRs);

            % for dynamical model, need continuous points
            usedDynInds{iBehv} = find(diff(behvNeurInds) == 1)+1;

            %also, as a control we can do subregions of each region, so
            %get the data for each subregion
            if labelTypeRunSubregions

                subRegionlabels = fineRegionCorrespondence{iBehv};
                for iSubRegion = 1:length(subRegionlabels)
                    [~, subRegionNeurInds, subRegionEMGs{iBehv}{iSubRegion}, subRegionFRs{iBehv}{iSubRegion}, subRegionFRsOrig{iBehv}{iSubRegion}] = getBehvData(regionAssignmentsFineNoBound, subRegionlabels(iSubRegion),...
                        binSize, useSmoothed, origDownsampEMGInd, sessionNames{iSess}, sessionArtifacts, normalizedEMG, normalizedFRs, allFRs);

                    usedDynSubRegionInds{iBehv}{iSubRegion} = find(diff(subRegionNeurInds) == 1)+1;
                end

                %Don't use subregions that only have a few amount of data points
                subRegionNPoints = cellfun(@(x) size(x,2), subRegionFRs{iBehv});
                subRegionNDynPoints = cellfun(@length, usedDynSubRegionInds{iBehv});
                usedSubregions{iBehv} = find(subRegionNPoints > subBehvMinPoints);
                if isempty(usedSubregions{iBehv})
                    subRegionLeastPoints(iBehv) = 0;
                    subRegionDynLeastPoints(iBehv) = 0;
                else
                    subRegionLeastPoints(iBehv) = min(subRegionNPoints(usedSubregions{iBehv}));
                    subRegionDynLeastPoints(iBehv) = min(subRegionNDynPoints(usedSubregions{iBehv}));
                end

            end

            % do shift controls to break behavior relationship with neur/emg
            for iShift = 1:nShifts
                %shift by at least 100 seconds
                shiftAmount(iShift) = randi(length(usedLabels)-100000/binSize*2)+100000/binSize;
                usedLabelsShifted = circshift(usedLabels,shiftAmount(iShift));

                % for shifts, only save the indices and extract the neural
                % activity on the fly, otherwise it will take too much
                % memory to save all the shift firing rates
                [~, behvNeurIndsShift{iBehv,iShift}, behvEMGsShift{iBehv,iShift}, ~, ~] = ...
                    getBehvData(usedLabelsShifted, usedLabelIds(iBehv), binSize, useSmoothed, ...
                    origDownsampEMGInd, sessionNames{iSess}, sessionArtifacts, normalizedEMG, normalizedFRs, allFRs);

                usedDynIndsShift{iBehv,iShift} = find(diff(behvNeurIndsShift{iBehv,iShift}) == 1)+1;
            end

        end

        goodSubregions{iSess,iLabelType} = usedSubregions;

        % for the some label types (e.g. classifier), there are some behaviors
        % where there are very few (or even no) points (mostly
        % jumping/jumpdown). Since I don't want to negatively influence
        % the firing rate calculations for all the other behaviors by
        % limiting them to such a small number of points, use a minimum
        % of minLabelPoints (default will use 20s of data),
        % and if any label types have less then those points, use all
        % of the points for that label type, and then save a variable
        % indicating that they have different points in the
        % calculations to inform future analyses and results.
        behvNPoints = cellfun(@(x) size(x,2), behvFRs);
        behvsLessThanMinPoints{iSess,iLabelType} = find(behvNPoints < minLabelPoints);
        leastPoints = min(behvNPoints(behvNPoints >= minLabelPoints));

        % also do it for points used in the dynamical system calculations
        % (but use behaviors based on full region inds)
        dynNPoints = cellfun(@length, usedDynInds);
        dynLeastPoints = min(dynNPoints(behvNPoints >= minLabelPoints));

        %calculate average firing rates for cells within each behv
        behvMeanFRsCell = cellfun(@(x) nanmean(x,2),behvFRsOrig,'UniformOutput',false);
        behvMeanFRs = cat(2,behvMeanFRsCell{:});

        %change to spks/sec
        behvMeanFRs = behvMeanFRs * 1000/binSize;

        %only use neurons above a certain firing rate within at least one behavior
        firingRateCutoff = 0.2;
        goodNeuronsAll = find(any(behvMeanFRs > firingRateCutoff,2));

        %more strict criterion we can use - use only neurons above a certain
        %firing rate across all behaviors
        reallyGoodNeuronsAll = find(all(behvMeanFRs > firingRateCutoff,2));

        %also divide into striatum and cortex
        goodNeuronsStr = intersect(goodNeuronsAll,1:length(striatumInds));
        goodNeuronsCtx = intersect(goodNeuronsAll,length(striatumInds)+1:size(allFRs,1));

        %     goodNeuronsStr = goodNeuronsStr(randperm(length(goodNeuronsStr),length(goodNeuronsCtx)));

        reallyGoodNeuronsStr = intersect(reallyGoodNeuronsAll,1:length(striatumInds));
        reallyGoodNeuronsCtx = intersect(reallyGoodNeuronsAll,length(striatumInds)+1:size(allFRs,1));

        % do PCA on all time points ----------------------------------------------------------------
        allBehvsFR = cat(2,behvFRs{:});
        allBehvsEMG = cat(2,behvEMGs{:});
        allBehvsInds = [0 cumsum(cellfun(@(x) size(x,2), behvFRs))];

        [pcaProjAllBehvsStr{iSess}, pcaTrajAllBehvsStr{iSess}, varExpAllBehvsStr{iSess}] = pca(allBehvsFR(goodNeuronsStr,:)');
        [pcaProjAllBehvsCtx{iSess}, pcaTrajAllBehvsCtx{iSess}, varExpAllBehvsCtx{iSess}] = pca(allBehvsFR(goodNeuronsCtx,:)');
        [pcaProjAllBehvsEmg{iSess}, pcaTrajAllBehvsEmg{iSess}, varExpAllBehvsEmg{iSess}] = pca(allBehvsEMG');

        %calculate dimensionality
        if runDimCalc
            [lbmleDimAllBehvs(1), paDimAllBehvs(1), cutoffDimAllBehvs(1)] = dimEst(allBehvsFR(goodNeuronsStr,:)', 5);
            [lbmleDimAllBehvs(2), paDimAllBehvs(2), cutoffDimAllBehvs(2)] = dimEst(allBehvsFR(goodNeuronsCtx,:)', 5);
            [lbmleDimAllBehvs(3), paDimAllBehvs(3), cutoffDimAllBehvs(3)] = dimEst(allBehvsEMG', 5);
        else
            paDimAllBehvs(1:3) = 5;
        end

        if runCCA && iLabelType == 1

            % do CCA on all time points (only need to do it once)
            [ccaNeurProjCombStr{iSess},ccaEMGProjCombStr{iSess},cannonCorrsCombStr(iSess,:),ccaNeurTrajCombStr{iSess},ccaEMGTrajCombStr{iSess}] = ...
                canoncorr(allBehvsFR(goodNeuronsStr,:)',allBehvsEMG(1:4,:)');
            [ccaNeurProjCombCtx{iSess},ccaEMGProjCombCtx{iSess},cannonCorrsCombCtx(iSess,:),ccaNeurTrajCombCtx{iSess},ccaEMGTrajCombCtx{iSess}] = ...
                canoncorr(allBehvsFR(goodNeuronsCtx,:)',allBehvsEMG(1:4,:)');

            %do CCA using only top PC's
            [ccaNeurProjCombLowDimStr{iSess},ccaEMGProjCombLowDimStr{iSess},cannonCorrsCombLowDimStr(iSess,:),...
                ccaNeurTrajCombLowDimStr{iSess},ccaEMGTrajCombLowDimStr{iSess}] = ...
                canoncorr(pcaTrajAllBehvsStr{iSess}(:,1:max(paDimAllBehvs(1),4)),allBehvsEMG(1:4,:)');
            [ccaNeurProjCombLowDimCtx{iSess},ccaEMGProjCombLowDimCtx{iSess},cannonCorrsCombLowDimCtx(iSess,:),...
                ccaNeurTrajCombLowDimCtx{iSess},ccaEMGTrajCombLowDimCtx{iSess}] = ...
                canoncorr(pcaTrajAllBehvsCtx{iSess}(:,1:max(paDimAllBehvs(2),4)),allBehvsEMG(1:4,:)');

            %shift EMG and neural data relative to each other for controls for CCA
            for iShift = 1:nShifts

                %shift by at least 5 seconds
                emgShiftAmount = randi(size(allBehvsEMG,2)-5000/binSize*2)+5000/binSize;

                shiftEMG = circshift(allBehvsEMG(1:4,:),emgShiftAmount,2)';
                [~,~,cannonCorrsCombStrShift(iSess,iShift,:)] = canoncorr(allBehvsFR(goodNeuronsStr,:)',shiftEMG);
                [~,~,cannonCorrsCombCtxShift(iSess,iShift,:)] = canoncorr(allBehvsFR(goodNeuronsCtx,:)',shiftEMG);

                [~,~,cannonCorrsCombLowDimStrShift(iSess,iShift,:)] = canoncorr(pcaTrajAllBehvsStr{iSess}(:,1:max(paDimAllBehvs(1),4)),shiftEMG);
                [~,~,cannonCorrsCombLowDimCtxShift(iSess,iShift,:)] = canoncorr(pcaTrajAllBehvsCtx{iSess}(:,1:max(paDimAllBehvs(2),4)),shiftEMG);

            end

        end

        %Free up some memory
        % clear allFRs
        % clear normalizedFRs
        % clear behvFRsOrig
        clear allBehvsFR allBehvsEMG

        %make pseudo-data, take random subset, of the top 50 PC's,
        %rotate it then project back up to
        %     for iShuff = 1:nShifts
        %
        %         for iBehv1 = 1:length(regionWatershedLabels)
        %             for iBehv2 = 1:length(regionWatershedLabels)
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
        %                 pcaNeurAligment10DimAll(iBehv1,iBehv2,iShuff) = calcPCAAlignment(psuedoData1,psuedoData2,psuedoProj1,psuedoProj2,10);
        %
        %             end
        %         end
        %
        %     end


        % now go through each region and run PCA, CCA, and LDS
        for iShift = 1:nShifts+1

            pcaTrajStr = {};
            pcaTrajCtx = {};
            pcaTrajEmg = {};
            pcaTrajSubStr = {};
            pcaTrajSubCtx = {};
            pcaTrajSubEmg = {};

            dynCtxFRs = {};
            dynStrFRs = {};
            dynCtxSubFRs = {};
            dynStrSubFRs = {};
            dynInds = {};
            dynSubInds = {};
            predictionsCtx = {};
            predictionsStr = {};
            predictionsSubCtx = {};
            predictionsSubStr = {};

            for iBehv = 1:nBehvs

                tic

                % since shifting uses the same code, fold shifts and the normal
                % unshifted run into the same code and use a for loop

                % first loop is unshifted
                if iShift == 1
                    isShift = false;
                    frsData = behvFRs{iBehv};
                    frsOrigData = behvFRsOrig{iBehv};
                    emgData = behvEMGs;

                    %downsample so same number of points in all regions (except for
                    %the behaviors where there are very few points
                    if leastPoints < size(behvFRs{iBehv},2)
                        timeIndsToUse{iSess,iLabelType}{iBehv} = randperm(size(behvFRs{iBehv},2),leastPoints);
                    else
                        timeIndsToUse{iSess,iLabelType}{iBehv} = 1:size(behvFRs{iBehv},2);
                    end

                    usedTimeInds = timeIndsToUse{iSess,iLabelType}{iBehv};

                    % for dynmical model systems too
                    if dynLeastPoints < length(usedDynInds{iBehv})
                        dynTimeIndsToUse{iSess,iLabelType}{iBehv} = randperm(length(usedDynInds{iBehv}),dynLeastPoints);
                    else
                        dynTimeIndsToUse{iSess,iLabelType}{iBehv} = 1:length(usedDynInds{iBehv});
                    end

                    usedDynTimeinds = dynTimeIndsToUse{iSess,iLabelType}{iBehv};
                    dynInds{iBehv} = usedDynInds{iBehv}(usedDynTimeinds);

                else

                    isShift = true;
                    frsData = normalizedFRs(:,behvNeurIndsShift{iBehv,iShift-1});
                    frsOrigData = allFRs(:,behvNeurIndsShift{iBehv,iShift-1});
                    emgData = behvEMGsShift(:,iShift -1);

                    if leastPoints < size(frsData,2)
                        shiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1} = randperm(size(frsData,2),leastPoints);
                    else
                        shiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1} = 1:size(frsData,2);
                    end

                    usedTimeInds = shiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1};

                    if dynLeastPoints < length(usedDynIndsShift{iBehv,iShift-1})
                        dynShiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1} = randperm(length(usedDynIndsShift{iBehv,iShift -1}),dynLeastPoints);
                    else
                        dynShiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1} = 1:length(usedDynIndsShift{iBehv,iShift -1});
                    end

                    usedDynTimeinds = dynShiftTimeIndsToUse{iSess,iLabelType}{iBehv,iShift-1};
                    dynInds{iBehv} = usedDynIndsShift{iBehv,iShift-1}(usedDynTimeinds);

                end

                if runPCA
                    % do PCA
                    [projStr, pcaTrajStr{iBehv}, varStr] = pca(frsData(goodNeuronsStr,usedTimeInds)');
                    [projCtx, pcaTrajCtx{iBehv}, varCtx] = pca(frsData(goodNeuronsCtx,usedTimeInds)');
                    [projEmg, pcaTrajEmg{iBehv}, varEmg] = pca(emgData{iBehv}(:,usedTimeInds)');

                    %calculate dimensionality
                    % to save time, just calculate the dimensionality of shift
                    % controls for just the first shift and use that for the rest
                    if runDimCalc
                        if iShift <= 2
                            [lbmleDimStr, paDimStr, cutoffDimStr] = dimEst(frsData(goodNeuronsStr,usedTimeInds)', 5);
                            [lbmleDimCtx, paDimCtx, cutoffDimCtx] = dimEst(frsData(goodNeuronsCtx,usedTimeInds)', 5);
                            [lbmleDimEmg, paDimEmg, cutoffDimEmg] = dimEst(emgData{iBehv}(:,usedTimeInds)', 5);
                        else
                            lbmleDimStr = lbmleDimShift{iSess,iLabelType}(iBehv,1,1);
                            lbmleDimCtx = lbmleDimShift{iSess,iLabelType}(iBehv,1,2);
                            lbmleDimEmg = lbmleDimShift{iSess,iLabelType}(iBehv,1,3);

                            paDimStr = paDimShift{iSess,iLabelType}(iBehv,1,1);
                            paDimCtx = paDimShift{iSess,iLabelType}(iBehv,1,2);
                            paDimEmg = paDimShift{iSess,iLabelType}(iBehv,1,3);

                            cutoffDimStr = cutoffDimShift{iSess,iLabelType}(iBehv,1,1);
                            cutoffDimCtx = cutoffDimShift{iSess,iLabelType}(iBehv,1,2);
                            cutoffDimEmg = cutoffDimShift{iSess,iLabelType}(iBehv,1,3);
                        end
                    else
                        lbmleDimStr = 5; paDimStr = 5; cutoffDimStr = 5;
                        lbmleDimCtx = 5; paDimCtx = 5; cutoffDimCtx = 5;
                        lbmleDimEmg = 5; paDimEmg = 5; cutoffDimEmg = 5;
                    end

                    % save to different variables
                    if isShift
                        pcaProjStrShift{iSess,iLabelType}{iBehv,iShift-1} = projStr;
                        pcaProjCtxShift{iSess,iLabelType}{iBehv,iShift-1} = projCtx;
                        pcaProjEmgShift{iSess,iLabelType}{iBehv,iShift-1} = projEmg;

                        if runDimCalc
                            lbmleDimShift{iSess,iLabelType}(iBehv,iShift-1,1:3) = [lbmleDimStr lbmleDimCtx lbmleDimEmg];
                            paDimShift{iSess,iLabelType}(iBehv,iShift-1,1:3) = [paDimStr paDimCtx paDimEmg];
                            cutoffDimShift{iSess,iLabelType}(iBehv,iShift,1:3) = [cutoffDimStr cutoffDimCtx cutoffDimEmg];
                        else
                            paDimShift{iSess,iLabelType}(iBehv,iShift-1,1:3) = 5;
                        end
                    else
                        pcaProjStr{iSess,iLabelType}{iBehv} = projStr;
                        varExpStr{iSess,iLabelType}{iBehv} = varStr;

                        pcaProjCtx{iSess,iLabelType}{iBehv} = projCtx;
                        varExpCtx{iSess,iLabelType}{iBehv} = varCtx;

                        pcaProjEmg{iSess,iLabelType}{iBehv} = projEmg;
                        varExpEmg{iSess,iLabelType}{iBehv} = varEmg;

                        if runDimCalc
                            lbmleDim{iSess,iLabelType}(iBehv,1:3) = [lbmleDimStr lbmleDimCtx lbmleDimEmg];
                            paDim{iSess,iLabelType}(iBehv,1:3) = [paDimStr paDimCtx paDimEmg];
                            cutoffDim{iSess,iLabelType}(iBehv,1:3) = [cutoffDimStr cutoffDimCtx cutoffDimEmg];
                        else
                            paDim{iSess,iLabelType}(iBehv,1:3) = 5;
                        end
                    end

                end % of run PCA block

                if runCCA

                    if isShift
                        dimToUse = paDimShift{iSess,iLabelType}(iBehv,iShift-1,1:2);
                    else
                        dimToUse = paDim{iSess,iLabelType}(iBehv,1:2);
                    end

                    %do CCA
                    [neurProjStr, emgProjStr, ccStr, ~, ~] = ...
                        canoncorr(frsData(goodNeuronsStr,usedTimeInds)', emgData{iBehv}(1:4,usedTimeInds)');
                    [neurProjCtx, emgProjCtx, ccCtx, ~, ~] = ...
                        canoncorr(frsData(goodNeuronsCtx,usedTimeInds)', emgData{iBehv}(1:4,usedTimeInds)');

                    %do CCA using only top PC's
                    [neurProjLowDStr, emgProjLowDStr, ccLowDStr, ~, ~] = ...
                        canoncorr(pcaTrajStr{iBehv}(:,1:max(dimToUse(1),4)), emgData{iBehv}(1:4,usedTimeInds)');
                    [neurProjLowDCtx, emgProjLowDCtx, ccLowDCtx, ~, ~] = ...
                        canoncorr(pcaTrajCtx{iBehv}(:,1:max(dimToUse(2),4)), emgData{iBehv}(1:4,usedTimeInds)');

                    % save to different variables
                    if isShift
                        ccaNeurProjStrShift{iSess,iLabelType}{iBehv,iShift-1} = neurProjStr;
                        ccaEMGProjStrShift{iSess,iLabelType}{iBehv,iShift-1} = emgProjStr;
                        cannonCorrsStrShift{iSess,iLabelType}(iBehv,iShift-1,:) = ccStr;

                        ccaNeurProjLowDimStrShift{iSess,iLabelType}{iBehv,iShift-1} = neurProjLowDStr;
                        ccaEMGProjLowDimStrShift{iSess,iLabelType}{iBehv,iShift-1} = emgProjLowDStr;
                        cannonCorrsLowDimStrShift{iSess,iLabelType}(iBehv,iShift-1,:) = ccLowDStr;

                        ccaNeurProjCtxShift{iSess,iLabelType}{iBehv,iShift-1} = neurProjCtx;
                        ccaEMGProjCtxShift{iSess,iLabelType}{iBehv,iShift-1} = emgProjCtx;
                        cannonCorrsCtxShift{iSess,iLabelType}(iBehv,iShift-1,:) = ccCtx;

                        ccaNeurProjLowDimCtxShift{iSess,iLabelType}{iBehv,iShift-1} = neurProjLowDCtx;
                        ccaEMGProjLowDimCtxShift{iSess,iLabelType}{iBehv,iShift-1} = emgProjLowDCtx;
                        cannonCorrsLowDimCtxShift{iSess,iLabelType}(iBehv,iShift-1,:) = ccLowDCtx;
                    else
                        ccaNeurProjStr{iSess,iLabelType}{iBehv} = neurProjStr;
                        ccaEMGProjStr{iSess,iLabelType}{iBehv} = emgProjStr;
                        cannonCorrsStr{iSess,iLabelType}(iBehv,:) = ccStr;

                        ccaNeurProjLowDimStr{iSess,iLabelType}{iBehv} = neurProjLowDStr;
                        ccaEMGProjLowDimStr{iSess,iLabelType}{iBehv} = emgProjLowDStr;
                        cannonCorrsLowDimStr{iSess,iLabelType}(iBehv,:) = ccLowDStr;

                        ccaNeurProjCtx{iSess,iLabelType}{iBehv} = neurProjCtx;
                        ccaEMGProjCtx{iSess,iLabelType}{iBehv} = emgProjCtx;
                        cannonCorrsCtx{iSess,iLabelType}(iBehv,:) = ccCtx;

                        ccaNeurProjLowDimCtx{iSess,iLabelType}{iBehv} = neurProjLowDCtx;
                        ccaEMGProjLowDimCtx{iSess,iLabelType}{iBehv} = emgProjLowDCtx;
                        cannonCorrsLowDimCtx{iSess,iLabelType}(iBehv,:) = ccLowDCtx;
                    end

                end % of run CCA block


                if runCCA && ~isShift
                    %shift EMG and neural data relative to each other for controls for CCA
                    for iCCAShift = 1:nShifts

                        %shift by at least 5 seconds
                        emgShiftAmount = randi(abs(size(emgData{iBehv},2)-5000/binSize*2))+5000/binSize;

                        shiftEMG = circshift(emgData{iBehv}(1:4,:), emgShiftAmount, 2)';
                        [~,~,cannonCorrsStrShift{iSess,iLabelType}(iBehv,iCCAShift,:)] = canoncorr(frsData(goodNeuronsStr,usedTimeInds)',shiftEMG(usedTimeInds,:));
                        [~,~,cannonCorrsCtxShift{iSess,iLabelType}(iBehv,iCCAShift,:)] = canoncorr(frsData(goodNeuronsCtx,usedTimeInds)',shiftEMG(usedTimeInds,:));

                        [~,~,cannonCorrsLowDimStrShift{iSess,iLabelType}(iBehv,iCCAShift,:)] = canoncorr(pcaTrajStr{iBehv}(:,1:max(dimToUse(1),4)),shiftEMG(usedTimeInds,:));
                        [~,~,cannonCorrsLowDimCtxShift{iSess,iLabelType}(iBehv,iCCAShift,:)] = canoncorr(pcaTrajCtx{iBehv}(:,1:max(dimToUse(2),4)),shiftEMG(usedTimeInds,:));

                    end
                end % of run CCA neur-musc shift block
                

                if runDyn

                    dynCtxFRs{iBehv} = frsOrigData(goodNeuronsCtx,:) + randn(length(goodNeuronsCtx),size(frsOrigData,2)) * dynNoiseStd;
                    pastDataCtx = dynCtxFRs{iBehv}(:,dynInds{iBehv}-1)';
                    currentDataCtx = dynCtxFRs{iBehv}(:,dynInds{iBehv})';

                    dynStrFRs{iBehv} = frsOrigData(goodNeuronsStr,:) + randn(length(goodNeuronsStr),size(frsOrigData,2)) * dynNoiseStd;
                    pastDataStr = dynStrFRs{iBehv}(:,dynInds{iBehv}-1)';
                    currentDataStr = dynStrFRs{iBehv}(:,dynInds{iBehv})';

                    testSizeFraction = 0.8;
                    ridgeKs = [1e4 1e2 1 1e-2 1e-4];

                    % Do dynamics modeling using PCA components rather than single neuron FRs
                    % [~, pcaTraj, vaf] = pca(frsData(goodNeuronsCtx,:)');
                    % usePCADims = find(cumsum(vaf)/sum(vaf)>0.9,1);
                    % behvCtxPCs = pcaTraj(usedDynInds,1:usePCADims);
                    % pastData = behvCtxPCs(1:end-1,:);
                    % currentData = behvCtxPCs(2:end,:);

                    [coeffsCtx, ridgeKCtx, perfCtx, predCtx] = ...
                        estimateLinearDynamics(pastDataCtx, currentDataCtx, testSizeFraction, ridgeKs, 5, dynNoId);
                    [coeffsStr, ridgeKStr, perfStr, predStr] = ...
                        estimateLinearDynamics(pastDataStr,currentDataStr,testSizeFraction,ridgeKs,5,dynNoId);

                    % save to different variables
                    if isShift
                        dynCoeffsCtxShift{iSess,iLabelType}{iBehv,iShift-1} = coeffsCtx;
                        dynBestRidgeKCtxShift{iSess,iLabelType}(iBehv,iShift-1) = ridgeKCtx;
                        dynPerformanceCtxShift{iSess,iLabelType}{iBehv,iShift-1} = perfCtx;

                        dynCoeffsStrShift{iSess,iLabelType}{iBehv,iShift-1} = coeffsStr;
                        dynBestRidgeKStrShift{iSess,iLabelType}(iBehv,iShift-1) = ridgeKStr;
                        dynPerformanceStrShift{iSess,iLabelType}{iBehv,iShift-1} = perfStr;
                    else
                        dynCoeffsCtx{iSess,iLabelType}{iBehv} = coeffsCtx;
                        dynBestRidgeKCtx{iSess,iLabelType}(iBehv) = ridgeKCtx;
                        dynPerformanceCtx{iSess,iLabelType}{iBehv} = perfCtx;

                        dynCoeffsStr{iSess,iLabelType}{iBehv} = coeffsStr;
                        dynBestRidgeKStr{iSess,iLabelType}(iBehv) = ridgeKStr;
                        dynPerformanceStr{iSess,iLabelType}{iBehv} = perfStr;
                    end

                end % of run dynamics block

                % also do it for each subregion individually
                if labelTypeRunSubregions && ~isShift
                    for iSubRegion = 1:length(usedSubregions{iBehv})

                        subTimeIndsToUse{iSess,iLabelType}{iBehv}{iSubRegion} = randperm(size(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)},2),subRegionLeastPoints(iBehv));
                        usedTimeInds = subTimeIndsToUse{iSess,iLabelType}{iBehv}{iSubRegion};

                        if runPCA
                            [pcaProjSubStr{iSess,iLabelType}{iBehv}{iSubRegion}, pcaTrajSubStr{iBehv}{iSubRegion}, varExpSubStr{iSess,iLabelType}{iBehv}{iSubRegion}] = ...
                                pca(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsStr,usedTimeInds)');
                            [pcaProjSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}, pcaTrajSubCtx{iBehv}{iSubRegion}, varExpSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}] = ...
                                pca(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsCtx,usedTimeInds)');
                            [pcaProjSubEmg{iSess,iLabelType}{iBehv}{iSubRegion}, pcaTrajSubEmg{iBehv}{iSubRegion}, varExpSubEmg{iSess,iLabelType}{iBehv}{iSubRegion}] = ...
                                pca(subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(:,usedTimeInds)');

                            %calculate dimensionality
                            if runDimCalc
                                [~, paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,1), cutoffDimSub{iSess,iLabelType}{iBehv}(iSubRegion,1)] = ...
                                    dimEst(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsStr,usedTimeInds)', 5);
                                [~, paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,2), cutoffDimSub{iSess,iLabelType}{iBehv}(iSubRegion,2)] = ...
                                    dimEst(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsCtx,usedTimeInds)', 5);
                                [~, paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,3), cutoffDimSub{iSess,iLabelType}{iBehv}(iSubRegion,3)] = ...
                                    dimEst(subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(:,usedTimeInds)', 5);
                            else
                                paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,1:3) = 5;
                            end
                        end

                        if runCCA
                            [ccaNeurProjSubStr{iSess,iLabelType}{iBehv}{iSubRegion}, ccaEMGProjSubStr{iSess,iLabelType}{iBehv}{iSubRegion}, cannonCorrsSubStr{iSess,iLabelType}{iBehv}(iSubRegion,:), ~, ~] = ...
                                canoncorr(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsStr,usedTimeInds)',...
                                subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(1:4,usedTimeInds)');
                            [ccaNeurProjSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}, ccaEMGProjSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}, cannonCorrsSubCtx{iSess,iLabelType}{iBehv}(iSubRegion,:), ~, ~] = ...
                                canoncorr(subRegionFRs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsCtx,usedTimeInds)',...
                                subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(1:4,usedTimeInds)');

                            %do CCA using only top PC's
                            [ccaNeurProjLowDimSubStr{iSess,iLabelType}{iBehv}{iSubRegion}, ccaEMGProjLowDimSubStr{iSess,iLabelType}{iBehv}{iSubRegion}, cannonCorrsLowDimSubStr{iSess,iLabelType}{iBehv}(iSubRegion,:), ~, ~] = ...
                                canoncorr(pcaTrajSubStr{iBehv}{iSubRegion}(:,1:max(paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,1),4)),...
                                subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(1:4,usedTimeInds)');
                            [ccaNeurProjLowDimSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}, ccaEMGProjLowDimSubCtx{iSess,iLabelType}{iBehv}{iSubRegion}, cannonCorrsLowDimSubCtx{iSess,iLabelType}{iBehv}(iSubRegion,:), ~, ~] = ...
                                canoncorr(pcaTrajSubCtx{iBehv}{iSubRegion}(:,1:max(paDimSub{iSess,iLabelType}{iBehv}(iSubRegion,2),4)),...
                                subRegionEMGs{iBehv}{usedSubregions{iBehv}(iSubRegion)}(1:4,usedTimeInds)');
                        end

                        if runDyn

                            subDynTimeIndsToUse{iSess,iLabelType}{iBehv}{iSubRegion} = randperm(length(usedDynSubRegionInds{iBehv}{usedSubregions{iBehv}(iSubRegion)}),subRegionDynLeastPoints(iBehv));
                            dynSubInds{iBehv}{iSubRegion} = usedDynSubRegionInds{iBehv}{usedSubregions{iBehv}(iSubRegion)}(subDynTimeIndsToUse{iSess,iLabelType}{iBehv}{iSubRegion});

                            dynCtxSubFRs{iBehv}{iSubRegion} = subRegionFRsOrig{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsCtx,:) + ...
                                randn(length(goodNeuronsCtx), size(subRegionFRsOrig{iBehv}{usedSubregions{iBehv}(iSubRegion)},2)) * dynNoiseStd;
                            pastDataCtx = dynCtxSubFRs{iBehv}{iSubRegion}(:,dynSubInds{iBehv}{iSubRegion}-1)';
                            currentDataCtx = dynCtxSubFRs{iBehv}{iSubRegion}(:,dynSubInds{iBehv}{iSubRegion})';

                            dynStrSubFRs{iBehv}{iSubRegion} = subRegionFRsOrig{iBehv}{usedSubregions{iBehv}(iSubRegion)}(goodNeuronsStr,:) + ...
                                randn(length(goodNeuronsStr),size(subRegionFRsOrig{iBehv}{usedSubregions{iBehv}(iSubRegion)},2)) * dynNoiseStd;
                            pastDataStr = dynStrSubFRs{iBehv}{iSubRegion}(:,dynSubInds{iBehv}{iSubRegion}-1)';
                            currentDataStr = dynStrSubFRs{iBehv}{iSubRegion}(:,dynSubInds{iBehv}{iSubRegion})';

                            testSizeFraction = 0.8;
                            ridgeKs = [1e6 1e4 1e2 1];

                            % Do dynamics modeling using PCA components rather than single neuron FRs
                            % [~, pcaTraj, vaf] = pca(frsData(goodNeuronsCtx,:)');
                            % usePCADims = find(cumsum(vaf)/sum(vaf)>0.9,1);
                            % behvCtxPCs = pcaTraj(usedDynInds,1:usePCADims);
                            % pastData = behvCtxPCs(1:end-1,:);
                            % currentData = behvCtxPCs(2:end,:);

                            [coeffsCtx, ridgeKCtx, perfCtx, predCtx] = ...
                                estimateLinearDynamics(pastDataCtx, currentDataCtx, testSizeFraction, ridgeKs, 5, dynNoId);
                            [coeffsStr, ridgeKStr, perfStr, predStr] = ...
                                estimateLinearDynamics(pastDataStr,currentDataStr,testSizeFraction,ridgeKs,5,dynNoId);

                            dynCoeffsSubCtx{iSess,iLabelType}{iBehv}{iSubRegion} = coeffsCtx;
                            dynBestRidgeKSubCtx{iSess,iLabelType}{iBehv}(iSubRegion) = ridgeKCtx;
                            dynPerformanceSubCtx{iSess,iLabelType}{iBehv}{iSubRegion} = perfCtx;

                            dynCoeffsSubStr{iSess,iLabelType}{iBehv}{iSubRegion} = coeffsStr;
                            dynBestRidgeKSubStr{iSess,iLabelType}{iBehv}(iSubRegion) = ridgeKStr;
                            dynPerformanceSubStr{iSess,iLabelType}{iBehv}{iSubRegion} = perfStr;

                        end

                    end % of subregions loop
                end % of subregions if block

                % % % do PCA and CCA with the behavior shift controls
                % % for iShift = 1:nShifts
                % %
                % %     if runPCA
                % %         % do PCA
                % %         shiftLeastPoints = min(cellfun(@(x) size(x,2),regionFRsShift(:,iShift)));
                % %         shiftTimeIndsToUse{iSess,iBehv,iShift} = randperm(size(regionFRsShift{iBehv,iShift},2),shiftLeastPoints);
                % %
                % %         [pcaProjStrShift{iSess,iBehv,iShift}, thisTrajStrShift, ~] = ...
                % %             pca(regionFRsShift{iBehv,iShift}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %         [pcaProjCtxShift{iSess,iBehv,iShift}, thisTrajCtxShift, ~] = ...
                % %             pca(regionFRsShift{iBehv,iShift}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %         [pcaProjEmgShift{iSess,iBehv,iShift}, ~, ~] = ...
                % %             pca(regionEMGsShift{iBehv,iShift}(:,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %
                % %         % to save time, just calculate the dimensionality of shift
                % %         % controls for just the first shift and use that for the rest
                % %         if iShift == 1
                % %             [lbmleDimShift(iBehv,iShift,1), paDimShift(iSess,iBehv,iShift,1), cutoffDimShift(iSess,iBehv,iShift,1)] = ...
                % %                 dimEst(regionFRsShift{iBehv,iShift}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iBehv,iShift})', 5);
                % %             [lbmleDimShift(iBehv,iShift,2), paDimShift(iSess,iBehv,iShift,2), cutoffDimShift(iSess,iBehv,iShift,2)] = ...
                % %                 dimEst(regionFRsShift{iBehv,iShift}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iBehv,iShift})', 5);
                % %             [lbmleDimShift(iBehv,iShift,3), paDimShift(iSess,iBehv,iShift,3), cutoffDimShift(iSess,iBehv,iShift,3)] = ...
                % %                 dimEst(regionEMGsShift{iBehv,iShift}(:,shiftTimeIndsToUse{iSess,iBehv,iShift})', 5);
                % %         end
                % %     end
                % %
                % %     if runCCA
                % %         % do CCA
                % %         [ccaNeurProjStrRegShift{iBehv,iShift},ccaEMGProjStrRegShift{iBehv,iShift},cannonCorrsStrRegShift(iSess,iBehv,iShift,:),~,~] = ...
                % %             canoncorr(regionFRsShift{iBehv,iShift}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iBehv,iShift})',...
                % %             regionEMGsShift{iBehv,iShift}(1:4,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %         [ccaNeurProjCtxRegShift{iBehv,iShift},ccaEMGProjCtxRegShift{iBehv,iShift},cannonCorrsCtxRegShift(iSess,iBehv,iShift,:),~,~] = ...
                % %             canoncorr(regionFRsShift{iBehv,iShift}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iBehv,iShift})',...
                % %             regionEMGsShift{iBehv,iShift}(1:4,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %
                % %         %                 %do CCA using only top PC's
                % %         %                 [ccaNeurProjLowDimStrRegShift{iBehv,iShift},ccaEMGProjLowDimStrRegShift{iBehv,iShift},cannonCorrsLowDimStrRegShift(iSess,iBehv,iShift,:),~,~] = ...
                % %         %                     canoncorr(thisTrajStrShift(:,1:paDimShift(iSess,iBehv,1,1)),...
                % %         %                     regionEMGsShift{iBehv,iShift}(1:4,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %         %                 [ccaNeurProjLowDimCtxRegShift{iBehv,iShift},ccaEMGProjLowDimCtxRegShift{iBehv,iShift},cannonCorrsLowDimCtxRegShift(iSess,iBehv,iShift,:),~,~] = ...
                % %         %                     canoncorr(thisTrajCtxShift(:,1:paDimShift(iSess,iBehv,1,2)),...
                % %         %                     regionEMGsShift{iBehv,iShift}(1:4,shiftTimeIndsToUse{iSess,iBehv,iShift})');
                % %     end
                % %
                % % end


                %         %do just correlations between neurons and EMG
                %         corrEMGAll{iBehv} = corr(regionFRs{iBehv}(reallyGoodNeuronsAll,timeIndsToUse{iBehv})',regionEMGs{iBehv}(:,timeIndsToUse{iBehv})');
                %         corrEMGStr{iBehv} = corr(regionFRs{iBehv}(reallyGoodNeuronsStr,timeIndsToUse{iBehv})',regionEMGs{iBehv}(:,timeIndsToUse{iBehv})');
                %         corrEMGCtx{iBehv} = corr(regionFRs{iBehv}(reallyGoodNeuronsCtx,timeIndsToUse{iBehv})',regionEMGs{iBehv}(:,timeIndsToUse{iBehv})');

                %do pairwise neuron correlations in addition to PCA
                %         corrMatrixStr{iBehv} = corr(regionFRs{iBehv}(goodNeuronsStr,timeIndsToUse{iSess,iLabelType}{iBehv})',regionFRs{iBehv}(goodNeuronsStr,timeIndsToUse{iSess,iLabelType}{iBehv})');
                %         corrMatrixStr{iBehv}(logical(diag(ones(length(corrMatrixStr{iBehv}),1),0))) = 0;
                %         corrMatrixCtx{iBehv} = corr(regionFRs{iBehv}(goodNeuronsCtx,timeIndsToUse{iSess,iLabelType}{iBehv})',regionFRs{iBehv}(goodNeuronsCtx,timeIndsToUse{iSess,iLabelType}{iBehv})');
                %         corrMatrixCtx{iBehv}(logical(diag(ones(length(corrMatrixCtx{iBehv}),1),0))) = 0;
                %         corrMatrixEmg{iBehv} = corr(regionEMGs{iBehv}(:,timeIndsToUse{iSess,iLabelType}{iBehv})',regionEMGs{iBehv}(:,timeIndsToUse{iSess,iLabelType}{iBehv})');
                %         corrMatrixEmg{iBehv}(logical(diag(ones(length(corrMatrixEmg{iBehv}),1),0))) = 0;

                disp(['Get region activity, Shift run ' num2str(iShift) ': ' num2str(toc)])

            end % of behvs loop

            alignmentDim = 10; %not sure what this is

            % now for each pair of regions, get alignment of PCA, CCA, and dynmaics
            plotExampleInfo.exampleSess = exampleSess;
            plotExampleInfo.exampleControl = exampleControl;
            plotExampleInfo.controlRun = iShift-1;
            plotExampleInfo.exampleBaseBehv = exampleBaseBehv;
            plotExampleInfo.exampleLowDivBehv= exampleLowDivBehv;
            plotExampleInfo.exampleHighDivBehv= exampleHighDivBehv;
            plotExampleInfo.exampleCCATimes = exampleCCATimes;
            

            inputVars.behvFRs = behvFRs;
            inputVars.behvEMGs = behvEMGs;
            inputVars.timeIndsToUse = timeIndsToUse{iSess,iLabelType};
            inputVars.goodNeuronsCtx = goodNeuronsCtx;
            inputVars.goodNeuronsStr = goodNeuronsStr;
            inputVars.nShuffs = nShifts;
            inputVars.isShift = isShift;
            inputVars.iSess = iSess;

            if isShift
                inputVars.runRotShuffs = false;
            else
                inputVars.runRotShuffs = true;
            end

            if runCCA
                if isShift
                    inputVars.ccaNeurProjCtx = ccaNeurProjCtxShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.ccaEMGProjCtx = ccaEMGProjCtxShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.cannonCorrsCtx = squeeze(cannonCorrsCtxShift{iSess,iLabelType}(:,iShift-1,:));
                    inputVars.ccaNeurProjStr = ccaNeurProjStrShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.ccaEMGProjStr = ccaEMGProjStrShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.cannonCorrsStr = squeeze(cannonCorrsStrShift{iSess,iLabelType}(:,iShift-1,:));
                else
                    inputVars.ccaNeurProjCtx = ccaNeurProjCtx{iSess,iLabelType};
                    inputVars.ccaEMGProjCtx = ccaEMGProjCtx{iSess,iLabelType};
                    inputVars.cannonCorrsCtx = cannonCorrsCtx{iSess,iLabelType};
                    inputVars.ccaNeurProjStr = ccaNeurProjStr{iSess,iLabelType};
                    inputVars.ccaEMGProjStr = ccaEMGProjStr{iSess,iLabelType};
                    inputVars.cannonCorrsStr = cannonCorrsStr{iSess,iLabelType};
                end
            end

            if runPCA
                inputVars.pcaTrajStr = pcaTrajStr;
                inputVars.pcaTrajCtx = pcaTrajCtx;
                inputVars.pcaTrajEmg = pcaTrajEmg;
                if isShift
                    inputVars.pcaProjCtx = pcaProjCtxShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.pcaProjStr = pcaProjStrShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.pcaProjEmg = pcaProjEmgShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.paDim = squeeze(paDimShift{iSess,iLabelType}(:,iShift-1,1:3));
                else
                    inputVars.pcaProjCtx = pcaProjCtx{iSess,iLabelType};
                    inputVars.pcaProjStr = pcaProjStr{iSess,iLabelType};
                    inputVars.pcaProjEmg = pcaProjEmg{iSess,iLabelType};
                    inputVars.paDim = paDim{iSess,iLabelType};
                end
            end

            if runDyn
                inputVars.dynInputDataCtx = dynCtxFRs;
                inputVars.dynInputDataStr = dynStrFRs;
                inputVars.dynNoId = dynNoId;
                inputVars.dynInds = dynInds;
                if isShift
                    inputVars.dynCoeffsCtx = dynCoeffsCtxShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.dynCoeffsStr = dynCoeffsStrShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.dynPerformanceCtx = dynPerformanceCtxShift{iSess,iLabelType}(:,iShift-1);
                    inputVars.dynPerformanceStr = dynPerformanceStrShift{iSess,iLabelType}(:,iShift-1);
                else
                    inputVars.dynCoeffsCtx = dynCoeffsCtx{iSess,iLabelType};
                    inputVars.dynCoeffsStr = dynCoeffsStr{iSess,iLabelType};
                    inputVars.dynPerformanceCtx = dynPerformanceCtx{iSess,iLabelType};
                    inputVars.dynPerformanceStr = dynPerformanceStr{iSess,iLabelType};
                end
            end

            outputVars = calcPairwiseSimilarity(inputVars,runCCA,runPCA,runDyn,plotExampleInfo);

            if runCCA
                if isShift
                    ccaPrinAngleShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaPrinAngleCtx;
                    ccaDiverShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDiverCtx;
                    ccaDropShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDropCtx;
                    ccaDropFracShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDropFracCtx;
                    ccaPrinAngleShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaPrinAngleStr;
                    ccaDiverShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDiverStr;
                    ccaDropShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDropStr;
                    ccaDropFracShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.ccaDropFracStr;
                else
                    ccaPrinAngleCtx{iLabelType}(iSess,:,:) = outputVars.ccaPrinAngleCtx;
                    ccaDiverCtx{iLabelType}(iSess,:,:) = outputVars.ccaDiverCtx;
                    ccaDropCtx{iLabelType}(iSess,:,:) = outputVars.ccaDropCtx;
                    ccaDropFracCtx{iLabelType}(iSess,:,:) = outputVars.ccaDropFracCtx;
                    ccaPrinAngleStr{iLabelType}(iSess,:,:) = outputVars.ccaPrinAngleStr;
                    ccaDiverStr{iLabelType}(iSess,:,:) = outputVars.ccaDiverStr;
                    ccaDropStr{iLabelType}(iSess,:,:) = outputVars.ccaDropStr;
                    ccaDropFracStr{iLabelType}(iSess,:,:) = outputVars.ccaDropFracStr;
                end
            end

            if runPCA
                if isShift
                    pcaDiverShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaDiverCtx;
                    pcaAlignShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAlignCtx;
                    pcaAngleShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAngleCtx;

                    pcaDiverShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaDiverStr;
                    pcaAlignShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAlignStr;
                    pcaAngleShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAngleStr;

                    pcaDiverShiftEmg{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaDiverEmg;
                    pcaAlignShiftEmg{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAlignEmg;
                    pcaAngleShiftEmg{iLabelType}(iSess,:,:,iShift-1) = outputVars.pcaAngleEmg;
                else
                    pcaDiverCtx{iLabelType}(iSess,:,:) = outputVars.pcaDiverCtx;
                    pcaAlignCtx{iLabelType}(iSess,:,:) = outputVars.pcaAlignCtx;
                    pcaAngleCtx{iLabelType}(iSess,:,:) = outputVars.pcaAngleCtx;

                    pcaDiverStr{iLabelType}(iSess,:,:) = outputVars.pcaDiverStr;
                    pcaAlignStr{iLabelType}(iSess,:,:) = outputVars.pcaAlignStr;
                    pcaAngleStr{iLabelType}(iSess,:,:) = outputVars.pcaAngleStr;

                    pcaDiverEmg{iLabelType}(iSess,:,:) = outputVars.pcaDiverEmg;
                    pcaAlignEmg{iLabelType}(iSess,:,:) = outputVars.pcaAlignEmg;
                    pcaAngleEmg{iLabelType}(iSess,:,:) = outputVars.pcaAngleEmg;

                    pcaDiverRotStr{iLabelType}(iSess,:,:,:) = outputVars.pcaDiverRotStr;
                    pcaDiverRotCtx{iLabelType}(iSess,:,:,:) = outputVars.pcaDiverRotCtx;
                    pcaDiverRotEmg{iLabelType}(iSess,:,:,:) = outputVars.pcaDiverRotEmg;
                    pcaAlignRotStr{iLabelType}(iSess,:,:,:) = outputVars.pcaAlignRotStr;
                    pcaAlignRotCtx{iLabelType}(iSess,:,:,:) = outputVars.pcaAlignRotCtx;
                    pcaAlignRotEmg{iLabelType}(iSess,:,:,:) = outputVars.pcaAlignRotEmg;
                    pcaAngleRotStr{iLabelType}(iSess,:,:,:) = outputVars.pcaAngleRotStr;
                    pcaAngleRotCtx{iLabelType}(iSess,:,:,:) = outputVars.pcaAngleRotCtx;
                    pcaAngleRotEmg{iLabelType}(iSess,:,:,:) = outputVars.pcaAngleRotEmg;
                end        
            end

            if runDyn
                if isShift
                    dynDropShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynDropCtx;
                    dynDropFracShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynDropFracCtx;
                    dynCCDropFracShiftCtx{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynCCDropFracCtx;
                    dynDropShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynDropStr;
                    dynDropFracShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynDropFracStr;
                    dynCCDropFracShiftStr{iLabelType}(iSess,:,:,iShift-1) = outputVars.dynCCDropFracStr;

                    dynNeursR2ChangeShiftCtx{iLabelType,iSess}(:,:,iShift-1) = outputVars.dynNeursR2ChangeCtx;
                    dynNeursR2ChangeShiftStr{iLabelType,iSess}(:,:,iShift-1) = outputVars.dynNeursR2ChangeStr;
                    dynNeursCCChangeShiftCtx{iLabelType,iSess}(:,:,iShift-1) = outputVars.dynNeursCCChangeCtx;
                    dynNeursCCChangeShiftStr{iLabelType,iSess}(:,:,iShift-1) = outputVars.dynNeursCCChangeStr;
                else
                    dynDropCtx{iLabelType}(iSess,:,:) = outputVars.dynDropCtx;
                    dynDropFracCtx{iLabelType}(iSess,:,:) = outputVars.dynDropFracCtx;
                    dynCCDropFracCtx{iLabelType}(iSess,:,:) = outputVars.dynCCDropFracCtx;
                    dynDropStr{iLabelType}(iSess,:,:) = outputVars.dynDropStr;
                    dynDropFracStr{iLabelType}(iSess,:,:) = outputVars.dynDropFracStr;
                    dynCCDropFracStr{iLabelType}(iSess,:,:) = outputVars.dynCCDropFracStr;

                    dynNeursR2ChangeCtx{iLabelType,iSess}(:,:) = outputVars.dynNeursR2ChangeCtx;
                    dynNeursR2ChangeStr{iLabelType,iSess}(:,:) = outputVars.dynNeursR2ChangeStr;
                    dynNeursCCChangeCtx{iLabelType,iSess}(:,:) = outputVars.dynNeursCCChangeCtx;
                    dynNeursCCChangeStr{iLabelType,iSess}(:,:) = outputVars.dynNeursCCChangeStr;
                end
            end

            % also calculate across subregions for each region
            for iBehv = 1:nBehvs

                if ~labelTypeRunSubregions || isempty(usedSubregions) ||isempty(usedSubregions{iBehv}) || isShift
                    continue
                end

                clear inputVars plotExampleInfo

                plotExampleInfo.exampleSess = nan;
                plotExampleInfo.exampleLowDivBehv = nan;
                plotExampleInfo.controlRun = nan;
                plotExampleInfo.exampleControl = nan;
                plotExampleInfo.exampleBaseBehv = nan;
                plotExampleInfo.exampleHighDivBehv = nan;
                plotExampleInfo.exampleCCATimes = nan;

                inputVars.behvFRs = subRegionFRs{iBehv}(usedSubregions{iBehv});
                inputVars.behvEMGs = subRegionEMGs{iBehv}(usedSubregions{iBehv});
                inputVars.timeIndsToUse = subTimeIndsToUse{iSess,iLabelType}{iBehv};
                inputVars.goodNeuronsCtx = goodNeuronsCtx;
                inputVars.goodNeuronsStr = goodNeuronsStr;
                inputVars.nShuffs = 0;
                inputVars.isShift = 0;
                inputVars.iSess = iSess;
                inputVars.runRotShuffs = false;

                if runCCA
                    inputVars.ccaNeurProjCtx = ccaNeurProjSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.ccaEMGProjCtx = ccaEMGProjSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.cannonCorrsCtx = cannonCorrsSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.ccaNeurProjStr = ccaNeurProjSubStr{iSess,iLabelType}{iBehv};
                    inputVars.ccaEMGProjStr = ccaEMGProjSubStr{iSess,iLabelType}{iBehv};
                    inputVars.cannonCorrsStr = cannonCorrsSubStr{iSess,iLabelType}{iBehv};
                end

                if runPCA
                    inputVars.pcaTrajStr = pcaTrajSubStr{iBehv};
                    inputVars.pcaTrajCtx = pcaTrajSubCtx{iBehv};
                    inputVars.pcaTrajEmg = pcaTrajSubEmg{iBehv};

                    inputVars.pcaProjCtx = pcaProjSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.pcaProjStr = pcaProjSubStr{iSess,iLabelType}{iBehv};
                    inputVars.pcaProjEmg = pcaProjSubEmg{iSess,iLabelType}{iBehv};
                    inputVars.paDim = paDimSub{iSess,iLabelType}{iBehv};
                end

                if runDyn
                    inputVars.dynInputDataCtx = dynCtxSubFRs{iBehv};
                    inputVars.dynInputDataStr = dynStrSubFRs{iBehv};
                    inputVars.dynNoId = dynNoId;
                    inputVars.dynInds = dynSubInds{iBehv};

                    inputVars.dynCoeffsCtx = dynCoeffsSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.dynCoeffsStr = dynCoeffsSubStr{iSess,iLabelType}{iBehv};
                    inputVars.dynPerformanceCtx = dynPerformanceSubCtx{iSess,iLabelType}{iBehv};
                    inputVars.dynPerformanceStr = dynPerformanceSubStr{iSess,iLabelType}{iBehv};
                end

                outputVars = calcPairwiseSimilarity(inputVars,runCCA,runPCA,runDyn,plotExampleInfo);

                if runCCA
                    ccaPrinAngleSubCtx{iSess,iBehv} = outputVars.ccaPrinAngleCtx;
                    ccaDiverSubCtx{iSess,iBehv} = outputVars.ccaDiverCtx;
                    ccaDropSubCtx{iSess,iBehv} = outputVars.ccaDropCtx;
                    ccaDropFracSubCtx{iSess,iBehv} = outputVars.ccaDropFracCtx;
                    ccaPrinAngleSubStr{iSess,iBehv} = outputVars.ccaPrinAngleStr;
                    ccaDiverSubStr{iSess,iBehv} = outputVars.ccaDiverStr;
                    ccaDropSubStr{iSess,iBehv} = outputVars.ccaDropStr;
                    ccaDropFracSubStr{iSess,iBehv} = outputVars.ccaDropFracStr;
                end

                if runPCA
                    pcaDiverSubCtx{iSess,iBehv} = outputVars.pcaDiverCtx;
                    pcaAlignSubCtx{iSess,iBehv} = outputVars.pcaAlignCtx;
                    pcaAngleSubCtx{iSess,iBehv} = outputVars.pcaAngleCtx;
                    pcaDiverSubStr{iSess,iBehv} = outputVars.pcaDiverStr;
                    pcaAlignSubStr{iSess,iBehv} = outputVars.pcaAlignStr;
                    pcaAngleSubStr{iSess,iBehv} = outputVars.pcaAngleStr;
                    pcaDiverSubEmg{iSess,iBehv} = outputVars.pcaDiverEmg;
                    pcaAlignSubEmg{iSess,iBehv} = outputVars.pcaAlignEmg;
                    pcaAngleSubEmg{iSess,iBehv} = outputVars.pcaAngleEmg;
                end

                if runDyn
                    dynDropSubCtx{iSess,iBehv} = outputVars.dynDropCtx;
                    dynDropFracSubCtx{iSess,iBehv} = outputVars.dynDropFracCtx;
                    dynCCDropFracSubCtx{iSess,iBehv} = outputVars.dynCCDropFracCtx;
                    dynDropSubStr{iSess,iBehv} = outputVars.dynDropStr;
                    dynDropFracSubStr{iSess,iBehv} = outputVars.dynDropFracStr;
                    dynCCDropFracSubStr{iSess,iBehv} = outputVars.dynCCDropFracStr;

                    dynNeursR2ChangeSubCtx{iSess,iBehv} = outputVars.dynNeursR2ChangeCtx;
                    dynNeursR2ChangeSubStr{iSess,iBehv} = outputVars.dynNeursR2ChangeStr;
                    dynNeursCCChangeSubCtx{iSess,iBehv} = outputVars.dynNeursCCChangeCtx;
                    dynNeursCCChangeSubStr{iSess,iBehv} = outputVars.dynNeursCCChangeStr;
                end

            end % of behvs loop for subregions

        end % of shifts loop

    end % of label types loop

end % of sessions loop


% for iSess = 1:3
% 
%     spreadCtx(iSess) = std(squareform(squeeze(pcaDiverCtx(iSess,:,:))));
%     spreadStr(iSess) = std(squareform(squeeze(pcaDiverStr(iSess,:,:))));
% 
%     for iBehv = 1:nBehvs
%         if length(pcaDiverSubCtx{iSess,iBehv}) >= 4
%             spreadSubCtx(iSess,iBehv) = std(squareform(squeeze(pcaDiverSubCtx{iSess,iBehv})));
%             spreadSubStr(iSess,iBehv) = std(squareform(squeeze(pcaDiverSubStr{iSess,iBehv})));
%         else
%             spreadSubCtx(iSess,iBehv) = nan;
%             spreadSubStr(iSess,iBehv) = nan;
%         end
%     end
% 
% end

% % figure;
% % plot(1:2,[spreadCtx; nanmean(spreadSubCtx,2)'],'o-','color',[0.5 0.5 0.5],'MarkerFaceColor','k')
% % hold on
% % plot(3:4,[spreadStr; nanmean(spreadSubStr,2)'],'o-','color',[0.5 0.5 0.5],'MarkerFaceColor','k')
% % line([0.8 1.2],repmat(mean(spreadCtx),1,2),'color','r')
% % line([1.8 2.2],repmat(mean(nanmean(spreadSubCtx,2)),1,2),'color','r')
% % line([2.8 3.2],repmat(mean(spreadStr),1,2),'color','r')
% % line([3.8 4.2],repmat(mean(nanmean(spreadSubStr,2)),1,2),'color','r')
% % set(gca,'xtick',1:4)
% % set(gca,'XTickLabel',{'Ctx inter-region','Ctx intra-region','Str inter-region','Str intra-region'})
% % set(gca,'TickDir','out')
% % set(gca,'FontSize',14)
% % ylabel('Hierarchy Spread (std)')
% % box off

commonSaveVars = {'binSize','useSmoothed','shiftAmount','subBehvMinPoints','behvsLessThanMinPoints','nShifts',...
    'minLabelPoints','labelTypeNames','goodSubregions','timeIndsToUse','shiftTimeIndsToUse'};

saveVarsPCA = {'pcaProjStr','varExpStr','pcaProjCtx','varExpCtx','pcaProjEmg','varExpEmg','paDim',
    'pcaProjStrShift','pcaProjCtxShift','pcaProjEmgShift','paDimShift','subTimeIndsToUse',...
    'pcaProjSubStr','varExpSubStr','pcaProjSubCtx','varExpSubCtx','pcaProjSubEmg','varExpSubEmg','paDimSub',...
    'pcaDiverStr','pcaAlignStr','pcaAngleStr','pcaDiverCtx','pcaAlignCtx','pcaAngleCtx','pcaDiverEmg','pcaAlignEmg','pcaAngleEmg',...
    'pcaDiverRotStr','pcaAlignRotStr','pcaAngleRotStr','pcaDiverRotCtx','pcaAlignRotCtx',...
    'pcaAngleRotCtx','pcaDiverRotEmg','pcaAlignRotEmg','pcaAngleRotEmg',...
    'pcaDiverShiftStr','pcaAlignShiftStr','pcaAngleShiftStr','pcaDiverShiftCtx','pcaAlignShiftCtx',...
    'pcaAngleShiftCtx','pcaDiverShiftEmg','pcaAlignShiftEmg','pcaAngleShiftEmg',...
    'pcaDiverSubStr','pcaAlignSubStr','pcaAngleSubStr','pcaDiverSubCtx','pcaAlignSubCtx',...
    'pcaAngleSubCtx','pcaDiverSubEmg','pcaAlignSubEmg','pcaAngleSubEmg'};

% saveVarsPCA = {'pcaProjAllRegionsStr','pcaTrajAllRegionsStr','varExpAllRegionsStr',...
%     'pcaProjAllRegionsCtx','pcaTrajAllRegionsCtx','varExpAllRegionsCtx',...
%     'pcaProjAllRegionsEmg','pcaTrajAllRegionsEmg','varExpAllRegionsEmg','timeIndsToUse','shiftTimeIndsToUse',...
%     'pcaProjStr','pcaTrajStr','varExpStr','pcaProjCtx','pcaTrajCtx','varExpCtx',...
%     'pcaProjEmg','pcaTrajEmg','varExpEmg','paDim','cutoffDim','pcaProjStrShift',...
%     'pcaProjCtxShift','pcaProjEmgShift','paDimShift','cutoffDimShift',...
%     'pcaDiverStr','pcaAlignStr','pcaAngleStr','pcaDiverCtx','pcaAlignCtx','pcaAngleCtx',...
%     'pcaDiverEmg','pcaAlignEmg','pcaAngleEmg','rotPaDim1Str','rotPaDim2Str',...
%     'pcaDiverRotStr','pcaAlignRotStr','pcaAngleRotStr','pcaDiverRotCtx','pcaAlignRotCtx','pcaAngleRotCtx',...
%     'pcaDiverRotEmg','pcaAlignRotEmg','pcaAngleRotEmg','pcaDiverShiftStr','pcaAlignShiftStr','pcaAngleShiftStr',...
%     'pcaDiverShiftCtx','pcaAlignShiftCtx','pcaAngleShiftCtx','pcaDiverShiftEmg','pcaAlignShiftEmg','pcaAngleShiftEmg'};


saveVarsDyn = {'dynTimeIndsToUse','dynShiftTimeIndsToUse','dynNoId','dynNoiseStd','ridgeKs',...
    'dynCoeffsCtx','dynBestRidgeKCtx','dynPerformanceCtx','dynCoeffsStr','dynBestRidgeKStr','dynPerformanceStr',...
    'dynCoeffsCtxShift','dynBestRidgeKCtxShift','dynPerformanceCtxShift','dynCoeffsStrShift','dynBestRidgeKStrShift','dynPerformanceStrShift',...
    'subDynTimeIndsToUse','dynCoeffsSubCtx','dynBestRidgeKSubCtx','dynPerformanceSubCtx','dynCoeffsSubStr','dynBestRidgeKSubStr','dynPerformanceSubStr',...
    'dynDropCtx','dynDropFracCtx','dynCCDropFracCtx','dynDropStr','dynDropFracStr','dynCCDropFracStr',...
    'dynNeursR2ChangeCtx','dynNeursR2ChangeStr','dynNeursCCChangeCtx','dynNeursCCChangeStr',...
    'dynDropShiftCtx','dynDropFracShiftCtx','dynCCDropFracShiftCtx','dynDropShiftStr','dynDropFracShiftStr','dynCCDropFracShiftStr',...
    'dynNeursR2ChangeShiftCtx','dynNeursR2ChangeShiftStr','dynNeursCCChangeShiftCtx','dynNeursCCChangeShiftStr',...
    'dynDropSubCtx','dynDropFracSubCtx','dynCCDropFracSubCtx','dynDropSubStr','dynDropFracSubStr','dynCCDropFracSubStr',...
    'dynNeursR2ChangeSubCtx','dynNeursR2ChangeSubStr','dynNeursCCChangeSubCtx','dynNeursCCChangeSubStr'};

if runPCA
    save('X:\David\AnalysesData\PCASubspaces.mat',[commonSaveVars saveVarsPCA],'-v7.3')
end

if runDyn
    save('X:\David\AnalysesData\LDSStability.mat',[commonSaveVars saveVarsDyn],'-v7.3')
end

saveVarsCCA = {'cannonCorrsCombStr','ccaNeurTrajCombCtx','cannonCorrsCombLowDimStr','cannonCorrsCombLowDimCtx',...
    'cannonCorrsCombStrShift','cannonCorrsCombCtxShift','cannonCorrsCombLowDimStrShift','cannonCorrsCombLowDimCtxShift',...
    'paDim','paDimShift','ccaNeurProjStr','ccaEMGProjStr','cannonCorrsStr','ccaNeurProjLowDimStr','ccaEMGProjLowDimStr','cannonCorrsLowDimStr',...
    'ccaNeurProjCtx','ccaEMGProjCtx','cannonCorrsCtx','ccaNeurProjLowDimCtx','ccaEMGProjLowDimCtx','cannonCorrsLowDimCtx',...
    'ccaNeurProjStrShift','ccaEMGProjStrShift','cannonCorrsStrShift','ccaNeurProjLowDimStrShift','ccaEMGProjLowDimStrShift','cannonCorrsLowDimStrShift',...
    'ccaNeurProjCtxShift','ccaEMGProjCtxShift','cannonCorrsCtxShift','ccaNeurProjLowDimCtxShift','ccaEMGProjLowDimCtxShift','cannonCorrsLowDimCtxShift',...
    'cannonCorrsStrShift','cannonCorrsCtxShift','cannonCorrsLowDimStrShift','cannonCorrsLowDimCtxShift',...
    'ccaNeurProjSubStr','ccaEMGProjSubStr','cannonCorrsSubStr','ccaNeurProjSubCtx','ccaEMGProjSubCtx','cannonCorrsSubCtx',...
    'ccaNeurProjLowDimSubStr','ccaEMGProjLowDimSubStr','cannonCorrsLowDimSubStr','ccaNeurProjLowDimSubCtx','ccaEMGProjLowDimSubCtx','cannonCorrsLowDimSubCtx'};

if runCCA
    save('X:\David\AnalysesData\CCAStability.mat',[commonSaveVars saveVarsCCA],'-v7.3')
end


% commonSaveVarsCCA = {'ccaNeurProjCombRegStr','ccaEMGProjCombRegStr','cannonCorrsCombRegStr','ccaNeurTrajCombRegStr',...
%     'ccaEMGTrajCombRegStr','ccaNeurProjCombRegCtx','ccaEMGProjCombRegCtx','cannonCorrsCombRegCtx','ccaNeurTrajCombRegCtx',...
%     'ccaEMGTrajCombRegCtx','ccaNeurProjCombRegLowDimStr','ccaEMGProjCombRegLowDimStr','cannonCorrsCombRegLowDimStr',...
%     'ccaNeurTrajCombRegLowDimStr','ccaEMGTrajCombRegLowDimStr','ccaNeurProjCombRegLowDimCtx','ccaEMGProjCombRegLowDimCtx',...
%     'cannonCorrsCombRegLowDimCtx','ccaNeurTrajCombRegLowDimCtx','ccaEMGTrajCombRegLowDimCtx','cannonCorrsCombRegStrShift',...
%     'cannonCorrsCombRegCtxShift','cannonCorrsCombRegLowDimStrShift','cannonCorrsCombRegLowDimCtxShift','ccaNeurProjStr',...
%     'ccaEMGProjStr','cannonCorrsStr','ccaNeurTrajStr','ccaEMGTrajStr','ccaNeurProjCtx','ccaEMGProjCtx','cannonCorrsCtx',...
%     'ccaNeurTrajCtx','ccaEMGTrajCtx','ccaNeurProjLowDimStr','ccaEMGProjLowDimStr','cannonCorrsLowDimStr','ccaNeurTrajLowDimStr',...
%     'ccaEMGTrajLowDimStr','ccaNeurProjLowDimCtx','ccaEMGProjLowDimCtx','cannonCorrsLowDimCtx','ccaNeurTrajLowDimCtx',...
%     'ccaEMGTrajLowDimCtx','timeIndsToUse','paDim','cutoffDim','cannonCorrsStrShift','cannonCorrsCtxShift',...
%     'cannonCorrsLowDimCtxShift','cannonCorrsLowDimStrShift','ccaPrinAngleCtx','ccaPrinAngleStr','ccaPrinAngleShiftCtx','ccaPrinAngleShiftStr',...
%     'ccaDiverCtx','ccaDiverStr','ccaDiverShiftCtx','ccaDiverShiftStr','ccaDropCtx','ccaDropStr','ccaDropShiftCtx','ccaDropShiftStr'};

switch lower(inputData)
    case 'humanannotated'
        %         save('PCASubspacesAnnotated','behvRegionLabels',commonSaveVarsPCA{:},'-v7.3');
        %         save('CCAAnalysisAnnotated','behvRegionLabels',commonSaveVarsCCA{:},'-v7.3')
    case 'umapregions'
        %         save('PCASubspaces','behvRegionLabels',commonSaveVarsPCA{:},'-v7.3');
        %         save('CCAAnalysis','behvRegionLabels',commonSaveVarsCCA{:},'-v7.3')
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


function [behvEMGInds, behvNeurInds, behvEMGs, behvFRs, behvFRsOrig] = getBehvData(behvAssignments, behvLabel,...
    binSize, useSmoothed, origDownsampEMGInd, sessionName, sessionArtifacts, normalizedEMG, normalizedFRs, allFRs)
%get the indices in the UMAP reduction for a behv

filePaths = getMouseDataNames(sessionName(1:4),sessionName,'CFA');
load(filePaths.VideoSyncFrames)

behvTimeInds = find(behvAssignments == behvLabel);
behvTimeInds(origDownsampEMGInd(behvTimeInds) > floor(frameEMGSamples{1}{end}(end)/20)) = [];

% Don't use artifact points
umapArtInds = sessionArtifacts.allArtUmapInds;
[~, behvArtInds] = intersect(behvTimeInds,umapArtInds);
behvTimeInds(behvArtInds) = [];

%get emg data for this behv
behvEMGInds = unique(round(origDownsampEMGInd(behvTimeInds)/binSize));
behvEMGInds(behvEMGInds==0) = [];

%get neural data for this behv
%sync umap(EMG) indices to neural indices
currentDir = pwd;
cd(filePaths.processedDataFolder)
behvNeurInds = round(NeurEMGSync(behvEMGInds*20*binSize,...
    frameEMGSamples, frameNeuropixelSamples, 'EMG')/(30*binSize));
cd(currentDir)

%don't use any time points past the video data
outOfBoundInds = find(behvNeurInds > floor(frameNeuropixelSamples{1}{end}(end)/(30*binSize))-1);

%also don't use points with nans or if the ind is 0
nanInds = find(isnan(behvNeurInds));
zeroInds = find(behvNeurInds==0);

%finally don't use artifact points (different from UMAP artifact points
%since I need to account for the gaussian smoothing)
if useSmoothed
    switch binSize
        case 100
            neurArtInds = find(ismember(behvNeurInds,sessionArtifacts.allArtNeurInds100ms50msSmooth));

        case 10
            neurArtInds = find(ismember(behvNeurInds,sessionArtifacts.allArtNeurInds10ms30msSmooth));

        case 1
            neurArtInds = find(ismember(behvNeurInds,sessionArtifacts.allArtNeurInds1ms10msSmooth));
    end
else
    neurArtInds = [];
end

%remove those bad time points
behvNeurInds([outOfBoundInds nanInds zeroInds neurArtInds]) = [];
behvEMGInds(:,[outOfBoundInds nanInds zeroInds neurArtInds]) = [];

behvFRs = normalizedFRs(:,behvNeurInds);
%save the un-normalized data as well for calculating firing rates
behvFRsOrig = allFRs(:,behvNeurInds);

behvEMGs = normalizedEMG(:,behvEMGInds);

%for removing with nan time points, don't use any neurons that are all nans
%for all points (which will cause us to remove all time points)
nonNanNeurons = find(~all(isnan(normalizedFRs),2));

% remove nan points
nanEMGInds = find(any(isnan(behvEMGs)));
nanFRInds = find(any(isnan(behvFRs(nonNanNeurons,:))));
behvEMGs(:,unique([nanEMGInds nanFRInds])) = [];
behvFRs(:,unique([nanEMGInds nanFRInds])) = [];
behvFRsOrig(:,unique([nanEMGInds nanFRInds])) = [];

behvEMGInds([nanEMGInds nanFRInds]) = [];
behvNeurInds([nanEMGInds nanFRInds]) = [];

end % of getBehvData function


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

end % of calcPCAAlignmentOld function


%


function [divMetric, elsayedAlign, prinAngle] = calcPCAAlignment(dat1,dat2,proj1,proj2,dim1,dim2,makePlots)

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

crossProjVars = diag(proj2(:,1:divDim1)'*cov(dat1)*proj2(:,1:divDim1));
selfProjVars = diag(proj1(:,1:divDim1)'*cov(dat1)*proj1(:,1:divDim1));
trajVars = var(traj1(:,1:divDim1));
nonZeroVarPCs = crossProjVars~=0 & selfProjVars~=0;
divergence1 = trajVars(nonZeroVarPCs)/sum(trajVars(nonZeroVarPCs)) * (abs(crossProjVars(nonZeroVarPCs) - selfProjVars(nonZeroVarPCs)) ./ (crossProjVars(nonZeroVarPCs) + selfProjVars(nonZeroVarPCs)));

crossProjVars = diag(proj1(:,1:divDim2)'*cov(dat2)*proj1(:,1:divDim2));
selfProjVars = diag(proj2(:,1:divDim2)'*cov(dat2)*proj2(:,1:divDim2));
trajVars = var(traj2(:,1:divDim2));
nonZeroVarPCs = crossProjVars~=0 & selfProjVars~=0;
divergence2 = trajVars(nonZeroVarPCs)/sum(trajVars(nonZeroVarPCs)) * (abs(crossProjVars(nonZeroVarPCs) - selfProjVars(nonZeroVarPCs)) ./ (crossProjVars(nonZeroVarPCs) + selfProjVars(nonZeroVarPCs)));

% dim1 = 10;
% dim2 = 10;
alignment1 = sum(trace(proj2(:,1:dim1)'*cov(dat1)*proj2(:,1:dim1))) / sum(trace(proj1(:,1:dim1)'*cov(dat1)*proj1(:,1:dim1)));
alignment2 = sum(trace(proj1(:,1:dim2)'*cov(dat2)*proj1(:,1:dim2))) / sum(trace(proj2(:,1:dim2)'*cov(dat2)*proj2(:,1:dim2)));

divMetric = mean([divergence1 divergence2]);
elsayedAlign = mean([alignment1 alignment2]);

[~,S,~] = svd(proj1(:,1:dim1)'*proj2(:,1:dim2));
prinAngle = acosd(S(1));

% prinAngle = subspace(proj1(:,1:dim1),proj2(:,1:dim2));

% make example plots
if makePlots
    plotColors = lines(6);
    figure;
    plot(diag(proj1(:,1:divDim1)'*cov(dat1)*proj1(:,1:divDim1))/sum(var(traj1(:,1:divDim1))),'color',plotColors(1,:))
    hold on
    plot(diag(proj2(:,1:divDim1)'*cov(dat1)*proj2(:,1:divDim1))/sum(var(traj1(:,1:divDim1))),'color',plotColors(2,:))
    box off
    set(gca,'XColor','k')
    set(gca,'YColor','k')
    set(gca,'linewidth',0.5)
    set(gca,'fontsize',6.5)
    set(gca,'TickDir','out')
    set(gca,'TickLength',[0.02 0.05])
    set(gca,'XTick',[0 25 50])
    set(gcf,'color','w')
    xlabel('Principal component')
    ylabel('Variance fraction')
    text(1,0.1,['Divergence = ' num2str(divergence1)])

end

end % of calcPCAAlignment function


%


function [divMetricRot, elsayedAlignRot, prinAngleRot, paDim1, paDim2] = randomRotationAlignment(pcaTraj1, pcaTraj2, pcaProj1, pcaProj2, dim1, dim2, makePlots)

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

[divMetricRot, elsayedAlignRot, prinAngleRot] = calcPCAAlignment(psuedoData1,psuedoData2,psuedoProj1,psuedoProj2,paDim1,paDim2,makePlots);

end % of randomRotationAlignment function


%


function [prinAngle, ccDiverg, ccDrop ccDropFrac] = calcCCAAlignmentMetrics(neur1,neur2,emg1,emg2,neurProj1,neurProj2,emgProj1,emgProj2,cc1,cc2)

%calc principle angle
%orthonormalize the projection matrices
neurProjOrtho1 = gson(neurProj1);
neurProjOrtho2 = gson(neurProj2);
[~, S, ~] = svd(neurProjOrtho1'*neurProjOrtho2);
Svalues = diag(S);
allPrinAngles = acosd(Svalues);

prinAngle = allPrinAngles(1);


%calc CC Divergence
centeredNeur1 = neur1 - mean(neur1,2);
centeredEMG1 = emg1 - mean(emg1,2);

trajCrossNeur1 = centeredNeur1'*neurProj2;
trajCrossEMG1 = centeredEMG1'*emgProj2;
crossCC1 = corr(trajCrossNeur1,trajCrossEMG1);
crossCC1 = abs(diag(crossCC1));
% cc1Ratio = sum(crossCC1)/sum(cc1);
cc1Ratio = (cc1/sum(cc1))'*(abs(crossCC1-cc1)./(crossCC1+cc1));
aveDrop1 = mean(cc1-crossCC1);
aveDropFrac1 = mean((cc1-crossCC1)./cc1);

centeredNeur2 = neur2 - mean(neur2,2);
centeredEMG2 = emg2 - mean(emg2,2);

trajCrossNeur2 = centeredNeur2'*neurProj1;
trajCrossEMG2 = centeredEMG2'*emgProj1;
crossCC2 = corr(trajCrossNeur2,trajCrossEMG2);
crossCC2 = abs(diag(crossCC2));
% cc2Ratio = sum(crossCC2)/sum(cc2);
cc2Ratio = (cc2/sum(cc2))'*(abs(crossCC2-cc2)./(crossCC2+cc2));
aveDrop2 = mean(cc2-crossCC2);
aveDropFrac2 = mean((cc2-crossCC2)./cc2);

ccDiverg = mean([cc1Ratio cc2Ratio]);
ccDrop = mean([aveDrop1 aveDrop2]);
ccDropFrac = mean([aveDropFrac1 aveDropFrac2]);

end % of calcCCAAlignmentMetrics function


%


function outputVars = calcPairwiseSimilarity(inputVars,runCCA,runPCA,runDyn,plotExampleInfo)

% process inputs
exampleSess = plotExampleInfo.exampleSess;
controlRun = plotExampleInfo.controlRun;
exampleControl = plotExampleInfo.exampleControl;
exampleBaseBehv = plotExampleInfo.exampleBaseBehv;
exampleLowDivBehv = plotExampleInfo.exampleLowDivBehv;
exampleHighDivBehv = plotExampleInfo.exampleHighDivBehv;
exampleCCATimes = plotExampleInfo.exampleCCATimes;

behvFRs = inputVars.behvFRs;
behvEMGs = inputVars.behvEMGs;
timeIndsToUse = inputVars.timeIndsToUse;
goodNeuronsCtx = inputVars.goodNeuronsCtx;
goodNeuronsStr = inputVars.goodNeuronsStr;
nShuffs = inputVars.nShuffs;
isShift = inputVars.isShift;
iSess = inputVars.iSess;
runRotShuffs = inputVars.runRotShuffs;

if runCCA
    ccaNeurProjCtx = inputVars.ccaNeurProjCtx;
    ccaEMGProjCtx = inputVars.ccaEMGProjCtx;
    cannonCorrsCtx = inputVars.cannonCorrsCtx;
    ccaNeurProjStr = inputVars.ccaNeurProjStr;
    ccaEMGProjStr = inputVars.ccaEMGProjStr;
    cannonCorrsStr = inputVars.cannonCorrsStr;
    % if inputVars.nShifts > 0
    %     ccaNeurProjCtxRegShift = inputVars.ccaNeurProjCtxRegShift;
    %     ccaEMGProjCtxRegShift = inputVars.ccaEMGProjCtxRegShift;
    %     cannonCorrsCtxRegShift = inputVars.cannonCorrsCtxRegShift;
    %     ccaNeurProjStrRegShift = inputVars.ccaNeurProjStrRegShift;
    %     ccaEMGProjStrRegShift = inputVars.ccaEMGProjStrRegShift;
    %     cannonCorrsStrRegShift = inputVars.cannonCorrsStrRegShift;
    % end
end

if runPCA
    pcaProjCtx = inputVars.pcaProjCtx;
    pcaProjStr = inputVars.pcaProjStr;
    pcaProjEmg = inputVars.pcaProjEmg;
    paDim = inputVars.paDim;
    pcaTrajStr = inputVars.pcaTrajStr;
    pcaTrajCtx = inputVars.pcaTrajCtx;
    pcaTrajEmg = inputVars.pcaTrajEmg;
    % if inputVars.nShifts > 0
    %     behvFRsShift = inputVars.behvFRsShift;
    %     behvEMGsShift = inputVars.behvEMGsShift;
    %     shiftTimeIndsToUse = inputVars.shiftTimeIndsToUse;
    %     pcaProjStrShift = inputVars.pcaProjStrShift;
    %     pcaProjCtxShift = inputVars.pcaProjCtxShift;
    %     pcaProjEmgShift = inputVars.pcaProjEmgShift;
    %     paDimShift = inputVars.paDimShift;
    % end
end

if runDyn
    dynInds = inputVars.dynInds;
    dynCoeffsCtx = inputVars.dynCoeffsCtx;
    dynCoeffsStr = inputVars.dynCoeffsStr;
    dynInputDataCtx = inputVars.dynInputDataCtx;
    dynInputDataStr = inputVars.dynInputDataStr;
    dynPerformanceCtx = inputVars.dynPerformanceCtx;
    dynPerformanceStr = inputVars.dynPerformanceStr;
end

nBehvs = length(behvFRs);
crossBehvAllTimes = false;

% go through each pair of behaviors
for iBehv1 = 1:nBehvs
    for iBehv2 = 1:nBehvs

        if iBehv1 < iBehv2 %since the metrics are symmetric, save time by just calculating upper triangle

            tic

            if iSess == exampleSess && iBehv1 == exampleBaseBehv && controlRun == exampleControl && ...
                    (iBehv2 == exampleLowDivBehv || iBehv2 == exampleHighDivBehv)
                makeExamplePlot = true;
            else
                makeExamplePlot = false;
            end

            if runCCA

                if crossBehvAllTimes
                    behv1UseTimes = 1:size(behvFRs{iBehv1},2);
                    behv2UseTimes = 1:size(behvFRs{iBehv2},2);
                else
                    behv1UseTimes = timeIndsToUse{iBehv1};
                    behv2UseTimes = timeIndsToUse{iBehv2};
                end

                [ccaPrinAngleCtx(iBehv1,iBehv2), ccaDiverCtx(iBehv1,iBehv2), ccaDropCtx(iBehv1,iBehv2), ccaDropFracCtx(iBehv1,iBehv2)] = calcCCAAlignmentMetrics(...
                    behvFRs{iBehv1}(goodNeuronsCtx,behv1UseTimes),behvFRs{iBehv2}(goodNeuronsCtx,behv2UseTimes),...
                    behvEMGs{iBehv1}(1:4,behv1UseTimes),behvEMGs{iBehv2}(1:4,behv2UseTimes),...
                    ccaNeurProjCtx{iBehv1},ccaNeurProjCtx{iBehv2},...
                    ccaEMGProjCtx{iBehv1},ccaEMGProjCtx{iBehv2},squeeze(cannonCorrsCtx(iBehv1,:))',squeeze(cannonCorrsCtx(iBehv2,:))');

                [ccaPrinAngleStr(iBehv1,iBehv2), ccaDiverStr(iBehv1,iBehv2), ccaDropStr(iBehv1,iBehv2), ccaDropFracStr(iBehv1,iBehv2)] = calcCCAAlignmentMetrics(...
                    behvFRs{iBehv1}(goodNeuronsStr,behv1UseTimes),behvFRs{iBehv2}(goodNeuronsStr,behv2UseTimes),...
                    behvEMGs{iBehv1}(1:4,behv1UseTimes),behvEMGs{iBehv2}(1:4,behv2UseTimes),...
                    ccaNeurProjStr{iBehv1},ccaNeurProjStr{iBehv2},...
                    ccaEMGProjStr{iBehv1},ccaEMGProjStr{iBehv2},squeeze(cannonCorrsStr(iBehv1,:))',squeeze(cannonCorrsStr(iBehv2,:))');


                % make example plots
                if makeExamplePlot

                    exampleNeur = behvFRs{iBehv1}(goodNeuronsCtx,:);
                    exampleEmg = behvEMGs{iBehv1}(1:4,:);
                    centeredNeur = exampleNeur - mean(exampleNeur,2);
                    centeredEMG = exampleEmg - mean(exampleEmg,2);

                    exampleNeurTraj = centeredNeur'*ccaNeurProjCtx{iBehv1};
                    exampleEmgTraj = centeredEMG'*ccaEMGProjCtx{iBehv1};

                    exampleCrossNeurTraj = centeredNeur'*ccaNeurProjCtx{iBehv2};
                    exampleCrossEmgTraj = centeredEMG'*ccaEMGProjCtx{iBehv2};

                    if ~exist('exampleCCAFigH')
                        exampleCCAFigH = figure('color','w');
                        tiledlayout(2,2)
                    else
                        figure(exampleCCAFigH)
                    end

                    if iBehv2 == exampleLowDivBehv
                        nexttile
                        plot(exampleEmgTraj(exampleCCATimes,1),'color',[0.5 0.5 0.5],'linewidth',0.75);
                        hold on
                        plot(exampleNeurTraj(exampleCCATimes,1),'color',lines(1),'LineWidth',0.75)
                        line([0 200],[-3 -3],'color','k','linewidth',2)
                        text(10,3,['CC = ' num2str(corr(exampleNeurTraj(:,1),exampleEmgTraj(:,1)))])
                        ylim([-5 5])
                        axis off

                        crossPlotColor = [1 0 0];
                    else
                        crossPlotColor = [1 0 1];
                    end

                    nexttile
                    plot(exampleCrossEmgTraj(exampleCCATimes,1),'color',[0.5 0.5 0.5],'linewidth',0.75);
                    hold on
                    plot(exampleCrossNeurTraj(exampleCCATimes,1))
                    line([0 200],[-3 -3],'color','k','linewidth',2)
                    text(10,3,['CC = ' num2str(corr(exampleCrossNeurTraj(:,1),exampleCrossEmgTraj(:,1)))])
                    ylim([-5 5])
                    axis off

                end

            end % of CCA block


            if runPCA

                [pcaDiverStr(iBehv1,iBehv2), pcaAlignStr(iBehv1,iBehv2), pcaAngleStr(iBehv1,iBehv2)] = ...
                    calcPCAAlignment(behvFRs{iBehv1}(goodNeuronsStr,timeIndsToUse{iBehv1})',...
                    behvFRs{iBehv2}(goodNeuronsStr,timeIndsToUse{iBehv2})',pcaProjStr{iBehv1},pcaProjStr{iBehv2},...
                    paDim(iBehv1,1),paDim(iBehv2,1),makeExamplePlot);

                [pcaDiverCtx(iBehv1,iBehv2), pcaAlignCtx(iBehv1,iBehv2), pcaAngleCtx(iBehv1,iBehv2)] = ...
                    calcPCAAlignment(behvFRs{iBehv1}(goodNeuronsCtx,timeIndsToUse{iBehv1})',...
                    behvFRs{iBehv2}(goodNeuronsCtx,timeIndsToUse{iBehv2})',pcaProjCtx{iBehv1},pcaProjCtx{iBehv2},...
                    paDim(iBehv1,2),paDim(iBehv2,2),makeExamplePlot);

                [pcaDiverEmg(iBehv1,iBehv2), pcaAlignEmg(iBehv1,iBehv2), pcaAngleEmg(iBehv1,iBehv2)] = ...
                    calcPCAAlignment(behvEMGs{iBehv1}(:,timeIndsToUse{iBehv1})',...
                    behvEMGs{iBehv2}(:,timeIndsToUse{iBehv2})',pcaProjEmg{iBehv1},pcaProjEmg{iBehv2},...
                    paDim(iBehv1,3),paDim(iBehv2,3),false);

            end % of PCA block

            if runDyn

                % first do cortex
                inputDataCtx1 = dynInputDataCtx{iBehv1}(:,dynInds{iBehv1}-1)';
                inputDataCtx2 = dynInputDataCtx{iBehv2}(:,dynInds{iBehv2}-1)';
                outputDataCtx1 = dynInputDataCtx{iBehv1}(:,dynInds{iBehv1})';
                outputDataCtx2 = dynInputDataCtx{iBehv2}(:,dynInds{iBehv2})';

                predictionsCtx1 = [];
                predictionsCtx2 = [];
                % for iChan = 1:size(inputDataCtx1,2)
                %     usedChans = 1:1:size(inputDataCtx1,2);
                %     if dynNoId
                %         usedChans(iChan) = [];
                %     end
                %     predictionsCtx1(:,iChan) = [ones(size(inputDataCtx1,1),1) inputDataCtx1(:,usedChans)-mean(inputDataCtx1(:,usedChans))] * arCoeffsCtx{iBehv2}(:,iChan);
                %     predictionsCtx2(:,iChan) = [ones(size(inputDataCtx2,1),1) inputDataCtx2(:,usedChans)-mean(inputDataCtx2(:,usedChans))] * arCoeffsCtx{iBehv1}(:,iChan);
                % end
                predictionsCtx1 = [ones(size(inputDataCtx1,1),1) inputDataCtx1-mean(inputDataCtx1)] * dynCoeffsCtx{iBehv2};
                predictionsCtx2 = [ones(size(inputDataCtx2,1),1) inputDataCtx2-mean(inputDataCtx2)] * dynCoeffsCtx{iBehv1};

                performance1 = calcPerformanceMetrics(predictionsCtx1',(outputDataCtx1-mean(outputDataCtx1))');
                performance2 = calcPerformanceMetrics(predictionsCtx2',(outputDataCtx2-mean(outputDataCtx2))');

                dynNeursR2ChangeCtx{iBehv1,iBehv2} = [performance1.R2 - dynPerformanceCtx{iBehv1}.R2, performance2.R2 - dynPerformanceCtx{iBehv2}.R2];
                
                dynPopR2ChangeCtx(iBehv1,iBehv2) = mean([performance1.R2Weighted - dynPerformanceCtx{iBehv1}.R2Weighted, ...
                    performance2.R2Weighted - dynPerformanceCtx{iBehv2}.R2Weighted]);
                dynPopR2FracChangeCtx(iBehv1,iBehv2) = mean([(performance1.R2Weighted - dynPerformanceCtx{iBehv1}.R2Weighted)/dynPerformanceCtx{iBehv1}.R2Weighted, ...
                    (performance2.R2Weighted - dynPerformanceCtx{iBehv2}.R2Weighted)/dynPerformanceCtx{iBehv2}.R2Weighted]);

                dynNeursCCChangeCtx{iBehv1,iBehv2} = [performance1.CC - dynPerformanceCtx{iBehv1}.CC, performance2.CC - dynPerformanceCtx{iBehv2}.CC];

                dynPopCCFracChangeCtx(iBehv1,iBehv2) = mean([mean(performance1.CC - dynPerformanceCtx{iBehv1}.CC)/mean(dynPerformanceCtx{iBehv1}.CC), ...
                    mean(performance2.CC - dynPerformanceCtx{iBehv2}.CC)/mean(dynPerformanceCtx{iBehv2}.CC)]);

                % next do striatum
                inputDataStr1 = dynInputDataStr{iBehv1}(:,dynInds{iBehv1}-1)';
                inputDataStr2 = dynInputDataStr{iBehv2}(:,dynInds{iBehv2}-1)';
                outputDataStr1 = dynInputDataStr{iBehv1}(:,dynInds{iBehv1})';
                outputDataStr2 = dynInputDataStr{iBehv2}(:,dynInds{iBehv2})';

                predictionsStr1 = [];
                predictionsStr2 = [];
                % for iChan = 1:size(inputDataStr1,2)
                %     usedChans = 1:1:size(inputDataStr1,2);
                %     if dynNoId
                %         usedChans(iChan) = [];
                %     end
                %     predictionsStr1(:,iChan) = [ones(size(inputDataStr1,1),1) inputDataStr1(:,usedChans)] * dynCoeffsStr{iBehv2}(:,iChan);
                %     predictionsStr2(:,iChan) = [ones(size(inputDataStr2,1),1) inputDataStr2(:,usedChans)] * dynCoeffsStr{iBehv1}(:,iChan);
                % end
                predictionsStr1 = [ones(size(inputDataStr1,1),1) inputDataStr1-mean(inputDataStr1)] * dynCoeffsStr{iBehv2};
                predictionsStr2 = [ones(size(inputDataStr2,1),1) inputDataStr2-mean(inputDataStr2)] * dynCoeffsStr{iBehv1};

                performance1 = calcPerformanceMetrics(predictionsStr1',(outputDataStr1-mean(outputDataStr1))');
                performance2 = calcPerformanceMetrics(predictionsStr2',(outputDataStr2-mean(outputDataStr2))');

                dynNeursR2ChangeStr{iBehv1,iBehv2} = [performance1.R2 - dynPerformanceStr{iBehv1}.R2, performance2.R2 - dynPerformanceStr{iBehv2}.R2];

                dynPopR2ChangeStr(iBehv1,iBehv2) = mean([performance1.R2Weighted - dynPerformanceStr{iBehv1}.R2Weighted, ...
                    performance2.R2Weighted - dynPerformanceStr{iBehv2}.R2Weighted]);
                dynPopR2FracChangeStr(iBehv1,iBehv2) = mean([(performance1.R2Weighted - dynPerformanceStr{iBehv1}.R2Weighted)/dynPerformanceStr{iBehv1}.R2Weighted, ...
                    (performance2.R2Weighted - dynPerformanceStr{iBehv2}.R2Weighted)/dynPerformanceStr{iBehv2}.R2Weighted]);

                dynNeursCCChangeStr{iBehv1,iBehv2} = [performance1.CC - dynPerformanceStr{iBehv1}.CC, performance2.CC - dynPerformanceStr{iBehv2}.CC];

                dynPopCCFracChangeStr(iBehv1,iBehv2) = mean([mean(performance1.CC - dynPerformanceStr{iBehv1}.CC)/mean(dynPerformanceStr{iBehv1}.CC), ...
                    mean(performance2.CC - dynPerformanceStr{iBehv2}.CC)/mean(dynPerformanceStr{iBehv2}.CC)]);

            end % of dynamics block

            %do shuff controls
            allRegDataStr = [behvFRs{iBehv1}(goodNeuronsStr,timeIndsToUse{iBehv1}) behvFRs{iBehv2}(goodNeuronsStr,timeIndsToUse{iBehv2})];
            allRegDataCtx = [behvFRs{iBehv1}(goodNeuronsCtx,timeIndsToUse{iBehv1}) behvFRs{iBehv2}(goodNeuronsCtx,timeIndsToUse{iBehv2})];
            allRegDataEmg = [behvEMGs{iBehv1}(:,timeIndsToUse{iBehv1}) behvEMGs{iBehv2}(:,timeIndsToUse{iBehv2})];

            if ~isShift  && runRotShuffs
                for iShuff = 1:nShuffs

                    % only make example plot for one of them
                    if iSess == exampleSess && iBehv1 == exampleBaseBehv && iShuff == exampleControl && ...
                            (iBehv2 == exampleLowDivBehv || iBehv2 == exampleHighDivBehv)
                        makeExamplePlotShuff = true;
                    else
                        makeExamplePlotShuff = false;
                    end

                    if runPCA

                        %for orth shuff, just do random rotation
                        individualProjsRot = 1;

                        if individualProjsRot
                            rotTrajStr1 = pcaTrajStr{iBehv1};
                            rotTrajStr2 = pcaTrajStr{iBehv2};
                            rotTrajCtx1 = pcaTrajCtx{iBehv1};
                            rotTrajCtx2 = pcaTrajCtx{iBehv2};
                            rotTrajEmg1 = pcaTrajEmg{iBehv1};
                            rotTrajEmg2 = pcaTrajEmg{iBehv2};

                            rotProjStr1 = pcaProjStr{iBehv1};
                            rotProjStr2 = pcaProjStr{iBehv2};
                            rotProjCtx1 = pcaProjCtx{iBehv1};
                            rotProjCtx2 = pcaProjCtx{iBehv2};
                            rotProjEmg1 = pcaProjEmg{iBehv1};
                            rotProjEmg2 = pcaProjEmg{iBehv2};
                        else
                            %                       rotTrajStr1 = pcaTrajAllBehvsStr(allBehvsInds(iBehv1)+1 : allBehvsInds(iBehv1+1),:);
                            %                       rotTrajStr2 = pcaTrajAllBehvsStr(allBehvsInds(iBehv2)+1 : allBehvsInds(iBehv2+1),:);
                            %                       rotTrajCtx1 = pcaTrajAllBehvsCtx(allBehvsInds(iBehv1)+1 : allBehvsInds(iBehv1+1),:);
                            %                       rotTrajCtx2 = pcaTrajAllBehvsCtx(allBehvsInds(iBehv2)+1 : allBehvsInds(iBehv2+1),:);
                            %                       rotTrajEmg1 = pcaTrajAllBehvsEmg(allBehvsInds(iBehv1)+1 : allBehvsInds(iBehv1+1),:);
                            %                       rotTrajEmg2 = pcaTrajAllBehvsEmg(allBehvsInds(iBehv2)+1 : allBehvsInds(iBehv2+1),:);

                            randInds1 = randperm(size(pcaTrajAllBehvsStr{iSess},1),leastPoints);
                            randInds2 = randperm(size(pcaTrajAllBehvsStr{iSess},1),leastPoints);

                            rotTrajStr1 = pcaTrajAllBehvsStr{iSess}(randInds1,:);
                            rotTrajStr2 = pcaTrajAllBehvsStr{iSess}(randInds2,:);
                            rotTrajCtx1 = pcaTrajAllBehvsCtx{iSess}(randInds1,:);
                            rotTrajCtx2 = pcaTrajAllBehvsCtx{iSess}(randInds2,:);
                            rotTrajEmg1 = pcaTrajAllBehvsEmg{iSess}(randInds1,:);
                            rotTrajEmg2 = pcaTrajAllBehvsEmg{iSess}(randInds2,:);

                            rotProjStr1 = pcaProjAllBehvsStr{iSess};
                            rotProjStr2 = pcaProjAllBehvsStr{iSess};
                            rotProjCtx1 = pcaProjAllBehvsCtx{iSess};
                            rotProjCtx2 = pcaProjAllBehvsCtx{iSess};
                            rotProjEmg1 = pcaProjAllBehvsEmg{iSess};
                            rotProjEmg2 = pcaProjAllBehvsEmg{iSess};

                        end

                        if iShuff == 1
                            %since calculating the pa dimensionality takes so
                            %long, just calculate it once for the first
                            %shuffle, since it doesn't change much across
                            %shuffles anyways
                            [pcaDiverRotStr(iBehv1,iBehv2,iShuff), pcaAlignRotStr(iBehv1,iBehv2,iShuff), pcaAngleRotStr(iBehv1,iBehv2,iShuff), rotPaDim1Str, rotPaDim2Str] = ...
                                randomRotationAlignment(rotTrajStr1, rotTrajStr2, rotProjStr1, rotProjStr2,[],[],makeExamplePlotShuff);

                            [pcaDiverRotCtx(iBehv1,iBehv2,iShuff), pcaAlignRotCtx(iBehv1,iBehv2,iShuff), pcaAngleRotCtx(iBehv1,iBehv2,iShuff), rotPaDim1Ctx, rotPaDim2Ctx] = ...
                                randomRotationAlignment(rotTrajCtx1, rotTrajCtx2, rotProjCtx1, rotProjCtx2,[],[],makeExamplePlotShuff);

                            [pcaDiverRotEmg(iBehv1,iBehv2,iShuff), pcaAlignRotEmg(iBehv1,iBehv2,iShuff), pcaAngleRotEmg(iBehv1,iBehv2,iShuff), rotPaDim1Emg, rotPaDim2Emg] = ...
                                randomRotationAlignment(rotTrajEmg1, rotTrajEmg2, rotProjEmg1, rotProjEmg2,[],[],false);
                        else
                            [pcaDiverRotStr(iBehv1,iBehv2,iShuff), pcaAlignRotStr(iBehv1,iBehv2,iShuff), pcaAngleRotStr(iBehv1,iBehv2,iShuff)] = ...
                                randomRotationAlignment(rotTrajStr1, rotTrajStr2, rotProjStr1, rotProjStr2,rotPaDim1Str,rotPaDim2Str,false);

                            [pcaDiverRotCtx(iBehv1,iBehv2,iShuff), pcaAlignRotCtx(iBehv1,iBehv2,iShuff), pcaAngleRotCtx(iBehv1,iBehv2,iShuff)] = ...
                                randomRotationAlignment(rotTrajCtx1, rotTrajCtx2, rotProjCtx1, rotProjCtx2,rotPaDim1Ctx,rotPaDim2Ctx,false);

                            [pcaDiverRotEmg(iBehv1,iBehv2,iShuff), pcaAlignRotEmg(iBehv1,iBehv2,iShuff), pcaAngleRotEmg(iBehv1,iBehv2,iShuff)] = ...
                                randomRotationAlignment(rotTrajEmg1, rotTrajEmg2, rotProjEmg1, rotProjEmg2,rotPaDim1Emg,rotPaDim2Emg,false);
                        end

                        % % % now for behavioral label shifts
                        % % if iShuff ~= 1
                        % %     makeExamplePlot = false;
                        % % end
                        % % [pcaDiverShiftStr(iBehv1,iBehv2,iShuff), pcaAlignShiftStr(iBehv1,iBehv2,iShuff), pcaAngleShiftStr(iBehv1,iBehv2,iShuff)] = ...
                        % %     calcPCAAlignment(behvFRsShift{iBehv1,iShuff}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iBehv1,iShuff})',...
                        % %     behvFRsShift{iBehv2,iShuff}(goodNeuronsStr,shiftTimeIndsToUse{iSess,iBehv2,iShuff})',pcaProjStrShift{iSess,iBehv1,iShuff},...
                        % %     pcaProjStrShift{iSess,iBehv2,iShuff},paDimShift(iSess,iBehv1,1,1),paDimShift(iSess,iBehv2,1,1),false);
                        % % 
                        % % [pcaDiverShiftCtx(iBehv1,iBehv2,iShuff), pcaAlignShiftCtx(iBehv1,iBehv2,iShuff), pcaAngleShiftCtx(iBehv1,iBehv2,iShuff)] = ...
                        % %     calcPCAAlignment(behvFRsShift{iBehv1,iShuff}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iBehv1,iShuff})',...
                        % %     behvFRsShift{iBehv2,iShuff}(goodNeuronsCtx,shiftTimeIndsToUse{iSess,iBehv2,iShuff})',pcaProjCtxShift{iSess,iBehv1,iShuff},...
                        % %     pcaProjCtxShift{iSess,iBehv2,iShuff},paDimShift(iSess,iBehv1,1,2),paDimShift(iSess,iBehv2,1,2),makeExamplePlot);
                        % % 
                        % % [pcaDiverShiftEmg(iBehv1,iBehv2,iShuff), pcaAlignShiftEmg(iBehv1,iBehv2,iShuff), pcaAngleShiftEmg(iBehv1,iBehv2,iShuff)] = ...
                        % %     calcPCAAlignment(behvEMGsShift{iBehv1,iShuff}(:,shiftTimeIndsToUse{iSess,iBehv1,iShuff})',...
                        % %     behvEMGsShift{iBehv2,iShuff}(:,shiftTimeIndsToUse{iSess,iBehv2,iShuff})',pcaProjEmgShift{iSess,iBehv1,iShuff},...
                        % %     pcaProjEmgShift{iSess,iBehv2,iShuff},paDimShift(iSess,iBehv1,1,3),paDimShift(iSess,iBehv2,1,3),false);

                    end % of PCA shuffs loop

                % do behavior label shifts for CCA analysis
                % % % if runCCA
                % % % 
                % % %     if crossBehvAllTimes
                % % %         behv1UseTimes = 1:size(behvFRsShift{iBehv1,iShuff},2);
                % % %         behv2UseTimes = 1:size(behvFRsShift{iBehv2,iShuff},2);
                % % %     else
                % % %         behv1UseTimes = shiftTimeIndsToUse{iSess,iBehv1,iShuff};
                % % %         behv2UseTimes = shiftTimeIndsToUse{iSess,iBehv2,iShuff};
                % % %     end
                % % % 
                % % %     [ccaPrinAngleShiftCtx(iSess,iBehv1,iBehv2,iShuff), ccaDiverShiftCtx(iSess,iBehv1,iBehv2,iShuff), ccaDropShiftCtx(iSess,iBehv1,iBehv2,iShuff)] = calcCCAAlignmentMetrics(...
                % % %         behvFRsShift{iBehv1,iShuff}(goodNeuronsCtx,behv1UseTimes),...
                % % %         behvFRsShift{iBehv2,iShuff}(goodNeuronsCtx,behv2UseTimes),...
                % % %         behvEMGsShift{iBehv1,iShuff}(1:4,behv1UseTimes),...
                % % %         behvEMGsShift{iBehv2,iShuff}(1:4,behv2UseTimes),...
                % % %         ccaNeurProjCtxRegShift{iBehv1,iShuff},ccaNeurProjCtxRegShift{iBehv2,iShuff},...
                % % %         ccaEMGProjCtxRegShift{iBehv1,iShuff},ccaEMGProjCtxRegShift{iBehv2,iShuff},...
                % % %         squeeze(cannonCorrsCtxRegShift(iSess,iBehv1,iShuff,:)),squeeze(cannonCorrsCtxRegShift(iSess,iBehv2,iShuff,:)));
                % % % 
                % % %     [ccaPrinAngleShiftStr(iSess,iBehv1,iBehv2,iShuff), ccaDiverShiftStr(iSess,iBehv1,iBehv2,iShuff), ccaDropShiftStr(iSess,iBehv1,iBehv2,iShuff)] = calcCCAAlignmentMetrics(...
                % % %         behvFRsShift{iBehv1,iShuff}(goodNeuronsStr,behv1UseTimes),...
                % % %         behvFRsShift{iBehv2,iShuff}(goodNeuronsStr,behv2UseTimes),...
                % % %         behvEMGsShift{iBehv1,iShuff}(1:4,behv1UseTimes),...
                % % %         behvEMGsShift{iBehv2,iShuff}(1:4,behv2UseTimes),...
                % % %         ccaNeurProjStrRegShift{iBehv1,iShuff},ccaNeurProjStrRegShift{iBehv2,iShuff},...
                % % %         ccaEMGProjStrRegShift{iBehv1,iShuff},ccaEMGProjStrRegShift{iBehv2,iShuff},...
                % % %         squeeze(cannonCorrsStrRegShift(iSess,iBehv1,iShuff,:)),squeeze(cannonCorrsStrRegShift(iSess,iBehv2,iShuff,:)));
                % % % 
                % % % end

                %                     %permutation shuff
                %                     shuffPerm = randperm(leastPoints*2);
                %                     [pcaProjShuff1Str, pcaTrajShuff1Str, varExpShuff1Str] = pca(allRegDataStr(:,shuffPerm(1:leastPoints))');
                %                     [pcaProjShuff2Str, pcaTrajShuff2Str, varExpShuff2Str] = pca(allRegDataStr(:,shuffPerm(leastPoints+1:end))');
                %
                %                     pcaNeurAligment10DimShuffStr(iBehv1,iBehv2,iShuff) = calcPCAAlignment(allRegDataStr(:,shuffPerm(1:leastPoints))',...
                %                         allRegDataStr(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Str,pcaProjShuff2Str,'auto');
                %
                %                     [pcaProjShuff1Ctx, pcaTrajShuff1Ctx, varExpShuff1Ctx] = pca(allRegDataCtx(:,shuffPerm(1:leastPoints))');
                %                     [pcaProjShuff2Ctx, pcaTrajShuff2Ctx, varExpShuff2Ctx] = pca(allRegDataCtx(:,shuffPerm(leastPoints+1:end))');
                %
                %                     pcaNeurAligment10DimShuffCtx(iBehv1,iBehv2,iShuff) = calcPCAAlignment(allRegDataCtx(:,shuffPerm(1:leastPoints))',...
                %                         allRegDataCtx(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Ctx,pcaProjShuff2Ctx,'auto');
                %
                %                     [pcaProjShuff1Emg, pcaTrajShuff1Emg, varExpShuff1Emg] = pca(allRegDataEmg(:,shuffPerm(1:leastPoints))');
                %                     [pcaProjShuff2Emg, pcaTrajShuff2Emg, varExpShuff2Emg] = pca(allRegDataEmg(:,shuffPerm(leastPoints+1:end))');
                %
                %                     pcaNeurAligment10DimShuffEmg(iBehv1,iBehv2,iShuff) = calcPCAAlignment(allRegDataEmg(:,shuffPerm(1:leastPoints))',...
                %                         allRegDataEmg(:,shuffPerm(leastPoints+1:end))',pcaProjShuff1Emg,pcaProjShuff2Emg,'auto');
                %
                %
                %                     [ccaNeurProjShuff1Ctx,ccaEMGProjShuff1Ctx,cannonCorrsShuff1Ctx,ccaNeurTrajShuff1Ctx,ccaEMGTrajShuff1Ctx] = ...
                %                         canoncorr(allRegDataCtx(:,shuffPerm(1:leastPoints))',allRegDataEmg(:,shuffPerm(1:leastPoints))');
                %                     [ccaNeurProjShuff2Ctx,ccaEMGProjShuff2Ctx,cannonCorrsShuff2Ctx,ccaNeurTrajShuff2Ctx,ccaEMGTrajShuff2Ctx] = ...
                %                         canoncorr(allRegDataCtx(:,shuffPerm(leastPoints+1:end))',allRegDataEmg(:,shuffPerm(leastPoints+1:end))');
                %
                %                     [ccaPrinAngleShuffCtx(iBehv1,iBehv2,iShuff), ccaAlignmentShuffCtx(iBehv1,iBehv2,iShuff)] = calcCCAAlignmentMetrics(...
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
                %                     [ccaPrinAngleShuffStr(iBehv1,iBehv2,iShuff), ccaAlignmentShuffStr(iBehv1,iBehv2,iShuff)] = calcCCAAlignmentMetrics(...
                %                         allRegDataStr(:,shuffPerm(1:leastPoints)),allRegDataStr(:,shuffPerm(leastPoints+1:end)),...
                %                         allRegDataEmg(:,shuffPerm(1:leastPoints)),allRegDataEmg(:,shuffPerm(leastPoints+1:end)),...
                %                         ccaNeurProjShuff1Str,ccaNeurProjShuff2Str,...
                %                         ccaEMGProjShuff1Str,ccaEMGProjShuff2Str,cannonCorrsShuff1Str,cannonCorrsShuff2Str);
                %
                %
                %                     reg1Proj = pcaProjCtx{iBehv1}(:);
                %                     reg1Proj = reg1Proj(randperm(length(reg1Proj),length(reg1Proj)));
                %                     reg1Proj = reshape(reg1Proj,size(pcaProjCtx{iBehv1},1),size(pcaProjCtx{iBehv1},2));
                %
                %                     reg2Proj = pcaProjCtx{iBehv2}(:);
                %                     reg2Proj = reg2Proj(randperm(length(reg2Proj),length(reg2Proj)));
                %                     reg2Proj = reshape(reg2Proj,size(pcaProjCtx{iBehv2},1),size(pcaProjCtx{iBehv2},2));
                %                     reg1Proj = pcaProjCtx{iBehv1}(:,randperm(size(pcaProjCtx{iBehv1},1),size(pcaProjCtx{iBehv1},1)));
                %                     reg2Proj = pcaProjCtx{iBehv2}(:,randperm(size(pcaProjCtx{iBehv1},1),size(pcaProjCtx{iBehv1},1)));
                %
                %                     pcaNeurAligment10DimCtxShuff(iBehv1,iBehv2,iShuff) = calcPCAAlignment(behvFRs{iBehv1}(goodNeuronsCtx,timeIndsToUse{iBehv1})',...
                %                         behvFRs{iBehv2}(goodNeuronsCtx,timeIndsToUse{iBehv2})',reg1Proj,reg2Proj,alignmentDim);
                end % of shuffs loop
            end % of shuffs loop if block

            %                 %do single neuron correlations
            %                 corrAlignmentStr(iBehv1,iBehv2) = corr(squareform(corrMatrixStr{iBehv1})',squareform(corrMatrixStr{iBehv2})');
            %                 corrAlignmentCtx(iBehv1,iBehv2) = corr(squareform(corrMatrixCtx{iBehv1})',squareform(corrMatrixCtx{iBehv2})');
            %                 corrAlignmentEmg(iBehv1,iBehv2) = corr(squareform(corrMatrixEmg{iBehv1})',squareform(corrMatrixEmg{iBehv2})');

        elseif iBehv1 == iBehv2

            if runPCA
                % same behaviors set to 0 divergence/1 alignemnt
                pcaDiverStr(iBehv1,iBehv2) = 0; pcaDiverCtx(iBehv1,iBehv2) = 0; pcaDiverEmg(iBehv1,iBehv2) = 0;
                pcaAlignStr(iBehv1,iBehv2) = 1; pcaAlignCtx(iBehv1,iBehv2) = 1; pcaAlignEmg(iBehv1,iBehv2) = 1;
                pcaAngleStr(iBehv1,iBehv2) = 0; pcaAngleCtx(iBehv1,iBehv2) = 0; pcaAngleEmg(iBehv1,iBehv2) = 0;

                pcaDiverRotStr(iBehv1,iBehv2,1:nShuffs) = 0; pcaDiverRotCtx(iBehv1,iBehv2,1:nShuffs) = 0; pcaDiverRotEmg(iBehv1,iBehv2,1:nShuffs) = 0;
                pcaAlignRotStr(iBehv1,iBehv2,1:nShuffs) = 1; pcaAlignRotCtx(iBehv1,iBehv2,1:nShuffs) = 1; pcaAlignRotEmg(iBehv1,iBehv2,1:nShuffs) = 1;
                pcaAngleRotStr(iBehv1,iBehv2,1:nShuffs) = 0; pcaAngleRotCtx(iBehv1,iBehv2,1:nShuffs) = 0; pcaAngleRotEmg(iBehv1,iBehv2,1:nShuffs) = 0;

                pcaDiverShiftStr(iBehv1,iBehv2,1:nShuffs) = 0; pcaDiverShiftCtx(iBehv1,iBehv2,1:nShuffs) = 0; pcaDiverShiftEmg(iBehv1,iBehv2,1:nShuffs) = 0;
                pcaAlignShiftStr(iBehv1,iBehv2,1:nShuffs) = 1; pcaAlignShiftCtx(iBehv1,iBehv2,1:nShuffs) = 1; pcaAlignShiftEmg(iBehv1,iBehv2,1:nShuffs) = 1;
                pcaAngleShiftStr(iBehv1,iBehv2,1:nShuffs) = 0; pcaAngleShiftCtx(iBehv1,iBehv2,1:nShuffs) = 0; pcaAngleShiftEmg(iBehv1,iBehv2,1:nShuffs) = 0;
            end

            if runCCA
                ccaPrinAngleCtx(iBehv1,iBehv2) = 0; ccaPrinAngleStr(iBehv1,iBehv2) = 0;
                ccaDiverCtx(iBehv1,iBehv2) = 0; ccaDiverStr(iBehv1,iBehv2) = 0;
                ccaDropCtx(iBehv1,iBehv2) = 0; ccaDropStr(iBehv1,iBehv2) = 0;
                ccaDropFracCtx(iBehv1,iBehv2) = 0; ccaDropFracStr(iBehv1,iBehv2) = 0;

                % ccaPrinAngleShiftCtx(iBehv1,iBehv2,1:nShuffs) = 0; ccaPrinAngleShiftStr(iBehv1,iBehv2,1:nShuffs) = 0;
                % ccaDiverShiftCtx(iBehv1,iBehv2,1:nShuffs) = 0; ccaDiverShiftStr(iBehv1,iBehv2,1:nShuffs) = 0;
                % ccaDropShiftCtx(iBehv1,iBehv2,1:nShuffs) = 0; ccaDropShiftStr(iBehv1,iBehv2,1:nShuffs) = 0;
            end

            if runDyn
                dynPopR2ChangeStr(iBehv1,iBehv2) = 0;
                dynPopR2ChangeCtx(iBehv1,iBehv2) = 0;
                dynPopR2FracChangeStr(iBehv1,iBehv2) = 0;
                dynPopR2FracChangeCtx(iBehv1,iBehv2) = 0;
                
                dynNeursCCChangeStr{iBehv1,iBehv2} = zeros(length(goodNeuronsStr),2);
                dynNeursCCChangeCtx{iBehv1,iBehv2} = zeros(length(goodNeuronsCtx),2);

                dynPopCCFracChangeStr(iBehv1,iBehv2) = 0;
                dynPopCCFracChangeCtx(iBehv1,iBehv2) = 0;

                dynNeursR2ChangeStr{iBehv1,iBehv2} = zeros(length(goodNeuronsStr),2);
                dynNeursR2ChangeCtx{iBehv1,iBehv2} = zeros(length(goodNeuronsCtx),2);
            end

        end

        disp(['Calc Alignment: ' num2str(toc)])

    end % of inner behavior loop
end % of outer behavior loop

%fill in lower triangle of the matrix
if runPCA
    pcaDiverStr = triu(squeeze(pcaDiverStr),1)' + triu(squeeze(pcaDiverStr));
    pcaDiverCtx = triu(squeeze(pcaDiverCtx),1)' + triu(squeeze(pcaDiverCtx));
    pcaDiverEmg = triu(squeeze(pcaDiverEmg),1)' + triu(squeeze(pcaDiverEmg));

    pcaAlignStr = triu(squeeze(pcaAlignStr),1)' + triu(squeeze(pcaAlignStr));
    pcaAlignCtx = triu(squeeze(pcaAlignCtx),1)' + triu(squeeze(pcaAlignCtx));
    pcaAlignEmg = triu(squeeze(pcaAlignEmg),1)' + triu(squeeze(pcaAlignEmg));

    pcaAngleStr = triu(squeeze(pcaAngleStr),1)' + triu(squeeze(pcaAngleStr));
    pcaAngleCtx = triu(squeeze(pcaAngleCtx),1)' + triu(squeeze(pcaAngleCtx));
    pcaAngleEmg = triu(squeeze(pcaAngleEmg),1)' + triu(squeeze(pcaAngleEmg));
end

if runCCA
    ccaPrinAngleStr = triu(squeeze(ccaPrinAngleStr),1)' + triu(squeeze(ccaPrinAngleStr));
    ccaPrinAngleCtx = triu(squeeze(ccaPrinAngleCtx),1)' + triu(squeeze(ccaPrinAngleCtx));

    ccaDiverStr = triu(squeeze(ccaDiverStr),1)' + triu(squeeze(ccaDiverStr));
    ccaDiverCtx = triu(squeeze(ccaDiverCtx),1)' + triu(squeeze(ccaDiverCtx));

    ccaDropStr = triu(squeeze(ccaDropStr),1)' + triu(squeeze(ccaDropStr));
    ccaDropCtx = triu(squeeze(ccaDropCtx),1)' + triu(squeeze(ccaDropCtx));

    ccaDropFracStr = triu(squeeze(ccaDropFracStr),1)' + triu(squeeze(ccaDropFracStr));
    ccaDropFracCtx = triu(squeeze(ccaDropFracCtx),1)' + triu(squeeze(ccaDropFracCtx));
end

if runDyn
    dynDropStr = triu(squeeze(dynPopR2ChangeStr),1)' + triu(squeeze(dynPopR2ChangeStr));
    dynDropCtx = triu(squeeze(dynPopR2ChangeCtx),1)' + triu(squeeze(dynPopR2ChangeCtx));

    dynDropFracStr = triu(squeeze(dynPopR2FracChangeStr),1)' + triu(squeeze(dynPopR2FracChangeStr));
    dynDropFracCtx = triu(squeeze(dynPopR2FracChangeCtx),1)' + triu(squeeze(dynPopR2FracChangeCtx));

    dynCCDropFracStr = triu(squeeze(dynPopCCFracChangeStr),1)' + triu(squeeze(dynPopCCFracChangeStr));
    dynCCDropFracCtx = triu(squeeze(dynPopCCFracChangeCtx),1)' + triu(squeeze(dynPopCCFracChangeCtx));

end

for iShuff = 1:nShuffs

    if runPCA
        pcaDiverRotStr(:,:,iShuff) = triu(squeeze(pcaDiverRotStr(:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotStr(:,:,iShuff)));
        pcaDiverRotCtx(:,:,iShuff) = triu(squeeze(pcaDiverRotCtx(:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotCtx(:,:,iShuff)));
        pcaDiverRotEmg(:,:,iShuff) = triu(squeeze(pcaDiverRotEmg(:,:,iShuff)),1)' + triu(squeeze(pcaDiverRotEmg(:,:,iShuff)));

        pcaAlignRotStr(:,:,iShuff) = triu(squeeze(pcaAlignRotStr(:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotStr(:,:,iShuff)));
        pcaAlignRotCtx(:,:,iShuff) = triu(squeeze(pcaAlignRotCtx(:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotCtx(:,:,iShuff)));
        pcaAlignRotEmg(:,:,iShuff) = triu(squeeze(pcaAlignRotEmg(:,:,iShuff)),1)' + triu(squeeze(pcaAlignRotEmg(:,:,iShuff)));

        pcaAngleRotStr(:,:,iShuff) = triu(squeeze(pcaAngleRotStr(:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotStr(:,:,iShuff)));
        pcaAngleRotCtx(:,:,iShuff) = triu(squeeze(pcaAngleRotCtx(:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotCtx(:,:,iShuff)));
        pcaAngleRotEmg(:,:,iShuff) = triu(squeeze(pcaAngleRotEmg(:,:,iShuff)),1)' + triu(squeeze(pcaAngleRotEmg(:,:,iShuff)));

        % pcaDiverShiftStr(:,:,iShuff) = triu(squeeze(pcaDiverShiftStr(:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftStr(:,:,iShuff)));
        % pcaDiverShiftCtx(:,:,iShuff) = triu(squeeze(pcaDiverShiftCtx(:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftCtx(:,:,iShuff)));
        % pcaDiverShiftEmg(:,:,iShuff) = triu(squeeze(pcaDiverShiftEmg(:,:,iShuff)),1)' + triu(squeeze(pcaDiverShiftEmg(:,:,iShuff)));
        % 
        % pcaAlignShiftStr(:,:,iShuff) = triu(squeeze(pcaAlignShiftStr(:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftStr(:,:,iShuff)));
        % pcaAlignShiftCtx(:,:,iShuff) = triu(squeeze(pcaAlignShiftCtx(:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftCtx(:,:,iShuff)));
        % pcaAlignShiftEmg(:,:,iShuff) = triu(squeeze(pcaAlignShiftEmg(:,:,iShuff)),1)' + triu(squeeze(pcaAlignShiftEmg(:,:,iShuff)));
        % 
        % pcaAngleShiftStr(:,:,iShuff) = triu(squeeze(pcaAngleShiftStr(:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftStr(:,:,iShuff)));
        % pcaAngleShiftCtx(:,:,iShuff) = triu(squeeze(pcaAngleShiftCtx(:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftCtx(:,:,iShuff)));
        % pcaAngleShiftEmg(:,:,iShuff) = triu(squeeze(pcaAngleShiftEmg(:,:,iShuff)),1)' + triu(squeeze(pcaAngleShiftEmg(:,:,iShuff)));
    end

    % if runCCA
    %     ccaPrinAngleShiftStr(:,:,iShuff) = triu(squeeze(ccaPrinAngleShiftStr(:,:,iShuff)),1)' + triu(squeeze(ccaPrinAngleShiftStr(:,:,iShuff)));
    %     ccaPrinAngleShiftCtx(:,:,iShuff) = triu(squeeze(ccaPrinAngleShiftCtx(:,:,iShuff)),1)' + triu(squeeze(ccaPrinAngleShiftCtx(:,:,iShuff)));
    % 
    %     ccaDiverShiftStr(:,:,iShuff) = triu(squeeze(ccaDiverShiftStr(:,:,iShuff)),1)' + triu(squeeze(ccaDiverShiftStr(:,:,iShuff)));
    %     ccaDiverShiftCtx(:,:,iShuff) = triu(squeeze(ccaDiverShiftCtx(:,:,iShuff)),1)' + triu(squeeze(ccaDiverShiftCtx(:,:,iShuff)));
    % 
    %     ccaDropShiftStr(:,:,iShuff) = triu(squeeze(ccaDropShiftStr(:,:,iShuff)),1)' + triu(squeeze(ccaDropShiftStr(:,:,iShuff)));
    %     ccaDropShiftCtx(:,:,iShuff) = triu(squeeze(ccaDropShiftCtx(:,:,iShuff)),1)' + triu(squeeze(ccaDropShiftCtx(:,:,iShuff)));
    % end

end

if runCCA
    outputVars.ccaPrinAngleCtx = ccaPrinAngleCtx;
    outputVars.ccaDiverCtx = ccaDiverCtx;
    outputVars.ccaDropCtx = ccaDropCtx;
    outputVars.ccaDropFracCtx = ccaDropFracCtx;
    outputVars.ccaPrinAngleStr = ccaPrinAngleStr;
    outputVars.ccaDiverStr = ccaDiverStr;
    outputVars.ccaDropStr = ccaDropStr;
    outputVars.ccaDropFracStr = ccaDropFracStr;
    % if inputVars.nShifts > 0
    %     outputVars.ccaPrinAngleShiftStr = ccaPrinAngleShiftStr;
    %     outputVars.ccaPrinAngleShiftCtx = ccaPrinAngleShiftCtx;
    %     outputVars.ccaDiverShiftStr = ccaDiverShiftStr;
    %     outputVars.ccaDiverShiftCtx = ccaDiverShiftCtx;
    %     outputVars.ccaDropShiftStr = ccaDropShiftStr;
    %     outputVars.ccaDropShiftCtx = ccaDropShiftCtx;
    % end
end

if runPCA
    outputVars.pcaDiverCtx = pcaDiverCtx;
    outputVars.pcaAlignCtx = pcaAlignCtx;
    outputVars.pcaAngleCtx = pcaAngleCtx;
    outputVars.pcaDiverStr = pcaDiverStr;
    outputVars.pcaAlignStr = pcaAlignStr;
    outputVars.pcaAngleStr = pcaAngleStr;
    outputVars.pcaDiverEmg = pcaDiverEmg;
    outputVars.pcaAlignEmg = pcaAlignEmg;
    outputVars.pcaAngleEmg = pcaAngleEmg;

    outputVars.pcaDiverRotStr = pcaDiverRotStr;
    outputVars.pcaDiverRotCtx = pcaDiverRotCtx;
    outputVars.pcaDiverRotEmg = pcaDiverRotEmg;
    outputVars.pcaAlignRotStr = pcaAlignRotStr;
    outputVars.pcaAlignRotCtx = pcaAlignRotCtx;
    outputVars.pcaAlignRotEmg = pcaAlignRotEmg;
    outputVars.pcaAngleRotStr = pcaAngleRotStr;
    outputVars.pcaAngleRotCtx = pcaAngleRotCtx;
    outputVars.pcaAngleRotEmg = pcaAngleRotEmg;

    % if inputVars.nShifts > 0
    %     outputVars.pcaDiverShiftStr = pcaDiverShiftStr;
    %     outputVars.pcaDiverShiftCtx = pcaDiverShiftCtx;
    %     outputVars.pcaDiverShiftEmg = pcaDiverShiftEmg;
    %     outputVars.pcaAlignShiftStr = pcaAlignShiftStr;
    %     outputVars.pcaAlignShiftCtx = pcaAlignShiftCtx;
    %     outputVars.pcaAlignShiftEmg = pcaAlignShiftEmg;
    %     outputVars.pcaAngleShiftStr = pcaAngleShiftStr;
    %     outputVars.pcaAngleShiftCtx = pcaAngleShiftCtx;
    %     outputVars.pcaAngleShiftEmg = pcaAngleShiftEmg;
    % end
end

if runDyn   
    outputVars.dynDropStr = dynDropStr;
    outputVars.dynDropCtx = dynDropCtx;
    outputVars.dynDropFracStr = dynDropFracStr;
    outputVars.dynDropFracCtx = dynDropFracCtx;
    outputVars.dynCCDropFracStr = dynCCDropFracStr;
    outputVars.dynCCDropFracCtx = dynCCDropFracCtx;

    outputVars.dynNeursR2ChangeStr = dynNeursR2ChangeStr;
    outputVars.dynNeursR2ChangeCtx = dynNeursR2ChangeCtx;
    outputVars.dynNeursCCChangeStr = dynNeursCCChangeStr;
    outputVars.dynNeursCCChangeCtx = dynNeursCCChangeCtx;
end


end % of calcPairwiseSimilarity function




%
