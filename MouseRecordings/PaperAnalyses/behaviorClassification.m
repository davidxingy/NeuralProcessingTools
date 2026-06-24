clear
close all

sessionNames = {'D020-062922-ArenaRecording', 'D024-111022-ArenaRecording', ...
    'D026-032923-ArenaRecording', 'D050-121825-ArenaRecording', 'D054-012126-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7 ...
    ];
behvRegionLabels = {'Climb Up','Climb Down','Jump Down/Across','Walk Flat/Grid','Rearing/Still','Grooming','Eating'};

classifierMethod = 'randomforrest';
nHist = 10;

lags = 0:2:40;
nShifts = 2;
nFolds = 5;
balanceClasses = true;

rng(2024)

for iSess = 1:length(sessionNames)

    behvAlignPerm = allBehvAlignPerms(iSess,:);
    filePaths = getMouseDataNames(sessionNames{iSess}(1:4),sessionNames{iSess},'CFA');
    sessionArtifacts = consolidateArtifactInds(sessionNames{iSess},'CFA');

    % load in spiking data
    binSize = 100;
    load(filePaths.VideoSyncFrames)

    % load in emg data
    load(filePaths.EMG1ms);
    % downsample by averaging
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


    if binSize == 10
        load(filePaths.NeuralFiringRates10msBins30msGauss,'allFRs','striatumInds','cortexInds');
    elseif binSize == 50
        load(filePaths.NeuralFiringRates50msBins50msGauss,'allFRs','striatumInds','cortexInds');
    elseif binSize == 100
        load(filePaths.NeuralFiringRates100msBins0msGauss,'allFRs','striatumInds','cortexInds');
    end

    % load UMAP projection
    load(filePaths.UMAPFile,'origDownsampEMGInd','regionAssignmentsFiltered','regionWatershedLabels','badEMGChans')

    %downsample the behavior labels
    %only keep bins where all of them belong to a label
    nDownSampBins = floor(length(regionAssignmentsFiltered)/binSize);
    regionAssignmentsFilteredReshape = reshape(regionAssignmentsFiltered(1:nDownSampBins*binSize),binSize,nDownSampBins);
    goodBins = find(all(regionAssignmentsFilteredReshape == mode(regionAssignmentsFilteredReshape)));

    % don't use bins with artifacts
    artInds = round(sessionArtifacts.allArtUmapInds/binSize);
    goodBins(ismember(goodBins,artInds)) = [];

    origDownsampEMGIndBinned = round(origDownsampEMGInd(goodBins*binSize-round(binSize/2))/binSize);
    regionAssignmentsFilteredBinned = regionAssignmentsFiltered(goodBins*binSize-round(binSize/2));

    %align neur data to reduction
    currentDir = pwd;
    cd(filePaths.processedDataFolder)
    reducNeurInds = round(NeurEMGSync(origDownsampEMGIndBinned*20*binSize, frameEMGSamples, frameNeuropixelSamples, 'EMG')/(30*binSize));
    cd(currentDir)

    % don't use any indices
    maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/(30*binSize))-1;
    reducIndsToUse = find(reducNeurInds < maxNeurSamples);
    origDownsampEMGIndBinned = origDownsampEMGIndBinned(reducIndsToUse);
    reducNeurInds = reducNeurInds(reducIndsToUse);

    %get the region classification for each point
    behvLabels = regionAssignmentsFilteredBinned(reducIndsToUse);

    % determine which neurons/emg chans to use
    goodNeuronsStr = find(nanmean(allFRs(1:length(striatumInds),:),2)*(1000/binSize)>0.1);
    usedStrInds = 1:length(striatumInds);
    usedStrInds = usedStrInds(goodNeuronsStr);

    goodNeuronsCtx = find(nanmean(allFRs(length(striatumInds)+1:end,:),2)*(1000/binSize)>0.1);
    usedCtxInds = length(striatumInds)+1:size(allFRs,1);
    usedCtxInds = usedCtxInds(goodNeuronsCtx);

    goodEmgs = setdiff(1:size(downsampEMG,1),badEMGChans);

    % add history to neural and emg data
    historyStrFRs = addHistory(allFRs(usedStrInds,:),nHist);
    historyCtxFRs = addHistory(allFRs(usedCtxInds,:),nHist);
    historyEmgs = addHistory(downsampEMG(goodEmgs,:),nHist);

    % do classifier analysis for both unshifted, and random time shifts (all in one for loop)
    for iShift = 2:nShifts+1

        % also if not shifts, do different time lags
        if iShift == 1
            nLags = length(lags);
            isShift = false;
        else
            nLags = 1;
            isShift = true;

            % get shift amount
            minShiftAmount = 100000/binSize;
            shiftAmount(iSess,iShift-1) = randperm(length(behvLabels)-2*minShiftAmount,1)+minShiftAmount;
        end

        for iLag = 1:nLags

            % get the data for this lag/shift
            if isShift
                usedEmgInds = origDownsampEMGIndBinned;
                usedFRInds = reducNeurInds;
                usedLabels = circshift(behvLabels,shiftAmount(iSess,iShift-1));

                usedStrFRs = historyStrFRs(:,usedFRInds);
                usedCtxFRs = historyCtxFRs(:,usedFRInds);
                usedEmgs = historyEmgs(:,usedEmgInds);

            else
                usedEmgInds = origDownsampEMGIndBinned - lags(iLag);
                usedFRInds = reducNeurInds - lags(iLag);

                badInds = unique([find(usedEmgInds < 1) find(usedFRInds < 1)]);
                usedLabels = behvLabels;
                usedEmgInds(badInds) = [];
                usedFRInds(badInds) = [];
                usedLabels(badInds) = [];

                usedStrFRs = historyStrFRs(:,usedFRInds);
                usedCtxFRs = historyCtxFRs(:,usedFRInds);
                usedEmgs = historyEmgs(:,usedEmgInds);
            end

            % don't use any points where there's nans in the data
            nanStrPoints = find(any(isnan(usedStrFRs),1));
            nanCtxPoints = find(any(isnan(usedCtxFRs),1));
            nanEmgPoints = find(any(isnan(usedEmgs),1));
            nanPoints = unique([nanStrPoints nanCtxPoints nanEmgPoints]);

            usedStrFRs(:,nanPoints) = [];
            usedCtxFRs(:,nanPoints) = [];
            usedEmgs(:,nanPoints) = [];
            usedLabels(nanPoints) = [];

            %down sample different behavioral classes to balance class sizes
            if balanceClasses

                behvLabelsIds = unique(usedLabels);
                leastPoints = min(histcounts(usedLabels));
                for iRegion = 1:length(behvLabelsIds)
                    regionInds = find(usedLabels==behvLabelsIds(iRegion));
                    regionIndsDownsamp{iRegion} = regionInds(sort(randperm(length(regionInds),leastPoints)));
                end

                usedLabels = usedLabels([regionIndsDownsamp{:}]);
                usedStrFRs = usedStrFRs(:,[regionIndsDownsamp{:}]);
                usedCtxFRs = usedCtxFRs(:,[regionIndsDownsamp{:}]);
                usedEmgs = usedEmgs(:,[regionIndsDownsamp{:}]);

            end

            %use cross validation
                    % cvBlocks = divideBlocks(1:length(usedLabels),nFolds);
            for iFold = 1:nFolds
                cvBlocks{iFold} = iFold:nFolds:length(usedLabels);
            end
            
            for iBlock = 1:length(cvBlocks)

                tic

                testInds = cvBlocks{iBlock};
                trainInds = [cvBlocks{setdiff(1:length(cvBlocks),iBlock)}];

                % % % % do PCA as feature extraction step if needed
                % % % nHistory = 1;
                % % % %             [strProj, strTraj, strVaf] = pca(addHistory(reducFRs(usedStrInds,:),nHistory)');
                % % % %             [ctxProj, ctxTraj, ctxVaf] = pca(addHistory(reducFRs(usedCtxInds,:),nHistory)');
                % % % %
                % % % %             [strProj, strTraj, strVaf] = pca(reducFRs(usedStrInds,:)');
                % % % %             [ctxProj, ctxTraj, ctxVaf] = pca(reducFRs(usedCtxInds,:)');
                % % % %
                % % % %             strFeatures = strTraj(:,1:find(cumsum(strVaf)/sum(strVaf)>0.9,1))';
                % % % %             ctxFeatures = ctxTraj(:,1:find(cumsum(ctxVaf)/sum(ctxVaf)>0.9,1))';
                % % % 
                % % % % do classification
                % % % strFeatures = addHistory(reducFRs(usedStrInds,:),nHistory);
                % % % ctxFeatures = addHistory(reducFRs(usedCtxInds,:),nHistory);
                % % % emgFeatures = addHistory(reducEMGs(goodEmgs,:),nHistory);

                % train
                switch lower(classifierMethod)
                    case 'lda'
                        strClassifier = fitcdiscr(usedStrFRs(:,trainInds)',usedLabels(trainInds));
                        ctxClassifier = fitcdiscr(usedCtxFRs(:,trainInds)',usedLabels(trainInds));
                        emgClassifier = fitcdiscr(usedEmgs(:,trainInds)',usedLabels(trainInds));
                    case 'knn'
                        strClassifier = fitcknn(usedStrFRs(:,trainInds)',usedLabels(trainInds),'NumNeighbors',20,'Standardize',0);
                        ctxClassifier = fitcknn(usedCtxFRs(:,trainInds)',usedLabels(trainInds),'NumNeighbors',20,'Standardize',0);
                        emgClassifier = fitcknn(usedEmgs(:,trainInds)',usedLabels(trainInds),'NumNeighbors',20,'Standardize',0);
                    case 'svm'
                        strClassifier = fitcecoc(usedStrFRs(:,trainInds)',usedLabels(trainInds));
                        ctxClassifier = fitcecoc(usedCtxFRs(:,trainInds)',usedLabels(trainInds));
                        emgClassifier = fitcecoc(usedEmgs(:,trainInds)',usedLabels(trainInds));
                    case 'randomforrest'
                        strClassifier = TreeBagger(100,usedStrFRs(:,trainInds)',usedLabels(trainInds),...
                            Method="classification", OOBPrediction="on");
                        ctxClassifier = TreeBagger(100,usedCtxFRs(:,trainInds)',usedLabels(trainInds),...
                            Method="classification", OOBPrediction="on");
                        emgClassifier = TreeBagger(100,usedEmgs(:,trainInds)',usedLabels(trainInds),...
                            Method="classification", OOBPrediction="on");
                end

                % test
                blockStrPred = str2double(string(predict(strClassifier,usedStrFRs(:,testInds)')));
                thisStrAcc = sum(blockStrPred == usedLabels(testInds)')/length(testInds);

                blockCtxPred = str2double(string(predict(ctxClassifier,usedCtxFRs(:,testInds)')));
                blockCtxAcc = sum(blockCtxPred == usedLabels(testInds)')/length(testInds);

                blockEmgPred = str2double(string(predict(emgClassifier,usedEmgs(:,testInds)')));
                blockEmgAcc = sum(blockEmgPred == usedLabels(testInds)')/length(testInds);

                blockTestLabels = usedLabels(testInds);

                if isShift
                    strPredLabelsShift{iSess,iShift-1,iBlock} = blockStrPred;
                    strBlockAccuracyShift(iSess,iShift-1,iBlock) = thisStrAcc;

                    ctxPredLabelsShift{iSess,iShift-1,iBlock} = blockCtxPred;
                    ctxBlockAccuracyShift(iSess,iShift-1,iBlock) = blockCtxAcc;

                    emgPredLabelsShift{iSess,iShift-1,iBlock} = blockEmgPred;
                    emgBlockAccuracyShift(iSess,iShift-1,iBlock) = blockEmgAcc;

                    testLabelsShift{iSess,iShift-1,iBlock} = blockTestLabels;
                else
                    strPredLabels{iSess,iLag,iBlock} = blockStrPred;
                    strBlockAccuracy(iSess,iLag,iBlock) = thisStrAcc;

                    ctxPredLabels{iSess,iLag,iBlock} = blockCtxPred;
                    ctxBlockAccuracy(iSess,iLag,iBlock) = blockCtxAcc;

                    emgPredLabels{iSess,iLag,iBlock} = blockEmgPred;
                    emgBlockAccuracy(iSess,iLag,iBlock) = blockEmgAcc;

                    testLabels{iSess,iLag,iBlock} = blockTestLabels;
                end

                disp(toc)
                
            end

            if isShift
                strAccuracyShift(iSess,iShift-1) = sum(cat(1,strPredLabelsShift{iSess,iShift-1,:}) == cat(2,testLabelsShift{iSess,iShift-1,:})') / ...
                    length(cat(2,testLabelsShift{iSess,iShift-1,:}));
                ctxAccuracyShift(iSess,iShift-1) = sum(cat(1,ctxPredLabelsShift{iSess,iShift-1,:}) == cat(2,testLabelsShift{iSess,iShift-1,:})') / ...
                    length(cat(2,testLabelsShift{iSess,iShift-1,:}));
                emgAccuracyShift(iSess,iShift-1) = sum(cat(1,emgPredLabelsShift{iSess,iShift-1,:}) == cat(2,testLabelsShift{iSess,iShift-1,:})') / ...
                    length(cat(2,testLabelsShift{iSess,iShift-1,:}));
            else
                strAccuracy(iSess,iLag) = sum(cat(1,strPredLabels{iSess,iLag,:}) == cat(2,testLabels{iSess,iLag,:})') / length(cat(2,testLabels{iSess,iLag,:}));
                ctxAccuracy(iSess,iLag) = sum(cat(1,ctxPredLabels{iSess,iLag,:}) == cat(2,testLabels{iSess,iLag,:})') / length(cat(2,testLabels{iSess,iLag,:}));
                emgAccuracy(iSess,iLag) = sum(cat(1,emgPredLabels{iSess,iLag,:}) == cat(2,testLabels{iSess,iLag,:})') / length(cat(2,testLabels{iSess,iLag,:}));
            end


        end % of lag loop

    end % of shift loop

end % of session loop

save('UMAPClassificationAnalysis_noSmooth','strAccuracy','ctxAccuracy','emgAccuracy','strAccuracyShift','ctxAccuracyShift',...
    'emgAccuracyShift','strPredLabels','ctxPredLabels','emgPredLabels','strBlockAccuracy','ctxBlockAccuracy','emgBlockAccuracy',...
    'testLabels','strPredLabelsShift','ctxPredLabelsShift','emgPredLabelsShift','strBlockAccuracyShift','ctxBlockAccuracyShift',...
    'emgBlockAccuracyShift','testLabelsShift','shiftAmount','lags')



%
