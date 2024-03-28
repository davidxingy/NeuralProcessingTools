clear
close all

sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    7 6 3 5 4 2 1; ...
    1 2 4 5 3 6 7 ...
    ];

behvRegionLabels = {'Climb Up','Climb Down','Jump Down/Across','Walk Flat/Grid','Rearing/Still','Grooming','Eating'};

classifierMethod = 'randomforrest';

for iSess = 1:length(sessionDirs)

    behvAlignPerm = allBehvAlignPerms(iSess,:);
    baseDir = sessionDirs{iSess};

    %load in spiking data
    % use 10ms bins
    downsampleFactor = 100;
    load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
    
    if downsampleFactor == 10
        load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs','striatumInds','cortexInds');
    elseif downsampleFactor == 50
        load(fullfile(baseDir,'ProcessedData','NeuralFiringRates50msBins50msGauss.mat'),'allFRs','striatumInds','cortexInds');
    elseif downsampleFactor == 100
        load(fullfile(baseDir,'ProcessedData','NeuralFiringRates100msBins50msGauss.mat'),'allFRs','striatumInds','cortexInds');
    end
    
    % load UMAP projection
    load(fullfile(baseDir,'ProcessedData','UMAP.mat'),'reduction','origDownsampEMGInd','gridXInds','gridYInds','regionAssignmentsFiltered','regionWatershedLabels')

    %downsample the behavior labels to 50ms
    %only keep 50ms bins where all 50 ms belong to a label
    nDownSampBins = floor(length(regionAssignmentsFiltered)/downsampleFactor);
    regionAssignmentsFilteredReshape = reshape(regionAssignmentsFiltered(1:nDownSampBins*downsampleFactor),downsampleFactor,nDownSampBins);
    goodBins = find(all(regionAssignmentsFilteredReshape == mode(regionAssignmentsFilteredReshape)));
    
    origDownsampEMGIndBinned = origDownsampEMGInd(goodBins*downsampleFactor-round(downsampleFactor/2));
    regionAssignmentsFilteredBinned = regionAssignmentsFiltered(goodBins*downsampleFactor-round(downsampleFactor/2));

    %align neur data to reduction
    currentDir = pwd;
    cd(fullfile(baseDir,'ProcessedData'))
    reducNeurInds = round(NeurEMGSync(origDownsampEMGIndBinned*20, frameEMGSamples, frameNeuropixelSamples, 'EMG')/(30*downsampleFactor));
    maxNeurSamples = round(frameNeuropixelSamples{1}{end}(end)/(30*downsampleFactor))-1;
    reducIndsToUse = find(reducNeurInds < maxNeurSamples);
    cd(currentDir)

    % get neural data corresponding to the UMAP time points
    reducFRs = allFRs(:,reducNeurInds(reducIndsToUse));
    clear allFRs

    % don't use any points where there's nans in the FR data
    nanReducPoints = any(isnan(reducFRs),1);
    reducIndsToUse(nanReducPoints) = [];

    reducFRs = reducFRs(:,~nanReducPoints);

    %get the region classification for each point
    behvLabels = regionAssignmentsFilteredBinned(reducIndsToUse);

    % do control by shifting behavior labels
    nShifts = 200;

    for iShift = 1:nShifts+1
        
        tic

        % first loop, just do actual classification without any shifting
        if iShift ~= 1
            minShiftAmount = 1000;
            shiftAmount(iSess,iShift-1) = randperm(length(behvLabels)-2*minShiftAmount,1)+minShiftAmount;
            behvLabels = circshift(behvLabels,shiftAmount(iSess,iShift-1));
        end

        %use 10-fold cross validation
        cvBlocks = divideBlocks(1:length(behvLabels),10);

%             for i = 1:10
%                 cvBlocks{i} = i:10:length(behvLabels);
%             end

        goodNeuronsStr = find(mean(reducFRs(1:length(striatumInds),:),2)*(1000/downsampleFactor)>0.2);
        usedStrInds = 1:length(striatumInds);
        usedStrInds = usedStrInds(goodNeuronsStr);

        goodNeuronsCtx = find(mean(reducFRs(length(striatumInds)+1:end,:),2)*(1000/downsampleFactor)>0.2);
        usedCtxInds = length(striatumInds)+1:size(reducFRs,1);
        usedCtxInds = usedCtxInds(goodNeuronsCtx);

        for iBlock = 1:length(cvBlocks)

            testInds = cvBlocks{iBlock};
            trainInds = [cvBlocks{setdiff(1:length(cvBlocks),iBlock)}];

            %         %do random forrest classification
            %         tic
            %         ctxRFClassifier = TreeBagger(10,reducFRs(length(striatumInds)+1:end,trainInds)',behvLabels(trainInds),...
            %             Method="classification",...
            %             OOBPrediction="on");
            %         trainTime(iSess,iBlock) = toc;
            %
            %         strRFPredLabels{iBlock} = str2double(string(predict(strRFClassifier,reducFRs(1:length(striatumInds),testInds)')));
            %         ctxRFPredLabels{iBlock} = str2double(string(predict(ctxRFClassifier,reducFRs(length(striatumInds)+1:end,testInds)')));
            %
            %         strRFBlockAccuracy(iSess,iBlock) = sum(strRFPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);
            %         ctxRFBlockAccuracy(iSess,iBlock) = sum(ctxRFPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);


            % do PCA as feature extraction step if needed
            nHistory = 2;
%             [strProj, strTraj, strVaf] = pca(addHistory(reducFRs(usedStrInds,:),nHistory)');
%             [ctxProj, ctxTraj, ctxVaf] = pca(addHistory(reducFRs(usedCtxInds,:),nHistory)');
% 
%             [strProj, strTraj, strVaf] = pca(reducFRs(usedStrInds,:)');
%             [ctxProj, ctxTraj, ctxVaf] = pca(reducFRs(usedCtxInds,:)');
% 
%             strFeatures = strTraj(:,1:find(cumsum(strVaf)/sum(strVaf)>0.9,1))';
%             ctxFeatures = ctxTraj(:,1:find(cumsum(ctxVaf)/sum(ctxVaf)>0.9,1))';

            % do classification
            strFeatures = addHistory(reducFRs(usedStrInds,:),nHistory);
            ctxFeatures = addHistory(reducFRs(usedCtxInds,:),nHistory);

            % train
            switch lower(classifierMethod)
                case 'lda'
                    strRFClassifier = fitcdiscr(strFeatures(:,trainInds)',behvLabels(trainInds));
                    ctxRFClassifier = fitcdiscr(ctxFeatures(:,trainInds)',behvLabels(trainInds));
                case 'knn'
                    strRFClassifier = fitcknn(strFeatures(:,trainInds)',behvLabels(trainInds),'NumNeighbors',20,'Standardize',0);
                    ctxRFClassifier = fitcknn(ctxFeatures(:,trainInds)',behvLabels(trainInds),'NumNeighbors',20,'Standardize',0);
                case 'svm'
                    strNBClassifier = fitcecoc(strFeatures(:,trainInds)',behvLabels(trainInds));
                    ctxNBClassifier = fitcecoc(ctxFeatures(:,trainInds)',behvLabels(trainInds));
                case 'randomforrest'
                    strRFClassifier = TreeBagger(100,strFeatures(:,trainInds)',behvLabels(trainInds),...
                        Method="classification", OOBPrediction="off");
                    ctxRFClassifier = TreeBagger(100,ctxFeatures(:,trainInds)',behvLabels(trainInds),...
                        Method="classification", OOBPrediction="off");
            end

            % test
            if iShift ==1

                strRFPredLabels{iBlock} = str2double(string(predict(strRFClassifier,strFeatures(:,testInds)')));
                strRFBlockAccuracy(iSess,iBlock) = sum(strRFPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);

                ctxRFPredLabels{iBlock} = str2double(string(predict(ctxRFClassifier,ctxFeatures(:,testInds)')));
                ctxRFBlockAccuracy(iSess,iBlock) = sum(ctxRFPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);

                testLabels{iBlock} = behvLabels(testInds);

            else

                strRFPredLabelsShift{iBlock,iShift-1} = str2double(string(predict(strRFClassifier,strFeatures(:,testInds)')));
                strRFBlockAccuracyShift(iSess,iBlock,iShift-1) = sum(strRFPredLabelsShift{iBlock,iShift-1} == behvLabels(testInds)')/length(testInds);

                ctxRFPredLabelsShift{iBlock,iShift-1} = str2double(string(predict(ctxRFClassifier,ctxFeatures(:,testInds)')));
                ctxRFBlockAccuracyShift(iSess,iBlock,iShift-1) = sum(ctxRFPredLabelsShift{iBlock,iShift-1} == behvLabels(testInds)')/length(testInds);

                testLabelsShift{iBlock,iShift-1} = behvLabels(testInds);

            end

            %         %Also use SVM
            %         %don't use low FR neurons
            %         goodNeuronsStr = find(mean(reducFRs(1:length(striatumInds),:),2)*100>0.5);
            %         goodNeuronsCtx = find(mean(reducFRs(length(striatumInds)+1:end,:),2)*100>0.5);
            %
            %         usedStrInds = 1:length(striatumInds);
            %         usedStrInds = usedStrInds(goodNeuronsStr);
            %         usedCtxInds = length(striatumInds)+1:size(reducFRs,1);
            %         usedCtxInds = usedCtxInds(goodNeuronsCtx);
            %
            %         strNBClassifier = fitcecoc(reducFRs(usedStrInds,trainInds)',behvLabels(trainInds));
            %         ctxNBClassifier = fitcecoc(reducFRs(usedCtxInds,trainInds)',behvLabels(trainInds));
            %
            %         strNBPredLabels{iBlock} = str2double(string(predict(strNBClassifier,reducFRs(usedStrInds,testInds)')));
            %         ctxNBPredLabels{iBlock} = str2double(string(predict(ctxNBClassifier,reducFRs(usedCtxInds,testInds)')));
            %
            %         strNBBlockAccuracy(iBlock) = sum(strNBPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);
            %         ctxNBBlockAccuracy(iBlock) = sum(ctxNBPredLabels{iBlock} == behvLabels(testInds)')/length(testInds);
            %

        end

        if iShift ==1
            strAccuracy(iSess) = sum(cat(1,strRFPredLabels{:}) == cat(2,testLabels{:})') / length(cat(2,testLabels{:}));
            ctxAccuracy(iSess) = sum(cat(1,ctxRFPredLabels{:}) == cat(2,testLabels{:})') / length(cat(2,testLabels{:}));
        else
            strAccuracyShift(iSess,iShift-1) = sum(cat(1,strRFPredLabelsShift{:,iShift-1}) == cat(2,testLabelsShift{:,iShift-1})') / length(cat(2,testLabelsShift{:,iShift-1}));
            ctxAccuracyShift(iSess,iShift-1) = sum(cat(1,ctxRFPredLabelsShift{:,iShift-1}) == cat(2,testLabelsShift{:,iShift-1})') / length(cat(2,testLabelsShift{:,iShift-1}));
        end

        disp(toc)

    end

end



% 
