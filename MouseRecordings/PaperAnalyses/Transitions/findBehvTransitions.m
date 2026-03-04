clear
close all

outputDir = 'X:\David\AnalysesData\DecodingInputData';

allAnimals = {'D020','D024','D026','D043','D047','D050','D054','D056'};

animalSessNames = {{'D020-062922-ArenaRecording'},... %D020
    {'D024-111022-ArenaRecording'},... %D024
    {'D026-032923-ArenaRecording'},... %D026
    {'D043-013125-ArenaRecording','D043-020625-ArenaRecording','D043-020525-ArenaRecording','D043-020425-ArenaRecording','D043-020325-ArenaRecording'},... %D043
    {'D047-090825-ArenaRecording','D047-090925-ArenaRecording','D047-091825-ArenaRecording','D047-091625-ArenaRecording'},... %D047
    {},... %D050
    {},... %D054
    {},... %D056
    };

animalBehvLabeledSessions = {[1],... %D020
    [1],... %D024
    [1],... %D026
    [1, 2],... %D043
    [1, 4],... %D047
    [],... %D050
    [],... %D054
    [],... %D056
    };

classifierTrainSessions = {{[1]};... %D020
    {[1]};... %D024
    {[1]};... %D026
    {[1, 4]};... %D043
    {[1, 4]};... %D047
    {[]};... %D050
    {[]};... %D054
    {[]};... %D056
    };

classifierPredSessions = {{[1]};... %D020
    {[1]};... %D024
    {[1]};... %D026
    {[1:5]};... %D043
    {[1:4]};... %D047
    {[]};... %D050
    {[]};... %D054
    {[]};... %D056
    };

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ... %D020
    1 2 3 4 5 6 7; ... %D024
    1 2 3 4 5 6 7; ... %D026
    1 2 3 4 5 6 7; ... %D043
    1 2 3 4 5 6 7; ... %D047
    1 2 3 4 5 6 7; ... %D050
    1 2 3 4 5 6 7; ... %D054
    1 2 3 4 5 6 7; ... %D056
    ];

classifierMethod = 'lda';
behvRegionFieldNames = {'climbup','climbdown','jumpdown','jumping','walkflat','walkgrid','rearing','still','grooming','eating'};
behvRegionLabels = {'Climb Up','Climb Down','Jump Down','Jump Across','Walk Flat','Walk Grid','Rear','Still','Groom','Eat'};

nCVFolds = 4;

for iAnimal = 4%1:length(allDirs)

    % load in data from all the sessions
    allFeatures = {};
    for iSess = 1:length(animalSessNames{iAnimal})

        dataNames = getMouseDataNames(allAnimals{iAnimal},animalSessNames{iAnimal}{iSess},'CFA');
        behvAlignPerm = allBehvAlignPerms(iAnimal,:);

        % determine if this is a training session (in which case it
        % should have labeled data)
        if any(iSess == cat(2,classifierTrainSessions{iAnimal}{:}))
            isTrainSess = true;
        else
            isTrainSess = false;
        end

        %load in behavior and behavior region labels
        load(dataNames.UMAPFile,'origDownsampEMGInd','freqData','reduction')

        if isTrainSess
            load(dataNames.UMAPFile,'regionAssignmentsFiltered','behvLabelsNoArt','analyzedBehaviors',...
                'regionWatershedLabels','regionAssignmentsFiltered','regionBehvAssignments')
        end

        %load in EMG
        load(dataNames.EMG1ms)

        % get EMG data corresponding to the UMAP time points
        reducEMGs = downsampEMG(:,origDownsampEMGInd);
        clear downsampEMG

        allFeatures{iSess} =  [freqData'; reducEMGs];
        clear freqData reducEMGs

        if isTrainSess
            % get behavior labels
            % Use data divded by human annotations
            manualLabeledInds{iSess} = find(behvLabelsNoArt ~= 0);
            unlabeledInds{iSess} = find(behvLabelsNoArt == 0);

            for iBehv = 1:length(behvRegionFieldNames)
                thisBehvLabel = find(contains(analyzedBehaviors,behvRegionFieldNames{iBehv}));
                behvInds{iSess,iBehv} = find(behvLabelsNoArt(manualLabeledInds{iSess})==thisBehvLabel);
            end

        end

    end %of sessions loop

    %go through and train models
    for iModel = 1:length(classifierTrainSessions{iAnimal})

        modelTrainSessions = classifierTrainSessions{iAnimal}{iModel};
        modelTrainedSessionNames = animalSessNames{iAnimal}(modelTrainSessions);

        % for this model, get the training data
        for iSess = 1:length(modelTrainSessions)

            sessTrainFeatures{iSess} = allFeatures{modelTrainSessions(iSess)}(:,manualLabeledInds{modelTrainSessions(iSess)});

            sessTrainLabels{iSess} = zeros(length(behvRegionFieldNames),length(manualLabeledInds{iSess}));
            for iBehv = 1:length(behvRegionFieldNames)
                sessTrainLabels{iSess}(iBehv,behvInds{modelTrainSessions(iSess),iBehv}) = 1;
            end
        end

        modelTrainFeatures = cat(2,sessTrainFeatures{:});
        modelTrainLabels = cat(2,sessTrainLabels{:});

        % for cross-fold validation, divide data into blocks and train
        % for each one
        for iFold = 1:nCVFolds

            cvFoldInds{iFold} = iFold:nCVFolds:size(modelTrainFeatures,2);

            testInds = cvFoldInds{iFold};
            trainInds = [cvFoldInds{setdiff(1:nCVFolds,iFold)}];

            % get a model for each behavior individually
            for ibehv = 1:length(behvRegionFieldNames)

                tic

                switch lower(classifierMethod)
                    case 'lda'
                        classifierModelCV{iFold,iBehv} = fitcdiscr(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'knn'
                        classifierModelCV{iFold,iBehv} = fitcknn(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)',...
                            'NumNeighbors',20,'Standardize',0);
                    case 'svm'
                        classifierModelCV{iFold,iBehv} = fitcecoc(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'naivebayes'
                        classifierModelCV{iFold,iBehv} = fitcnb(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'randomforrest'
                        classifierModelCV{iFold,iBehv} = TreeBagger(2000,modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)',...
                            Method="classification", OOBPrediction="off");
                end

                % predict on test data
                [outputs,confScores] = predict(classifierModelCV{iFold,iBehv},modelTrainFeatures(:,testInds)');
                cvPredLabels{iFold}(iBehv,:) = str2double(string(outputs));
                cvPredProbs{iFold}(iBehv,:) = confScores;
                cvAccuracy{iFold}(iBehv,:) = sum(cvPredLabels{iFold}(iBehv,:) == modelTrainLabels(iBehv,testInds)')/length(testInds);

                disp(['Animal ' allAnimals{iAnimal} ', Model ' num2str(iModel) ', ' behvRegionFieldNames{iBehv}, ...
                    ', Fold ' num2str(iFold) ' time: ' num2str(toc)]);

            end %of behaviors loop

        end %of cv folds loop

        % Now just get one model trained on all the training data
        classifierModel = {};
        for ibehv = 1:length(behvRegionFieldNames)

            switch lower(classifierMethod)
                case 'lda'
                    classifierModel{iBehv} = fitcdiscr(modelTrainFeatures',modelTrainLabels(iBehv,:)');
                case 'knn'
                    classifierModel{iBehv} = fitcknn(modelTrainFeatures',modelTrainLabels(iBehv,:)',...
                        'NumNeighbors',20,'Standardize',0);
                case 'svm'
                    classifierModel{iBehv} = fitcecoc(modelTrainFeatures',modelTrainLabels(iBehv,:)');
                case 'naivebayes'
                    classifierModel{iBehv} = fitcnb(modelTrainFeatures',modelTrainLabels(iBehv,:)');
                case 'randomforrest'
                    classifierModel{iBehv} = TreeBagger(2000,modelTrainFeatures',modelTrainLabels(iBehv,:)',...
                        Method="classification", OOBPrediction="off");
            end

        end %of behaviors loop

        % use them to predict behavior classes for the desired sessions
        modelPredSessions = classifierPredSessions{iAnimal}{iModel};
        for iSess = 1:length(modelPredSessions)

            predLabels = [];
            predProbs = [];
            for iBehv = 1:length(behvRegionFieldNames)
                [outputs,confScores] = predict(classifierModel{iBehv},allFeatures{modelPredSessions(iSess)}');
                predLabels(iBehv,:) = str2double(string(outputs));
                predProbs(iBehv,:) = confScores;
            end

            % now save back to the session's ProcessedData folder
            dataNames = getMouseDataNames(allAnimals{iAnimal},animalSessNames{iAnimal}{modelPredSessions(iSess)},'CFA');
            save(fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers'),'classifierModelCV','cvPredLabels','cvPredProbs',...
                'cvAccuracy','classifierModel','predLabels','predProbs','classifierMethod','behvRegionLabels','modelTrainedSessionNames')

        end %of pred sessions loop


    end %of model loop

end %of animal loop


%
