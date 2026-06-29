clear
close all

allAnimals = {'D020','D024','D026','D043','D047','D050','D054','D056'};

animalSessNames = {{'D020-062922-ArenaRecording'},... %D020
    {'D024-111022-ArenaRecording'},... %D024
    {'D026-032923-ArenaRecording'},... %D026
    {'D043-013125-ArenaRecording','D043-020625-ArenaRecording','D043-020525-ArenaRecording','D043-020425-ArenaRecording','D043-020325-ArenaRecording','D043-021825-ArenaRecording'},... %D043
    {'D047-090825-ArenaRecording','D047-090925-ArenaRecording','D047-091825-ArenaRecording','D047-091625-ArenaRecording'},... %D047
    {'D050-120925-ArenaRecording','D050-121825-ArenaRecording','D050-120625-ArenaRecording','D050-120525-ArenaRecording'},... %D050
    {'D054-011626-ArenaRecording','D054-012126-ArenaRecording','D054-011326-ArenaRecording','D054-012426-ArenaRecording'},... %D054
    {'D056-012726-ArenaRecording','D056-020926-ArenaRecording','D056-012926-ArenaRecording','D056-013126-ArenaRecording'},... %D056
    };

animalBehvLabeledSessions = {[1],... %D020
    [1],... %D024
    [1],... %D026
    [1, 2],... %D043
    [1, 4],... %D047
    [1, 2, 4],... %D050 (should do two models, 1+4 predict 1 3 and 4, and 2 predict 2)
    [2, 3, 4],... %D054 (should do three models, 3 predict 1 and 3, 2 predict 2, and 4 predict 4)
    [1, 2, 4],... %D056 (should do two models, 1+4 predict 1 3 and 4, and 2 predict 2)
    };

classifierTrainSessions = {{[1]};... %D020
    {[1]};... %D024
    {[1]};... %D026
    {[1, 2]};... %D043
    {[1, 4]};... %D047
    {[2]};... %D050
    {[3]};... %D054
    {[2]};... %D056
    };

classifierPredSessions = {{[1]};... %D020
    {[1]};... %D024
    {[1]};... %D026
    {[6]};... %D043
    {[1:4]};... %D047
    {[2]};... %D050
    {[1, 3]};... %D054
    {[2]};... %D056
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
behvRegionFieldNames = {'climbdown','climbup','eating','grooming','jumpdown','jumping','rearing','still','walkflat','walkgrid'};
behvRegionLabels = {'Climb Down','Climb Up','Eat','Groom','Jump Down','Jump Across','Rear','Still','Walk Flat','Walk Grid'};

nCVFolds = 4;

for iAnimal = 1:length(animalSessNames)

    % % % load in data from all the sessions
    % % allFeatures = {};
    % % for iSess = 1:length(animalSessNames{iAnimal})
    % % 
    % %     dataNames = getMouseDataNames(allAnimals{iAnimal},animalSessNames{iAnimal}{iSess},'CFA');
    % %     behvAlignPerm = allBehvAlignPerms(iAnimal,:);
    % % 
    % %     % determine if this is a training session (in which case it
    % %     % should have labeled data)
    % %     if any(iSess == cat(2,classifierTrainSessions{iAnimal}{:}))
    % %         isTrainSess = true;
    % %     else
    % %         isTrainSess = false;
    % %     end
    % % 
    % %     %load in behavior and behavior region labels
    % %     load(dataNames.UMAPFile,'origDownsampEMGInd','freqData','reduction')
    % % 
    % %     if isTrainSess
    % %         load(dataNames.UMAPFile,'regionAssignmentsFiltered','behvLabelsNoArt','analyzedBehaviors',...
    % %             'regionWatershedLabels','regionAssignmentsFiltered','regionBehvAssignments')
    % %     end
    % % 
    % %     %load in EMG
    % %     load(dataNames.EMG1ms)
    % % 
    % %     % get EMG data corresponding to the UMAP time points
    % %     reducEMGs = downsampEMG(:,origDownsampEMGInd);
    % %     clear downsampEMG
    % % 
    % %     allFeatures{iSess} =  [freqData'; reducEMGs];
    % %     clear freqData reducEMGs
    % % 
    % %     if isTrainSess
    % %         % get behavior labels
    % %         % Use data divded by human annotations
    % %         manualLabeledInds{iSess} = find(behvLabelsNoArt ~= 0);
    % %         unlabeledInds{iSess} = find(behvLabelsNoArt == 0);
    % % 
    % %         for iBehv = 1:length(behvRegionFieldNames)
    % %             thisBehvLabel = find(contains(analyzedBehaviors,behvRegionFieldNames{iBehv}));
    % %             behvInds{iSess,iBehv} = find(behvLabelsNoArt(manualLabeledInds{iSess})==thisBehvLabel);
    % %         end
    % % 
    % %     end
    % % 
    % % end %of sessions loop

    %go through and train models
    for iModel = 1:length(classifierTrainSessions{iAnimal})

        modelTrainSessions = classifierTrainSessions{iAnimal}{iModel};
        modelTrainedSessionNames = animalSessNames{iAnimal}(modelTrainSessions);

        % for this model, get the training data
        for iSess = 1:length(modelTrainSessions)

            dataNames = getMouseDataNames(allAnimals{iAnimal},modelTrainedSessionNames{iSess},'CFA');
            behvAlignPerm = allBehvAlignPerms(iAnimal,:);

            %load in behavior and behavior region labels
            load(dataNames.UMAPFile,'freqData','reduction','regionAssignmentsFiltered','behvLabelsNoArt','analyzedBehaviors',...
                'origDownsampEMGInd','regionWatershedLabels','regionAssignmentsFiltered','regionBehvAssignments')

            %load in EMG
            load(dataNames.EMG1ms)

            % get EMG data corresponding to the UMAP time points
            reducEMGs = downsampEMG(:,origDownsampEMGInd);
            clear downsampEMG

            allFeatures =  [freqData'; reducEMGs];
            clear freqData reducEMGs

            % remove artifact points (some sessions have bad EMG based on UMAP)
            badPoints = find(isnan(regionAssignmentsFiltered));

            allFeatures(:,badPoints) = [];
            behvLabelsNoArt(badPoints) = [];

            % get behavior labels
            % Use data divded by human annotations
            manualLabeledInds = find(behvLabelsNoArt ~= 0);

            for iBehv = 1:length(behvRegionFieldNames)
                thisBehvLabel = find(contains(analyzedBehaviors,behvRegionFieldNames{iBehv}));
                behvInds{iBehv} = find(behvLabelsNoArt(manualLabeledInds)==thisBehvLabel);
            end

            sessTrainFeatures{iSess} = allFeatures(:,manualLabeledInds);

            sessTrainLabels{iSess} = zeros(length(behvRegionFieldNames),length(manualLabeledInds));
            for iBehv = 1:length(behvRegionFieldNames)
                sessTrainLabels{iSess}(iBehv,behvInds{iBehv}) = 1;
            end

        end % of train sessions loop

        % combine all training sessions data
        modelTrainFeatures = cat(2,sessTrainFeatures{:});
        modelTrainLabels = cat(2,sessTrainLabels{:});

        clear sessTrainFeatures sessTrainLabels allFeatures

        % for cross-fold validation, divide data into blocks and train
        % for each one
        cvPredLabels = {}; cvPredProbs = {}; cvAccuracy = {};
        for iFold = 1:nCVFolds
            cvFoldInds{iFold} = iFold:nCVFolds:size(modelTrainFeatures,2);
        end
        for iFold = 1:nCVFolds

            testInds = cvFoldInds{iFold};
            trainInds = [cvFoldInds{setdiff(1:nCVFolds,iFold)}];

            % get a model for each behavior individually
            for iBehv = 1:length(behvRegionFieldNames)

                tic

                switch lower(classifierMethod)
                    case 'lda'
                        classifierModelCV = fitcdiscr(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'knn'
                        classifierModelCV = fitcknn(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)',...
                            'NumNeighbors',20,'Standardize',0);
                    case 'svm'
                        classifierModelCV = fitcecoc(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'naivebayes'
                        classifierModelCV = fitcnb(modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)');
                    case 'randomforrest'
                        classifierModelCV = TreeBagger(2000,modelTrainFeatures(:,trainInds)',modelTrainLabels(iBehv,trainInds)',...
                            Method="classification", OOBPrediction="off");
                end

                % predict on test data
                [outputs,confScores] = predict(classifierModelCV,modelTrainFeatures(:,testInds)');
                cvPredLabels{iFold}(iBehv,:) = str2double(string(outputs));
                cvPredProbs{iFold}(iBehv,:,:) = confScores;
                cvAccuracy{iFold}(iBehv,:) = sum(cvPredLabels{iFold}(iBehv,:) == modelTrainLabels(iBehv,testInds))/length(testInds);

                disp(['Animal ' allAnimals{iAnimal} ', Model ' num2str(iModel) ', ' behvRegionFieldNames{iBehv}, ...
                    ', Fold ' num2str(iFold) ' time: ' num2str(toc)]);

            end %of behaviors loop

        end %of cv folds loop

        clear classifierModelCV

        % Now just get one model trained on all the training data
        classifierModel = {};
        for iBehv = 1:length(behvRegionFieldNames)

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

        clear modelTrainFeatures modelTrainLabels

        % because the model variable is so huge, should save it first, then
        % load it for each behavior to predict later
        dataNames = getMouseDataNames(allAnimals{iAnimal},modelTrainedSessionNames{1},'CFA');
        modelSavedFile = fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers.mat');
        save(modelSavedFile,'cvPredLabels','cvPredProbs','cvAccuracy','classifierModel',...
            'classifierMethod','behvRegionLabels','modelTrainedSessionNames','-v7.3')

        clear classifierModel

        % use them to predict behavior classes for the desired sessions
        modelPredSessions = classifierPredSessions{iAnimal}{iModel};
        modelPredSessionNames = animalSessNames{iAnimal}(modelPredSessions);

        for iSess = 1:length(modelPredSessions)

            % load in the data 
            dataNames = getMouseDataNames(allAnimals{iAnimal},modelPredSessionNames{iSess},'CFA');

            %load in EMG and frequency transform
            load(dataNames.EMG1ms)
            load(dataNames.UMAPFile,'origDownsampEMGInd','freqData','reduction','regionAssignmentsFiltered')

            reducEMGs = downsampEMG(:,origDownsampEMGInd);
            clear downsampEMG

            allFeatures =  [freqData'; reducEMGs];
            clear freqData reducEMGs

            predLabels = [];
            predProbs = [];
            for iBehv = 1:length(behvRegionFieldNames)

                % now load in the model and just use it for the current
                % behavior to save memory
                load(modelSavedFile,'classifierModel')
                thisModel = classifierModel{iBehv};
                clear classifierModel

                [outputs,confScores] = predict(thisModel,allFeatures');
                predLabels(iBehv,:) = str2double(string(outputs));
                predProbs(iBehv,:,:) = confScores;
            end
            clear allFeatures

            % remove UMAP artifact points
            badPoints = find(isnan(regionAssignmentsFiltered));
            predLabels(:,badPoints) = nan;
            predProbs(:,badPoints,:) = nan;

            clear thisModel

            % save predictions
            if exist(fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers.mat'), 'file')
                save(fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers'),'modelSavedFile',...
                    'predLabels','predProbs','classifierMethod','behvRegionLabels','modelTrainedSessionNames','-append')
            else
                save(fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers.mat'),'modelSavedFile',...
                    'predLabels','predProbs','classifierMethod','behvRegionLabels','modelTrainedSessionNames','-v7.3')
            end

        end %of pred sessions loop


    end %of model loop

end %of animal loop


%
