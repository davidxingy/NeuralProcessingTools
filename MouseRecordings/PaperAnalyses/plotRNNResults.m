clear
close all


sessionNames = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

baseRnnFolder = 'X:\David\AnalysesData\pythonRNN\pythonRNN';

for iSess = 2:2%length(sessionNames)
   
    resultsDir = fullfile(baseRnnFolder,'BehaviorAveraged_CV_StriatumShiftInputs');

    behvNames = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};
    muscleNames = {'Right Biceps','Right Triceps','Right ECR','Right PL'};


    usedBehvs = {'Walk - Right Paw Strike','Climb - Right Paw Off','Eat - Right Paw Move','Groom - Right Paw Lift','Rear - Right Paw Lift','Jump Down - Right Paw Strike'};
    usedBehvsLabels = {'Walk','Climb Up','Eat','Groom','Rear','Jump Down'};


    % load in data
    NPYFiles = string(ls(resultsDir));
    NPYFiles = NPYFiles(3:end);

    for iFile = 1:length(NPYFiles)
        [~, varName,varExt] = fileparts(strip(NPYFiles{iFile}));
        if strcmpi(varExt,'.npy')
            resultsStruct.(varName) = readNPY(fullfile(resultsDir,NPYFiles{iFile}));
        end
    end

    load(fullfile(resultsDir,'..','data','limbEventData5ms.mat'),'behvLabels','strData','ctxData','emgData','testTrials','trainTrials')
    load(fullfile(resultsDir,'savedTrainTrials.mat'))

    behvLabels = behvLabels(testTrials{1});
    % make averages if single trial
    if size(resultsStruct.modelTestInputs,1) > 6
        for iBehv = 1:7
            modelInputsBehvAve(iBehv,:,:) = mean(resultsStruct.modelInputs(behvLabels==iBehv,:,:),1);
            modelOutputsBehvAve(iBehv,:,:) = mean(resultsStruct.modelOutputs(behvLabels==iBehv,:,:),1);
            targetOutputsBehvAve(iBehv,:,:) = mean(resultsStruct.targetOutputs(behvLabels==iBehv,:,:),1);
        end
    end



end



% 
