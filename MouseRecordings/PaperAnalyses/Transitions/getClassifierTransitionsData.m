clear
close all

allAnimals = {'D020','D024','D026','D043','D047','D050','D054','D056'};

animalSessNames = {{'D020-062922-ArenaRecording'},... %D020
    {'D024-111022-ArenaRecording'},... %D024
    {'D026-032923-ArenaRecording'},... %D026
    {'D043-013125-ArenaRecording','D043-020625-ArenaRecording','D043-020525-ArenaRecording','D043-020425-ArenaRecording','D043-020325-ArenaRecording'},... %D043
    {'D047-090825-ArenaRecording','D047-090925-ArenaRecording','D047-091825-ArenaRecording','D047-091625-ArenaRecording'},... %D047
    {'D050-120925-ArenaRecording','D050-121825-ArenaRecording','D050-120625-ArenaRecording','D050-120525-ArenaRecording'},... %D050
    {'D054-011626-ArenaRecording','D054-012126-ArenaRecording','D054-011326-ArenaRecording','D054-012426-ArenaRecording'},... %D054
    {'D056-012726-ArenaRecording','D056-020926-ArenaRecording','D056-012926-ArenaRecording','D056-013126-ArenaRecording'},... %D056
    };

% The behaviors in the classifier output aren't in the same order as the
% ones in the UMAP file, so realign them
% behvAlignPerm = [2 1 10 9 3 4 7 8 5 6]; Changed this in classifier code so the behaviors align
behvAlignPerm = 1:10;

% parameters to use
scoreCutoff = 0.9;
minBoutLength = 50;
initMinBaseline = 200;
initMinLength = 200;
termMinBaseline = 200;
termMinLength = 200;

preTransPoints = 200;
postTransPoints = 200;

nShifts = 100;

for iAnimal = 1:length(allAnimals)

    for iSess = 1:length(animalSessNames{iAnimal})

        % get session file paths
        thisSessName = animalSessNames{iAnimal}{iSess};
        sessFilepaths = getMouseDataNames(allAnimals{iAnimal},thisSessName,'CFA');

        % load in session data
        load(sessFilepaths.neuronDataStruct)
        load(sessFilepaths.NeuralFiringRates10msBins30msGauss,'allFRs','striatumInds','cortexInds')
        load(sessFilepaths.UMAPFile,'origDownsampEMGInd')
        load(sessFilepaths.VideoSyncFrames)
        load(sessFilepaths.EMGSingleBehvClassifiers,'predLabels','predProbs','classifierMethod','behvRegionLabels')
        load(sessFilepaths.UMAPOverlayNeurons)

        nBehvs = length(behvRegionLabels);
        assert(size(origDownsampEMGInd,2)==size(predProbs,2),'Classifier Outputs not same number of data points as UMAP points!')
        classifierLabels= zeros(1,size(predProbs,1));
        classifierBehvs = behvRegionLabels(behvAlignPerm);

        predLabels = predLabels(behvAlignPerm,:);
        predProbs = predProbs(behvAlignPerm,:,:);

        %get transitions
        initFRs = {};
        termFRs = {};
        initFRsShift = {};
        termFRsShift = {};
        for iBehv = 1:nBehvs

            behvScores = predProbs(iBehv,:,2);
            nonBehvScores = predProbs(setdiff(1:nBehvs,iBehv),:,2);

            behvPoints = behvScores > scoreCutoff & ~any(nonBehvScores > scoreCutoff,1);

            % remove instances that are too short
            behvPointsStarts = find(behvPoints(2:end) & ~behvPoints(1:end-1))+1;
            behvPointsEnds = find(~behvPoints(2:end) & behvPoints(1:end-1));

            if behvPoints(1) == 1
                behvPointsStarts = [1, behvPointsStarts];
            end
            if behvPoints(end) == 1
                behvPointsEnds = [behvPointsEnds, length(behvPoints)];
            end

            behvBoutLengths = behvPointsEnds - behvPointsStarts;
            shortLengths = find(behvBoutLengths < minBoutLength);

            for iBout = 1:length(shortLengths)
                behvPoints(behvPointsStarts(shortLengths(iBout)):behvPointsEnds(shortLengths(iBout))) = 0;
            end

            classifierLabels(behvPoints) = iBehv;

            % find initiations and terminations
            behvPointsStarts = find(behvPoints(2:end) & ~behvPoints(1:end-1))+1;
            behvPointsEnds = find(~behvPoints(2:end) & behvPoints(1:end-1));
            if behvPoints(1) == 1
                behvPointsStarts = [1, behvPointsStarts];
            end
            if behvPoints(end) == 1
                behvPointsEnds = [behvPointsEnds, length(behvPoints)];
            end

            % only count an initiation if it is at least a certain amount
            % of time away from the previous termination and if the bout is
            % of a minimun length  (to avoid noisy
            % fluctuations from the classifier output)
            goodInits = behvPointsStarts( behvPointsStarts - [1, behvPointsEnds(1:end-1)] > initMinBaseline & ...
                behvPointsEnds - behvPointsStarts > initMinLength );
            goodTerms = behvPointsEnds( behvPointsEnds - behvPointsStarts > termMinLength & ...
                [behvPointsStarts(2:end), length(behvPoints)] - behvPointsEnds > termMinBaseline );

            % get neural activity
            currentDir = pwd;
            cd(sessFilepaths.processedDataFolder)
            initNeurInds = round(NeurEMGSync(origDownsampEMGInd(goodInits)*20,...
                frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
            termNeurInds = round(NeurEMGSync(origDownsampEMGInd(goodTerms)*20,...
                frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
            cd(currentDir)

            % don't use inits or terms that are outside of the session or
            % are nan
            initNeurInds(isnan(initNeurInds)) = [];
            termNeurInds(isnan(termNeurInds)) = [];

            initNeurInds(initNeurInds >= floor(frameNeuropixelSamples{1}{end}(end)/300)-postTransPoints) = [];
            termNeurInds(termNeurInds >= floor(frameNeuropixelSamples{1}{end}(end)/300)-postTransPoints) = [];

            for iInit = 1:length(initNeurInds)
                initFRs{iBehv}(:,:,iInit) = allFRs(:,initNeurInds(iInit)-preTransPoints:initNeurInds(iInit)+postTransPoints);
            end

            for iTerm = 1:length(termNeurInds)
                termFRs{iBehv}(:,:,iTerm) = allFRs(:,termNeurInds(iTerm)-preTransPoints:termNeurInds(iTerm)+postTransPoints);
            end

            % get shifts
            shiftAmount = [];
            for iShift = 1:nShifts

                shiftAmount(iShift) = 0;
                while shiftAmount(iShift) < 3000 | shiftAmount(iShift) > size(allFRs,2) - 3000
                    shiftAmount(iShift) = randi(size(allFRs,2),1);
                end

                % just add the shift and wrap around those that go through the end of the session
                initNeurIndsShift = initNeurInds + shiftAmount(iShift);
                initNeurIndsShift(initNeurIndsShift >= size(allFRs,2)) = initNeurIndsShift(initNeurIndsShift >= size(allFRs,2)) - size(allFRs,2);
                initNeurIndsShift(initNeurIndsShift <= preTransPoints) = initNeurIndsShift(initNeurIndsShift <= preTransPoints) + preTransPoints + 1;
                initNeurIndsShift(initNeurIndsShift >= size(allFRs,2) - postTransPoints) = initNeurIndsShift(initNeurIndsShift >= size(allFRs,2) - postTransPoints) - postTransPoints - 1;

                termNeurIndsShift = termNeurInds + shiftAmount(iShift);
                termNeurIndsShift(termNeurIndsShift >= size(allFRs,2)) = termNeurIndsShift(termNeurIndsShift >= size(allFRs,2)) - size(allFRs,2);
                termNeurIndsShift(termNeurIndsShift <= preTransPoints) = termNeurIndsShift(termNeurIndsShift <= preTransPoints) + preTransPoints + 1;
                termNeurIndsShift(termNeurIndsShift >= size(allFRs,2) - postTransPoints) = termNeurIndsShift(termNeurIndsShift >= size(allFRs,2) - postTransPoints) - postTransPoints - 1;

                for iInit = 1:length(initNeurIndsShift)
                    data = allFRs(:,initNeurIndsShift(iInit)-preTransPoints:initNeurIndsShift(iInit)+postTransPoints);
                    initFRsShift{iBehv}(:,:,iInit,iShift) = data;
                end

                for iTerm = 1:length(termNeurIndsShift)
                    data = allFRs(:,termNeurIndsShift(iTerm)-preTransPoints:termNeurIndsShift(iTerm)+postTransPoints);
                    termFRsShift{iBehv}(:,:,iTerm,iShift) = data;
                end

            end % of shifts

%             % also got behavior specificities (using same approach as with the 7
%             % behavior zones)
%             currentDir = pwd;
%             cd(fullfile(sessions{iSess},'ProcessedData'))
%             behvNeurInds = round(NeurEMGSync(origDownsampEMGInd(find(behvPoints))*20,...
%                 frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
%             cd(currentDir)
%             behvNeurInds(isnan(behvNeurInds)) = [];
%             behvNeurInds(behvNeurInds > size(allFRs,2)) = [];
%             behvAveFRs{iSess}(:,iBehv) = nanmean(allFRs(:,behvNeurInds),2);
% 
%             for iShift = 1:nShifts
%                 shiftAmount = 0;
%                 while shiftAmount < 6000 | shiftAmount > size(allFRs,2) - 6000
%                     shiftAmount = randi(size(allFRs,2),1);
%                 end
%                 shiftBehvs = circshift(behvPoints,shiftAmount);
%                 currentDir = pwd;
%                 cd(fullfile(sessions{iSess},'ProcessedData'))
%                 behvNeurShiftInds = round(NeurEMGSync(origDownsampEMGInd(find(shiftBehvs))*20,...
%                     frameEMGSamples, frameNeuropixelSamples, 'EMG')/300);
%                 cd(currentDir)
%                 behvNeurShiftInds(isnan(behvNeurShiftInds)) = [];
%                 behvNeurShiftInds(behvNeurShiftInds > size(allFRs,2)) = [];
%                 behvAveFRsShift{iSess}(:,iBehv,iShift) = nanmean(allFRs(:,behvNeurShiftInds),2);
%             end

        end %behaviors

%         shiftsMean = mean(behvAveFRsShift{iSess},3);
%         shiftsStd = std(behvAveFRsShift{iSess},[],3);
%         modulation.str{iSess} = behvAveFRs{iSess}(1:length(striatumInds),:) >= shiftsMean(1:length(striatumInds),:) + 3*shiftsStd(1:length(striatumInds),:);
%         modulation.ctx{iSess} = behvAveFRs{iSess}(length(striatumInds)+1:end,:) >= shiftsMean(length(striatumInds)+1:end,:) + 3*shiftsStd(length(striatumInds)+1:end,:);
%         for iBehv = 1:nBehvs
%             regionSpec.str{iSess}(:,iBehv) = modulation.str{iSess}(:,iBehv) & ~any(modulation.str{iSess}(:,setdiff(1:nBehvs,iBehv)),2);
%             regionSpec.ctx{iSess}(:,iBehv) = modulation.ctx{iSess}(:,iBehv) & ~any(modulation.ctx{iSess}(:,setdiff(1:nBehvs,iBehv)),2);
%         end

        save(sessFilepaths.UMAPFile,'classifierLabels','classifierBehvs','-append')
        save(sessFilepaths.EMGClassifierInits,'initFRs','termFRs','initFRsShift','termFRsShift','shiftAmount','goodInits','goodTerms','-v7.3')


    end % of sessions loop

end  % of animals loop



%
