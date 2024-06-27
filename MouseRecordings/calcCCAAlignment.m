clear

baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording';

load('NeuralFiringRates1msBins10msGauss.mat','cortexInds','striatumInds','rejectedNeurons')
load(fullfile(baseDir,'ProcessedData','UMAP'),'analyzedBehaviors','behvLabelsNoArt',...
    'origDownsampEMGInd','regionAssignmentsNoBound','regionBehvAssignments','regionWatershedLabels')
load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'))
load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))

% load('emgDecodingM1.mat','cvFoldsIndsAll')

% load(fullfile(baseDir,'ProcessedData','EpochedData1ms.mat'))
% behaviorNames = fieldnames(behavioralData);

% look at different lags
lags = -100:10:100; %in bin sizes, which for now is 10s of ms
lags = 0;


% remove low-fr neurons
meanFRLimit = 0.0001;
for iBehv = 1:length(behaviorNames)
    allFRs{iBehv} = behavioralData.(behaviorNames{iBehv}).allBoutFRs(regionInds,:);
end
allFRsMean = nanmean(cat(2,allFRs{:}),2);
highRegionInds = regionInds(allFRsMean > meanFRLimit);

goodEMGChans = 1:6;

nPerms = 200;

behaviorNames(9) = [];

for iLag = 1:length(lags)
    
    for iBehv1 = 1:length(behaviorNames)

        for iBehv2 = 1:length(behaviorNames)

            [allNeurs1, allEMG1, boutNeur1, boutEMG1] = extractAndProcessData(behavioralData, behaviorNames{iBehv1}, highRegionInds, goodEMGChans, 0);

            [allNeurs2, allEMG2, boutNeur2, boutEMG2] = extractAndProcessData(behavioralData, behaviorNames{iBehv2}, highRegionInds, goodEMGChans, 0);

            if iBehv1 ~= iBehv2

                %do PCA on neural data first
                %         [projPCA{iBehv,iLag}, trajPCA{iBehv,iLag}, ~, ~, pcaVAF{iBehv,iLag}] = pca(allNeurs');
                %CCA

                %do random permutation shuffling of the two
                bothNeurs = [allNeurs1 allNeurs2];
                bothEMGs = [allEMG1 allEMG2];

                for iPerm = 1:nPerms

                    permInds = randperm(size(bothNeurs,2));
                    split1Inds = permInds(1:size(allNeurs1,2));
                    split2Inds = permInds(size(allNeurs1,2)+1:end);

                    permNeur1 = bothNeurs(:,split1Inds);
                    permEMG1 = bothEMGs(:,split1Inds);

                    permNeur2 = bothNeurs(:,split2Inds);
                    permEMG2 = bothEMGs(:,split2Inds);

                    %do CCA
                    [projNeurPerm1,projEMGPerm1,cannonCorrsPerm1,trajNeurPerm1,trajEMGPerm1] = canoncorr(permNeur1',permEMG1');
                    [projNeurPerm2,projEMGPerm2,cannonCorrsPerm2,trajNeurPerm2,trajEMGPerm2] = canoncorr(permNeur2',permEMG2');

                    %get alignment metrics
                    [prinAngle ccAlignment] = calcAlignmentMetrics(permNeur1,permNeur2,permEMG1,permEMG2,...
                        projNeurPerm1,projNeurPerm2,projEMGPerm1,projEMGPerm2,cannonCorrsPerm1,cannonCorrsPerm2);

                    permPrinAngles(iBehv1,iBehv2,iPerm) = prinAngle;
                    permCCAlignment(iBehv1,iBehv2,iPerm) = ccAlignment;

                end

                %do CCA without shuffling
                [projNeur1,projEMG1,cannonCorrs1,trajNeur1,trajEMG1] = canoncorr(allNeurs1',allEMG1');
                [projNeur2,projEMG2,cannonCorrs2,trajNeur2,trajEMG2] = canoncorr(allNeurs2',allEMG2');

                %get alignment metrics
                [prinAngle ccAlignment] = calcAlignmentMetrics(allNeurs1,allNeurs2,allEMG1,allEMG2,...
                    projNeur1,projNeur2,projEMG1,projEMG2,cannonCorrs1,cannonCorrs2);

                crossPrinAngles(iBehv1,iBehv2) = prinAngle;
                crossCCAlignment(iBehv1,iBehv2) = ccAlignment;

                %         neurCCA{iBehv,iLag} = allNeurs';
                %         emgCCA{iBehv,iLag} = allEMG';

            else
                [allProjNeurs{iBehv1},allProjEMGs{iBehv1},allCannonCorrs{iBehv1},allTrajNeurs{iBehv1},allTrajEMGs{iBehv1}] = ...
                    canoncorr(allNeurs1',allEMG1');
            end

%             %decoding
%             nHist = 10;
%             histSkip = 10;
%             nFolds = 10;
% 
%             % add history
%             neurHist1 = {};
%             neurHist2 = {};
%             emgHist1 = {};
%             emgHist2 = {};
%             for iTrial = 1:length(boutNeur1)
%                 neurHist1{iTrial} = addHistorySkip(boutNeur1{iTrial}, nHist, histSkip);
%                 emgHist1{iTrial} = boutEMG1{iTrial};
%                 neurHist1{iTrial}(:,1:nHist*histSkip) = [];
%                 emgHist1{iTrial}(:,1:nHist*histSkip,:) = [];
%             end
% 
%             for iTrial = 1:length(boutNeur2)
%                 neurHist2{iTrial} = addHistorySkip(boutNeur2{iTrial}, nHist, histSkip);
%                 emgHist2{iTrial} = boutEMG2{iTrial};
%                 neurHist2{iTrial}(:,1:nHist*histSkip) = [];
%                 emgHist2{iTrial}(:,1:nHist*histSkip,:) = [];
%             end
% 
%             allNeurHist1 = cat(2,neurHist1{:});
%             allEMGHist1 = cat(2,emgHist1{:});
% 
%             nanInds = find(any(isnan(allEMGHist1)));
%             nanInds = unique([nanInds find(any(isnan(allNeurHist1)))]);
% 
%             allNeurHist1(:,nanInds)=[];
%             allEMGHist1(:,nanInds)=[];
% 
%             allNeurHist1 = allNeurHist1(:,1:29000);
%             allEMGHist1 = allEMGHist1(:,1:29000);
% 
%             allNeurHist2 = cat(2,neurHist2{:});
%             allEMGHist2 = cat(2,emgHist2{:});
% 
%             nanInds = find(any(isnan(allEMGHist2)));
%             nanInds = unique([nanInds find(any(isnan(allNeurHist2)))]);
% 
%             allNeurHist2(:,nanInds)=[];
%             allEMGHist2(:,nanInds)=[];
% 
%             allNeurHist2 = allNeurHist2(:,1:29000);
%             allEMGHist2 = allEMGHist2(:,1:29000);
% 
%             if iBehv1 == iBehv2
%                 % do cross vaidation on individual behaviors
%                 [decodingPerformanceCV{iBehv1}, cvFoldsInds{iBehv1}] = decodingCrossValidation(allNeurHist1,allEMGHist1, nFolds);
% 
%             elseif iBehv1 < iBehv2
% %                 % do cross validation on concatenated behaviors
% %                 [decodingPerformanceCrossCV{iBehv1,iBehv2}, cvFoldsIndsCross{iBehv1,iBehv2}] = decodingCrossValidation([allNeurHist1 allNeurHist2], [allEMGHist1 allEMGHist2], 2);
% 
%                 %train on one then decode on the other
% %                 decodingPerformanceGen{iBehv1,iBehv2} = decodingGeneralization(allNeurHist1,allEMGHist1,allNeurHist2,allEMGHist2);
% %                 decodingPerformanceGen{iBehv2,iBehv1} = decodingGeneralization(allNeurHist2,allEMGHist2,allNeurHist1,allEMGHist1);
% 
% 
%                 [decodingPerformanceGen{iBehv1,iBehv2}, decodingPerformanceGenControl{iBehv1,iBehv2}, testFoldInds{iBehv1,iBehv2}] = decodingCombinedGeneralization...
%                     (allNeurHist1,allEMGHist1,allNeurHist2,allEMGHist2,nFolds,behaviorNames{iBehv1},behaviorNames{iBehv2});
%                 
%                 [decodingPerformanceGen{iBehv2,iBehv1}, decodingPerformanceGenControl{iBehv2,iBehv1}, testFoldInds{iBehv2,iBehv1}] = decodingCombinedGeneralization...
%                     (allNeurHist2,allEMGHist2,allNeurHist1,allEMGHist1,nFolds,behaviorNames{iBehv2},behaviorNames{iBehv1});
% 
%             end


        end

        
    end
    
end

metricUsedMuscs = 1:3;

for iBehv1 = 1:length(behaviorNames)
        for iBehv2 = 1:length(behaviorNames)

            if iBehv1 == iBehv2
                continue
            end
            selfCVPerformance = decodingPerformanceCV{iBehv2}.R2;
            crossPerformance = decodingPerformanceGen{iBehv1,iBehv2}.R2;

            meanR2diff(iBehv1,iBehv2) = mean(decodingPerformanceGenControl{iBehv1,iBehv2}.R2 - decodingPerformanceGen{iBehv1,iBehv2}.R2);
            meanR2perct(iBehv1,iBehv2) = mean((decodingPerformanceGenControl{iBehv1,iBehv2}.R2 - decodingPerformanceGen{iBehv1,iBehv2}.R2)./decodingPerformanceCV{iBehv2}.R2);

            meanCCdiff(iBehv1,iBehv2) = mean(decodingPerformanceGenControl{iBehv1,iBehv2}.CC(metricUsedMuscs) - decodingPerformanceGen{iBehv1,iBehv2}.CC(metricUsedMuscs));
            meanCCperct(iBehv1,iBehv2) = mean((decodingPerformanceGenControl{iBehv1,iBehv2}.CC(metricUsedMuscs) - decodingPerformanceGen{iBehv1,iBehv2}.CC(metricUsedMuscs))./decodingPerformanceGenControl{iBehv1,iBehv2}.CC(metricUsedMuscs));

            overfittingMetric(iBehv1,iBehv2) = mean(abs(decodingPerformanceGenControl{iBehv1,iBehv2}.R2));
        end
end

% cluster
cgObj = clustergram(meanCCdiff);
clustOrdering = cellfun(@(x) str2num(x),cgObj.RowLabels);

imagesc(meanCCdiff(clustOrdering,clustOrdering))
set(gca,'XTickLabel',behaviorNames(clustOrdering))
set(gca,'YTickLabel',behaviorNames(clustOrdering))

colorbarH = colorbar;
colorbarH.Ticks = 0:5:70;

% by cc alignment
figure
% cgObj = clustergram(ccAlignment);
% clustOrdering = cellfun(@(x) str2num(x),cgObj.RowLabels);

imagesc(ccAlignment(clustOrdering,clustOrdering))
set(gca,'XTickLabel',behaviorNames(clustOrdering))
set(gca,'YTickLabel',behaviorNames(clustOrdering))

colorbarH = colorbar;
colorbarH.Ticks = 0:0.05:0.45;



function [allNeurs, allEMG, boutNeur, boutEMG] = extractAndProcessData(epochedData, behaviorName, neurChans, emgChans, lag)

        neurData = epochedData.(behaviorName).boutFRs;
        emgData = epochedData.(behaviorName).boutEMGs;

        boutNeur = {};
        boutEMG = {};
        for iBout = 1:length(neurData)
            
%             thisBoutEMG = downsample(emgData{iBout}',10)';
%             thisBoutEMG = thisBoutEMG(:,2:end);
%             thisBoutNeur = neurData{iBout}(:,1:size(emg10ms,2));
            
            thisBoutNeur = neurData{iBout}(neurChans,:);
            thisBoutEMG = emgData{iBout}(emgChans,:);

            if lag < 0
                boutNeur{iBout} = thisBoutNeur(:,-1*lag+1:end);
                boutEMG{iBout} = thisBoutEMG(:,1:end+lag);
            elseif lag > 0
                boutNeur{iBout} = thisBoutNeur(:,1:end-lag);
                boutEMG{iBout} = thisBoutEMG(:,lag+1:end);
            else
                boutNeur{iBout} = thisBoutNeur;
                boutEMG{iBout} = thisBoutEMG;
            end
            
        end
        
        allNeurs = cat(2,boutNeur{:});
        allEMG = cat(2,boutEMG{:});
        
        nanInds = find(any(isnan(allEMG)));
        nanInds = unique([nanInds find(any(isnan(allNeurs)))]);

        allNeurs(:,nanInds)=[];
        allEMG(:,nanInds)=[];

        allNeurs(allNeurs<0.001) = 0;

%         allNeurs = allNeurs(:,1:37000);
%         allEMG = allEMG(:,1:37000);

end


function [prinAngle ccAlignment] = calcAlignmentMetrics(neur1,neur2,emg1,emg2,neurProj1,neurProj2,emgProj1,emgProj2,cc1,cc2)

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
cc1Ratio = sum(crossCC1)/sum(cc1);


centeredNeur2 = neur2 - mean(neur2,2);
centeredEMG2 = emg2 - mean(emg2,2);

trajCrossNeur2 = centeredNeur2'*neurProj1;
trajCrossEMG2 = centeredEMG2'*emgProj1;
crossCC2 = corr(trajCrossNeur2,trajCrossEMG2);
crossCC2 = abs(diag(crossCC2));
cc2Ratio = sum(crossCC2)/sum(cc2);

ccAlignment = mean([cc1Ratio cc2Ratio]);

end


function [decodingPerformance, cvFoldsInds] = decodingCrossValidation(allNeurHist,allEMGHist,nFolds)

% divide into test and training sets
% cvFoldsInds = cvFoldsIndsAll(iBehv,:);
cvFoldsInds = divideBlocks(randperm((size(allNeurHist,2))), nFolds);
% cvFoldsIndsAll(iBehv,:) = cvFoldsInds;

for iFold = 1:nFolds

    testInds = cvFoldsInds{iFold};
    testInds(testInds>size(allNeurHist,2)) = [];
    trainInds = setdiff(1:size(allNeurHist,2), testInds);
    trainInds(trainInds>size(allNeurHist,2)) = [];

    %get data and mean center
    trainNeur = allNeurHist(:,trainInds);
    trainNeur = trainNeur(:,:)-mean(trainNeur,2);
    trainEMG = allEMGHist(:,trainInds);
    trainEMGMeans = mean(trainEMG,2);
    trainEMG = trainEMG(:,:)-trainEMGMeans;
    
    
    %do training only on half of the data
    trainNeurHalf = allNeurHist(:,trainInds(1:round(length(trainInds)/2)));
    trainNeurHalf = trainNeurHalf(:,:)-mean(trainNeurHalf,2);
    trainEMGHalf = allEMGHist(:,trainInds);
    trainEMGHalf = trainEMGHalf(:,:)-mean(trainEMGHalf,2);

    testNeurFold = allNeurHist(:,testInds);
    testNeurFold = testNeurFold(:,:)-mean(testNeurFold,2);
    testNeur = testNeurFold;
    testEMG = allEMGHist(:,testInds);
    realEMG{iFold} = testEMG(:,:)-mean(testEMG,2);

    %get weights
    trainWeights = trainEMG/trainNeur;

    %decode
    estEMG{iFold} = trainWeights*testNeurFold;

end

%get performance metrics
%get metrics for all steps concatenated
decodingPerformance = calcPerformanceMetrics(cat(2,estEMG{:}), cat(2,realEMG{:}));

end



function [decodingPerformanceCrossBehv, decodingPerformanceOneBehv, testFoldInds] = decodingCombinedGeneralization...
    (allNeurHist1,allEMGHist1,allNeurHist2,allEMGHist2,nFolds, behv1Name, behv2Name)

%first do first behavior
foldInds1 = divideBlocks(randperm((size(allNeurHist1,2))), nFolds);
testFoldInds = foldInds1;

for iFold = 1:length(foldInds1)
    
    %train on half within behavior, half other behavior
    testInds1 = foldInds1{iFold};
    trainInds1 = setdiff(1:size(allNeurHist1,2), testInds1);
    
    numTrainPoints = min(floor(length(trainInds1))/2,size(allNeurHist2,2));
    trainNeurCrossBehv = [allNeurHist1(:,trainInds1(1:numTrainPoints)) allNeurHist2(:,1:numTrainPoints)];
    trainNeurCrossBehv = trainNeurCrossBehv(:,:)-mean(trainNeurCrossBehv,2);
    trainNeurOneBehv = allNeurHist1(:,trainInds1(1:numTrainPoints*2));
    trainNeurOneBehv = trainNeurOneBehv(:,:)-mean(trainNeurOneBehv,2);
    
    trainEMGCrossBehv = [allEMGHist1(:,trainInds1(1:numTrainPoints)) allEMGHist2(:,1:numTrainPoints)];
    trainEMGCrossBehv = trainEMGCrossBehv-mean(trainEMGCrossBehv,2);
    trainEMGOneBehv = allEMGHist1(:,trainInds1(1:numTrainPoints*2));
    trainEMGOneBehv = trainEMGOneBehv(:,:)-mean(trainEMGOneBehv,2);
    
    %get weights
    trainWeightsCrossBehv = trainEMGCrossBehv/trainNeurCrossBehv;
    trainWeightsOneBehv = trainEMGOneBehv/trainNeurOneBehv;

    %decode
    testNeurFold = allNeurHist1(:,testInds1);
    testNeurFold = testNeurFold-mean(testNeurFold,2);
    testEMG = allEMGHist1(:,testInds1);
    realEMG{iFold} = testEMG-mean(testEMG,2);
    
    estEMGCrossBehv{iFold} = trainWeightsCrossBehv*testNeurFold;
    estEMGOneBehv{iFold} = trainWeightsOneBehv*testNeurFold;
    
end

%get performance metrics
decodingPerformanceCrossBehv = calcPerformanceMetrics(cat(2,estEMGCrossBehv{:}), cat(2,realEMG{:}));
decodingPerformanceOneBehv = calcPerformanceMetrics(cat(2,estEMGOneBehv{:}), cat(2,realEMG{:}));

% make plots
estEMGCrossBehvAll = cat(2,estEMGCrossBehv{:});
estEMGOneBehvAll = cat(2,estEMGOneBehv{:});
realEMGAll = cat(2,realEMG{:});
foldIndsAll = cat(2,foldInds1{:});

% change back to real emg order
invFoldInds = findInvPermInds(foldIndsAll);
estEMGCrossBehvPlot = estEMGCrossBehvAll(:,invFoldInds);
estEMGOneBehvPlot = estEMGOneBehvAll(:,invFoldInds);
realEMGPlot = realEMGAll(:,invFoldInds);

figH = figure('color','w');

tileHandle = tiledlayout(size(realEMGPlot,1),1);
tileHandle.TileSpacing = 'compact';
tileHandle.Padding = 'compact';

for iChan = 1:size(realEMGPlot,1)
    nexttile
    plot(realEMGPlot(iChan,:),'k','LineWidth',2)
    hold on
    plot(estEMGOneBehvPlot(iChan,:))
    plot(estEMGCrossBehvPlot(iChan,:))

    if iChan == 1
        legend('Actual','Control Estimated','Cross-Behv Estimated')
        title(['Decoding ' behv1Name ' from ' behv2Name])
    end
end

savefig(figH,fullfile('Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording\ProcessedData\EMGCCA\CrossDecodingPlots',...
    ['Decoding_' behv1Name '_from_' behv2Name]))

close(figH)

end



function decodingPerformance = decodingGeneralization(allNeurHistTrain,allEMGHistTrain,allNeurHistTest,allEMGHistTest)

%get data and mean center
trainNeur = allNeurHistTrain-mean(allNeurHistTrain,2);
trainEMGMeans = mean(allEMGHistTrain,2);
trainEMG = allEMGHistTrain-trainEMGMeans;

testNeur = allNeurHistTest-mean(allNeurHistTest,2);
realEMG = allEMGHistTest-mean(allEMGHistTest,2);

%get weights
trainWeights = trainEMG/trainNeur;

%decode
estEMG = trainWeights*testNeur;

%get performance metrics
%get metrics for all steps concatenated
decodingPerformance = calcPerformanceMetrics(estEMG, realEMG);

end


function dataHist = addHistorySkip(data, nHist, skip)

    for iHist = 1:nHist
        
        shiftAmount = iHist*skip;
        shiftData{iHist} = [zeros(size(data,1), shiftAmount) data(:, 1:end-shiftAmount)];
        
    end

    dataHist = [data; cat(1,shiftData{:})];

end


% 
