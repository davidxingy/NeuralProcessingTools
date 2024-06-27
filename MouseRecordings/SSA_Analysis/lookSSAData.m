clear
close all

% do analysis for each of the datasets
recordingSessions = {
    'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording', ...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording', ...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording' ...
};

behvAlignPerm = [
    1 2 3 4 5 6 7; ...
    7 6 3 5 4 2 1; ...
    1 2 4 5 3 6 7 ...
    ];


for iSess = 1:length(recordingSessions)

    baseDir = recordingSessions{iSess};
    load(fullfile(baseDir,'ProcessedData','SSA','allTimePointsSSA.mat'))
    ssaUMAP = load(fullfile(baseDir,'ProcessedData','UMAPFRs','SSAStrRegionProps.mat'),'regionAveSigs','regionAveSigsShuff');
%     pcaUMAP = load(fullfile(baseDir,'ProcessedData','UMAPFRs','PCACtxRegionProps.mat'),'regionAveSigs','regionAveSigsShuff');
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'allFRs')

    nanInds = find(any(isnan(allFRs)));
    allFRs(:,nanInds) = [];
    ssaInputData = allFRs(goodNeurons{2},:)';
    ssaInputData = ssaInputData - mean(ssaInputData);
    clear allFRs
        
    %PCA
    [pcaWeights, pcaTrajs] = pca(ssaInputData);

    ssaTrajVars{iSess} = var(ssaResults{2}.trajs);
    pcaVars{iSess} = var(pcaTrajs);
    origDataVars{iSess} = var(ssaInputData);

    %calculate multi-behavioral specificity
    ssaAveSigNorm = abs(ssaUMAP.regionAveSigs(:,behvAlignPerm(iSess,:)))./sum(abs(ssaUMAP.regionAveSigs(:,behvAlignPerm(iSess,:))),2);
    ssaMultiSpec{iSess} = sqrt(squeeze(sum(ssaAveSigNorm.^2,2)));

    ssaAveSigShuffNorm = abs(ssaUMAP.regionAveSigsShuff(:,behvAlignPerm(iSess,:),:))./sum(abs(ssaUMAP.regionAveSigsShuff(:,behvAlignPerm(iSess,:),:)),2);
    ssaMultiSpecShuff{iSess} = sqrt(squeeze(sum(ssaAveSigShuffNorm.^2,2)));

%     pcaAveSigNorm = abs(pcaUMAP.regionAveSigs(1:20,behvAlignPerm(iSess,:)))./sum(abs(pcaUMAP.regionAveSigs(1:20,behvAlignPerm(iSess,:))),2);
%     pcaMultiSpec{iSess} = sqrt(squeeze(sum(pcaAveSigNorm.^2,2)));

%     pcaAveSigShuffNorm = abs(pcaUMAP.regionAveSigsShuff(1:20,behvAlignPerm(iSess,:),:))./sum(abs(pcaUMAP.regionAveSigsShuff(1:20,behvAlignPerm(iSess,:),:)),2);
%     pcaMultiSpecShuff{iSess} = sqrt(squeeze(sum(pcaAveSigShuffNorm.^2,2)));

    %calc CDF for the multibehavioral selectivity
    multiSpecCdfVals = 0:0.01:1;
    ssaMultiSpecCdfFreq{iSess} = calcCDF(ssaMultiSpec{iSess},multiSpecCdfVals);
%     pcaMultiSpecCdfFreq{iSess} = calcCDF(pcaMultiSpec{iSess},multiSpecCdfVals);

end

catSessSSASpec = calcCDF(cat(1,ssaMultiSpec{:}),multiSpecCdfVals);
% catSessPCASpec = calcCDF(cat(1,pcaMultiSpec{:}),multiSpecCdfVals);




function freq = calcCDF(inputData,range)

% remove nans
inputData(isnan(inputData)) = [];

for i = 1:length(range)

    freq(i) = sum(inputData <= range(i)) / length(inputData);

end

end



% 
