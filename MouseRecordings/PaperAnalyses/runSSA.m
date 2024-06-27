clear
close all

% baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording';
baseDir = 'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording';
% baseDir = 'X:\David\ArenaRecordings\D026-032923-ArenaRecording';
nSSADims = 20;
timeBinSize = 10;


pathToEnvExe = 'C:\Users\david\Anaconda3\envs\SSA\python.exe';

if timeBinSize == 10
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'))
elseif timeBinSize == 1
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates1msBins10msGauss.mat'))
end

% remove nan points
nanInds = find(any(isnan(allFRs)));

allFRs(:,nanInds) = [];
cortexFRs(:,nanInds) = [];
striatumFRs(:,nanInds) = [];

% remove neurons with too low of a firing rate
meanFRs = nanmean(allFRs,2);
goodNeurons{3} = find(any(meanFRs*(1000/timeBinSize)>0.2,2));
goodNeurons{1} = intersect(goodNeurons{3},striatumInds);
goodNeurons{2} = intersect(goodNeurons{3},cortexInds);


% run on all three regions
regionNames = {'Striatum','Cortex','StriatumAndCortex'};

for iArea = 1:3

    nSSADims = length(goodNeurons{iArea});
    data = allFRs(goodNeurons{iArea},:)';
    dataCent{iArea} = data - mean(data);
    ssaResults{iArea} = callSSA(dataCent{iArea}, [], nSSADims, pathToEnvExe);

end
