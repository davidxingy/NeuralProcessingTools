function [behavioralData, allNeurInds, allEMGInds] = splitIntoAnnotatedBehaviors(baseRecordingFolder,neuralDataFile, analyzedBehaviors, badEMGChans)

load(fullfile(baseRecordingFolder,'ProcessedData','VideoSyncFrames.mat'))
load(fullfile(baseRecordingFolder,'ProcessedData','neuronDataStruct.mat'))
load(fullfile(baseRecordingFolder,'ProcessedData',neuralDataFile),'noNanFRs')
load(fullfile(baseRecordingFolder,'Neuropixels','artifactTimestamps.mat'))
load(fullfile(baseRecordingFolder,'ProcessedData','BehaviorAnnotations','BehaviorLabels.mat'))


% load and concatenate all EMG data
allEMG = [];

allProcessedFiles = string(ls(fullfile(baseRecordingFolder,'processedData')));
processedEMGFiles = allProcessedFiles(contains(allProcessedFiles,'_ProcessedEMG_Block'));

% make sure to load the EMG data parts in order
for iPart = 1:length(processedEMGFiles)
    fileNameParts = split(processedEMGFiles{iPart}(1:end-4),'_');
    baseName = fileNameParts{1};
    partNum(iPart) = str2double(fileNameParts{end}(5:end));
end

for iPart = 1:length(processedEMGFiles)

    fileInd = find(partNum == iPart);
    load(fullfile(baseRecordingFolder,'processedData',processedEMGFiles{fileInd}))

    allEMG = [allEMG processedEMG];

end

clear processedEMG
clear filteredEMG

% also load in metaData
load(fullfile(baseRecordingFolder,'ProcessedData',[baseName '_ProcessedEMG_MetaData']))


% downsample and align to neural data
neuralEMGOffset = round(frameNeuropixelSamples{1}{1}(1)/30 - frameEMGSamples{1}{1}(1)/20);
downsampEMG = downsample(allEMG',20)';

maxSamples = frameNeuropixelSamples{1}{end}(end);

artifacts = histcounts(artifactTS,1:30:maxSamples);
artifactNeurBins = find(convGauss(artifacts, 1, 10));

% get datapoints which has neural artifacts
artifactEMGInds = find(histcounts(...
    (artifactTS-frameNeuropixelSamples{1}{1}(1))/30*20+frameEMGSamples{1}{1}(1),1:20:maxSamples));

% get datapoints which has emg artifacts
artifactEMGInds = unique([artifactEMGInds find(histcounts(removedInds,1:20:size(allEMG,2)))]);

clear allEMG

for iBehv = 1:length(analyzedBehaviors)
    
    iSection = 1;
    boutFRs = {};
    boutEMGs = {};
    behav = analyzedBehaviors{iBehv};
    
    behavInd = find(strcmp(string(behaviors),behav));
    for iVid = 1:length(behaviorFrames{behavInd})
        
        frameInds = behaviorFrames{behavInd}{iVid};
        
        if isempty(frameInds)
            continue
        end
        
        %divide into bouts
        boutTransitions = find(diff(frameInds)>1);
        boutStarts = [1 boutTransitions+1];
        boutEnds = [boutTransitions length(frameInds)];
        
        for iBout = 1:length(boutStarts)
            boutNeuralStart = round(frameNeuropixelSamples{1}{iVid}(frameInds(boutStarts(iBout)))/30);
            boutNeuralEnd = round(frameNeuropixelSamples{1}{iVid}(frameInds(boutEnds(iBout)))/30);
            boutNeuralInds = boutNeuralStart:boutNeuralEnd;
            [~, badNeurInds] = intersect(boutNeuralInds,artifactNeurBins);
%             boutNeuralInds(badBoutInds)=[];
            
            boutEMGStart = round(frameEMGSamples{1}{iVid}(frameInds(boutStarts(iBout)))/20);
            boutEMGEnd = round(frameEMGSamples{1}{iVid}(frameInds(boutEnds(iBout)))/20);
            
            boutSegDiff = boutEMGEnd-boutEMGStart - (boutNeuralEnd-boutNeuralStart);
            if boutSegDiff == 0
                boutEMGInds = boutEMGStart:boutEMGEnd;
            elseif abs(boutSegDiff) <= 1
                boutEMGInds = boutEMGStart:boutEMGStart+boutNeuralEnd-boutNeuralStart;
            else
                error('Greater than 1 ind difference in neural and EMG bout segment lengths!')
            end
            
            [~, badEMGInds] = intersect(boutEMGInds,artifactEMGInds);
%             boutEMGInds(badBoutInds)=[];
            
            thisBoutFR = noNanFRs(:,boutNeuralInds);
            thisBoutFR(:,badNeurInds) = nan;
            boutFRs{iSection} = thisBoutFR;
            allNeurInds{iBehv,iSection} = boutNeuralInds;

            thisBoutEMG = downsampEMG(:,boutEMGInds);
            thisBoutEMG(:,badEMGInds) = nan;
            boutEMGs{iSection} = thisBoutEMG;
            allEMGInds{iBehv,iSection} = boutEMGInds;
            
            iSection = iSection+1;
        end
        
        allBoutFRs = [boutFRs{:}];
        allBoutEMGs = [boutEMGs{:}];
        
        behavioralData.(behav).boutFRs = boutFRs;
        behavioralData.(behav).allBoutFRs = allBoutFRs;
                
        behavioralData.(behav).boutEMGs = boutEMGs;
        behavioralData.(behav).allBoutEMGs = allBoutEMGs;
        
    end
end


% 
