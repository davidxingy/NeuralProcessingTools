function [behavioralData, allNeurInds, allEMGInds, unusedBouts] = splitIntoAnnotatedBehaviors(baseRecordingFolder, neuralDataFile, analyzedBehaviors, badEMGChans, binSize, gaussStd)

if isempty(neuralDataFile)
    onlyEMG = true;
    allNeurInds = [];
else
    onlyEMG = false;
end

load(fullfile(baseRecordingFolder,'ProcessedData','VideoSyncFrames.mat'))
load(fullfile(baseRecordingFolder,'ProcessedData','BehaviorAnnotations','BehaviorLabels.mat'))

if ~onlyEMG
    load(fullfile(baseRecordingFolder,'ProcessedData','neuronDataStruct.mat'))
    load(fullfile(baseRecordingFolder,'ProcessedData',neuralDataFile),'noNanFRs')
    load(fullfile(baseRecordingFolder,'Neuropixels','artifactTimestamps.mat'))
end

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
downsampEMG = downsample(allEMG',(20*binSize))';

if ~onlyEMG
    neuralEMGOffset = round(frameNeuropixelSamples{1}{1}(1)/(30*binSize) - frameEMGSamples{1}{1}(1)/(20*binSize));
    maxSamples = frameNeuropixelSamples{1}{end}(end);
    artifacts = histcounts(artifactTS,1:(30*binSize):maxSamples);
    artifactNeurBins = find(convGauss(artifacts, 1, gaussStd));
    % get datapoints which has neural artifacts
    artifactEMGInds = find(histcounts(...
        (artifactTS-frameNeuropixelSamples{1}{1}(1))/3*2+frameEMGSamples{1}{1}(1),1:(20*binSize):maxSamples));
end

% get datapoints which has emg artifacts
if onlyEMG
    artifactEMGInds = find(histcounts(removedInds,1:(20*binSize):size(allEMG,2)));
else
    artifactEMGInds = unique([artifactEMGInds find(histcounts(removedInds,1:(20*binSize):size(allEMG,2)))]);
end

clear allEMG
downsampEMG(badEMGChans,:) = [];

for iBehv = 1:length(analyzedBehaviors)
    
    iSection = 1;
    iAllSection = 1;
    boutFRs = {};
    boutEMGs = {};
    behav = analyzedBehaviors{iBehv};
    unusedBouts{iBehv} = [];

    
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
            
            if ~onlyEMG
                
                boutNeuralStart = round(frameNeuropixelSamples{1}{iVid}(frameInds(boutStarts(iBout)))/(30*binSize));
                boutNeuralEnd = round(frameNeuropixelSamples{1}{iVid}(frameInds(boutEnds(iBout)))/(30*binSize));
                boutNeuralInds = boutNeuralStart:boutNeuralEnd;
                [~, badNeurInds] = intersect(boutNeuralInds,artifactNeurBins);
                %             boutNeuralInds(badBoutInds)=[];
                
                %check to make sure the data is valid (there are some datasets
                %where there was dropped data packets from spikeGLX)
                if isnan(boutNeuralStart) || isnan(boutNeuralEnd)
                    unusedBouts{iBehv}(iAllSection) = 1;
                    iAllSection = iAllSection + 1;
                    continue
                else
                    unusedBouts{iBehv}(iAllSection) = 0;
                    iAllSection = iAllSection + 1;
                end
                
            end
            
            boutEMGStart = round(frameEMGSamples{1}{iVid}(frameInds(boutStarts(iBout)))/(20*binSize));
            boutEMGEnd = round(frameEMGSamples{1}{iVid}(frameInds(boutEnds(iBout)))/(20*binSize));
            
            if ~onlyEMG
                
                boutSegDiff = boutEMGEnd-boutEMGStart - (boutNeuralEnd-boutNeuralStart);
                if boutSegDiff == 0
                    boutEMGInds = boutEMGStart:boutEMGEnd;
                elseif abs(boutSegDiff) <= 1
                    boutEMGInds = boutEMGStart:boutEMGStart+boutNeuralEnd-boutNeuralStart;
                else
                    error('Greater than 1 ind difference in neural and EMG bout segment lengths!')
                end
                
            else
                boutEMGInds = boutEMGStart:boutEMGEnd;
            end
            
            [~, badEMGInds] = intersect(boutEMGInds,artifactEMGInds);
%             boutEMGInds(badBoutInds)=[];
            
            if ~onlyEMG
                thisBoutFR = noNanFRs(:,boutNeuralInds);
                thisBoutFR(:,badNeurInds) = nan;
                boutFRs{iSection} = thisBoutFR;
                allNeurInds{iBehv,iSection} = boutNeuralInds;
            end
            
            thisBoutEMG = downsampEMG(:,boutEMGInds);
            thisBoutEMG(:,badEMGInds) = nan;
            boutEMGs{iSection} = thisBoutEMG;
            allEMGInds{iBehv,iSection} = boutEMGInds;
            
            iSection = iSection+1;
        end
        
        if ~onlyEMG
            allBoutFRs = [boutFRs{:}];
            behavioralData.(behav).boutFRs = boutFRs;
            behavioralData.(behav).allBoutFRs = allBoutFRs;
        end
        
        allBoutEMGs = [boutEMGs{:}];
        behavioralData.(behav).boutEMGs = boutEMGs;
        behavioralData.(behav).allBoutEMGs = allBoutEMGs;
        
    end
end


% 
