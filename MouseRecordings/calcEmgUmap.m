function calcEmgUmap(baseDir, recordingSessionName, projDownSamp, nUMapNeighbors, varargin)
% calcUMAP(baseDir,recordingSessionName, badEMGChans, makePlots)
% function to perform UMAP dimensionality reduction using EMG inputs via
% the Berman method (Berman 2014), e.g. doing wavelet decomponsition. Does
% it for 1ms time bins.
% 
% Inputs:
% baseDir -                 String containing the filepath to the directory of
%                           the recording session
% 
% recordingSessionName -    String containing the name of recording session
% 
% projDownSamp -            Double indicating the amount to downsample the
%                           time series by before finding the
%                           projection. Used to avoid "streakiness",
%                           typically 50.
% 
% nUMapNeighbors -          Double indicating the number of neighbors when
%                           calculating UMAP projection, typically 100
% 
% badEMGChans -             Optional. 1xN array, indicating if any of the
%                           EMG channels are bad and shouldn't be used.
%                           Default empty.
% 
% subsetToUse -             Optional. 1xN array of logicals. If there is a
%                           subset of (1ms) time points to run the UMAP on.
%                           N should be the same as the total timepoints
%                           that is being used
% 
% useManualAnnotations -    Optional. True/False indicating whether to
%                           calculated (train) the original projection
%                           using only timepoints that have manual
%                           annotations. Default true
% 
% nUMAPDims -               Optional. Dimension of the UMAP Projection.
%                           Default 2
% 
% David Xing 2/20/2023

% parse optional input parameters
narginchk(3, 8)

if length(varargin)>=1
    if isempty(varargin{1})
        badEMGChans = [];
    else
        badEMGChans = varargin{1};
    end
else
    badEMGChans = [];
end

if length(varargin)>=2
    if isempty(varargin{2})
        subsetToUse = [];
    else
        subsetToUse = varargin{2};
    end
else
    subsetToUse = [];
end

if length(varargin)>=3
    if isempty(varargin{3})
        useManualAnnotations = true;
    else
        useManualAnnotations = varargin{3};
    end
else
    useManualAnnotations = true;
end

if length(varargin)>=4
    if isempty(varargin{4})
        nUMAPDims = 2;
    else
        nUMAPDims = varargin{4};
    end
else
    nUMAPDims = 2;
end


% load downsampled EMG Data
load(fullfile(baseDir,'ProcessedData','EMG1ms.mat'))
load(fullfile(baseDir,'ProcessedData',[recordingSessionName '_ProcessedEMG_MetaData.mat']))
nDownsampPoints = size(downsampEMG,2);

% % % % load video and neuropixel sync
% load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
% % % 
% % % % load in neuropixel artifact times
% % % load(fullfile(baseDir,'Neuropixels','artifactTimestamps.mat'))
% % % maxSamples = frameNeuropixelSamples{1}{end}(end);
% % % 
% % % % get datapoints which has neural lick artifacts and also don't use datapoints
% % % % where the neural data was removed (e.g. D024 where packets were dropped)
% % % neuroRemovedInds = getRejectedDataInds(frameNeuropixelSamples);
% % % badEMGInds = NeurEMGSync([lickArtifactTS neuroRemovedInds], frameEMGSamples, frameNeuropixelSamples, 'Neuropixel');
% % % % get it in downsampled 1ms bins
% % % artifactEMGInds = find(histcounts(badEMGInds,1:20:maxSamples));
% % % % artifactEMGInds = find(histcounts(...
% % % %     (lickArtifactTS-frameNeuropixelSamples{1}{1}(1))/30*20+frameEMGSamples{1}{1}(1),1:20:maxSamples));
% % % artifactEMGInds(artifactEMGInds > size(downsampEMG,2)) = [];
% % % % also don't use bins which have EMG artifacts
% % % artifactEMGInds = unique([artifactEMGInds find(histcounts(removedInds,1:20:size(downsampEMG,2)*20))]);
artifactEMGInds = find(histcounts(removedInds,1:20:size(downsampEMG,2)*20));

% load in behavaior labels
load(fullfile(baseDir,'ProcessedData','BehaviorAnnotations','BehaviorLabels.mat'))
analyzedBehaviors = behaviors([1:10]); %{'grooming','eating','walkgrid','walkflat','rearing','climbup','climbdown','still','jumping','jumpdown'};

% load in manually annotated data to use as the supervised UMAP projection
% load('EpochedData1ms.mat')

allData = downsampEMG;
badChans = badEMGChans;
allData(badChans,:) = [];
dataChannels = 1:size(allData,1);
% frameSyncs = frameEMGSamples;
artifactBins = artifactEMGInds;

% Berman frequency decomposition method
waveletParams.omega0 = 5; %default used by Berman
waveletParams.numPeriods = 50;
waveletParams.samplingFreq = 1000;
waveletParams.maxF = 20;
waveletParams.minF = 0.5;

% remove artifact inds and also run frequency decomposition for each
% continuous-nonartifact segment to avoid transitions contaminating the
% wavelet decompositon)

% load in artifacts from NP
load(fullfile(baseDir,'Neuropixels','artifactTimestamps.mat'),'lickArtifactTS')
load(fullfile(baseDir,'ProcessedData','VideoSyncFrames.mat'))
lickEMGInds = round(NeurEMGSync(lickArtifactTS,frameEMGSamples,frameNeuropixelSamples,'Neuropixel')/20);
lickEMGInds(isnan(lickEMGInds)) = [];

removedIndsDownSamp = unique([artifactBins lickEMGInds]);

if removedIndsDownSamp(1) == 1 || removedIndsDownSamp(end) == size(allData,2)
%     error('Beginning or last sample is artifact, need to add code to handle this case')
end
artifactTransitions = [find(diff(removedIndsDownSamp)~=1)];
artifactSegStarts = [removedIndsDownSamp(1) removedIndsDownSamp(artifactTransitions+1)];
artifactSegEnds = [removedIndsDownSamp(artifactTransitions) removedIndsDownSamp(end)];

% Don't use any segments that are too short and may not get a good
% frequency calculation, these are refered to as bad segs
minSegLength = 1000;

for iSeg = 1:length(artifactSegStarts)+1

    if iSeg == 1
        segData = allData(:,1:artifactSegStarts(1)-1);
        if size(segData,2) < minSegLength || removedIndsDownSamp(1) == 1
            freqDataSeg{iSeg} = [];
            badSegIndsCell{iSeg} = 1:artifactSegStarts(1)-1;
        else
            [freqDataSeg{iSeg}, f] = findWavelets(segData', size(allData,1), waveletParams);
        end

    elseif iSeg == length(artifactSegStarts)+1
        segData = allData(:,artifactSegEnds(iSeg-1)+1:end);
        if size(segData,2) < minSegLength || removedIndsDownSamp(end) == size(allData,2)
            freqDataSeg{iSeg} = [];
            badSegIndsCell{iSeg} = artifactSegEnds(iSeg-1)+1:size(allData,2);
        else
            [freqDataSeg{iSeg}, ~] = findWavelets(segData', size(allData,1), waveletParams);
        end

    else
        segData = allData(:,artifactSegEnds(iSeg-1)+1:artifactSegStarts(iSeg)-1);
        if size(segData,2) < minSegLength
            freqDataSeg{iSeg} = [];
            badSegIndsCell{iSeg} = artifactSegEnds(iSeg-1)+1:artifactSegStarts(iSeg)-1;
        else
            [freqDataSeg{iSeg}, ~] = findWavelets(segData', size(allData,1), waveletParams);
        end

    end

end

freqData = cat(1,freqDataSeg{:});
badSegInds = cat(2,badSegIndsCell{:});
clear freqDataSeg badSegIndsCell

% also save the original inds before removing the artifacts and bad segs
origDownsampEMGInd = 1:size(allData,2);
origDownsampEMGInd([removedIndsDownSamp badSegInds]) = [];

% remove the artifacts and the bad segs
allData(:,[removedIndsDownSamp badSegInds]) = [];

clear allData
clear downsampEMG

% if makePlots
%     figure
%     plot(1,1)
%     hold on
%     imagesc(freqData(:,1:waveletParams.numPeriods)')
%     plot(allData(1,:)/100+1)
%     ylim([0.5 waveletParams.numPeriods+0.5])
%     set(gca,'YTick',1:waveletParams.numPeriods)
%     set(gca,'YTickLabel',f)
% end

% map behavior labels onto the time points
[~, ~, boutEMGInds] = splitIntoAnnotatedBehaviors(baseDir,[], analyzedBehaviors, [], 1, 10);

% assign behavioral labels to each of the points
behvLabels = zeros(1,nDownsampPoints); %use downsampEMG which was before the artifact inds were removed
for iBehv = 1:length(analyzedBehaviors)
    
    allBehvInds = cat(2,boutEMGInds{iBehv,:});
    behvLabels(allBehvInds) = iBehv;
    
end

% remove artifact inds from the behavioral labels
behvLabelsNoArt = behvLabels;
behvLabelsNoArt([removedIndsDownSamp badSegInds]) = [];

annotatedBehvLabels = find(behvLabelsNoArt~=0);

if isempty(subsetToUse)
    indsToUse = 1:size(freqData,1);
else
    if length(subsetToUse) ~= size(freqData,1)
        error('subsetToUse must be the same size as the input data to the UMAP!')
    end
    indsToUse = find(subsetToUse);
end

%generally found these values to be a good amount
% nUMapNeighbors = 100; 
% projDownSamp = 50;

% use K-L divergence as the distance metric a la Berman
distKL = @(x,y) log(y./sum(y,2))*(x/sum(x))'*-1+repmat((x/sum(x))*log((x/sum(x))'),size(y,1),1);

% Do UMAP
if useManualAnnotations
    % run initial UMAP projection using only manually annotated points
    inputFeatures = freqData(intersect(indsToUse,[annotatedBehvLabels]),:);
    [reduction,umap,clusterIdentifiers,extras] = run_umap(inputFeatures(1:projDownSamp:end,:),...
        'n_components',nUMAPDims,'n_neighbors',nUMapNeighbors,'save_template_file', 'projUMAPTemplate.mat','Distance',distKL);
else
    % run UMAP using all points
    inputFeatures = freqData(indsToUse,:);
    [reduction,umap,clusterIdentifiers,extras] = run_umap(inputFeatures(1:projDownSamp:end,:),...
        'n_components',nUMAPDims,'n_neighbors',nUMapNeighbors,'save_template_file', 'projUMAPTemplate.mat','Distance',distKL);
end

save('UMAP_MixedTemplate','reduction','umap','clusterIdentifiers','extras','freqData','analyzedBehaviors','behvLabels','behvLabelsNoArt','origDownsampEMGInd','indsToUse','-v7.3')
% divide into chunks to do the mapping since it's way too many data points
% to do the projection all at once
nChunks = 50;
reduction = [];
umapProps = {};
clusterIdentifiers = [];
umapExtras = {};

for iChunk = 1:nChunks
    [reduction_chunk,umap_chunk,clusterIdentifiers_chunk,extras_chunk] = run_umap(freqData(iChunk:nChunks:end,:),'template_file','projUMAPTemplate.mat');
    reduction(iChunk:nChunks:size(freqData,1),:) = reduction_chunk;
    clusterIdentifiers(iChunk:nChunks:size(freqData,1)) = clusterIdentifiers_chunk;
    umapProps{iChunk} = umap_chunk;
    umapExtras{iChunk} = extras_chunk;
end

% save projection
save('UMAP_Mixed','reduction','umapProps','clusterIdentifiers','umapExtras','freqData','analyzedBehaviors','behvLabels','behvLabelsNoArt','origDownsampEMGInd','-v7.3')


% 
