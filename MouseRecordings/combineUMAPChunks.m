function reduction = combineUMAPChunks(baseDir,umapInputName,chunkBaseName,nChunks)
%COMBINEUMAPCHUNKS Summary of this function goes here
%   Detailed explanation goes here

load(fullfile(baseDir,umapInputName));

reduction = [];
for iChunk = 1:nChunks

    load(fullfile(baseDir,'UMAPChunks',[chunkBaseName num2str(iChunk)]),'reduction_chunk')
    reduction(iChunk:nChunks:length(origDownsampEMGInd),:) = reduction_chunk;

end

save(fullfile(baseDir,'UMAPMixed.mat'),'reduction','freqData','origDownsampEMGInd','analyzedBehaviors','behvLabels','behvLabelsNoArt','-v7.3')


