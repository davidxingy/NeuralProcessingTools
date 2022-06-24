function syncInds = syncNeuropixels(rawDataFilename,vidSyncFile, pulseDuration, nPulsesPerTrain)

syncThresh = 50;

[baseDir,filename] = fileparts(rawDataFilename);

% load in metadata
filesInBaseDir = string(ls(baseDir));

metaFilename = filesInBaseDir(contains(filesInBaseDir,'ap.meta'));
if size(metaFilename) ~= 1
    warning('Unable to find unique meta file!')
    syncInds = [];
    return
end

meta = ReadMeta(strtrim(metaFilename{1}),baseDir);

% total channel count
nChans = str2double(meta.nSavedChans);

% total number of samples
totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));

% load sync channel data in segments since can't load in all at once
syncSignal = zeros(1,totSamps);
segLength = 1000000;
nSegs = ceil(totSamps/segLength);
readFid = fopen(rawDataFilename, 'rb');
for iSeg = 1:nSegs
    
    if iSeg == nSegs
        samples2Load = totSamps - (nSegs-1)*segLength;
    else
        samples2Load = segLength;
    end
    
    dataArray = fread(readFid, [nChans, samples2Load], 'int16=>double');
    syncSignal(1+segLength*(iSeg-1):segLength*(iSeg-1)+samples2Load) = dataArray(end,:);
    
end

% get timing of rising edges
risingEdges = find(syncSignal(2:end)>=syncThresh & syncSignal(1:end-1)<syncThresh)+1;

vidSyncVars = load(vidSyncFile);

if ~isfield(vidSyncVars,'risingEdgeFrames') || ~isfield(vidSyncVars,'ledValues')
    
    warning('Video sync file did not contain the ''risingEdgeFrames'' and  ''ledValues'' variables!')
    syncInds = [];
    return 
    
end

syncInds = syncFLIR(risingEdges,vidSyncVars.risingEdgeFrames,pulseDuration,nPulsesPerTrain,30000,...
    cellfun(@(x) cellfun(@length, x),vidSyncVars.ledValues,'un',0));




% 

