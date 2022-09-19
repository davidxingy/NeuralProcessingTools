function removeArtifactFromBin(binDIR,binFilename,newBinFilename)
% binDIR = 'Z:\David\ArenaRecordings\NeuropixelsTest\D012-022422-ArenaRecording\Neuropixels';
% binFilename = 'D012-022422-ArenaRecording_g0_t0.imec0.ap.bin';
% newBinFilename = 'D012-022422-ArenaRecording_g0_t0_RemovedArtifact.imec0.ap.bin';

savePlots = true;

% load in metadata
meta = ReadMeta(binFilename, binDIR);

% total channel count
nChans = str2double(meta.nSavedChans);

% total number of samples
totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));

% variable holding list of all detected artifact timestamps
artifactTS = [];
smallArtifactTS  = [];
lickArtifactTS = [];

% save plots
if savePlots
    plotHandles = [];
    currentH = 1;
    %make directory for figures
    if ~exist(fullfile(binDIR,'ArtifactRemovalPlots'),'dir')
        mkdir(fullfile(binDIR,'ArtifactRemovalPlots'))
    end
end


% load in data, in segments of 1 000 000 samples  (to avoid running out of
% memory)
segLength = 1000000;
nSegs = ceil(totSamps/segLength);
baselineLength = 1000000;

%open file for reading
readFid = fopen(fullfile(binDIR, binFilename), 'rb');
writeFid = fopen(fullfile(binDIR, newBinFilename), 'wb');

for iSeg = 1:nSegs
    
    if iSeg == nSegs
        samples2Load = totSamps - (nSegs-1)*segLength;
    else
        samples2Load = segLength;
    end

%     fseek(readFid, segLength * 2 * nChans * 100, 'bof');
    dataArray = fread(readFid, [nChans, samples2Load], 'int16=>double');
    
    %plot pre-removal data
    if savePlots
        if iSeg == 2
            currentH = 2;
        end
        plotHandles(currentH) = figure('Color','w','Units','normalized','OuterPosition',[0 0 1 1],'Visible','off');
        plot(dataArray(1:10,:)','color',[0.3 0.3 0.3],'LineWidth',2.5)
        
    end
    
    %detect artifact samples
    [allArtInds smallArtInds lickArtInds] = artifactDetection(dataArray);
    
    %get data to replace the artifact with
    if iSeg == 1
        
        if any(allArtInds<30000)
            warning('Baseline has artifact, aborting')
            return
        end
        
        baseline = dataArray(1:end-1,1:baselineLength);

    end
    
    %add to list
    artifactTS = [artifactTS allArtInds+segLength*(iSeg-1)];
    smallArtifactTS = [smallArtifactTS smallArtInds+segLength*(iSeg-1)];
    lickArtifactTS = [lickArtifactTS lickArtInds+segLength*(iSeg-1)];
    
    %now, remove from the current block, the detected artifacts in the
    %current block, and the artifacts bleeding over from the previous block
    currentBlockInds = allArtInds(allArtInds > 0 & allArtInds <= segLength);
    
    if length(currentBlockInds)>baselineLength
        warning('More artifact indices than baseline, aborting')
        return
    end
    
    dataArray(1:end-1,currentBlockInds) = baseline(:,1:length(currentBlockInds));
    if iSeg ~=1
        dataArray(1:end-1,nextBlockInds) = baseline(:,1:length(nextBlockInds));
    end
    
    %add samples less than 0 to previous block
    prevBlockInds = allArtInds(allArtInds<=0) + segLength;
    
    %save samples after the block end for next block
    nextBlockInds = allArtInds(allArtInds>segLength) - segLength;
    
    if iSeg >= 2
        %now, for the previous block, remove the artifacts that were bleeding
        %over from this block
        prevBlockArray(1:end-1,prevBlockInds) = baseline(:,1:length(prevBlockInds));
        prevBlockArray = int16(prevBlockArray);
        %add final processed data without artifact to the output plot
        if savePlots
                currentH = saveOutputFig(plotHandles, prevBlockArray(1:10,:)', currentH, ...
                    fullfile(binDIR,'ArtifactRemovalPlots',[binFilename '_Block' num2str(iSeg-1) '.png']));
        end
        
        %and save to bin file
        writtenCount = fwrite(writeFid,prevBlockArray,'int16');
        if writtenCount ~= numel(prevBlockArray)
            error('Not all datapoints written to bin file! @Seg %u', iSeg);
        end
        
    end
    
    %current block becomes previous block
    prevBlockArray = int16(dataArray);
    
    %if it's the last block, then ignore any holdovers, and just save this
    %last block
    if iSeg == nSegs
        
        %save currrent data block to bin file
        writtenCount = fwrite(writeFid,prevBlockArray,'int16');
        if writtenCount ~= numel(prevBlockArray)
            error('Not all datapoints written to bin file! @Seg %u', iSeg);
        end
        
        %and save final figure
        if savePlots
            currentH = saveOutputFig(plotHandles, prevBlockArray(1:10,:)', plotHandles,...
                fullfile(binDIR,'ArtifactRemovalPlots',[binFilename '_Block' num2str(iSeg) '.png']));
        end

    end
    
end

% close file
fclose(readFid);
fclose(writeFid);
close all

% save list of all artifact timestamps
save(fullfile(binDIR, 'artifactTimestamps.mat'),'artifactTS','smallArtifactTS','lickArtifactTS')

end


function [allArtInds smallArtInds lickArtInds] = artifactDetection(dataArray)
% function for finding indices of artifacts in the data

% criterion is more than 80% of all channels experience a large increase in
% signal amplitude over the last 3 samples
% also, any voltage levels past 300 is considered a lick artifact
% numbers based on emperical observations of D012-022422-ArenaRecording 
smallArtThresh = 150;
smallArtDiffSamples = 3;
smallArtFractionChans = 0.8;

lickArtDiffThresh = 750;
lickArtDiffSamples = 3;
lickArtDiffFractionChans = 0.8;
lickArtAbsThresh = 400;
lickArtAbsFractionChans = 0.25; %For such a large value, need less channels to cross it to be considered an artifact

% first small artifact
changes = abs(dataArray(:,1:end-smallArtDiffSamples) - dataArray(:,smallArtDiffSamples+1:end));
smallArts = [zeros(1, smallArtDiffSamples) sum(changes >= smallArtThresh,1) >= smallArtFractionChans * size(dataArray,1)];
% next lick artifacts
changes = abs(dataArray(:,1:end-lickArtDiffSamples) - dataArray(:,lickArtDiffSamples+1:end));
lickArts = [zeros(1, lickArtDiffSamples) sum(changes >= lickArtDiffThresh,1) >= lickArtDiffFractionChans * size(dataArray,1)];
lickArts = lickArts | sum(abs(dataArray) >= lickArtAbsThresh,1) >= lickArtAbsFractionChans * size(dataArray,1);

% add a buffer zone around artifact for removal as well
smallArtPreBuffer = 20;
smallArtPostBuffer = 50;
lickArtPreBuffer = 3000;
lickArtPostBuffer = 10000;

% first small artifacts
smallBuffInds = getBufferInds(smallArts,smallArtPreBuffer,smallArtPostBuffer);

% next lick artifacts
lickBuffInds = getBufferInds(lickArts,lickArtPreBuffer,lickArtPostBuffer);

smallArtInds = unique([find(smallArts) smallBuffInds]);
lickArtInds = unique([find(lickArts) lickBuffInds]);
allArtInds = unique([smallArtInds lickArtInds]);


end



function buffInds = getBufferInds(artBins,nPreBuffer,nPostBuffer)
% helper function to get buffer inds surrounding artifacts

% first, pre-artifact buffer
artStarts = find(artBins(2:end)==1 & artBins(1:end-1)==0)+1;
preBuffInds = zeros(length(artStarts),nPreBuffer);
for iPreStart = 1:length(artStarts)
    
    preBuffInds(iPreStart,:) = artStarts(iPreStart)-nPreBuffer:artStarts(iPreStart)-1;
    
end

% next, post-artifact buffer
artStops = find(artBins(2:end)==0 & artBins(1:end-1)==1);
postBuffInds = zeros(length(artStops),nPostBuffer);
for iPostStop = 1:length(artStops)
    
    postBuffInds(iPostStop,:) = artStops(iPostStop)+1:artStops(iPostStop)+nPostBuffer;
    
end

buffInds = [preBuffInds(:); postBuffInds(:)]';

end


function currentH = saveOutputFig(plotHandles, plotData, currentH, filename)
% helper function to save output plot
if currentH == 1
    currentH = 2;
else
    currentH = 1;
end

figure(plotHandles(currentH));
hold on;
plot(plotData,'LineWidth',1.8);
box off
saveas(gcf,filename,'png')
close(plotHandles(currentH));

end


% 
