function [ledValues, risingEdgeFrames, chunkFiles] = extractVidSync(vidDir, vidBaseName, thresh)

if isempty(thresh)
    LEDThresh = 150;
else
    LEDThresh = thresh;
end

% get all file with the basename
allFilenames = string(ls(vidDir));
vidFilenames = allFilenames(contains(allFilenames,vidBaseName));

% divide into chunks and files
chunkValues = unique(extract(vidFilenames,'Vid' + digitsPattern));
for iChunk = 1:length(chunkValues)
    
    %get the files related to this chunk
    unsortedChunkFiles = vidFilenames(contains(vidFilenames,chunkValues{iChunk}));
    
    %make sure they are in correct order (based on last numbering in the filename)
    fileNums = extract(unsortedChunkFiles,digitsPattern);
    if size(fileNums,2) == 1
        fileNums = fileNums';
    end
    [~, sortInd] = sort(str2double(fileNums(:,end)));
    chunkFiles{iChunk} = unsortedChunkFiles(sortInd);
    
end

% now get LED sync signal for each file

for iChunk = 1:length(chunkFiles)
    
    for iFile = 1:length(chunkFiles{iChunk})
        
        %load in the video
        v = VideoReader([vidDir '\' chunkFiles{iChunk}{iFile}]);
        
        %load in first frame 
        firstFrame = readFrame(v);
        
        %get ROI for the LED
        if iFile==1 && iChunk==1
            ledROI = GetLEDExtractionRegions(double(rgb2gray(firstFrame)),1,[1 size(firstFrame,2)] ,[1 size(firstFrame,1)]);
            close all;
        end
        
        %get LED intensity
        [ledStatus, intensities] = ExtractLEDs([vidDir '\' chunkFiles{iChunk}{iFile}],...
            1, 200, 'gray', 0, ledROI);
        ledValues{iChunk}{iFile} = intensities.mean;
        
        %also threshold at 150
        risingEdgeFrames{iChunk}{iFile} = find(intensities.mean(2:end) >= LEDThresh & intensities.mean(1:end-1) < LEDThresh) + 1;
        
    end
    
end


% 
