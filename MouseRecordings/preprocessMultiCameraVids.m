function preprocessMultiCameraVids(vidDir,baseFilename)
%PREPROCESSMULTICAMERAVIDS Summary of this function goes here
%   Detailed explanation goes here

allFiles = string(ls(vidDir));
nonCalibrationFiles = allFiles(~cellfun(@(x) contains(x,'CameraCalibration'),string(allFiles)));

frameInfoFiles = nonCalibrationFiles(cellfun(@(x) contains(x,'FrameNums'),string(nonCalibrationFiles)));

nCameras = length(frameInfoFiles);
readerH = VideoReader(fullfile(vidDir,[baseFilename '_Camera0-0000.avi']));
sessionFps = readerH.FrameRate;

for iCamera = 1:nCameras    
    frameInfo{iCamera} = dlmread(fullfile(vidDir,[baseFilename '-Camera' num2str(iCamera-1) 'FrameNums.txt']));
    droppedFrames{iCamera} = find(diff(frameInfo{iCamera}(:,1))>1);
    
    cameraVideoFiles = strtrim(nonCalibrationFiles(cellfun(@(x) contains(x,['_Camera' num2str(iCamera-1)]),string(nonCalibrationFiles))));
    currentWriteFrame = 1;

    for iFile = 1:length(cameraVideoFiles)
        
        tic
        
        readerH = VideoReader(fullfile(vidDir,cameraVideoFiles{iFile}));
        newFileName = [cameraVideoFiles{iFile}(1:end-4) '_FrameFilled.avi'];
        writerH = VideoWriter(newFileName);
        writerH.FrameRate = readerH.FrameRate;
        
        open(writerH)
        currentReadFrame = 1;
        
        % if multiple files, need to offset current frame by the number of
        % frames in preceding files
        if iFile == 1
            nPrevFrames(1) = 0;
        end
        nPrevFrames(iFile+1) = readerH.NumFrames;
            
        while readerH.hasFrame
           
            frame = readFrame(readerH);
            writeVideo(writerH,frame);
            
            % if at a skipped frame, repeat the same frame
            if any(currentReadFrame + sum(nPrevFrames(1:iFile)) == droppedFrames{iCamera})
                writeVideo(writerH,frame);
                currentWriteFrame = currentWriteFrame + 1;
            end
        
            currentReadFrame = currentReadFrame + 1;
            currentWriteFrame = currentWriteFrame + 1;
            
        end
        
        close(writerH)
        
        % double check that we went through the expected number of frames
        if readerH.NumFrames ~= currentReadFrame-1
            warning(['Frame counts off for camera ' num2str(iCamera-1) ', file ' num2str(iFile-1)])
        end
        
        disp(['Camera ' num2str(iCamera-1) ', file ' num2str(iFile-1) ' complete, time: ' num2str(toc)])
        
    end %files
    
    % double check that we wrote the expected number of frames
    if length(frameInfo{iCamera})+length(droppedFrames{iCamera}) ~= currentWriteFrame-1
        warning(['Total written frame counts off for camera ' num2str(iCamera-1)])
    end
    
end %cameras
    
save(fullfile(vidDir,'droppedFrames'),'droppedFrames');

end %main function

