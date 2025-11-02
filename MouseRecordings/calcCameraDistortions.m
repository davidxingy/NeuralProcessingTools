function distortionParams = calcCameraDistortions(baseDir,baseFilename,squareSize)
% distortionParams = calcCameraDistortions(baseDir,baseFilename,squareSize)
% Uses matlab function estimateCameraParameters() to get the distortion
% (not specifically fisheye distortion, but wide angle distortion)
% parameters of the camera lenses based on a calibration session with
% a checkerboard, with the size of the board squares indicated by squareSize

nCameras = 4;

for iCam = 1:nCameras

    calVid = VideoReader(fullfile(baseDir,[baseFilename '_Camera' num2str(iCam-1) '-0000.avi']));
    
    % read frames
    allFrames = uint8([]);
    vidNFrames = calVid.NumFrames;

    % only use every 4 frames to save time and memory
    usedFrames = 1:4:vidNFrames;

    tic
    
    for iFrame = 1:length(usedFrames)

        frame = read(calVid,usedFrames(iFrame));
        allFrames(:,:,1,iFrame) = uint8(frame(:,:,1));

    end
    
    disp(['Loaded in video data from Camera ' num2str(iCam) ', time: ' num2str(toc)])

    tic
    
    [imagePoints,boardSize,usedFrames] = detectCheckerboardPoints(allFrames);%,'HighDistortion',false,'MinCornerMetric',0.05);
    worldPoints = generateCheckerboardPoints(boardSize,squareSize);
    
    disp(['Computed checkerboard points from Camera ' num2str(iCam) ', time: ' num2str(toc)])
    
    tic
    
    imageSize = [size(allFrames,1) size(allFrames,2)];
    
    %Using more points takes a lot more time (increases nonlinearly), so
    %only use around 100 frames.
    usedPoints = randperm(size(imagePoints,3),100);
    distortionParams{iCam} = estimateCameraParameters(imagePoints(:,:,usedPoints),worldPoints,'ImageSize',imageSize);

    disp(['Computed distortion parameters from Camera ' num2str(iCam) ', time: ' num2str(toc)])
    
end


end %of function


% 
