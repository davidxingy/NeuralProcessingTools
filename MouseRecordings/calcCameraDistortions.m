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
    
    [imagePoints,boardSize,usedFrames] = detectCheckerboardPoints(allFrames,'MinCornerMetric',0.5,'PartialDetections',false);%,'HighDistortion',false,'MinCornerMetric',0.05);
    worldPoints = generateCheckerboardPoints(boardSize,squareSize);
    
    % sometimes the checkboard detection isn't great, so do some filtering
    % to remove those instances
    usedFramesInds = find(usedFrames);
    boardWidth = squeeze(max(imagePoints(:,1,:),[],1)-min(imagePoints(:,1,:),[],1));
    boardHeight = squeeze(max(imagePoints(:,2,:),[],1)-min(imagePoints(:,2,:),[],1));
    boardRatio = boardWidth./boardHeight;
    minHeight = 100;
    minWidth = 100;
    minRatio = 0.6;
    maxRatio = 3;
    
    badPoints = boardWidth<minWidth | boardHeight<minHeight | boardRatio<minRatio | boardRatio>maxRatio;
    usedFrames(usedFramesInds(badPoints)) = 0;
    usedFramesInds(badPoints) = [];
    imagePoints(:,:,badPoints) = [];
    
    
    disp(['Computed checkerboard points from Camera ' num2str(iCam) ', time: ' num2str(toc)])
    disp(['Used ' num2str(length(find(usedFrames))) ' frames'])
    
    tic
    
    imageSize = [size(allFrames,1) size(allFrames,2)];
    
    %Using more points takes a lot more time (increases nonlinearly), so
    %only use around 1000 frames.
    usedPoints = randperm(size(imagePoints,3),min(size(imagePoints,3),1000));
    distortionParams{iCam} = estimateCameraParameters(imagePoints(:,:,usedPoints),worldPoints,'ImageSize',imageSize);

    disp(['Computed distortion parameters from Camera ' num2str(iCam) ', time: ' num2str(toc)])
    disp(['Reprojection error: ' num2str(distortionParams{iCam}.MeanReprojectionError)])
    
end


end %of function


% 
