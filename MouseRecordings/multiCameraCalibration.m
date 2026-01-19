% function caliParams = multiCameraCalibration(baseDir,baseFilename,lensCaliDir,squareSize)

load(fullfile(lensCaliDir,'DistortionParams.mat'))

% get calibration videos
dirFiles = dir(baseDir);
calibrationVids = dirFiles(contains(string({dirFiles.name}),baseFilename) & ...
    contains(string({dirFiles.name}),'.avi'));

% check that it matches the number of cameras from the lens calibration
% folder
assert(length(distortionParams) == length(calibrationVids),...
    'Number of detected calibration vids and number of lens calibrations not equal!');

% go through each video and the the final frame (which should be the one
% that contains the intended view of the checkerboard to use for
% calibration)
for iCam = 1:length(calibrationVids)
    
    vReader = VideoReader(calibrationVids(iCam).name);
    frame = read(vReader,inf);
    
    % undistort
    [frameUndist,distParamsPost{iCam}] = undistortImage(frame,distortionParams{iCam},OutputView="same");
    
    processedFrame = imadjust(frame,...
        stretchlim(frame,0.01),[]);
    
    [detectedPoints{iCam},detectedBoardSize{iCam}] = detectCheckerboardPoints(frameUndist,'MinCornerMetric',0.25);
    
end

% mostPoints = 
worldPoints = generateCheckerboardPoints(boardSize,squareSize);

% use matlab's built in calibration function
params = estimateMultiCameraParameters(imagePoints,worldPoints,virtualIntrinsics);


% 
