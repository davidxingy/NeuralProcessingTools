function distortionParams = calcCameraDistortions(baseDir,baseFilename,squareSize)
%CALCCAMERADISTORTIONS Summary of this function goes here
%   Detailed explanation goes here

nCameras = 4;

for iCam = 1:nCameras

    calVid = VideoReader(fullfile(baseDir,[baseFilename '_Camera' num2str(iCam) '-0000.avi']));
    
    % read frames
    allFrames = uint8([]);
    vidNFrames = calVid.NumFrames;

    % only use every 4 frames to save time and memory
    usedFrames = 1:4:vidNFrames;

    for iFrame = 1:length(usedFrames)

        frame = read(calVid,usedFrames(iFrame));
        allFrames(:,:,1,iFrame) = uint8(frame(:,:,1));

    end

    [imagePoints,boardSize,usedFrames] = detectCheckerboardPoints(allFrames);%,'HighDistortion',false,'MinCornerMetric',0.05);
    worldPoints = generateCheckerboardPoints(boardSize,squareSize);

    imageSize = [size(allFrames,1) size(allFrames,2)];
    distortionParams{iCam} = estimateFisheyeParameters(imagePoints,worldPoints,imageSize);

end


end %of function


% 
