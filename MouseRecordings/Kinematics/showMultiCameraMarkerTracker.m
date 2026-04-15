function [plotH, functionVars]= showMultiCameraMarkerTracker(currentFrame,initiatePlot,functionVarsIn,functionInitParams)
% [plotH, functionVars]= showMultiCameraMarkerTracker(currentFrame,initiatePlot,functionVarsIn,argsIn)
% Function to work with my Matlab MarkerTracker to show all frames of all cameras during a multi-camera arena recording. 
% Called with the plotMarkerTrakerPhysData.m function
% 
% Inputs:
% 
% currentFrame      - Standard input for MarkerTracker synced functions,
%                     used by Markertracker to indicate the current frame
% 
% initiatePlot      - Standard input for MarkerTracker synced functions,
%                     used by Markertracker to indicate whether the call to
%                     this function is for initializing the setup
% 
% functionVarsIn    - Standard input for MarkerTracker synced functions,
%                     used by Markertracker to pass variables to this
%                     function specified by functionVars
% 
% functionInitParams - Standard input for MarkerTracker synced functions,
%                     used to pass function specific parameters to this
%                     plotting function. Parameters for this function are
%                     described below
% 

if initiatePlot

    % get plot parameters
    baseDir = functionInitParams.baseDir; %folder with video files
    cameraInds = functionInitParams.cameraInds; %for all the cameras
    primaryCamInd = functionInitParams.primaryCamInd; %Which camera is being used by the main marker tracker gui (use index start by 1 for this) 
    primaryCamFileNum = functionInitParams.primaryCamFileNum; %Which file from that camera is being used by the tracker gui (use index start by 1 for this)
    baseFileName = functionInitParams.baseFileName; %Common base for the name of all the vid files
    flipVert = functionInitParams.flipVert; %for each of the cameras, determine whether to vertically flip the image

    plotH = figure('Color','w','Visible','on','units','pixels','OuterPosition',[100, 50, 1000, 1000]);
    tiledlayout(length(cameraInds)-1,1,'Padding','tight','TileSpacing','none')

    allVidFiles = dir(baseDir);

    % Assuming that video files are .avi and that we're only using "_FrameFilled" files
    allVidFiles = allVidFiles(contains(string({allVidFiles.name}),'.avi'));
    allVidFiles = allVidFiles(contains(string({allVidFiles.name}),'_FrameFilled'));

    iSecondaryCam = 0;
    for iCamera = 1:length(cameraInds)

        %Camera file format should be: "[recording name]_Camera[camera number in 1 digit]-[file number in 4 digits]_FrameFilled.avi"
        %e.g. 'D050-120525-ArenaRecording_Camera2-0000_FrameFilled.avi'

        camFilenameInds = contains(string({allVidFiles.name}),[baseFileName '_Camera' num2str(cameraInds(iCamera)) '-' ]);
        camFileNames = string({allVidFiles(camFilenameInds).name});

        if isempty(camFileNames)
            error(['No video files found for camera ' num2str(cameraInds(iCamera))]);
        end

        %now go through each of the files corresponding to this camera
        for iFile = 1:length(camFileNames)

            % make video file reader object
            vidReaders{iCamera,iFile} = VideoReader(fullfile(baseDir,camFileNames{iFile}));
            fileNumFrames(iCamera,iFile) = vidReaders{iCamera,iFile}.NumFrames;
            
        end

        %make plot axes for secondary cameras
        if iCamera ~= primaryCamInd
            iSecondaryCam = iSecondaryCam + 1;
            nexttile
            initialFrame = vidReaders{iCamera,1}.readFrame;
            if flipVert(iCamera)
                initialFrame = fliplr(flipud(initialFrame));
            end
            axesHandles{iSecondaryCam} = imagesc(initialFrame);
            axis off
            flipSecondary(iSecondaryCam) = flipVert(iCamera);
        end

    end

    % double check that all the cameras have the same number of total frames
%     if length(unique(sum(fileNumFrames,2))) ~= 1
%         error('Cameras have different number of total frames!')
%     end

    functionVars{1} = vidReaders;
    functionVars{2} = fileNumFrames;
    functionVars{3} = axesHandles;
    functionVars{4} = primaryCamInd;
    functionVars{5} = primaryCamFileNum;
    functionVars{6} = flipSecondary;

else

    % go through each secondary camera
    vidReaders = functionVarsIn{1};
    fileNumFrames = functionVarsIn{2};
    axesHandles = functionVarsIn{3};
    primaryCamInd = functionVarsIn{4};
    primaryCamFileNum = functionVarsIn{5};
    flipSecondary = functionVarsIn{6};

    secondaryCameras = 1:length(vidReaders);
    secondaryCameras(primaryCamInd) = [];

    % see what the current frame is in context of all video files combined
    frameSums = [0 cumsum(fileNumFrames(primaryCamInd,:))];
    realCurrentFrame = frameSums(primaryCamFileNum) + currentFrame;

    for iCamera = 1:length(secondaryCameras)

        %see which file to use
        frameCutoffs = [0 cumsum(fileNumFrames(secondaryCameras(iCamera),:))]; 
        
        fileToUse = find(frameCutoffs >= realCurrentFrame,1)-1;

        fileFrame = realCurrentFrame - frameCutoffs(fileToUse);

        %read the frame
        frame = read(vidReaders{secondaryCameras(iCamera),fileToUse},fileFrame);

        if flipSecondary(iCamera)
            frame = fliplr(flipud(frame));
        end

        %plot 
        axesHandles{iCamera}.CData = frame;

    end

end



% 
