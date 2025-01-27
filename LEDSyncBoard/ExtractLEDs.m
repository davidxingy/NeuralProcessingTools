function [ledStatus, intensities]=extractLEDs(filename, nRegions, thresh, colorChannel, saveOutputVideo, regions)
% [ledStatus, intensities, regions, usedThreshs]=extractLED(filename, nRegions, thresh, colorChannel, saveOutputVideo, prevDefRegions)
%
% Function which extracts from a video containing LEDs the frames where
% each LED was turned on or off.
%
%
% Outputs:
% ledStatus - NxM array of bools with N video frames by M LEDs/ROIs indicating
% whether each LED was on or not for each frame
%
% intensities - NxM array of doubles which gives the average intensity of
% each ROI for each frame
%
% regions - cell array of ROIs which correspond to each LED. Each cell
% contains the x,y coordinates, and width/height (in pixels) of the rectangle
% defining the ROI.
%
% usedThreshs - the thresholds used for each ROI (useful if user
% got thresholds automatically)
%
%
% Inputs:
% filename - string containing the file (include the extension) to be used (can be a
% .mat file or a video file)
%
% nRegions - number of LEDs/ROIs to extract
%
% thresh - the ROI average intensity threshold to count the LED as on. Can
% be one value to use for all ROIs, or can be vector of values for each ROI.
% Can use 'auto' for the function to automatically calculate the best
% threshold to use for each ROI.
%
% colorChannel - which color channel to use for the intensity calculation.
% Can be 'red', 'blue', 'green', or 'grey'
%
% saveOutputVideo - bool indicating if you want to save a video with the
% frames drawing region boxes whenever the LED status was detected to be
% on. Useful for verification but will take alot longer.
%
% regions - A cell array of rectangle coordinates defining the regions of interest
% (ROIs) which will be used to get the LEDs. Used the function
% GetLEDExtractionRegions() for a quick way to get these regions
% 
% David Xing
% Last Updated 3/2/2018

% if its a mat file, then just load in the file
isVideo = true;
if(strcmp(filename(end-3:end),'.mat'))
    frames=load(filename);
    varName=fieldnames(frames);
    nFrames=size(frames.(varName{1}),3);
    isVideo = false;
    frameRate = 30; %default set to 30FPS
    minValue = min(min(min(frames.(varName{1}))));
    maxValue = max(max(max(frames.(varName{1}))));
else
    % read in video
    try
        vidReader = VideoReader(filename);
    catch
        error(['Cannot open ' filename '!'])
    end
    nFrames = round(vidReader.Duration*vidReader.FrameRate);
    frameRate = vidReader.FrameRate;
end

% initialize output
ledStatus=zeros(nFrames,nRegions);

if (isVideo)
    firstFrame = readFrame(vidReader);
else
    firstFrame = squeeze(frames.(varName{1})(:,:,1));
end

% thickness (in pixels) to use for drawing region boxes
boxThickness=1;

% If the user did not input ROIs, ask user to define ROIs from first frame
%Otherwise, just display the first frame with the predefined ROI inputs
if (length(regions)~=nRegions)
    error('regions must be of length nRegions!');
end

firstFrameRects=firstFrame;
for iROI=1:nRegions
    rect=round(regions{iROI});
    
    %draw rectangles around the regions
    firstFrameRects=AddRect(firstFrameRects,rect, boxThickness);
end

%draw frame
imshow(firstFrameRects,[min(min(min(firstFrameRects))), min(max(max(firstFrameRects)))]);

%display numbers
for iROI=1:nRegions
    rect=round(regions{iROI});
    text(rect(1)+rect(3)/2,rect(2)+rect(4)/2, num2str(iROI),'color',[1 1 1]);
end


% Now that we have all our regions, go through all frames and get the
% average intensity in each region
intensities.mean=zeros(nFrames, nRegions);
intensities.std=zeros(nFrames, nRegions);

% first frame
if (isVideo)
    firstFrameColor=GetFrameColorChannel(firstFrame,colorChannel);
else
    firstFrameColor = firstFrame;
end

for iRegion=1:nRegions
    ROI=firstFrameColor(regions{iRegion}(2):regions{iRegion}(2)+regions{iRegion}(4),...
        regions{iRegion}(1):regions{iRegion}(1)+regions{iRegion}(3));
    intensities.mean(1,iRegion)=mean(ROI(:));
    intensities.std(1,iRegion)=std(ROI(:));
end

clear firstFrame firstFrameColor;

% all remaining frames
iFrame=2;
while true
    if (isVideo)
        if(~hasFrame(vidReader))
            break
        end
        frame=readFrame(vidReader);
        frameColor=GetFrameColorChannel(frame,colorChannel);
    else
        if(iFrame > nFrames)
            break
        end
        frameColor = squeeze(frames.(varName{1})(:,:,iFrame));
    end
    
    
    for iRegion=1:nRegions
        ROI=frameColor(regions{iRegion}(2):regions{iRegion}(2)+regions{iRegion}(4),...
            regions{iRegion}(1):regions{iRegion}(1)+regions{iRegion}(3));
        intensities.mean(iFrame,iRegion)=mean(ROI(:));
        intensities.std(iFrame,iRegion)=std(ROI(:));
    end
    
%     disp(iFrame)
    iFrame = iFrame+1;
end

% if user wants to auto-calculate thresholds, look at the distribution of
% average intensities in each region and define a threshold based on
% mixture of gaussians model
if (strcmpi(class(thresh),'char') && strcmpi(thresh,'auto'))
    
    for iROI=1:nRegions
        %fit GM model
        inputs = [intensities.mean(:,iROI), intensities.std(:,iROI)];
        try
            GMModel = fitgmdist(inputs,2);
        catch
            %mixture model failed, most likely due to there being only 1
            %distrbution (always on or always off), use average threshhold
            %from other ROIs in that case
            warning(['Region ' num2str(iROI) ' failed to find GMModel, most likely LED is always off or always on, '...
                'will use average threshold of other ROIs as threshold value for this region.']);
            thresholds(iROI)=NaN;
            continue;
        end
        
        %save model
        thresholds{iROI}=GMModel;
        
        %cluster the frames by the model
        inputsClusters = cluster(GMModel,inputs)-1;
        
        %cluster 1 should always be the one with higher std (LED on).
        if(mean(inputs(inputsClusters==0,2)) > mean(inputs(inputsClusters==1,2)))
            inputsClusters = ~inputsClusters;
        end
        
        %finally save to classification output variable
        ledStatus(:,iROI) = inputsClusters;

    end
    
else
    
    %use defined threshold
    if (length(thresh)==1)
        thresholds=repmat(thresh,1,nRegions);
    elseif (length(thresh)~=nRegions)
        error('Number of thresholds must be equal to number of ROIs or be equal to 1!');
    else
        thresholds=thresh;  
    end
    
end

% Now threshold for all intensities
ledStatus = intensities.mean > repmat(thresholds,nFrames,1);
usedThreshs = thresholds;

% Finally, write output video file if user requested
if (saveOutputVideo)
    
    %name output video
    vidWriter=VideoWriter([filename(1:end-4) '_processed'],'MPEG-4');
    vidWriter.FrameRate=frameRate;
    open(vidWriter);
    
    %load each frame of the video (again) and add boxes if LED is detected
    %to be on, then save the frame
    if (isVideo)
        vidReader.CurrentTime=0;
    end
    
    iFrame=1;
    while (iFrame<=nFrames)
        
        if(isVideo)
            frame=readFrame(vidReader);
        else
            frame = squeeze(frames.(varName{1})(:,:,iFrame));
        end
        
        for iROI=1:nRegions
            
            if (ledStatus(iFrame,iROI))
                %draw rectangles around the regions
                rect=round(regions{iROI});
                frame=AddRect(frame, rect, boxThickness);
            end
            
        end
        
        writeVideo(vidWriter, uint8(frame));
        
        disp(iFrame);
        iFrame = iFrame+1;
        
    end
    
    close(vidWriter)
    
end

function outFrame=GetFrameColorChannel(inFrame,channel)
% from rgb image, get the specified channel

switch channel
    case 'red'
        outFrame=double(squeeze(inFrame(:,:,1)));
    case 'blue'
        outFrame=double(squeeze(inFrame(:,:,3)));
    case 'green'
        outFrame=double(squeeze(inFrame(:,:,2)));
    case {'gray', 'grey'}
        outFrame=double(rgb2gray(inFrame));
end



function outFrame=AddRect(inFrame, rect, thichness)
% add a white rectangle to image

outFrame = inFrame;

outFrame(rect(2):rect(2)+rect(4),rect(1):rect(1)+thichness-1,:)=255;
outFrame(rect(2):rect(2)+thichness-1,rect(1):rect(1)+rect(3),:)=255;
outFrame(rect(2):rect(2)+rect(4),(rect(1)+rect(3)-thichness+1):(rect(1)+rect(3)),:)=255;
outFrame((rect(2)+rect(4)-thichness+1):(rect(2)+rect(4)),rect(1):rect(1)+rect(3),:)=255;


%
