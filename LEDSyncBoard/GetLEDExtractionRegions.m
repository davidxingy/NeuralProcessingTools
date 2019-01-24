function regions = GetLEDExtractionRegions(image,nRegions,xlims,ylims)
% regions = GetLEDExtractionRegions(image,nRegions,xlims,ylims)
% 
% A function which acts as a crude GUI to get a cell array of rectangles 
% ([xmin ymin width height]) for each of the regions of interest (ROIs)
% used in ExtractLEDs(). User selects each region from a figure.
% 
% Outputs:
%   regions - a 1 x nRegions cell array containing the extract regions
%   
% Inputs:
%   image - the n x m x 1 image to show the user to get the regions from
% 
%   nRegions - the number of regions to be extracted
% 
%   xlims - a 1 x 2 array containing the lower and upper x limits of the
%           inputted image to zoom in on (in pixels).
% 
%   ylims - a 1 x 2 array containing the lower and upper y limits of the
%           inputted image to zoom in on (in pixels).
% 
% David Xing
% Last Updated 3/2/2018

% thickness to display the rectangles
boxThickness = 1;

% initialize memory
regions=cell(1,nRegions);

% crop to the sub region of the image that we want to display
subImage=image(ylims(1):ylims(2),xlims(1):xlims(2));
% add an "undo" region at the top of the image
imageWithUndo = [zeros(15, xlims(2)-xlims(1)+1)+max(max(subImage)); subImage];

% display the image
imshow(imageWithUndo,[min(min(subImage)), max(max(subImage))]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
text(round(size(subImage,2)/3),5,'UNDO')

imageRects=image;

%use getrect to let user select all regions
iROI=1;
while iROI<=nRegions
    
    rect=round(getrect());
    
    if (rect(2)<=15)
        if (iROI==1)
            continue;
        end
        iROI = iROI-1;
        imageRects=RemoveRect(imageRects, image, regions{iROI}, boxThickness);
    else
        
        %change the x and y cornor to be in reference to the whole image, not
        %the cropped and undo-added image
        rect(1)= rect(1)+xlims(1)-1;
        rect(2)= rect(2)+ylims(1)-1-15;
        
        %save to regions variable
        regions{iROI}=rect;
        
        %draw the selected rectangle on the frame so use can see what has
        %already been defined (use thickness of 2)
        imageRects=AddRect(imageRects, rect, boxThickness);

        iROI = iROI + 1;
    end
    
    %Display the updated image
    subImage=imageRects(ylims(1):ylims(2),xlims(1):xlims(2));
    imageWithUndo = [zeros(15, xlims(2)-xlims(1)+1)+max(max(subImage)); subImage];
    imshow(imageWithUndo,[min(min(imageRects)), max(max(imageRects))]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    text(round(size(subImage,2)/3),5,'UNDO')
    
            
    %display numbers of all the previously defined regions
    for iPrevROI=1:iROI-1
        text(regions{iPrevROI}(1)+regions{iPrevROI}(3)/2-xlims(1),...
            regions{iPrevROI}(2)+regions{iPrevROI}(4)/2-ylims(1)+15,...
            num2str(iPrevROI),'color',[1 1 1]);
    end
    
end


function outFrame=AddRect(inFrame, rect, thichness)
% add a white rectangle to image

outFrame = inFrame;

outFrame(rect(2):rect(2)+rect(4),rect(1):rect(1)+thichness-1,:)=255;
outFrame(rect(2):rect(2)+thichness-1,rect(1):rect(1)+rect(3),:)=255;
outFrame(rect(2):rect(2)+rect(4),(rect(1)+rect(3)-thichness+1):(rect(1)+rect(3)),:)=255;
outFrame((rect(2)+rect(4)-thichness+1):(rect(2)+rect(4)),rect(1):rect(1)+rect(3),:)=255;

function outFrame=RemoveRect(inFrame, originalFrame, rect, thichness)
% remove the white rectangle from image

outFrame = inFrame;

outFrame(rect(2):rect(2)+rect(4),rect(1):rect(1)+thichness-1,:)=...
    originalFrame(rect(2):rect(2)+rect(4),rect(1):rect(1)+thichness-1,:);
outFrame(rect(2):rect(2)+thichness-1,rect(1):rect(1)+rect(3),:)=...
    originalFrame(rect(2):rect(2)+thichness-1,rect(1):rect(1)+rect(3),:);
outFrame(rect(2):rect(2)+rect(4),(rect(1)+rect(3)-thichness+1):(rect(1)+rect(3)),:)=...
    originalFrame(rect(2):rect(2)+rect(4),(rect(1)+rect(3)-thichness+1):(rect(1)+rect(3)),:);
outFrame((rect(2)+rect(4)-thichness+1):(rect(2)+rect(4)),rect(1):rect(1)+rect(3),:)=...
    originalFrame((rect(2)+rect(4)-thichness+1):(rect(2)+rect(4)),rect(1):rect(1)+rect(3),:);



%

