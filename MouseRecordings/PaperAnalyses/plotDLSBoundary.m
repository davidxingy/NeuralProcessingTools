clear

striatalDepthCutoffs = [2600 2400 2400];
eatBoundaries = [1080 900 860]; %from the tip

sessionDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

% first load in allen atlas
% directory of reference atlas files
annotation_volume_location = 'X:\David\HistologyImages\AllanReferenceMap\annotation_volume_10um_by_index.npy';
structure_tree_location = 'X:\David\HistologyImages\AllanReferenceMap\structure_tree_safe_2017.csv';
template_volume_location = 'X:\David\HistologyImages\AllanReferenceMap\template_volume_10um.npy';

av = readNPY(annotation_volume_location);
st = loadStructureTree(structure_tree_location);
tv = readNPY(template_volume_location);

ccfBregma = [540 0 570]; %[AP, DV, ML]

% plot brain
figure
sliceImage = squeeze(tv(ccfBregma(1),:,:));

% the atlas plots background as black, change to white
leftBackgroundBound = 230;
rightBackgroundBound = 900;
topBackgroundBound = 150;
botBackgroundBound = 650;
backgroundThresh = 40;
mask = sliceImage;
mask(topBackgroundBound:botBackgroundBound,leftBackgroundBound:rightBackgroundBound) = backgroundThresh*2;
sliceImage(mask<backgroundThresh) = max(max(sliceImage));
imagesc(sliceImage)
colormap(gray)
hold on
axis off
set(gcf,'color','w')

plotColors = lines(7);

% plot the probe tract and eating boundary location for each animal
for iAnimal = 1:length(sessionDirs)
    load(fullfile(sessionDirs{iAnimal},'ProcessedData','histologyCoords'))

    % align white matter tract
    wmtLabel = 1199;
    wmtTableInd = find(borders_table.avIndex==1199);
    wmtBorder = mean([borders_table.upperBorder(wmtTableInd) borders_table.lowerBorder(wmtTableInd)]);

    % plot probe
    plot3([insertionSite(3)+probeTraj(3)*-50 insertionSite(3)+probeTraj(3)*400*shrinkage_factor], [insertionSite(2)+probeTraj(2)*-50 insertionSite(2)+probeTraj(2)*400*shrinkage_factor],...
        [insertionSite(1)+probeTraj(1)*-50 insertionSite(1)+probeTraj(1)*400*shrinkage_factor],'color',[1 1 1 0.5],'LineWidth',1.5)

    % get boundary
    wmtEphys = (4000-striatalDepthCutoffs(iAnimal))*shrinkage_factor;
    offset(iAnimal) = wmtBorder - wmtEphys;

    eatBoundaryDepth = (4000 - eatBoundaries(iAnimal))*shrinkage_factor;% + offset;

    boundaryX = [insertionSite(3)+probeTraj(3)*(eatBoundaryDepth*0.1) - 25, insertionSite(3)+probeTraj(3)*(eatBoundaryDepth*0.1) + 25];
    boundaryY = repmat(insertionSite(2)+probeTraj(2)*(eatBoundaryDepth*0.1),1,2);
    boundaryZ = repmat(insertionSite(1)+probeTraj(1)*(eatBoundaryDepth*0.1),1,2);

    plot3(boundaryX,boundaryY,boundaryZ,':','color',plotColors(iAnimal,:),'LineWidth',2)

end

% scale bar
line([100 200], [800 800],'color','k','linewidth',2)

% 

