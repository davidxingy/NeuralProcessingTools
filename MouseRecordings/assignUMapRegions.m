    function [regionAssignments, regionAssignmentsNoBound, regionAssignmentsFiltered] = assignUMapRegions(umapFile, densityGaussStd, nGridPoints, nRegions, jumpMinDuration,templateUMAPFile)
% [regionAssignments, regionAssignmentsNoBound, regionAssignmentsFiltered] = assignUMapRegions(umapFile, densityGaussStd, nGridPoints, nRegions, jumpMinDuration)
% 
% Function to cluster UMap embedding space into distinct regions using
% watershed and manual assignment of clusters. Will also do some cleaning
% up of the time-series of the assigned regions. Will also save the region
% assignments to the same file that contains the umap results.
% 
% Inputs:
% umapFile          - File which contains the UMap results (must contain
%                     'reduction' variable)
% 
% densityGaussStd   - When getting density plot, the std of the gaussian to
%                     convlve with (default 0.6)
% 
% nGridPoints       - Number of grid points for the density plot (default 1001)
% 
% nRegions          - Total number of regions/clusters to use
% 
% jumpMinDuration   - When cleaning up the time-series of the cluster
%                     labels, the duration of a jump to a different cluster
%                     which will be filtered out (in ms).
% 
% templateUMAPFile  - If the region boundaries have already been assigned
%                     by a previous session, can use the boundaries from
%                     that and will just divy up the data in the current
%                     session (so no user input required).
% 
% David Xing, last updated 2/10/2023

% set defaults
if isempty(densityGaussStd)
    densityGaussStd = 0.55;
end

if isempty(nGridPoints)
    nGridPoints = 2001;
end

if isempty(nRegions)
    nRegions = 3;
end

if isempty(jumpMinDuration)
    jumpMinDuration = 1;
end

load(umapFile,'reduction','origDownsampEMGInd')

if isempty(templateUMAPFile) || nargin < 6

    load(umapFile,'analyzedBehaviors','behvLabelsNoArt')

    origDownsampEMGInd = origDownsampEMGInd(behvLabelsNoArt~=0);

    % get density and watershed
    pointsToUse = reduction;
    pointsToUse(behvLabelsNoArt == 0,:) = [];
    behvLabelsNoBack = behvLabelsNoArt;
    behvLabelsNoBack(behvLabelsNoArt == 0) = [];

    reducLimsX1 = min(min(reduction(:,1)))*1.5-2;
    reducLimsX2 = max(max(reduction(:,1)))*1.5+2;
    reducLimsY1 = min(min(reduction(:,2)))*1.5-2;
    reducLimsY2 = max(max(reduction(:,2)))*1.5+2;

    [gridXInds,gridYInds, density] = findPointDensity(pointsToUse,densityGaussStd,nGridPoints,[reducLimsX1 reducLimsX2 reducLimsY1 reducLimsY2]);
    densityNormalized = (max(max(density))-density)/max(max(density));
    densityNormalizedProcessed = imhmin(densityNormalized,0.03);
    watershedRegions = watershed(densityNormalizedProcessed);

    figH = figure;
    axH = axes('Parent', figH);

    % plot the density
    imagesc(gridXInds,gridYInds,density)
    [boundaryYs,boundaryXs] = find(watershedRegions==0);
    hold on
    % also plot the labeled points
    colormap = turbo(length(analyzedBehaviors));
    for iBehv = 1:length(analyzedBehaviors)
        plotData = reduction(behvLabelsNoArt==iBehv,:);
        plot(plotData(1:1:end,1),plotData(1:1:end,2),'.','color',colormap(iBehv,:),'MarkerSize',0.1)
    end
    % plot watershed boundaries
    plot(gridXInds(boundaryXs),gridYInds(boundaryYs),'k.')
    set(gca,'ydir','normal')


    for iRegion = 1:nRegions

        iRegionSelection = 1;
        contGettingPoints = 1;
        disp(['Select points for Region ' num2str(iRegion)])

        regionWatershedLabels{iRegion} = [];
        while contGettingPoints

            selectRect = true;
            if ~selectRect
                selectedPoint = ginput(1);
                [~, xPoint] = min(abs(selectedPoint(1)-gridXInds));
                [~, yPoint] = min(abs(selectedPoint(2)-gridYInds));
                selectedRegionLabels = watershedRegions(yPoint,xPoint);
                regionWatershedLabels{iRegion}(iRegionSelection) = selectedRegionLabels;
            else
                selectedRect = drawrectangle();
                rectPoints = [selectedRect.Vertices(1,1), selectedRect.Vertices(3,1), ...
                    selectedRect.Vertices(1,2), selectedRect.Vertices(3,2),];
                [~, xPointStart] = min(abs(rectPoints(1)-gridXInds));
                [~, xPointStop] = min(abs(rectPoints(2)-gridXInds));
                [~, yPointStart] = min(abs(rectPoints(3)-gridYInds));
                [~, yPointStop] = min(abs(rectPoints(4)-gridYInds));
                delete(selectedRect)
                selectedRegionLabels = unique(watershedRegions(min([yPointStart yPointStop]):max([yPointStart yPointStop]),...
                    min([xPointStart xPointStop]):max([xPointStart xPointStop])));
                selectedRegionLabels(selectedRegionLabels==0) = [];
                regionWatershedLabels{iRegion} = unique([regionWatershedLabels{iRegion} selectedRegionLabels']);
            end
            iRegionSelection = iRegionSelection + 1;

            %mask out the region that was selected (so it's easier for the user
            %to see what has already been selected)

            if iRegionSelection == 2
                mask = zeros(nGridPoints,nGridPoints);
                for iSubRegion = 1:length(selectedRegionLabels)
                    mask(watershedRegions == selectedRegionLabels(iSubRegion)) = 1;
                end
                maskH = imagesc(gridXInds, gridYInds, watershedRegions/max(max(watershedRegions))/1000, 'AlphaData', mask);
            else
                for iSubRegion = 1:length(selectedRegionLabels)
                    mask(watershedRegions == selectedRegionLabels(iSubRegion)) = 1;
                end
                set(maskH,'AlphaData',mask);
            end

            morePointsInput = input('Select more points? [Y/N]: ','s');
            if strcmpi(morePointsInput,'y')
                contGettingPoints = 1;
            elseif strcmpi(morePointsInput,'n')
                contGettingPoints = 0;
            else
                continue
            end

        end

    end

    % now remap watershed labels to merge the areas that were selected to be
    % part of the same region

    largestLabel = max(max(watershedRegions));
    for iRegion = 1:nRegions

        newRegionLabel = largestLabel + iRegion;
        for iSubregion = 1:length(regionWatershedLabels{iRegion})
            watershedRegions(watershedRegions == regionWatershedLabels{iRegion}(iSubregion)) = newRegionLabel;
        end
        regionWatershedLabels{iRegion} = newRegionLabel;

    end
    regionWatershedLabels = [regionWatershedLabels{:}];

    % also remove any boundary points that were separating areas within the
    % same region
    % those boundries now are surrounded by the same label

    for iBoundry = 1:length(boundaryXs)

        boundaryAreaLabels = unique(watershedRegions(max(boundaryYs(iBoundry)-1,1):min(boundaryYs(iBoundry)+1,nGridPoints),...
            max(boundaryXs(iBoundry)-1,1):min(boundaryXs(iBoundry)+1,nGridPoints)));

        if ~any(boundaryAreaLabels==0)
            error('Something wrong when removing intra-region boundaries from watershed')
        end

        boundaryAreaLabels(boundaryAreaLabels==0) = [];

        if length(boundaryAreaLabels) == 1
            watershedRegions(boundaryYs(iBoundry),boundaryXs(iBoundry)) = boundaryAreaLabels;
        end

    end

else
    %just load the watershed region boundaries from template file
    load(templateUMAPFile,'gridXInds','gridYInds','watershedRegions','regionWatershedLabels','regionBehvAssignments')
    nRegions = length(regionWatershedLabels);
    nGridPoints = length(gridXInds);
end

% assign regions to each time point in reduction by finding the closest
% grid point to each of the points
for iPoint = 1:size(reduction,1)
    [~, minXInds(iPoint)] = min(abs(reduction(iPoint,1) - gridXInds));
    [~, minYInds(iPoint)] = min(abs(reduction(iPoint,2) - gridYInds));
end

regionSpaceinds = sub2ind([size(watershedRegions,1),size(watershedRegions,2)],minYInds,minXInds);
regionAssignments = double(watershedRegions(regionSpaceinds));

% go through each region and see if we can associate particular behaviors
% to each of them
if isempty(templateUMAPFile) || nargin < 6

    % first find the number of labeled points for each behavior within each
    % subregion
    for iReg = 1:length(regionWatershedLabels)

        %get the labels at each of the subregions
        regionInds = find(regionAssignments == regionWatershedLabels(iReg));
        regionLabels{iReg} = behvLabelsNoArt(regionInds);
        regionLabels{iReg}(regionLabels{iReg} == 0) = [];

        %save the number of labels of each behavior at each subregion
        for iBehv = 1:length(analyzedBehaviors)
            numBehvPoints(iBehv,iReg) = sum(regionLabels{iReg} == iBehv);
        end

    end

    % for each behavior, assign to region based on a threshold
    percentageInRegionThresh = 70;
    regionBehvAssignments = cell(1,length(regionWatershedLabels));
    for iBehv = 1:length(analyzedBehaviors)
        behvRegionPercentages = numBehvPoints(iBehv,:)/sum(numBehvPoints(iBehv,:))*100;
        associatedRegions = find(behvRegionPercentages > percentageInRegionThresh);

        for iReg = 1:length(associatedRegions)
            regionBehvAssignments{associatedRegions(iReg)} = [regionBehvAssignments{associatedRegions(iReg)} iBehv];
        end
    end

end

figure;
hold on;
colormap = turbo(length(regionWatershedLabels));

backgroundMarkerSize = 0.5;
behvMarkerSize = 2;

for iRegion = 1:length(regionWatershedLabels)
   plotH(iRegion) = plot(reduction(regionAssignments==regionWatershedLabels(iRegion),1),...
        reduction(regionAssignments==regionWatershedLabels(iRegion),2),'.','color',colormap(iRegion,:),'MarkerSize',behvMarkerSize);
end

% first, fix boundry points
boundryInds = find(regionAssignments==0);
regionAssignmentsNoBound = regionAssignments;
for iInd = 1:length(boundryInds)

%     %for boundry points that from and back to the same region, just assign
%     %to that region
%     prevAssignment = find(regionAssignments(1:boundryInds(iInd))~=0,1,'last');
%     nextAssignment = find(regionAssignments(boundryInds(iInd):end)~=0,1,'first')+boundryInds(iInd)-1;
%     if regionAssignments(prevAssignment) == regionAssignments(nextAssignment)
%         regionAssignmentsNoBound(boundryInds(iInd)) = regionAssignments(prevAssignment);
%     end

    %Assign to the closest region
    distX = reduction(:,1)-reduction(boundryInds(iInd),1);
    distX(boundryInds) = inf;
    distY = reduction(:,2)-reduction(boundryInds(iInd),2);
    distY(boundryInds) = inf;
    [~, closestPoint] = min(sqrt(distX.^2+distY.^2));
    closestRegion = regionAssignments(closestPoint);
    regionAssignmentsNoBound(boundryInds(iInd)) = closestRegion;

end

% next, do some "cleaning up" where if there is a very short jump to a
% different region, remove it
regionAssignmentsFiltered = regionAssignmentsNoBound;
allTransitions = find(regionAssignmentsNoBound(1:end-1)~=regionAssignmentsNoBound(2:end));

shortJumps = find(diff(allTransitions) <= jumpMinDuration);
for iJump = 1:length(shortJumps)

%     if iJump ~= length(shortJumps) && allTransitions(shortJumps(iJump+1)) == allTransitions(shortJumps(iJump)+1)
%         continue
%     end
% 
%     if iJump ~= 1 && allTransitions(shortJumps(iJump-1)) == allTransitions(shortJumps(iJump)-1)
%         continue
%     end

    if regionAssignmentsFiltered(allTransitions(shortJumps(iJump))) == regionAssignmentsFiltered(allTransitions(shortJumps(iJump)+1)+1) 
        regionAssignmentsFiltered(allTransitions(shortJumps(iJump)):allTransitions(shortJumps(iJump)+1)) = regionAssignmentsNoBound(allTransitions(shortJumps(iJump)));
    end
end

% save
save(umapFile,'regionAssignments','regionAssignmentsNoBound','regionAssignmentsFiltered','regionBehvAssignments',...
    'gridXInds','gridYInds','watershedRegions','regionWatershedLabels','-append')


% 
