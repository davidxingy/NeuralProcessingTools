function [axH, clusterLabelPerm] = customDendrogram(linkageValues,varargin)
% [axH clusterLabelPerm] = customDendrogram(linkageValues,[xValues],[axH],[plotProps])
% Plot dendrogram from linkage output. Have more control over plotting than dendrogram()

% original number of nodes
nOrigNodes = size(linkageValues,1) + 1;

% get labels of clusters
allClusterLabels = linkageValues(:,1:2);
sortedClusterLabels = sort(allClusterLabels(:));

% get the default ordering of the elements based on clustering
tmpDendH = figure('Visible','off');
[~,~,clusterLabelPerm] = dendrogram(linkageValues);
close(tmpDendH);
clusterLabelPermInv = findInvPermInds(clusterLabelPerm);

if length(varargin)>=1
    if isempty(varargin{1})
        xValues = clusterLabelPermInv;
    else
        xValues = varargin{1};
    end
else
    xValues = clusterLabelPerm;
end

if length(varargin)>=2
    if isempty(varargin{2})
        figure;
        axH = axes();
    else
        axH = varargin{2};
    end
else
    figure;
    axH = axes();
end

if length(varargin)>=3
    if isempty(varargin{3})
        plotProps = {};
    else
        plotProps = varargin{3};
    end
else
    plotProps = {};
end

if length(varargin)>=3
    if isempty(varargin{3})
        plotProps = {};
    else
        plotProps = varargin{3};
    end
else
    plotProps = {};
end

assert(length(xValues) == nOrigNodes);


% get the x and y positions of the lines
for iLab = 1:length(sortedClusterLabels)

    if iLab <= nOrigNodes
        %for the original nodes
        allXValues(iLab) = xValues(iLab);
        allYValues(iLab) = 0;
    else
        %for the clusters
        allXValues(iLab) = mean([allXValues(linkageValues(iLab - nOrigNodes,1)) allXValues(linkageValues(iLab - nOrigNodes,2))]);
        allYValues(iLab) = linkageValues(iLab - nOrigNodes,3);
    end

end
allYValues(length(sortedClusterLabels)+1) = linkageValues(end,3);

% plot actual lines for the dendrogram
hold on

for iDist = 1:size(linkageValues)

    lineHs = [];
    %make horizontal line
    lineHs(1) = line([allXValues(linkageValues(iDist,1)) allXValues(linkageValues(iDist,2))],repmat(linkageValues(iDist,3),1,2));
    lineHs(2) = line(repmat(allXValues(linkageValues(iDist,1)),1,2),[allYValues(linkageValues(iDist,1)) allYValues(iDist+nOrigNodes)]);
    %make vertical lines
    lineHs(3) = line(repmat(allXValues(linkageValues(iDist,2)),1,2),[allYValues(linkageValues(iDist,2)) allYValues(iDist+nOrigNodes)]);

    %set line properties
    for iLine = 1:length(lineHs)
        for iPlotProp = 1:floor(length(plotProps)/2)
            set(lineHs(iLine),plotProps{iPlotProp*2-1},plotProps{iPlotProp*2});
        end
    end

end

set(gca,'XTick',1:nOrigNodes)
set(gca,'XTickLabel',clusterLabelPerm)




%
