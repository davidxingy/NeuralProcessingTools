function [barH] = barScatterPlot(data,varargin)
% plotH = barScatterPlot(data,plotErrors,plotDots,jitter,drawLines,noBars)
% 
% Inputs:
% data          - MxN cell array containing the data to plot, where M is
%                 the number of catagories (i.e. the number of bar groups)
%                 and N is the number of series( i.e. the number of bars in
%                 each bar group). Each cell contains a vector containing
%                 the individual samples for each bar to plot
% 
% plotErrors    - String, indicates what to use for the error bars. Options
%                 are 'std', for standard deviation, or 'sem' for standard
%                 error of the mean, or 'none' for no error bars. Default 
%                 'sem'.
% 
% plotDots      - Boolean array. Same size as the data cell (MxN), which for
%                 each set of samples just indicates whether to plot the
%                 individual samples as a scatter plot overlaid on the bar
%                 plots or not. Default false.
% 
% jitter        - Same size as data, For plotting individual samples,
%                 specify the amount to offset each sample by (in terms of
%                 fraction relative to the bar width). Default 0
% 
% drawLines     - For drawing individual lines for each sample across series, 
%                 for each bar group. A Nx2 array, specifying the starting 
%                 and ending indices of the series to connect with lines.
%                 Default empty array (no lines drawn).
% 
% noBars        - If instead of drawing a whole error bar rectangle, just
%                 want to draw a line instead. Default false.
% 

% parse inputs
if length(varargin)>=1
    if isempty(varargin{1})
        plotErrors = 'sem';
    else
        plotErrors = varargin{1};
    end
else
    plotErrors = 'sem';
end

if length(varargin)>=2
    if isempty(varargin{2})
        plotDots = zeros(size(data));
    else
        plotDots = varargin{2};
    end
else
    plotDots = zeros(size(data));
end

if length(varargin)>=3
    if isempty(varargin{3})
        jitter = cellfun(@(x) zeros(length(x),1), data, 'un', 0);
    else
        jitter = varargin{3};
    end
else
    jitter = cellfun(@(x) zeros(length(x),1), data, 'un', 0);
end

if length(varargin)>=4
    if isempty(varargin{4})
        drawLines = [];
    else
        drawLines = varargin{4};
    end
else
    drawLines = [];
end

if length(varargin)>=5
    if isempty(varargin{5})
        noBars = false;
    else
        noBars = varargin{5};
    end
else
    noBars = false;
end


plotColors = lines(size(data,2));
scatterSize = 4;

% Dim 1 is the number of bar groups, dim 2 is the number of bars in each group
nCategories = size(data,1);
nSeries = size(data,2);
nSamples = cellfun(@(x) length(x), data);

plotXBase = 1 : (nSeries+1) : (nCategories*(nSeries+1));

% get stats
dataMeans = cellfun(@mean, data);
dataStds = cellfun(@std, data);
dataSems = dataStds./sqrt(nSamples);

% make figure
hold on;

for iSeries = 1:nSeries

    % plot bars or lines as the means
    if ~noBars
        barH(iSeries) = bar(plotXBase+iSeries-1,dataMeans(:,iSeries),0.33333,'EdgeColor','none','FaceColor',plotColors(iSeries,:)*1.05,'FaceAlpha',0.5);
    else
        line([plotXBase+iSeries-1.5; plotXBase+iSeries-0.5],repmat(dataMeans(:,iSeries)',2,1),'linewidth',3,'color',plotColors(iSeries,:))
    end

    %make error bars
    switch lower(plotErrors)
        case 'sem'
            errorbar(plotXBase+iSeries-1,dataMeans(:,iSeries),dataSems(:,iSeries),'LineStyle','none','LineWidth',1.5,'Color','k','CapSize',8)
        case 'std'
            errorbar(plotXBase+iSeries-1,dataMeans(:,iSeries),dataStds(:,iSeries),'LineStyle','none','LineWidth',1.5,'Color','k','CapSize',8)
        case 'none'

    end

    %Make scatter plot if requested
    for iCategory = 1:nCategories
        if plotDots(iCategory,iSeries)
            plot(plotXBase(iCategory)+iSeries-1+jitter{iCategory,iSeries},data{iCategory,iSeries},'o','MarkerFaceColor',plotColors(iSeries,:)*0.9,'MarkerEdgeColor','none','MarkerSize',scatterSize)
        end
    end

end

%draw lines between categories within each group
for iLine = 1:size(drawLines,1)
    for iCategory = 1:nCategories
        line(repmat([plotXBase(iCategory)+drawLines(iLine,1)-1; plotXBase(iCategory)+drawLines(iLine,2)-1],1,nSamples(iCategory,1)),...
            [data{iCategory,drawLines(iLine,1)}; data{iCategory,drawLines(iLine,2)}],'linewidth',2,'color',[0.5 0.5 0.5])
    end
end


% make pretty
box off
set(gca,'XColor','k')
set(gca,'YColor','k')
set(gca,'FontSize',14)
set(gca,'linewidth',2)
set(gca,'TickDir','out')
set(gca,'XTick',plotXBase+mean(1:nSeries)-1)
set(gcf,'color','w')


% 
