function fig_h=plotWithDiffs(x,plotData,varargin)
% function fig_h=plotWithDiffs([x],plotData,[plotParameters])
% 
% function to plot data like with normal plot(), but it also plots the
% difference (along the the columns, aka the first dimension) in a subplot
% underneath. Links the two axes together.
% 
% Inputs:
% x -   the x-axis data. Can just input an empty array and will
%       automatically use the array indices, like in normal plot().
% 
% plotData -    The y-axis data. If it's a matrix, will plot each of the
%               columns separately. Differences will also be calculated
%               along the columns.
% plotParameters -  Anything that goes into plot() to specify plot type,
%                   e.g. Color, LineWidth, '.-', ect.
% 
% Outputs:
% fig_h -   Handle to the created figure with the plots

% make the axes
fig_h=figure();
ax(1)=subplot(2,1,1);
ax(2)=subplot(2,1,2);

if size(plotData,1)==1 && size(plotData,2)>1
    plotData=plotData';
end

% plot normal data
axes(ax(1))
if isempty(x)
    plot(plotData,varargin{:})
else
    plot(x,plotData,varargin{:})
end

% plot differences
axes(ax(2))

diffs=[NaN(1,size(plotData,2)); diff(double(plotData),1,1)];
if isempty(x)
    plot(diffs,varargin{:})
else
    plot(x,diffs,varargin{:})
end

% link the x axis of the two
linkaxes(ax,'x')