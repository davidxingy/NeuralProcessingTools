function rasterplot(spikes,inputtype,plottype,ax,times,yPos,varargin)
% fighandle=rasterplot(spikes,inputtype,plottype,ax)
% 
% Creates a raster plot from the spikes of multiple neurons. For inputtype
% 'trains', spikes is a NxT matrix with N neurons and T time bins, with
% each time point being a 0 or 1 indicating whether there was a spike or
% not.  For inputtype 'times', spikes is a cell array with N cells, with
% each cell containing the timing of the spikes for each neuron.
% 
% To use lines to mark each spike, use '|' as the plottype input. To
% use dots, input the normal matlab string to indicate the type of marker
% to use (eg '*' for stars, 'o' for circles) as plottype.
% 
% David Xing 6/10/2015

if nargin==3
    ax=gca;
    times=1:size(spikes,2);
elseif nargin==4
    times=1:size(spikes,2);
elseif ~strcmpi(class(ax),'matlab.graphics.axis.Axes')
    varargin=[ax varargin];
    ax=gca;
end

if ~exist('yPos')
    yPos = 1:length(spikes);
end

if nargin == 7

    plotProps = varargin{1};

else
    
    plotProps = [];

end


axes(ax);
switch lower(inputtype)
    case 'trains'
        numneurons=size(spikes,1);
        
        for whichneuron=1:numneurons
            hold on
            
            spiketimes=times(find(spikes(whichneuron,:)==1));
            numspikes=length(spiketimes);
            
            switch plottype
                case '|'
                    line(repmat(spiketimes,2,1),repmat([yPos(whichneuron)-0.4; yPos(whichneuron)+0.4],1,numspikes),'color','k','linewidth',1.5)
                otherwise
                    plot(spiketimes,repmat(yPos(whichneuron),1,numspikes),[plottype 'k'],'markersize',7)
            end
            
        end
        hold off
        ylim([0 numneurons+1])
        
    case 'times'
        numneurons=length(spikes);
        for whichneuron=1:numneurons
            numspikes=length(spikes{whichneuron});
            hold on
            
            switch plottype
                case '|'
                    line([spikes{whichneuron}';spikes{whichneuron}'],...
                        repmat([yPos(whichneuron)-0.4; yPos(whichneuron)+0.4],1,numspikes),'color','k','linewidth',3);
                otherwise
                    plotH = scatter(spikes{whichneuron},repmat(yPos(whichneuron),1,numspikes),[plottype 'k'],'sizedata', 20);
                    for iPlotProp = 1:length(plotProps)/2

                        set(plotH,plotProps{iPlotProp*2-1},plotProps{iPlotProp*2})

                    end

            end
            
        end
        hold off
        ylim([min(yPos)-1 max(yPos)+1])
end


% 
