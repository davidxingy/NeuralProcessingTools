function [plotH, functionVars]= plotUMAPWithMarkerTracker(currentFrame,initiatePlot,functionVarsIn)
% [plotH, functionVars]= plotUMAP(currentFrame,initiatePlot,functionVarsIn)
% function to plot uMAP (or any two scatter plot) with respect to video
% recording

if initiatePlot

    load('UMAPTest.mat','reduction','behvLabels','analyzedBehaviors')
    load('VideoSyncFrames.mat')

    behvLabelsDown = behvLabels(1:5:end);

    figH = figure('Color','w','Visible','on','units','pixels','OuterPosition',[100, 50, 1000, 1000]);
    hold on;
    colormap = turbo(length(analyzedBehaviors));
    colormap(9,:) = colormap(9,:)*1.2;
    
    backgroundMarkerSize = 0.5;
    behvMarkerSize = 2;
    trailLength = 5;

    plotH(length(analyzedBehaviors)+1) = plot(reduction(behvLabelsDown==0,1),reduction(behvLabelsDown==0,2),'.','color',[0.7 0.7 0.7],'MarkerSize',backgroundMarkerSize);

    for iBehv = 1:length(analyzedBehaviors)
        plotH(iBehv) = plot(reduction(behvLabelsDown==iBehv,1),reduction(behvLabelsDown==iBehv,2),'.','color',colormap(iBehv,:),'MarkerSize',behvMarkerSize);
    end

    analyzedBehaviors = {'grooming','eating','walkgrid','walkflat','rearing','climbup','climbdown','still','jumping','jumpdown'};

    legendNames = [analyzedBehaviors 'background'];
    legendH = legend(plotH,legendNames,'Box','off','FontSize',14);
    for iLabel = 1:length(legendH.String)
        legendH.String{iLabel} = ['\color[rgb]{' num2str(plotH(iLabel).Color) '} ' legendH.String{iLabel}];
    end

    %plot current marker
    if behvLabelsDown(1) == 0
        frameMarkerColor = [0.7 0.7 0.7];
    else
        frameMarkerColor = colormap(behvLabelsDown(1),:);
    end
    functionVars{1} = plot(reduction(1,1),reduction(1,2),'o','MarkerEdgeColor','k','MarkerFaceColor',frameMarkerColor,'MarkerSize',5);

    %plot trailing tail
    functionVars{2} = plot(reduction(1:trailLength,1),reduction(1:trailLength,2),'color','k');

    axis off
    set(gcf,'Color','w')
    set(gca,'LineWidth',1.5)
    set(gca,'FontSize',12)
    set(gca,'TickDir','out')

    functionVars{3} = frameEMGSamples;
    functionVars{4} = trailLength;
    functionVars{5} = reduction;
    functionVars{6} = behvLabelsDown;

else

    iVid = 1;
    %get the current frame
    rawEMGInd = functionVarsIn{3}{1}{iVid}(currentFrame);
    uMAPEMGInd = round(rawEMGInd/20/5);
    
    functionVarsIn{1}.XData = functionVarsIn{5}(uMAPEMGInd,1);
    functionVarsIn{1}.YData = functionVarsIn{5}(uMAPEMGInd,2);

    colormap = turbo(10);
    colormap(9,:) = colormap(9,:)*1.2;

    if functionVarsIn{6}(uMAPEMGInd) == 0
        functionVarsIn{1}.MarkerFaceColor = [0.7 0.7 0.7];
    else
        functionVarsIn{1}.MarkerFaceColor = colormap(functionVarsIn{6}(uMAPEMGInd),:);
    end

    if uMAPEMGInd > functionVarsIn{4}
        functionVarsIn{2}.XData = functionVarsIn{5}(uMAPEMGInd-functionVarsIn{4}:uMAPEMGInd,1);
        functionVarsIn{2}.YData = functionVarsIn{5}(uMAPEMGInd-functionVarsIn{4}:uMAPEMGInd,2);
    end

end



% 
