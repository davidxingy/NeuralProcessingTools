clear
close all


sessionNames = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
    'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
    'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

for iSess = 1:length(sessionNames)
   
    load(fullfile(sessionNames{iSess},'ProcessedData','EpochedData10msLick.mat'));

    %compare EMG autocorrelation
    for iMuscle = 1:8
        [acfFlat(:,iMuscle), lags] = autocorr(behavioralData.walkflat.allBoutEMGs(iMuscle,:),'NumLags',50);

        for iBout = 1:length(behavioralData.walkgrid.boutEMGs)
            [acfGrid(:,iMuscle), lags] = autocorr(behavioralData.walkgrid.allBoutEMGs(iMuscle,:),'NumLags',50);
        end
    end

    figure;
    shadedErrorBar(lags*10,mean(acfFlat'),std(acfFlat')/sqrt(size(acfFlat,2)),'lineprops','r')
    hold on
    shadedErrorBar(lags*10,mean(acfGrid'),std(acfGrid')/sqrt(size(acfGrid,2)),'lineprops','b')
    xlabel('Lag (ms)')
    ylabel('Autocorrelation')
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    set(gca,'tickdir','out')
    set(gcf,'color','w')
    legend('Walk Flat','Walk Grid','box','off')

    %compare average activity of neurons
    meanFRsFlat = nanmean(behavioralData.walkflat.allBoutFRs,2);
    meanFRsGrid = nanmean(behavioralData.walkgrid.allBoutFRs,2);

    %fit linear regression
    linFit = fitlm(meanFRsFlat, meanFRsGrid);
    fitR2(iSess) = linFit.Rsquared.Ordinary;
    fitSlope(iSess) = linFit.Coefficients.Estimate(2);
    fitIntercept(iSess) = linFit.Coefficients.Estimate(1);
    cc(iSess) = corr(meanFRsFlat, meanFRsGrid);

    figure;
    plot(meanFRsFlat*100,meanFRsGrid*100,'.','MarkerSize',10)
    xlims = get(gca,'XLim');
    ylims = get(gca,'YLim');
    box off
    line(xlims, xlims*fitSlope(iSess)+fitIntercept(iSess),'linestyle','--','color','r')
    text(5,ylims(2)*0.9,{['R^2 = ' num2str(fitR2(iSess))],['Slope = ' num2str(fitSlope(iSess))]},'fontsize',14)
    ylim(ylims)
    xlabel('Walk flat ave FRs (spks/sec)')
    ylabel('Walk grid ave FRs (spks/sec)')
    set(gca,'linewidth',2)
    set(gca,'fontsize',14)
    set(gca,'tickdir','out')
    set(gcf,'color','w')

    %next look at subspace alignment
    goodNeurons = find(meanFRsFlat*100 > 0.2 | meanFRsGrid*100 > 0.2);

    flatNanInds = any(isnan(behavioralData.walkflat.allBoutFRs(goodNeurons,:)));
    flatPCAInput = behavioralData.walkflat.allBoutFRs(goodNeurons,~flatNanInds)';
    flatPCAInput = flatPCAInput - mean(flatPCAInput);
    [pcaWeightsFlat,pcaTrajsFlat] = pca(flatPCAInput);

    gridNanInds = any(isnan(behavioralData.walkgrid.allBoutFRs(goodNeurons,:)));
    gridPCAInput = behavioralData.walkgrid.allBoutFRs(goodNeurons,~gridNanInds)';
    gridPCAInput = gridPCAInput - mean(gridPCAInput);
    [pcaWeightsGrid,pcaTrajsGrid] = pca(gridPCAInput);

    pcaTrajsFlatCross = flatPCAInput*pcaWeightsGrid;
    pcaTrajsGridCross = gridPCAInput*pcaWeightsFlat;

    figure;
    tiledlayout(2,1)
    nexttile
    plot(var(pcaTrajsFlat),'.-','MarkerSize',10)
    hold on
    plot(var(pcaTrajsFlatCross),'.-','MarkerSize',10)
    ylabel('Variance')
    title('Flat Wallking PCA Var')
    legend('Projected into own subspace','Projected into walk grid subspace','box','off')
    box off
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',13)
    set(gca,'tickdir','out')

    nexttile
    plot(var(pcaTrajsGrid),'.-','MarkerSize',10)
    hold on
    plot(var(pcaTrajsGridCross),'.-','MarkerSize',10)
    xlabel('PC')
    ylabel('Variance')
    title('Grid Wallking PCA Var')
    legend('Projected into own subspace','Projected into walk flat subspace','box','off')
    box off

    set(gcf,'color','w')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',13)
    set(gca,'tickdir','out')

end



% 
