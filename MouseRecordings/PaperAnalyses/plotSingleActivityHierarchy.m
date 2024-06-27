clear

allDirs = {'Z:\David\ArenaRecordings\NeuropixelsTest\D020-062922-ArenaRecording',...
            'Z:\David\ArenaRecordings\NeuropixelsTest\D024-111022-ArenaRecording',...
            'X:\David\ArenaRecordings\D026-032923-ArenaRecording'};

allBehvAlignPerms = [
    1 2 3 4 5 6 7; ...
    1 2 3 4 5 6 7; ...
    1 2 4 5 3 6 7 ...
    ];

allAnimalLabels = {'D020','D024','D026'};
inputData = 'humanannotated';

plotColors = lines(7);
nShifts = 100;

for iAnimal = 1:length(allDirs)

    baseDir = allDirs{iAnimal};
    load(fullfile(baseDir,'ProcessedData','NeuralFiringRates10msBins30msGauss.mat'),'cortexFRs','striatumFRs','cortexInds','striatumInds')

    maxCtxFRs = max(cortexFRs,[],2)*100;
    maxStrFRs = max(striatumFRs,[],2)*100;

    meanCtxFRs = nanmean(cortexFRs,2)*100;
    stdCtxFRs = nanstd(cortexFRs,[],2)*100;

    meanStrFRs = nanmean(striatumFRs,2)*100;
    stdStrFRs = nanstd(striatumFRs,[],2)*100;

    switch lower(inputData)

        case 'umapregions'

            behvRegionLabels = {'Climb Up','Climb Down','Misc/Jump','Walk','Misc/Rear/Still','Groom','Eat'};
            load(fullfile(baseDir,'ProcessedData','UMAPFRs','NeuronRegionProps.mat'))
            load(fullfile(baseDir,'ProcessedData','UMAP.mat'), 'regionBehvAssignments','analyzedBehaviors')
            regionAveSigs = regionAveSigs(:,allBehvAlignPerms(iAnimal,:))/10;
            regionAveSigsShuff = regionAveSigsShuff(:,allBehvAlignPerms(iAnimal,:),1:nShifts,:)/10;
            
            exampleAnimal = 3;
            exampleBaseBehv = 1;
            exampleCrossBehv1 = 2;
            exampleCrossBehv2 = 7;

        case 'humanannotated'

            behvRegionFieldNames = {'climbup','climbdown','jumpdown','jumping','walkflat','walkgrid','rearing','still','grooming','eating'};
            behvRegionLabels = {'Climb Up','Climb Down','Jump Down','Jump Across','Walk Flat','Walk Grid','Rear','Still','Groom','Eat'};

            load(fullfile(baseDir,'ProcessedData','EpochedData10ms.mat'))
            sessStructNames = fieldnames(behavioralData);

            %get average firing rate across behaviors
            regionAveSigs = [];
            regionAveSigsShuff = [];
            for iBehv = 1:length(behvRegionFieldNames)

                behvStructInd(iBehv) = find(strcmpi(sessStructNames,behvRegionFieldNames{iBehv}));
                regionAveSigs(:,iBehv) = nanmean(behavioralData.(sessStructNames{behvStructInd(iBehv)}).allBoutFRs,2)*100;

                %also get information for doing the shift controls
                behvBoutFRs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutFRs;
                behvBoutEMGs{iBehv} = behavioralData.(behvRegionFieldNames{iBehv}).boutEMGs;

                %get the time of each bout, to put in order for getting the shift controls
                behvBoutStartInds{iBehv} = cellfun(@(x) x(1),allNeurInds(behvStructInd(iBehv),~cellfun(@isempty,allNeurInds(behvStructInd(iBehv),:))));
                behvBoutLabels{iBehv} = cellfun(@(x) repmat(iBehv,1,size(x,2)), behvBoutFRs{iBehv}, 'un', 0);

            end

            % do shift controls to break behavior relationship with neur/emg
            % first combine all behaviors (sorting bouts by time)
            [~, sortPerm] = sort([behvBoutStartInds{:}]);
            allBoutFRs = [behvBoutFRs{:}];
            allBoutFRs = allBoutFRs(sortPerm);
            allBoutEMGs = [behvBoutEMGs{:}];
            allBoutEMGs = allBoutEMGs(sortPerm);
            allBoutLabels = [behvBoutLabels{:}];
            allBoutLabels = allBoutLabels(sortPerm);

            allBoutFRsCat = cat(2,allBoutFRs{:});
            allBoutEMGsCat = cat(2,allBoutEMGs{:});
            allBoutLabelsCat = cat(2,allBoutLabels{:});

            % don't use points with nans
            catFRNans = find(any(isnan(allBoutFRsCat)));
            catEMGNans = find(any(isnan(allBoutEMGsCat)));

            allBoutFRsCat(:,unique([catFRNans catEMGNans])) = [];
            allBoutEMGsCat(:,unique([catFRNans catEMGNans])) = [];
            allBoutLabelsCat(:,unique([catFRNans catEMGNans])) = [];

            % do shift then divide into separate behaviors again
            for iShift = 1:nShifts

                %shift by at least 60 seconds
                shiftAmount(iShift) = randi(length(allBoutLabelsCat)-6000*2)+6000;
                labelsShift = circshift(allBoutLabelsCat,shiftAmount(iShift));

                % get FRs and EMGs based on new shifted labels
                for iBehv = 1:length(behvRegionFieldNames)
                    shiftBehvInds = find(labelsShift == iBehv);
                    regionFRsShift{iBehv,iShift} = allBoutFRsCat(:,shiftBehvInds);
                    regionEMGsShift{iBehv,iShift} = allBoutEMGsCat(:,shiftBehvInds);
                    regionAveSigsShuff(:,iBehv,iShift) = nanmean(regionFRsShift{iBehv,iShift},2)*100;
                end

            end

            exampleAnimal = 3;
            exampleBaseBehv = 1;
            exampleCrossBehv1 = 2;
            exampleCrossBehv2 = 10;

    end


    for iRegion1 = 1:size(regionAveSigs,2)
        for iRegion2 = 1:size(regionAveSigs,2)

            % get both R2 of linear fit as well as pearson correlation
            region1CtxFRs = regionAveSigs(length(striatumInds)+1:end,iRegion1)./stdCtxFRs;
            region2CtxFRs = regionAveSigs(length(striatumInds)+1:end,iRegion2)./stdCtxFRs;

            nanFRNeurons = find(isnan(region1CtxFRs) | isnan(region2CtxFRs));
            region1CtxFRs(nanFRNeurons) = [];
            region2CtxFRs(nanFRNeurons) = [];

            linReg = fitlm(region1CtxFRs,region2CtxFRs);
            r2FitCtx(iRegion1,iRegion2,iAnimal) = linReg.Rsquared.Ordinary;
            corrFitCtx(iRegion1,iRegion2,iAnimal) = corr(region1CtxFRs,region2CtxFRs);


            region1StrFRs = regionAveSigs(1:length(striatumInds),iRegion1)./stdStrFRs;
            region2StrFRs = regionAveSigs(1:length(striatumInds),iRegion2)./stdStrFRs;

            nanFRNeurons = find(isnan(region1StrFRs) | isnan(region2StrFRs));
            region1StrFRs(nanFRNeurons) = [];
            region2StrFRs(nanFRNeurons) = [];

            linReg = fitlm(region1StrFRs,region2StrFRs);
            r2FitStr(iRegion1,iRegion2,iAnimal) = linReg.Rsquared.Ordinary;
            corrFitStr(iRegion1,iRegion2,iAnimal) = corr(region1StrFRs,region2StrFRs);

            if iRegion1 == iRegion2
                r2FitCtx(iRegion1,iRegion2,iAnimal) = 1;
                r2FitStr(iRegion1,iRegion2,iAnimal) = 1;
                corrFitCtx(iRegion1,iRegion2,iAnimal) = 1;
                corrFitStr(iRegion1,iRegion2,iAnimal) = 1;
            end

            % next do it for two controls: randomly permuting the two
            % classes, as well as randomly shifting behavior labels
            % ******TODO: permute control (first need to get individual
            % FRs for each region)
            
            % next do behavior shifts controls
            for iShift = 1:nShifts

                if iRegion1 == iRegion2
                    r2FitCtxShift(iRegion1,iRegion2,iAnimal,iShift) = 1;
                    r2FitStrShift(iRegion1,iRegion2,iAnimal,iShift) = 1;
                    corrFitCtxShift(iRegion1,iRegion2,iAnimal,iShift) = 1;
                    corrFitStrShift(iRegion1,iRegion2,iAnimal,iShift) = 1;

                    continue
                end

                % get both R2 of linear fit as well as pearson correlation
                region1CtxFRs = regionAveSigsShuff(length(striatumInds)+1:end,iRegion1,iShift)./stdCtxFRs;
                region2CtxFRs = regionAveSigsShuff(length(striatumInds)+1:end,iRegion2,iShift)./stdCtxFRs;

                nanFRNeurons = find(isnan(region1CtxFRs) | isnan(region2CtxFRs));
                region1CtxFRs(nanFRNeurons) = [];
                region2CtxFRs(nanFRNeurons) = [];

                linReg = fitlm(region1CtxFRs,region2CtxFRs);
                r2FitCtxShift(iRegion1,iRegion2,iAnimal,iShift) = linReg.Rsquared.Ordinary;
                corrFitCtxShift(iRegion1,iRegion2,iAnimal,iShift) = corr(region1CtxFRs,region2CtxFRs);


                region1StrFRs = regionAveSigsShuff(1:length(striatumInds),iRegion1,iShift)./stdStrFRs;
                region2StrFRs = regionAveSigsShuff(1:length(striatumInds),iRegion2,iShift)./stdStrFRs;

                nanFRNeurons = find(isnan(region1StrFRs) | isnan(region2StrFRs));
                region1StrFRs(nanFRNeurons) = [];
                region2StrFRs(nanFRNeurons) = [];

                linReg = fitlm(region1StrFRs,region2StrFRs);
                r2FitStrShift(iRegion1,iRegion2,iAnimal,iShift) = linReg.Rsquared.Ordinary;
                corrFitStrShift(iRegion1,iRegion2,iAnimal,iShift) = corr(region1StrFRs,region2StrFRs);

            end
            

            %make example plot
            if iAnimal == exampleAnimal && iRegion1 == exampleBaseBehv && (iRegion2 == exampleCrossBehv1 || iRegion2 == exampleCrossBehv2)

                if iRegion2 == exampleCrossBehv1
                    exampleCorrsH = figure;
                    tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
                    nexttile
                    lineColor = [1 0 0];
                elseif iRegion2 == exampleCrossBehv2
                    figure(exampleCorrsH)
                    nexttile
                    lineColor = [0.7 0 0.8];
                end

                plot(regionAveSigs(length(striatumInds)+1:end,iRegion1)./stdCtxFRs,...
                    regionAveSigs(length(striatumInds)+1:end,iRegion2)./stdCtxFRs,...
                    'o','MarkerSize',4,'MarkerEdgeColor','none','MarkerFaceColor',lines(1));

                xAxRange = get(gca,'XLim');
                xAxRange = [-0.05 2];
                yAxRange = [-0.05 2];

                hold on;
                line(xAxRange,xAxRange*linReg.Coefficients.Estimate(2)+linReg.Coefficients.Estimate(1),...
                    'linestyle','--','linewidth',1.5,'color',lineColor)

                box off
                set(gca,'FontSize',14)
                set(gca,'linewidth',1.5)
                set(gca,'TickDir','out')
                set(gcf,'Color','w')
                set(gca,'XColor',[0 0 0])
                set(gca,'YColor',[0 0 0])
                xlabel([behvRegionLabels{iRegion1} ' normalized FR'])
                ylabel([behvRegionLabels{iRegion2} ' normalized FR'])
                xlim(xAxRange)
                ylim(yAxRange)
                
                text(0.1,1.9,['Corr = ' num2str(corrFitCtx(iRegion1,iRegion2,iAnimal))],'color',lineColor,'FontSize',14)

            end

        end
    end

    % get hierarchy measurement
    hierarchyMeanCtx(iAnimal) = mean(squareform(1-corrFitCtx(:,:,iAnimal)));
    hierarchyMeanStr(iAnimal) = mean(squareform(1-corrFitStr(:,:,iAnimal)));

    hierarchyStdCtx(iAnimal) = std(squareform(1-corrFitCtx(:,:,iAnimal)));
    hierarchyStdStr(iAnimal) = std(squareform(1-corrFitStr(:,:,iAnimal)));

    % get for shifts too
    for iShift = 1:nShifts
        hierarchyMeanCtxShift(iAnimal,iShift) = mean(squareform(1-corrFitCtxShift(:,:,iAnimal,iShift)));
        hierarchyMeanStrShift(iAnimal,iShift) = mean(squareform(1-corrFitStrShift(:,:,iAnimal)));

        hierarchyStdCtxShift(iAnimal,iShift) = std(squareform(1-corrFitCtxShift(:,:,iAnimal,iShift)));
        hierarchyStdStrShift(iAnimal,iShift) = std(squareform(1-corrFitStrShift(:,:,iAnimal,iShift)));
    end

    % plot example hiearchy
    if iAnimal == exampleAnimal

        % first do cortex
        figure('Units','normalized','OuterPosition',[0.1 0.1 0.445 0.7])
        tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
        nexttile
        plotTitle = 'Cortex';
        imagesc(corrFitCtx(:,:,iAnimal))
        set(gca,'XTick',1:length(behvRegionLabels))
        set(gca,'YTick',1:length(behvRegionLabels))
        set(gca,'XTickLabelRotation',30)
        set(gca,'XTickLabel',behvRegionLabels)
        set(gca,'YTickLabel',behvRegionLabels)
        set(gcf,'Color','w')
        set(gca,'LineWidth',1.5)
        set(gca,'FontSize',13)
        set(gca,'TickLength',[0 0])
        cH = colorbar;
        caxis([0.4 1])
        cH.Label.String = 'Correlation';
        title(plotTitle)

        % next do striatum
        nexttile
        plotTitle = 'Striatum';
        imagesc(corrFitStr(:,:,iAnimal))
        set(gca,'XTick',1:length(behvRegionLabels))
        set(gca,'YTick',1:length(behvRegionLabels))
        set(gca,'XTickLabelRotation',30)
        set(gca,'XTickLabel',behvRegionLabels)
        set(gca,'YTickLabel',behvRegionLabels)
        set(gcf,'Color','w')
        set(gca,'LineWidth',1.5)
        set(gca,'FontSize',13)
        set(gca,'TickLength',[0 0])
        cH = colorbar;
        caxis([0.4 1])
        cH.Label.String = 'Correlation';
        title(plotTitle)

        % plot dendrograms
        % first cortex
        nexttile
        plotMatrix = 1-corrFitCtx(:,:,iAnimal);
        [~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});
        % add an example shift control dendrogram
        plotMatrix = 1-corrFitCtxShift(:,:,iAnimal,1);
        customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(3,:) 0.7]})

        set(gca,'TickDir','out')
        set(gca,'LineWidth',2)
        set(gca,'XColor',[0 0 0])
        set(gca,'YColor',[0 0 0])
        set(gca,'FontSize',14)
        ylabel('1 - Correlation')
        set(gca,'XTickLabel',behvRegionLabels(dendPerm))
        set(gca,'XTickLabelRotation',30)
        set(gca,'TickDir','out')
        ylim([0 0.4])
        axH = gca;

        % next striatum
        nexttile
        plotMatrix = 1-corrFitStr(:,:,iAnimal);
        [~, dendPerm] = customDendrogram(linkage(squareform(plotMatrix)),[],gca,{'Linewidth',2,'Color','k'});
        % add an example shift control dendrogram
        plotMatrix = 1-corrFitStrShift(:,:,iAnimal,12);
        customDendrogram(linkage(squareform(plotMatrix)),findInvPermInds(dendPerm)+0.15,gca,{'Linewidth',2,'color',[plotColors(3,:) 0.7]})

        set(gca,'TickDir','out')
        set(gca,'LineWidth',2)
        set(gca,'XColor',[0 0 0])
        set(gca,'YColor',[0 0 0])
        set(gca,'FontSize',14)
        ylabel('1 - Correlation')
        set(gca,'XTickLabel',behvRegionLabels(dendPerm))
        set(gca,'XTickLabelRotation',30)
        set(gca,'TickDir','out')
        ylim([0 0.4])
        axH = gca;

    end

end
    
% now make the summary plots
figure;

plotJitter = [-0.2 0 0.2];
barScatterPlot({hierarchyStdCtx;hierarchyStdStr},'none',[1; 1],repmat({plotJitter},2,1),[],true);
hold on

% also plot controls
for iAnimal = 1:length(allDirs)

    controlJitter = randn(1,nShifts)/100;
    
    % first do cortex
    plotH = scatter(controlJitter+plotJitter(iAnimal)+1,hierarchyStdCtxShift(iAnimal,:),10,plotColors(3,:),'filled');
    alpha(plotH,0.2)
    % then striatum
    plotH = scatter(controlJitter+plotJitter(iAnimal)+3,hierarchyStdStrShift(iAnimal,:),10,plotColors(3,:),'filled');
    alpha(plotH,0.2)

end

xlim([0 4])
set(gca,'XTickLabel',{'Cortex','Striatum'})
ylabel('Hierarchy spread (s.d.)')



% 
