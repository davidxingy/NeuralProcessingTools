function metrics=calcPerformanceMetrics(estimates, actuals, varargin)
% [metrics]=calcPerformanceMetrics(estimates, actuals, [dim], [ignoreNaN])
% 
% This function will calculate measures of decoder performance (currently
% implements R2, MSE, CC, and VAF) from the decoded and actual signals. Can
% input an array, and will calculate the metrics for each row or column, or
% a cell of arrays. If there are NaN values, can specify whether that will
% automatically return a NaN, or to calculate the metrics excluding the NaN
% values.
% 
% Inputs:
% estimates -   the decoded signal. Can be any NxM array of numbers, or, a
%               cell array of of NxM arrays (in which case, the metrics will
%               be calculated for each row/column of each cell). The
%               outputs will also be cell array of the metrics in this
%               case.
% 
% actuals -     the real signal which the estmated signal will be compared
%               against. Must be the same size as estimate.
% 
% dim -         Optional. If estimate and actual are arrays (or cells of arrays)
%               then this input specifies the dimension to calculate the metrics
%               along. Default 2 (calculate for each row).
% 
% ignoreNaN -   If this is set to true, then if either estimate or actual
%               has a NaN value as an observation, it will ignore that
%               observaion when calculating the metrics. If set to false, 
%               it just won't calculate the metrics and return NaN for the
%               metrics. Default true.
% 
% outlierLvl -  If this is not empty, then the calculation will remove
%               time points in the estimates that may be outliers (it'll
%               remove anything above the Xth percentile, where X is the
%               value inputted in outlierLvl) before calculating the
%               perforance metrics. Default no removal
% 
% weightChans - If this is not empty, then the calculation for weighted R2s
%               will only use these channels. If this input is not used,
%               then it will use all the channels. Default empty (use all chans)
% 
% Outputs:
% metrics -     A struct containing all the calulated metrics. Currently
%               contains the following fields:
% 
%   R2:         Coeeficient of determination, calculated according to the
%               equation:
%               1 - (residual sum of squares) / (total sum of squares)
% 
%   MSE:        Mean squared error, calculated as:
%               sum(square(residuals)) / npoints
% 
%   CC:         The Pearson correlation coefficient, calculated as:
%               covariance(estimate, actual) / (std(estimate) * std(actual))
% 
%   VAF:        Variance accounted for, calculated as the fraction of the
%               estimate variance over the real variance:
%               var(estimate)/var(actual)
% 
%   WeightedR2: For multi channel calculations, this finds a R2 across all
%               channels, which does a weighted mean of the R2s using the
%               varaiance of the real activity in the channel as the weight
% 
% David Xing, last updated 7/26/2024

% 10/17/2018 - Added VAF output
% 
% 7/26/2024 - Added outlier removal option
% 
% 8/5/2024 - Added Weighted R2 option

% parse inputs, set defaults
narginchk(2,6);
% dimension
if nargin>2
    if(isempty(varargin{1}))
        dim = 2;
    else
        dim = varargin{1};
    end
else
    dim = 2;
end

% whether or not to ignore nans
if nargin>3
    if(isempty(varargin{2}))
        ignoreNaN = true;
    else
        ignoreNaN= varargin{2};
    end
else
    ignoreNaN=true;
end

% ourlier removal
if nargin>4
    if(isempty(varargin{3}))
        removeOutliers = false;
        outlierLvl = 100;
    else
        removeOutliers = true;
        outlierLvl = varargin{3};
    end
else
    removeOutliers = false;
    outlierLvl = 100;
end

% weighted average channels
if nargin>5
    if(isempty(varargin{4}))
        weightChans = [];
    else
        weightChans = varargin{4};
    end
else
    weightChans = [];
end

% cell to hold all the outputs since it's easier to iterate across, will convert to struct later
% 1  - R2
% 2  - MSE
% 3  - CC
% 4  - VAF
% 5  - Weighted R2
metricsCell={{},{},{},{},{}};

% change to cell array if not already one
isCell=true;
if ~strcmpi(class(estimates),'cell')
    estimates={estimates};
    actuals={actuals};
    isCell=false;
end
    
% now calculate for each cell
for iCell=1:length(estimates)
    
    %initialize with zeros
    if dim==1
        nSamples=size(actuals{iCell},1);
        nTrials=size(actuals{iCell},2);
        
        for iMetric=1:length(metricsCell)
            metricsCell{iMetric}{iCell}=zeros(1,nTrials);
        end
        
    elseif dim==2
        nSamples=size(actuals{iCell},2);
        nTrials=size(actuals{iCell},1);
        
        for iMetric=1:length(metricsCell)
            metricsCell{iMetric}{iCell}=zeros(nTrials,1);
        end
    end
    
    if isempty(weightChans)
        weightChans = 1:nTrials;
    end
        
    %do for each row/column
    for iTrial=1:nTrials
        if dim==1
            trialEst=estimates{iCell}(:,iTrial);
            trialReal=actuals{iCell}(:,iTrial);
        elseif dim==2
            trialEst=estimates{iCell}(iTrial,:);
            trialReal=actuals{iCell}(iTrial,:);
        end
        
         %remove nans if we're told to ignore nans
         nanValues=union(find(isnan(trialEst)), find(isnan(trialReal)));
         if ignoreNaN
             trialEst(nanValues)=[];
             trialReal(nanValues)=[];
         elseif ~empty(nanValues)
             %just output NaN for the whole trial
             for iMetric=1:length(metricsCell)
                 metricsCell{iMetric}{iCell}(iTrial)=NaN;
             end
             continue;
         end

         %remove outliers if desired
         if removeOutliers
             outlierThresh = prctile(abs(trialEst),outlierLvl);
             outlierInds = find(abs(trialEst) >= outlierThresh);
             trialEst(outlierInds) = [];
             trialReal(outlierInds) = [];
         end
         
         %R2:
         R2 = 1-sum((trialEst-trialReal).^2)/sum((trialReal-mean(trialReal)).^2);
         metricsCell{1}{iCell}(iTrial)=R2;
         
         %MSE:
         MSE = sum((trialEst-trialReal).^2)/length(trialEst);
         metricsCell{2}{iCell}(iTrial)=MSE;
         
         %CC:
         covarMat = cov(trialEst,trialReal);
         CC = covarMat(1,2)/(sqrt(covarMat(1,1))*sqrt(covarMat(2,2)));
         metricsCell{3}{iCell}(iTrial)=CC;
         
         %VAF:
         VAF = var(trialEst)/var(trialReal);
         metricsCell{4}{iCell}(iTrial)=VAF;
         
         %Get the variance in the real for the weighted R2 metric:
         chanVar(iTrial) = var(trialReal);
         
    end
    
    %calculate weighted R2
    weightedR2 = chanVar(weightChans)*metricsCell{1}{iCell}(weightChans)/sum(chanVar(weightChans));
    metricsCell{5}{iCell} = weightedR2;
    
end

% Convert to output struct. Also, change back to array from cell if the variables were given as arrays
if (~isCell)
    metrics.R2=metricsCell{1}{1};
    metrics.MSE=metricsCell{2}{1};
    metrics.CC=metricsCell{3}{1};
    metrics.VAF=metricsCell{4}{1};
    metrics.R2Weighted=metricsCell{5}{1};
else
    metrics.R2=metricsCell{1};
    metrics.MSE=metricsCell{2};
    metrics.CC=metricsCell{3};
    metrics.VAF=metricsCell{4};
    metrics.R2Weighted = metricsCell{5};
end


% 
