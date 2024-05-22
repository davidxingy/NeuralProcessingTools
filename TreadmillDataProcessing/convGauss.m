function smoothedFRs = convGauss(binnedSpiketimes, binSize, convStd, causal)
% smoothedFRs = convGauss(binnedSpiketimes, binSize, convStd)
% 
% function to smooth binned spike counts to a firing rate by convolving
% with a gaussian. The size of the gaussian kernal will be from around -3
% standard devations to +3 standard deviations.
% 
% Inputs:
% binnedSpiketimes -    MxN matrix, containing the spike counts. Each
%                       row will be considered a channel, and the smoothing
%                       will be applied to each channel.
% 
% binSize -             Number indicating how large each bin is (in ms).
% 
% convStd -             How large a standard deviation for the gaussian
%                       kernal will be (in ms).
% 
% Outputs:
% smoothedFRs -         MxN matrix of the smoothed spike counts
% 
% David Xing 3/9/2020


% go from around -3 stds to 3 stds
gaussX = round(-3*convStd):binSize:round(3*convStd);
gaussY=normpdf(gaussX,0,convStd)*binSize;

if causal
    gaussY(1:floor(length(gaussY)/2)) = 0;
end

for iChan = 1:size(binnedSpiketimes,1)
    smoothedFRs(iChan,:)=conv(binnedSpiketimes(iChan,:), gaussY, 'same');
end