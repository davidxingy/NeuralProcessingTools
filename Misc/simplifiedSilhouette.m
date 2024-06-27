function clustScores = simplifiedSilhouette(data,labels)
% s = simplifiedSilhouette(data,labels)
% 
% function to calculate simplified silhouette scores were the distances to
% cluster centroids are used rather than to every point, which makes it
% much faster to compute

uniqueLabels = unique(labels);

%calculate centroids (using median instead of means)
for iLab = 1:length(uniqueLabels)
    cents(iLab,:) = median(data(labels == uniqueLabels(iLab),:));
end

% calculate the coefficients
for iLab = 1:length(uniqueLabels)

    intraClustDists = sqrt(sum((data(labels == uniqueLabels(iLab),:) - cents(iLab,:)).^2,2))';

    allOtherClustDists = [];
    for iLab2 = 1:length(uniqueLabels)
        if iLab ~= iLab2
            allOtherClustDists(iLab2,:) = sqrt(sum((data(labels == uniqueLabels(iLab),:) - cents(iLab2,:)).^2,2));
        else
            allOtherClustDists(iLab2,:) = zeros(1,sum(labels == uniqueLabels(iLab)));
        end
    end
    interClustDists = max(allOtherClustDists);

    clustScores{iLab} = (interClustDists - intraClustDists)./max([interClustDists; intraClustDists]);

end