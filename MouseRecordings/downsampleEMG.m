function downsampEMG = downsampleEMG(baseDir, downsampAmount)
% downsampledEMG = downsampleEMG(baseDir, downsampAmount);
% 
% Get the full emg recording duration and downsample and save in a separate
% file.

% load and concatenate all EMG data
allEMG = [];

allProcessedFiles = string(ls(fullfile(baseDir,'processedData')));
processedEMGFiles = allProcessedFiles(contains(allProcessedFiles,'_ProcessedEMG_Block'));

% make sure to load the EMG data parts in order
for iPart = 1:length(processedEMGFiles)
    fileNameParts = split(processedEMGFiles{iPart}(1:end-4),'_');
    baseName = fileNameParts{1};
    partNum(iPart) = str2double(fileNameParts{end}(5:end));
end

for iPart = 1:length(processedEMGFiles)

    fileInd = find(partNum == iPart);
    load(fullfile(baseDir,'processedData',processedEMGFiles{fileInd}))

    allEMG = [allEMG processedEMG];

end
clear processedEMG
clear filteredEMG

% downsample
downsampEMG = decimate(allEMG',downsampAmount)';