function dataNames = getMouseDataNames(mouseID,baseSessionName,probeRegion)
% dataNames = getMouseDataNames(mouseID,baseSessionName,probeRegion)
% Function to get the paths and file names (e.g. raw NP file, or video
% files) for the specified animal and session name. Since the animal
% recordings spanned many years where the type of implants and naming
% conventions have been changing, this function serves as a map to get the
% different processed data files associated with each animal

validMice = {'D020','D024','D026','D036','D040','D041','D043','D047','D050','D054','D056'};

% check the mouse ID is valid
assert(any(strcmp(validMice,mouseID)),'mouseID input does not match any of the valid mice!')

% D020 and D024 are in the Miri local server (Z drive)
% Rest of the animals are on resfiles (X drive)
switch mouseID
    case {'D020', 'D024'}
        baseServerFolder = 'Z:\David\ArenaRecordings\NeuropixelsTest';
    otherwise
        baseServerFolder = 'X:\David\ArenaRecordings';
end

% D020, D024, and D026 are NP 1.0's, the rest are 2.0's
switch mouseID
    case {'D020', 'D024', 'D026'}
        npProbeType = 1;
    otherwise
        npProbeType = 2;
end

% D020, D024, D026, and D043 have only one probe implanted, the rest have 2
switch mouseID
    case {'D020', 'D024', 'D026', 'D043'}
        nImplantedProbes = 1;
    otherwise
        nImplantedProbes = 2;
end

% check that if it is a two probe animal, the probe region input is
% specified, and it can only be 'CFA' or 'tjMC'
if nImplantedProbes == 2
    assert(~isempty(probeRegion),'For this animal, must specify which probe region to get!');
    assert(any(strcmpi(string({'CFA','tjMC'}),probeRegion)),'For this animal, must specify which probe region to get!');
end

% for the two probe animals, currently the first probe (_imec0) corresponds
% to CFA/striatum probe, and the second probe (_imec1) corresponds to the
% tjMC probe, but I'm still specifying it here in case future recordings
% have a different correspondence
if nImplantedProbes == 2
    switch probeRegion
        case 'CFA'
            probeID = 0;
        case 'tjMC'
            probeID = 1;
    end
else
    probeID = nan;
end


% add higher level folders first
dataNames.videoFolder = fullfile(baseServerFolder,baseSessionName,'Video');
dataNames.processedDataFolder = fullfile(baseServerFolder,baseSessionName,'ProcessedData');
dataNames.Intan = fullfile(baseServerFolder,baseSessionName,'Intan');
if nImplantedProbes == 1
    dataNames.npDataFolder = fullfile(baseServerFolder,baseSessionName,'Neuropixels');
else
    dataNames.npDataFolder = fullfile(baseServerFolder,baseSessionName,'Neuropixels',[baseSessionName '_g0_imec' num2str(probeID)]);
end


% add neuron data structures and firing rates
if nImplantedProbes == 1
    appendSuffix = '.mat';
elseif nImplantedProbes == 2
    appendSuffix = ['_' probeRegion '.mat'];
end
dataNames.neuronDataStruct = fullfile(dataNames.processedDataFolder,['neuronDataStruct' appendSuffix]);

% the current possible firing rate files that I could've generated:
firingRateFileNames = {'NeuralFiringRates1msBins0msGauss','NeuralFiringRates1msBins10msGauss',...
    'NeuralFiringRates5msBins30msGauss','NeuralFiringRates5msBins0msGauss',...
    'NeuralFiringRates10msBins30msGauss','NeuralFiringRates10msBins0msGauss',...
    'NeuralFiringRates100msBins50msGauss','NeuralFiringRates100msBins0msGauss'};

for iFRType = 1:length(firingRateFileNames)
    dataNames.(firingRateFileNames{iFRType}) = fullfile(dataNames.processedDataFolder,[firingRateFileNames{iFRType} appendSuffix]);
end


% add Downsampled EMG data
dataNames.EMG1ms = fullfile(dataNames.processedDataFolder,'EMG1ms.mat');


% add EMG meta data file
dataNames.emgMetaData = fullfile(dataNames.processedDataFolder,[baseSessionName '_ProcessedEMG_MetaData.mat']);


% add Histology coordinates file
dataNames.histologyCoords = fullfile(dataNames.processedDataFolder,['histologyCoords' appendSuffix]);


% add neural-EMG sync script file
dataNames.NeurEMGSync = fullfile(dataNames.processedDataFolder,'NeurEMGSync.m');


% add UMAP file
dataNames.UMAPFile = fullfile(dataNames.processedDataFolder,'UMAP.mat');


% add sync indices file
dataNames.VideoSyncFrames = fullfile(dataNames.processedDataFolder,'VideoSyncFrames.mat');


% get UMAP activity overlay analysis output files
if nImplantedProbes == 1
    appendSuffix = '';
elseif nImplantedProbes == 2
    appendSuffix = ['_' probeRegion];
end
dataNames.NeuronRegionProps = fullfile(dataNames.processedDataFolder,['UMAPFRs' appendSuffix],'NeuronRegionProps.mat');
dataNames.MuscleRegionProps = fullfile(dataNames.processedDataFolder,['UMAPFRs' appendSuffix],'MuscleRegionProps.mat');


% add manual behavior annotation data
dataNames.BehaviorAnnotationLabels = fullfile(dataNames.processedDataFolder,'BehaviorAnnotations','BehaviorLabels.mat');


% add NP data file (post-local median sub, post-artifact removal). Some of
% the earlier animals had weird naming conventions
switch mouseID
    case {'D020'}
        dataNames.ArtRemovedNPFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0_RemovedArtifact.imec0.ap.bin']);
    case {'D024','D026'}
        dataNames.ArtRemovedNPFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0_RemovedArtifact_LocalMedianSubtr.imec0.ap.bin']);
    otherwise
        if nImplantedProbes == 1
            dataNames.ArtRemovedNPFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0_LocalMedianSubtr_OptoArtRemoved.imec0.ap.bin']);
        else
            dataNames.ArtRemovedNPFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0_LocalMedianSubtr_OptoArtRemoved.imec' num2str(probeID) '.ap.bin']);
        end
end

% add the meta file for it too
if nImplantedProbes == 1
    dataNames.npMetaFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0.imec0.ap.meta']);
else
    dataNames.ArtRemovedNPFile = fullfile(dataNames.npDataFolder,[baseSessionName '_g0_t0.imec' num2str(probeID) '.ap.meta']);
end
% ****There's some special cases where the meta file isn't valid because
% the recording was terminated prematurely. These include D043-013125-ArenaRecording


% get NP artifact timestamps
dataNames.npArtifactTimestamps = fullfile(dataNames.npDataFolder,'artifactTimestamps.mat');

% also, for some sessions, there were issues with baseline drift hitting
% the amplifier rails on the NP recordings, these files have a separate
% "flatlineDetections" file containing those artfiacts
dataNames.flatlineDetections = fullfile(dataNames.npDataFolder,'flatlineDetections.mat');

% add EMG-activity based single behavior classification analysis file
dataNames.EMGSingleBehvClassifiers = fullfile(dataNames.processedDataFolder,'EMGSingleBehvClassifiers.mat');


% and the behavior initations/termination triggered activity based on the
% classificaiton analysis
dataNames.EMGClassifierInits = fullfile(dataNames.processedDataFolder,'EMGClassifierInits.mat');


% add data file for EMG-activation triggered trials (+ associated controls)
if nImplantedProbes == 1
    appendSuffix = '';
elseif nImplantedProbes == 2
    appendSuffix = ['_' probeRegion];
end
dataNames.EMGTrigInitiations = fullfile(dataNames.processedDataFolder,['EMGTrigInitiations' appendSuffix '.mat']);

% add data file for kinematics from Lightning Pose
dataNames.litPoseKinematics = fullfile(dataNames.processedDataFolder,'litPoseKinematics.mat');

% add data file for Arena spatial positions and coordinate transforms
dataNames.arenaSpatialCoords = fullfile(dataNames.processedDataFolder,'arenaSpatialCoords.mat');

end % of function



% 


