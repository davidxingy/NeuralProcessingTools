clear;

baseDIR = 'X:\David\ArenaRecordings\DT019-121924-OptoTest\BottomMostShanks';
binFile = 'DT019-121924-OptoTest_BottomMostShanks_g0_t0_LocalMedianSubtr_OptoArtRemoved.imec0.ap.bin';
metaFile = 'DT019-121924-OptoTest_BottomMostShanks_g0_t0.imec0.ap.meta';

meta = ReadMeta(metaFile, baseDIR);
nChans = str2double(meta.nSavedChans);
totSamps = round(str2double(meta.fileTimeSecs) * str2double(meta.imSampRate));

% nChans = 385;
% totSamps = 69736615;

sampleRate = 30000;

readFid = fopen(fullfile(baseDIR, binFile), 'rb');

load(fullfile(baseDIR,'syncSignal.mat'))

checkRunNumber = 1;
runNTrials = 50;

optoThresh = 40;
optoOns = find(syncSignal(2:end) > optoThresh & syncSignal(1:end-1) <= optoThresh);
optoOffs = find(syncSignal(2:end) < optoThresh & syncSignal(1:end-1) >= optoThresh); 
pulseNums = (checkRunNumber-1)*runNTrials+1:checkRunNumber*runNTrials;

for iPulse = 1:length(pulseNums)
    
    %load in 50 ms before and 100 ms after each opto pulse
    nPreSamples = 30*50;
    nPostSamples = 30*100;
    startReadSample = optoOns(pulseNums(iPulse)) - nPreSamples;
    fseek(readFid, startReadSample * 2 * nChans, 'bof');
    rawData = fread(readFid, [nChans, nPreSamples+nPostSamples], 'int16=>double');
    
    pulseSigs(:,:,iPulse) = rawData;
    
end

pulseStd = movstd(pulseSigs,100,0,2);

plotChanRadius = 10;
plotChanCenter = 100;

figure;
plot(squeeze(pulseSigs(plotChanCenter,:,:)) + (1:runNTrials)*20)
hold on
plot(pulseSigs(end,:,1),'k','LineWidth',2)

figure;
plot(squeeze(pulseSigs(plotChanCenter-plotChanRadius:plotChanCenter+plotChanRadius,:,10))' + (1:plotChanRadius*2+1)*20)
hold on
plot(pulseSigs(end,:,1),'k','LineWidth',2)


figure;
plot(mean(mean(pulseStd,1),3))
hold on
plot(pulseSigs(end,:,1)>0,'k','LineWidth',2)

figure;
plot(squeeze(pulseStd(plotChanCenter-plotChanRadius:plotChanCenter+plotChanRadius,:,10))' + (1:plotChanRadius*2+1)*20)
hold on
plot(pulseSigs(end,:,1),'k','LineWidth',2)


% 
