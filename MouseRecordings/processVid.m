clear; close all;

vidDir = 'C:\Users\Arena\Videos\';
vidNames = {'D002-081221-EMGTest-08122021144913-000'};
vidNSegs = [1];

emgDir = 'C:\Users\Arena\Desktop\D002\WakingUp_210812_142741\';
load([emgDir 'ProcessedEMG'])
load([emgDir 'VidLedSyncs'])

channelNames = {'PL','ECR','Bi','Tri'};
channelPos = [75 175 275 375];
channelColors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560];
envScale = [6 6 8 3];
channelHeight = 100;

% match syncs
ledSyncInds = find(ledIntensity{2}(2:end)>200 & ledIntensity{2}(1:end-1)<200)+1;
samplesPerFrame = 20000/30;
syncInds = [syncInds{:}];
iTrain = 1;
syncIndIntan{1} = [];
syncIndVid{1} = [];
for iPulse = 1:length(syncInds)-1

    syncIndIntan{iTrain} = [syncIndIntan{iTrain} syncInds(iPulse)];
%     syncIndVid{iTrain} = [syncIndVid{iTrain} ledSyncInds(iPulse)];
    
    if syncInds(iPulse+1)-syncInds(iPulse) > 30000
        iTrain = iTrain+1;
        syncIndIntan{iTrain} = [];
%         syncIndVid{iTrain} = [];
    end
    
end

% just use first pulse for now
frame2SampleParams = [ledSyncInds(1) syncIndIntan{2}(1)];

for iVid = 1:length(vidNames)
    for iSeg = 1:vidNSegs(iVid)
        
        vIn = VideoReader([vidDir vidNames{iVid} num2str(iSeg-1) '.avi']);
        vOut = VideoWriter([vidDir vidNames{iVid} num2str(iSeg-1) '_EMG.avi']);
        open(vOut);

        frameInd = 1;
        while hasFrame(vIn)
            
            frame = readFrame(vIn);
            %get corresponding emg sample number
            sampleInd = round((frameInd - frame2SampleParams(1))*samplesPerFrame + frame2SampleParams(2));
            
            %plot 2 second of data
            emgSig = filteredEMG(:,(sampleInd-40000):sampleInd)*-4;
            emgEnv = processedEMG(:,(sampleInd-40000):sampleInd)*-4;
            
            %plot channels
            h = figure('Visible','off');
            imshow(frame);
            hold on;
            
            for iChan = 1:4
                text(100,channelPos(iChan),channelNames{iChan},'FontSize',20,'Color','w');
                plot(linspace(300,1900,size(emgSig,2)), emgSig(iChan,:)+channelPos(iChan), 'color', channelColors(iChan,:), 'LineWidth',0.75)
                plot(linspace(300,1900,size(emgEnv,2)), emgEnv(iChan,:)*envScale(iChan)+channelPos(iChan), 'color', channelColors(iChan,:), 'LineWidth',1.5)
            end
            hold off
            
            outFrame = getframe(h);
            writeVideo(vOut,outFrame);
            
            close(h)
            
            frameInd = frameInd+1;
            disp(frameInd)
            
        end
        
        close(vOut);

    end
        
end