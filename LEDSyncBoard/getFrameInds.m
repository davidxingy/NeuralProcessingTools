function [framePulseStart, framePulseMean, badCounterSeg, badCounterSegInfo, badGridSeg, badGridSegInfo] = getFrameInds(startFrame, ledStatus)

% post-processing of LED synching after extracint LED on/offs
prevCounterValue=0;
stopFrame=size(ledStatus,1);

%initialize memory
framePulseStart=zeros(1,stopFrame-startFrame);
framePulseMean=zeros(1,stopFrame-startFrame);
badCounterSeg=zeros(1,stopFrame-startFrame);
badGridSeg=zeros(1,stopFrame-startFrame);

% the way I defined my regions was correct:
% 
%   18 28 38 ... 88
% 8 17 27 37 ... 87 96
% .  .  .  .      .  .
% .  .  .  .      .  .
% .  .  .  .      .  .      Actual LED movement
% 2  .  .  .      .  .
% 1 10 20 30 ... 80 89
%   9  19 29 ... 79
% 

% mappedInd=8:-1:1; %first column
% 
% for iCol=1:8 %columns 2-9
%     mappedInd=[mappedInd, ((iCol-1)*10+18):-1:((iCol-1)*10)+9];
% end
% mappedInd=[mappedInd  88:-1:81]; %last column

mappedInd=1:96;

cycleSize = 97; %there was 1 blank at the end of each cycle

%Now go through each frame and determine the pulse number based on the LEDs
for iFrame=1:stopFrame-startFrame
    
    %First, get the binary clock counter value=============================
    clockLit=ledStatus(iFrame+startFrame-1,97:end);
    counterValue=bi2de(clockLit);
    
    %if counter value changed, it should increment by 1
    if (counterValue~=prevCounterValue)
        
        if (counterValue-1~=prevCounterValue && counterValue-2~=prevCounterValue)
            % Uh-oh! The counter value isn't the same as the previous
            % value, or 1 more than the previous value.
            % sometimes counter LEDs bleed through to the next frame, eg:
            % previous count:               1 1 1 0 1 0 1 ...
            % current count (should be):    0 0 0 1 1 0 1 ...
            % current count (actuallt is):  1 1 1 1 1 0 1 ... <-lower bits bleed
            
            %find the smallest bit that is not on in the previous counter
            ceilBit=find(de2bi(prevCounterValue,16)==0,1);
            
            %if it is the smallest bit, then it's not the bleed problem,
            %something else is wrong, notify user
            if (ceilBit==1)
                warning(['Inconsistent counter value; At frame ' num2str(iFrame+startFrame-1)])
                badCounterSeg(iFrame)=1;
                disp(['Prev clock: ' num2str(ledStatus(iFrame+startFrame-2,97:end))])
                disp(['This clock: ' num2str(ledStatus(iFrame+startFrame-1,97:end))])
                badCounterSegInfo{iFrame}=[ledStatus(iFrame+startFrame-2,97:end); ledStatus(iFrame+startFrame-1,97:end)];
            else
                %see if removing the bleed fixes it
                modClockLit=clockLit;
                modClockLit(1:ceilBit-1)=0;
                
                if(bi2de(modClockLit)-1~=prevCounterValue)
                    %something else is the problem, notify user
                    warning(['Inconsistent counter value; At frame ' num2str(iFrame+startFrame-1)])
                    disp(['Prev clock: ' num2str(ledStatus(iFrame+startFrame-2,97:end))])
                    disp(['This clock: ' num2str(ledStatus(iFrame+startFrame-1,97:end))])
                    badCounterSegInfo{iFrame}=[ledStatus(iFrame+startFrame-2,97:end); ledStatus(iFrame+startFrame-1,97:end)];
                    badCounterSeg(iFrame)=1;
                else
                    %did increment by one when accounted for bleed, we're good
                    counterValue = bi2de(modClockLit);
                end
            end
            
        end
        
    end
    prevCounterValue = counterValue;
    
    %Next, get the LED grid values=========================================
    gridLit=sort(mappedInd(logical(ledStatus(iFrame+startFrame-1,1:96))));
    
    startLED=min(gridLit);%+floor((min(gridLit)-1)/8); %account for blanks
    stopLED=max(gridLit);%+floor((max(gridLit)-1)/8);
    
    check{iFrame}=diff(gridLit);
    
    if isempty(gridLit)
        %uh oh, no LEDs were on!
        badGridSeg(iFrame)=1;
        framePulseStart(iFrame)=NaN;
        framePulseMean(iFrame)=NaN;
        continue;
    end
    
    %see if it is continuous strip of LEDs (sometimes, one of the LEDs in 
    %the strip wasn't detected as on, this means the difference is at most 
    %2, we'll still count this as ok)
            
    %Maybe the LEDs lapped around
    circGrid=toeplitz([gridLit(1) fliplr(gridLit(2:end))], gridLit);
    
    if(all(diff(gridLit,1,2)==1 | diff(gridLit,1,2)==2) || length(gridLit)==1)
        
        %we're good
        framePulseStart(iFrame)=startLED+counterValue*cycleSize;
        framePulseMean(iFrame)=(startLED+stopLED)/2+counterValue*cycleSize;
        
    elseif (any(all(diff(circGrid,1,2)==1 | diff(circGrid,1,2)==2 | ...
            diff(circGrid,1,2)==-95 | diff(circGrid,1,2)==-94,2)))
        
        %LEDs lapped around 
        
        %I'm going to assume that the counter is therefore at the higher
        %value (there was counter bleeding)
        
        gridLit(gridLit>44)=gridLit(gridLit>44)-88;
        startLED=min(gridLit)+floor((min(gridLit)-1)/8); %account for blanks
        stopLED=max(gridLit)+floor((max(gridLit)-1)/8);
        
        framePulseStart(iFrame)=startLED+counterValue*cycleSize;
        framePulseMean(iFrame)=(startLED+stopLED)/2+counterValue*cycleSize;
    else
        
        %ok, some random LEDs probably got erroneously classified as on
        warning(['Bad grid values; At frame ' num2str(iFrame+startFrame-1)])
        disp(['Grid values: ' num2str(gridLit)])
        badGridSeg(iFrame)=1;
        badGridSegInfo{1,iFrame}=gridLit;
        
        %we'll just take the longest "continuous" (ie up to 1 missing led
        %in the strip) strip, eg:
        % x x o x x o o o x o x x x x o o x x ... (x - off, o - on)
        %     ^     |_______|         ^ ^
        % error  this is the strip    error
        
        gridDiffs=diff(gridLit);

%         if (gridDiffs(1)-gridDiffs(end)==-87 || gridDiffs(1)-gridDiffs(end)==-86)
%             %wrap around
%             gridDiffs=[gridDiffs gridDiffs(1)-gridDiffs(end)+88];
%         end
            
        gridDiffs(gridDiffs==2)=1; %again ,its ok to skip 1 in the strip
                
        contPositions=find(gridDiffs==1); %contPositions is the location in gridLit of where the next values are continuous
        if isempty(contPositions)
            %no two continuous LEDs, nothing we can do...
            framePulseStart(iFrame)=NaN;
            framePulseMean(iFrame)=NaN;
            
        elseif (length(contPositions)==1)
            %there's only a single two-LED long continuous strip
            gridLit=gridLit(contPositions:contPositions+1);
            startLED=min(gridLit)+floor((min(gridLit)-1)/8); %account for blanks
            stopLED=max(gridLit)+floor((max(gridLit)-1)/8);
            
            disp(['Fixed grid values: ' num2str(gridLit)])
            badGridSegInfo{2,iFrame}=gridLit;
            
            framePulseStart(iFrame)=startLED+counterValue*cycleSize;
            framePulseMean(iFrame)=(startLED+stopLED)/2+counterValue*cycleSize;
            
        else
            %find the longest strip (in contPositions)
            iStrip=1;
            iStripStart=1;
            iStripLength=0;
            stripStarts=[];
            stripLengths=[];
            for iPos=1:length(contPositions)-1
                
                if (iPos==length(contPositions)-1) %hit the end
                    if (contPositions(iPos)+1==contPositions(iPos+1))
                        %
                        stripStarts(iStrip)=iStripStart;
                        stripLengths(iStrip)=iStripLength+2;
                    else
                        stripStarts(iStrip)=iStripStart;
                        stripLengths(iStrip)=iStripLength+1;
                        
                        stripStarts(iStrip+1)=length(contPositions);
                        stripLengths(iStrip+1)=1;
                    end
                elseif (contPositions(iPos)+1==contPositions(iPos+1))
                    iStripLength = iStripLength + 1;
                else
                    stripStarts(iStrip)=iStripStart;
                    stripLengths(iStrip)=iStripLength+1;
                    iStrip = iStrip+1;
                    
                    iStripStart=iPos+1;
                    iStripLength=0;
                end
            end
            
            if(sum(max(stripLengths)==stripLengths)>1)
                %more than one strip with the biggest strip length, can't
                %decide which strip to use
                framePulseStart(iFrame)=NaN;
                framePulseMean(iFrame)=NaN;
            else
                %use the biggest strip
                [~, iMaxStrip]=max(stripLengths);
                startPos=contPositions(stripStarts(iMaxStrip));
                endPos=contPositions(stripStarts(iMaxStrip)+stripLengths(iMaxStrip)-1);
                
                gridLit=gridLit(startPos:endPos+1);
                
                startLED=min(gridLit)+floor((min(gridLit)-1)/8); %account for blanks
                stopLED=max(gridLit)+floor((max(gridLit)-1)/8);
                
                disp(['Fixed grid values: ' num2str(gridLit)])
                badGridSegInfo{2,iFrame}=gridLit;
                
                framePulseStart(iFrame)=startLED+counterValue*cycleSize;
                framePulseMean(iFrame)=(startLED+stopLED)/2+counterValue*cycleSize;
            end
            
        end
        
    end
    
end


%
