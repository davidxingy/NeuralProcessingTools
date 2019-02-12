function [processed_data, dropped_sig, dropped_sig_intol, artifact_inds, artifacts_intol, flatLineInds, flatlines_intol, noiseLevel]=artifactRejection(chan_data, maxNumSamples)
% [processed_data, dropped_sig, dropped_sig_intol, artifact_inds, artifacts_intol, flatLineInds, flatlines_intol]=artifact_rejection(chan_data, [maxNumSamples])
%
% chan_data is the NxM matrix contianing the M-sample signal for N channels
% maxNumSamples is an optional input parameter. For very large input data
% matrices, the function might run out of memory. Set the max number of
% samples to run the algorithm on and the code will automatically divide
% the data into blocks of size NxmaxNumSamples and run the algorithm on
% each of them (and merging the results at the end).
%
% Returns the cleaned up signal for all changes as processed_data
% Returns dropped_sig as the indices where the signal from transmitter was dropped
% (flatlines in all channels).
% Returns the cell array artifact_inds as the indices where it found
% artifacts in the signal (each cell is for each channel).
% Returns the cell array flatLineInds as the indices where it found
% flatlines (10 or more constant values) in each of the channels.
%
% In addition, returns the indices where it disregards periods of
% dropped/artifact/flatline if it was less than 2 ms, and also counts data
% between drops/artifacts/flatlines as part of the drop/artifact/flatline
% if it is less than 30 samples. These are the _intol outputs
%
% Finally, the program can calculate the noise level of the signal in 1000-
% sample blocks by finding the standard deviation of the 98% percentile
% data (so no outliers). This is returned in the noiseLevel output.
% 
% If you just want to get the dropped_sig (i.e. for quickly determining what %
% of the signal was dropped), just specify the first 2 or 3 outputs when
% calling the function.
% 
% 
% Finally, sometimes the signal, when it is very close to the top and
% bottom rails (limits), it will "flip" and jump to the other rail. The
% processed_data will remove these jumps
% 
% David Xing
% Last updated 1/10/2018


% first divide data into blocks if neccesary.
isDivided = false;
if (nargin == 2) && (size(chan_data,2) > maxNumSamples)
    
    %More samples than the max number of samples
    isDivided = true;
    
    %split into blocks
    clear dataBlocks
    nBlocks=ceil(size(chan_data,2)/maxNumSamples);
    blockStartInds = zeros(1,nBlocks);
    blockEndInds = zeros(1,nBlocks);
    
    for iBlock = 1:nBlocks
        blockStartInds(iBlock) = (iBlock-1)*maxNumSamples+1;
        blockEndInds(iBlock) = min(iBlock*maxNumSamples, size(chan_data ,2));
        dataBlocks{iBlock} = chan_data(:, blockStartInds(iBlock):blockEndInds(iBlock));
    end
    
    fprintf('Split data into %i blocks \n',nBlocks);
    clear chan_data;
else
    dataBlocks = {chan_data};
    clear chan_data;
end

blockEndsInJump=false(1,size(dataBlocks{1},1)); %since we could be dividing blocks right in the middle of a jump

for iBlock = 1:length(dataBlocks)
    
    data = dataBlocks{iBlock};
    if (isDivided)
        fprintf('Running on Block %i',iBlock);
    end
    
    % dropped signal analysis (bad data across all channels)=======================================
    
    % define when signal is dropped when for 4 consecutive data points the signal is changing by 1
    % or less (sometimes in some channels the signal looks constant but for some reason
    % it moves up or down one digital value [based on
    % 201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001, never more
    % than one]) in all 96 channels (based on
    % 201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001, there aren't
    % any data points where between 44 and 95 channels satisfy this criteron,
    % so I think it's a good criteron for classfying when signal drops vs when
    % by chance all electrodes happen to have the same signal for 4 data
    % points)
    
    fprintf('\n Finding dropped signals...')
    chan_diffs=(abs(data(:,1:end-4)-data(:,2:end-3)))<=1&...
        (abs(data(:,2:end-3)-data(:,3:end-2)))<=1&...
        (abs(data(:,3:end-2)-data(:,4:end-1)))<=1&...
        (abs(data(:,4:end-1)-data(:,5:end)))<=1;
    
    dropped_sig=int32(find(sum(chan_diffs)==96));
    clear chan_diffs
    
    % Define what is tolerable for data analysis (could be different depending
    % on if using spikes or LFP) in terms of lost signal
    
    dropped_limit=60; %2 ms or less of bad data will be ok for now
    dropped_spacing=30; %if theres less than 30 samples between bad samples, count all of it as bad
    
    % Now get the intolerable time periods
    [dropped_sig_intol]=find_intolerable(dropped_sig, dropped_spacing, dropped_limit,size(data,2));
    
    %if the data is blocked, the indices are shifted by the previous blocks sizes
    if (isDivided)
        allDroppedSig{iBlock} = dropped_sig+blockStartInds(iBlock)-1;
        allDroppedSigIntol{iBlock} = dropped_sig_intol+blockStartInds(iBlock)-1;
    end
    
    % User just wanted to get dropped sig
    if nargout==2 || nargout==3
        processed_data{iBlock}=data;
        continue
    end
    
    %First fix any issues with voltages hitting the limits and flipping to opposite rails==========================================
    
    fprintf('\n Cleaning up rail flips')
    
    processed_data{iBlock}=int32(data);

    jumpThresh=1.5e4; %these jumps are very large
    fixRailJumpByPairs=false;
    
    for iChan=1:size(data,1)
        fprintf('.')
        
        %find jump times
        jumpUpInds=sort(find(diff(processed_data{iBlock}(iChan,:))>jumpThresh))+1;
        jumpDownInds=sort(find(diff(processed_data{iBlock}(iChan,:))<-1*jumpThresh))+1;
        
        if fixRailJumpByPairs
            %see if the previous block ended in the middle of a segment jump
            if iBlock~=1
                if blockEndsInJump(iChan)
                    %if it did, next check if there is a jump between the
                    %beginning of this block and the end of the last block
                    if abs(processed_data{iBlock}(iChan,1)-processed_data{iBlock-1}(iChan,end))>jumpThresh
                        %we are still in the middle of a segment jump:
                        %`'`'`'         |<-block divider  `'`'`'`'
                        %      ,.,.,..,.|,.,.,,.,.,.,.,.,.,
                        %   Last Block  |    Current block
                        
                        %we need to move everything from the first point up to
                        %the first jump in this block up/down
                        firstJumpInd=min([jumpUpInds(1) jumpDownInds(1)]);
                        jumpAmount=diff(processed_data{iBlock}(iChan,firstJumpInd-1:firstJumpInd));
                        processed_data{iBlock}(iChan,1:firstJumpInd-1)=...
                            processed_data{iBlock}(iChan,1:firstJumpInd-1)+jumpAmount;
                        
                        %and remove the first jump from the list of jumps
                        jumpUpInds(jumpUpInds==firstJumpInd)=[];
                        jumpDownInds(jumpDownInds==firstJumpInd)=[];
                        
                    else
                        %the segment jump just ended with the last point of the
                        %previous block (no jump because we had already moved
                        %everything up in the last block):
                        %`'`'`'              |`'`'`'`''`'`'`'`'`'
                        %      ,,.,.,.,.,.,.,| <-block divider
                        %    Last Block      |    Current block
                        
                        %don't need to do anything
                    end
                    
                else
                    %Check if there is a jump between the
                    %beginning of this block and the end of the last block
                    if abs(processed_data{iBlock}(iChan,1)-processed_data{iBlock-1}(iChan,end))>jumpThresh
                        %if there is a jump, the the very first point is the
                        %start of a segement shifting:
                        %`'`'`'`'`'`'`'`'| <-block divider     `'`'`'`'`'
                        %                |.,.,.,.,.,.,,.,.,.,.,
                        %   Last Block   |    Current block
                        
                        %we need to add the very first index to the list of
                        %jumps
                        if processed_data{iBlock}(iChan,1)>processed_data{iBlock-1}(iChan,end)
                            jumpUpInds=[1 jumpUpInds];
                        else
                            jumpDownInds=[1 jumpDownInds];
                        end
                        
                    else
                        %no segments shifts
                        %`'`'`'`'`'`'`'`'|`'`'`'`'`'`'`'`'`'`'`'`
                        %                | <-block divider
                        %   Last Block   |    Current block
                        
                        %don't need to do anything
                    end
                    
                end
            end
            
            %jump times should always be in pairs, or at most 1 more of one of
            %them
            if abs(length(jumpUpInds)-length(jumpDownInds))>1
                warning('Number of upwards rail jumps does not equal the number of downwards rail jumps!')
            end
            
            %no jumps, don't need to do anything
            if isempty(jumpUpInds) && isempty(jumpDownInds)
                continue
            end
            
            %now, get pairs of up/down or down/up
            if isempty(jumpUpInds)
                %no up jumps, should be only one down jump
                initialJumpInds=jumpDownInds;
                returningJumpInds=[];
            elseif isempty(jumpDownInds)
                %no down jumps, should be only one up jump
                initialJumpInds=jumpUpInds;
                returningJumpInds=[];
            else
                
                allJumpInds=sort([jumpUpInds jumpDownInds]);
                
                isPair=false;
                initialJumpInds=[];
                returningJumpInds=[];
                for iJump=1:length(allJumpInds)
                    
                    if isPair
                        %second part of pair, go to next jump to find new pair
                        isPair=false;
                        continue
                    end
                    
                    %if we're at the last jump, and it isn't part of a pair,
                    %then we've ended the block in the middle of a jump
                    if iJump==length(allJumpInds)
                        initialJumpInds=[initialJumpInds allJumpInds(iJump)];
                        
                        continue
                    end
                    
                    %the next closest jump should be a jump in the opposite
                    %direction to complete the pair
                    if any(allJumpInds(iJump)==jumpUpInds)
                        if all(allJumpInds(iJump+1)~=jumpDownInds)
                            warning('An upwards Jump isn''t followed by a downwards jump!')
                            continue
                        end
                    else
                        if all(allJumpInds(iJump+1)~=jumpUpInds)
                            warning('A downwards Jump isn''t followed by am upwards jump!')
                            continue
                        end
                    end
                    
                    initialJumpInds=[initialJumpInds allJumpInds(iJump)];
                    returningJumpInds=[returningJumpInds allJumpInds(iJump+1)];
                    isPair=true;
                end
                
            end
            
            %now go and fix all the jumped segments
            [processed_data{iBlock}(iChan,:), blockEndsInJump(iChan)]=fixJumpSegments(...
                processed_data{iBlock}(iChan,:), initialJumpInds, returningJumpInds);
            
        else

            %do the fast and easy way of fixing jumps (just shift
            %everything after the jump back,) however this method may
            %result in some offsets
            
            %Shift everything in this block to the last point of the last
            %block, since things may have been offset (and the offset could
            %be small, so it might slip through the jump threshold)
            if iBlock~=1
                initialJumpAmount=int32(data(iChan,1))-processed_data{iBlock-1}(iChan,end);
                if initialJumpAmount>0
                    jumpUpInds=[1 jumpUpInds];
                elseif initialJumpAmount<0
                    jumpDownInds=[1 jumpDownInds];
                end
            end
            
            %now go through each jump and remove it by shifting everything after
            %the jump by the jump amount
            for iJump=1:length(jumpUpInds)
                
                %see how big the jump is
                if jumpUpInds(iJump)==1
                    jumpAmount=initialJumpAmount;
                else
                    jumpAmount=diff(int32(data(iChan,jumpUpInds(iJump)-1:jumpUpInds(iJump))));
                end
                
                %shift everything after the jump by the opposite amount
                processed_data{iBlock}(iChan,jumpUpInds(iJump):end)=processed_data{iBlock}(iChan,jumpUpInds(iJump):end)-jumpAmount;
                
            end
            
            %same for downward jumps
            for iJump=1:length(jumpDownInds)
                
                if jumpDownInds(iJump)==1
                    jumpAmount=initialJumpAmount;
                else
                    jumpAmount=diff(int32(data(iChan,jumpDownInds(iJump)-1:jumpDownInds(iJump))));
                end
                
                processed_data{iBlock}(iChan,jumpDownInds(iJump):end)=processed_data{iBlock}(iChan,jumpDownInds(iJump):end)-jumpAmount;
                
            end
            
        end
    
    end
    
    
    % Receiver artifact analysis (different for each channel)=======================================
    
    % Looked at distribution of the difference between subsequent time points
    % (based on 201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001, chan 60)
    % and found that it drops off at around 200 digital values. I used to
    % use 200 as a hard-coded cutoff but I realized that the amount that
    % the signal changes between samples can vary across recordings and
    % even within the same recording! I'm therefore going to divide the
    % data into windows of 30 000 samples and defining the cutoff as the
    % 99.9th percentile of the data for each window.
    
    fprintf('\n Finding artifacts')
    
%     cutoff=200;
    
    artifact_inds={};
    false_artifact_inds={};
    artifacts_intol={};
        
    % go through and do the analysis for each channel
    for whichchan=1:size(data,1)
        
        fprintf('.')
        
        %split into windows of at least 1000 samples (up to 1999 if
        %there's lots of samples left over)
        windowSize=1000;
        nWindows=floor(size(data,2)/windowSize);
        
        for iWindow=1:nWindows
        
            %get the window of data
            startInd=(iWindow-1)*windowSize+1;
            if iWindow==nWindows
                endInd=size(data,2);
            else
                endInd=iWindow*windowSize;
            end
            signal=processed_data{iBlock}(whichchan,startInd:endInd);
            diffs=diff(signal);

            %get the cutoff (7x standard deviation, where std is calculated
            %from the 98th percetile of data to remove outliers).
            percentileCutoff=prctile(abs(diffs),98);
            cutoff=7*std(double(diffs(diffs<percentileCutoff & diffs>-1*percentileCutoff)));
            
            %also calculate the overall noise level (std) of the actual
            %signal, not just the diffs
            percentileCutoff=prctile(abs(signal),98);
            noiseLevel{iBlock}(whichchan,iWindow)=...
                std(double(signal(signal<percentileCutoff & signal>-1*percentileCutoff)));
            
            
            %find time points where signal jumps bigger than the cutoff
            drops=find(diffs<=-cutoff);
            rises=find(diffs>=cutoff);
            
            %Looking at chan 60 of the
            %201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001 data, using
            %200 as a cutoff, seems to find natural looking jumps which
            %don't seem to be artifact. Artifacts should only happen when there is
            %a big jump, followed by a few constant values (no change after the
            %jump), and then a return to previous values before the jump.
            %(In some artifacts, there could be a series of jumps in the same direction
            %before going to constant value/jumping back)
            %If it stays at a constant value for too long (>10, based on the distribution
            %of constant lengths in
            %201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001 chan 60)
            %after a jump, then consider it an artifact, even if it doesn't jump
            %back in the opposite direction.
            
            %go through each jump and see if it looks like what I describe above
            %first look at drops
            [isdrop, drop_inds, non_drop_inds, drop_waveforms, non_drop_waveforms]=...
                confirm_artifact(signal, diffs, drops, rises, whichchan, false);
            
            %then look at rises
            [isrise, rise_inds, non_rise_inds, rise_waveforms, non_rise_waveforms]=...
                confirm_artifact(signal, diffs, rises, drops, whichchan, false);
            
            %combine artifacts derrived from large drops and rises (add 1 due to
            %offset by one from diff, and add the window offset)
            artifact_inds{whichchan}{iWindow}=uint32(unique(sort([drop_inds+1 rise_inds+1])')'+windowSize*(iWindow-1));
            
            %false artifacts are ones that didn't get rejected by both runs
            nonartifact_inds=unique(sort([non_drop_inds+1 non_rise_inds+1]))+windowSize*(iWindow-1);
            false_artifact_inds{whichchan}{iWindow}=setdiff(nonartifact_inds,artifact_inds{whichchan}{iWindow});
            if isempty(false_artifact_inds{whichchan}{iWindow})
                false_artifact_inds{whichchan}{iWindow}=[];
            end
            
            
            %Now, as with the signal dropout, want to find intolerable artifact
            %time periods
            artifact_limit=60; %2 ms or less of bad data will be ok for now
            artifact_spacing=30; %if theres less than 30 samples between bad samples, count all of it as bad
            
            [artifacts_intol{whichchan}{iWindow}]=...
                find_intolerable(artifact_inds{whichchan}{iWindow}, artifact_spacing, artifact_limit,length(signal));
%             plot(signal)
%             hold on
%             plot(artifact_inds{whichchan}{iWindow}-windowSize*(iWindow-1),signal(artifact_inds{whichchan}{iWindow}-windowSize*(iWindow-1)),'*');
        end
        
        %combine all the windows
        artifact_inds{whichchan}=[artifact_inds{whichchan}{:}];
        false_artifact_inds{whichchan}=[false_artifact_inds{whichchan}{:}];
        artifacts_intol{whichchan}=uint32([artifacts_intol{whichchan}{:}]);
    end
    
    % Flatline analysis (different for each channel, unsure what the origin is)=========================
    
    % Looked at distribution of the number of constant consecutive points,
    % decided that any flatlines with length greater than 10 should be
    % considered a flatline artifact, and not physiologial signal (based on
    % 201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001, chan 60).
    % Won't alter flatline signals, since not much to do to alleviate them, but
    % will let the user be aware of where these happen.
    
    fprintf('\n Finding flatlines')
    
    for whichchan=1:size(data,1)
        fprintf('.')
        
        %find flatlines
        signal=data(whichchan,:);
        flatLineInds{whichchan}=findFlatlines(signal, 11);
        
        artifact_limit=60; %2 ms or less of bad data will be ok for now
        artifact_spacing=30; %if theres less than 30 samples between bad samples, count all of it as bad
        % Now get the intolerable time periods
        [flatlines_intol{whichchan}]=...
            find_intolerable(flatLineInds{whichchan}, artifact_spacing, artifact_limit,length(signal));
    end
    
    
    %Now, try to alleviate artifacts with linear interpolation==========================================
    
    fprintf('\n Cleaning up signal')
        
    for whichchan=1:size(data,1)
        fprintf('.')
        
        % first combine flatlines and switching artifacts
        artifacts_merge=...
            unique([flatLineInds{whichchan} dropped_sig artifact_inds{whichchan}]);
        
        % if artifacts are closely spaced together (less than 20 samples) consider
        % it all artifact (for when there's bunch of switching)
        artifacts_all=find_intolerable(artifacts_merge, 20, 0, size(data,2));
        
        % if there are no artifacts, don't need to interpolate
        if (isempty(artifacts_all))
            continue
        end
        
        %find start and end indices of artifact periods
        bad_jump_inds=find(diff(artifacts_all)>1);
        
        interp_starts=[artifacts_all(1) artifacts_all(bad_jump_inds+1)]-1;
        interp_starts(interp_starts==0)=1; interp_starts=unique(interp_starts);
        interp_ends=[artifacts_all(bad_jump_inds) artifacts_all(end)]+1;
        interp_ends(interp_ends>size(data,2))=size(data,2);
        
        assert(length(interp_starts)==length(interp_ends))
        
        %go through each artifact period and linearly interpolate
        for whichartifact=1:length(interp_starts)
            beg_val=double(processed_data{iBlock}(whichchan,interp_starts(whichartifact)));
            end_val=double(processed_data{iBlock}(whichchan,interp_ends(whichartifact)));
            numsamples=interp_ends(whichartifact)-interp_starts(whichartifact)+1;
            lin_fill=linspace(beg_val,end_val,numsamples);
            
            processed_data{iBlock}(whichchan,interp_starts(whichartifact):interp_ends(whichartifact))=lin_fill;
        end
    end
    fprintf('\n')
    
    %if the data is blocked, the indices are shifted by the previous blocks sizes
    if (isDivided)
        allArtifactInds(iBlock,:) = cellfun(@(x) x+blockStartInds(iBlock)-1, artifact_inds, 'un', 0);
        allArtifactIndsIntol(iBlock,:) = cellfun(@(x) x+blockStartInds(iBlock)-1, artifacts_intol, 'un', 0);
        allFlatlineInds(iBlock,:) = cellfun(@(x) x+blockStartInds(iBlock)-1, flatLineInds, 'un', 0);
        allFlatlineIndsIntol(iBlock,:) = cellfun(@(x) x+blockStartInds(iBlock)-1, flatlines_intol, 'un', 0);
    end
    
end

processed_data = cat(2, processed_data{:});

if exist('noiseLevel')
    noiseLevel = cat(2, noiseLevel{:});
end

if (isDivided)
    %merge dropped sigs
    dropped_sig = cat(2, allDroppedSig{:});
    clear allDroppedSig;
    dropped_sig_intol = cat(2,allDroppedSigIntol{:});
    clear allDroppedSigIntol;
    
    if nargout==2 || nargout==3
        return;
    end
    
    %merge artifact and flatline inds
    for iChan = 1:size(allArtifactInds,2)
        artifact_inds{iChan} = cat(2, allArtifactInds{:,iChan});
        artifacts_intol{iChan} = cat(2, allArtifactIndsIntol{:,iChan});
        flatLineInds{iChan} = cat(2,allFlatlineInds{:,iChan});
        flatlines_intol{iChan} = cat(2,allFlatlineIndsIntol{:,iChan});
    end
end


%end of main function


function [intolerable_samples]=find_intolerable(badsamples, merge_spacing, tolerance, nSig)

% transform good samples in between bad samples to bad based on merge_spacing
bad_sig=zeros(1,nSig);
bad_sig(badsamples)=1;

intolerable_samples=uint32([]);

if isempty(badsamples)
    return
end

for whichbadsample=2:length(badsamples)
    if (badsamples(whichbadsample)-badsamples(whichbadsample-1)<=merge_spacing)
        bad_sig(badsamples(whichbadsample-1):badsamples(whichbadsample))=1;
    end
end

if all(bad_sig)
    %the whole trial is all bad
    intolerable_samples=find(bad_sig);
    return
elseif ~any(bad_sig)
    %no bad periods
    return
end

% find bad period starts and stops
bad_periods=find(bad_sig(2:end)~=bad_sig(1:end-1));

% check that trial doesn't start or stop as bad (since that will mess with
% my algorithm for finding time blocks)
if (bad_sig(bad_periods(1))==1)
    bad_periods=[1 bad_periods];
elseif (bad_sig(bad_periods(end))==1)
    bad_periods=[bad_periods length(bad_sig)];
end

% go through and remove bad periods if their length is within tolerance
intolerable_sig=bad_sig;
for whichperiod=1:(length(bad_periods)/2)
    if (bad_periods(2*whichperiod)-bad_periods(2*whichperiod-1)<=tolerance)
        intolerable_sig(bad_periods(2*whichperiod-1):bad_periods(2*whichperiod))=0;
    end
end

% output new indices
intolerable_samples=find(intolerable_sig);





function [isartifact, artifact_inds, nonartifact_inds, artifact_waveforms, non_artifact_waveforms]=...
    confirm_artifact(signal, diffs, jump, jump_opposite, chan, plot_waveforms)

isartifact=[];
artifact_inds=uint32([]);
artifact_waveforms={};
non_artifact_waveforms={};


%make figures if want to plot
if (plot_waveforms)
    artifact_fig=figure; title(['Chan: ' num2str(chan) ', Artifact Waveforms']); hold on;
    nonartifact_fig=figure; title(['Chan: ' num2str(chan) ', Non-artifact Waveforms']); hold on;
end

artifact_inds_cell={};
nonartifact_inds_cell={};

for whichjump=1:length(jump)
    
    %The jump is very large (>500 based on histogram of diffs) <---*Not doing this anymore!!!*
    %
    %OR
    %
    %it is followed or preceded by another jump in the opposite
    %direction, with only constant values (+/- 1 change in digital value)
    %or same direction jumps in between.
    %
    %OR
    %
    %It is followed by or preceded by a long series of constants (>10,
    %based on the distribution of constant lengths in
    %201612201054-Starbuck_Treadmill-Array1480_Right-Trial00001 chan 60)
    
    next_ind=1;
    isartifact(whichjump)=false;
    while true
        
        if (next_ind>=11)
            %it is followed by more than 30 constants/jumps, consider
            %it an artifact (reason in main function comments)
            isartifact(whichjump)=true;
        end
        
        if (jump(whichjump)+next_ind==length(diffs) || jump(whichjump)==length(diffs))
            %we've reached the end
            break
            
        elseif (abs(diffs(jump(whichjump)+next_ind))<=1)
            %becomes constant value, keep going
            next_ind=next_ind+1;
            
        elseif (any(jump(whichjump)+next_ind==jump))
            %or is another jump in the same direction, keep going
            next_ind=next_ind+1;
            
        elseif (any(jump(whichjump)+next_ind==jump_opposite))
            %next point is a jump in the opposite direction, consider it artifact
            isartifact(whichjump)=true;
            next_ind=next_ind+1;
            
        else %next point is not a jump or a constant, end
            break
        end
        
    end
    
    if isartifact(whichjump)
        %if its an artifact, the end index of the artifact is until it goes
        %to a non-constant and non-jump
        start_ind=jump(whichjump);
        end_ind=jump(whichjump)+next_ind-2;
    end
    
%     %When a jump is very large, it has to be bad, find the closest jump in
%     %the opposite direction and say all datapoints in between are artifacts
%     if (~isartifact(whichjump) && abs(diffs(jump(whichjump)))>=500)  %If it's already been classified as artifact, don't need to do this
%         
%         isartifact(whichjump)=true;
%         [~,closeOppositeJump]=min(abs(jump_opposite-jump(whichjump)));
%         
%         %save ranges of the artifact
%         if (jump_opposite(closeOppositeJump)>jump(whichjump))
%             start_ind=jump(whichjump);
%             end_ind=jump_opposite(closeOppositeJump)-1;
%         else
%             start_ind=jump_opposite(closeOppositeJump);
%             end_ind=jump(whichjump)-1;
%         end
%         
%         %sometimes there are loooong periods of constants and the voltage
%         %changes alot at the end, without having an opposite jump in the beginning. In these
%         %cases, don't used the closest opposite jump, just ignore for now,
%         %as I will catch these with the long contants check
%         if (abs(start_ind-end_ind)>=1000)
%             isartifact(whichjump)=false;
%         end
%         
%     end
    
    %if it is an artifact, some jumps have trains of constants before it,
    %so look backwards to get rid of any constants
    if (isartifact(whichjump))
        back_ind=start_ind;
        
        while back_ind>=2 %stop if hit the beginning of the file
            if (abs(signal(back_ind)-signal(back_ind-1))<=1)
                %there's some constants
                back_ind=back_ind-1;
            else
                %end of constants
                break
            end
            
            if (back_ind==1)
                %reached begining of file
                break
            end
        end
        
    end
    
    %save the time points defined as artifact if it is an artifact
    if isartifact(whichjump)
        artifact_inds_cell{whichjump}=back_ind:end_ind;
        
        %and save artifact waveform
        artifact_waveforms{whichjump}=signal(back_ind:end_ind)-...
            signal(jump(whichjump));
        
        %plot if desired
        if (plot_waveforms)
            figure(artifact_fig); plot(artifact_waveforms{whichjump})
        end
    else
        nonartifact_inds_cell{whichjump}=jump(whichjump):jump(whichjump)+10;
        
        %save waveform of non-artifacts too (5 points before and after jump)
        non_artifact_waveforms{whichjump}=signal(max(jump(whichjump)-5,1):min(jump(whichjump)+next_ind+4,length(signal)))-...
            signal(jump(whichjump));
        
        %plot if desired
        if (plot_waveforms)
            figure(nonartifact_fig); plot(non_artifact_waveforms{whichjump})
        end
    end
    
end

artifact_inds=cat(2,artifact_inds_cell{:});
nonartifact_inds=cat(2,nonartifact_inds_cell{:});



function [flatInds, flatLengths]=getFlatLengths(signal)
flatLength=0;
flatLengths=[]; flatInds=uint32([]);
for iSignal=1:length(signal)-1
    if (signal(iSignal)==signal(iSignal+1))
        flatLength=flatLength+1;
    else
        flatLengths=[flatLengths flatLength];
        flatInds=[flatInds iSignal-flatLength];
        flatLength=0;
    end
end



function flatLineInds=findFlatlines(signal, minLength)
for iDiff=1:minLength+1
    shifts(:,iDiff)=signal(iDiff:end-minLength+iDiff-1);
end

diffs=diff(shifts,1,2);
clear shifts;

foundInds=find(sum(abs(diffs),2)==0);
init_length=length(foundInds)*10;
if isempty(foundInds)
    flatLineInds=[];
    return;
end
clear diffs;

flat_starts_inds=find(diff(foundInds)>1)+1;

flat_starts=[foundInds(1); foundInds(flat_starts_inds)];
flatLineInds=zeros(1,init_length);
clear foundInds

counter=1;
for iFlatline=1:length(flat_starts)
    
    nextPoint=0;
    while true
        if (flat_starts(iFlatline)+nextPoint)==length(signal)
            flatLineInds(counter:counter+nextPoint)=...
                flat_starts(iFlatline):flat_starts(iFlatline)+nextPoint;
            counter=counter+nextPoint+1;
            break
        elseif signal(flat_starts(iFlatline)+nextPoint)~=signal(flat_starts(iFlatline)+nextPoint+1)
            flatLineInds(counter:counter+nextPoint)=...
                flat_starts(iFlatline):flat_starts(iFlatline)+nextPoint;
            counter=counter+nextPoint+1;
            break
        else
            nextPoint=nextPoint+1;
        end
    end
end

if counter<=init_length
    flatLineInds(counter:end)=[];
end

flatLineInds=unique(sort(flatLineInds));


function [fixedData, endInJump]=fixJumpSegments(data2fix, initialJumpInds, returningJumpInds)
% function that fixes segements of data that have jumped up or down by
% shifting those segments back up or down

fixedData=data2fix;

%Go through each segment (i.e. each pair of initial-returning jumps)
for iJump=1:length(initialJumpInds)
    
    if length(returningJumpInds)<iJump
        %if no more jumps, then we have ended the block in the middle of 
        %the segment. Shift everything until the last point.
        jumpAmount=diff(data2fix(initialJumpInds(iJump)-1:initialJumpInds(iJump)));
        fixedData(initialJumpInds(iJump):end)=fixedData(initialJumpInds(iJump):end)-jumpAmount;
        endInJump=true;
        
    else
        %if there's more jumps then the next closest jump should be the
        %returning jump
        nextClosestJump=min([returningJumpInds(returningJumpInds>initialJumpInds(iJump)),...
            initialJumpInds(initialJumpInds>initialJumpInds(iJump))]);
        if isempty(nextClosestJump) || ~any(returningJumpInds==nextClosestJump)
            warning('An initialil Jump isn''t followed by a returning jump (in the opposite direction)!')
            continue
        end
        
        %now shift the segment (between the jumps) back
        jumpAmount=diff(data2fix(returningJumpInds(iJump)-1:returningJumpInds(iJump)));
        fixedData(initialJumpInds(iJump):returningJumpInds(iJump)-1)=...
            fixedData(initialJumpInds(iJump):returningJumpInds(iJump)-1)+jumpAmount;
        endInJump=false;
    end
    
end
            
            
            
%
