function [dMat, traceOut] = Widefield_analogToDesign(traceIn, stdThresh, trialCnt, sourceRate, targRate, motorIdx, gaussShift)
% code to create a peri-event design matrix based on an analog trace. Trace
% should be continous and will be reshaped into a trial structure to create
% a design matrix for every trial individually.
% Inputs:   traceIn = analog trace.
%           stdThresh = threshold for event detection in SDUs
%           trialCnt = number of trials in the dataset. used to reshape analog trace.
%           sourceRate = sampling rate of analog trace
%           targRate = sampling rate of design matrix.
%           motorIdx = index for peri-event matrix.
%           gaussShift = variable for subsampling in case model is used with gaussian convolution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traceIn = (traceIn - prctile(traceIn,1))./ nanstd(traceIn); %minimum values are at 0, signal in standard deviation units
traceOut = traceIn; %return normalized analog trace
traceIn = traceIn > stdThresh; %take activity above supplied threshold as indicator for event occurence
traceIn = diff([0; traceIn]) == 1; %find event onsets
traceIn = reshape(traceIn,[],trialCnt); %reshape to trial structure
frames = size(traceIn,1) / (sourceRate/targRate); %second dimension is trials so first should be frames per trial when taking differences in sampling rate into account

dMat = cell(1,trialCnt);
for iTrials = 1:trialCnt
    
    trace = logical(histcounts(find(traceIn(:,iTrials)), 0: sourceRate/targRate : (sourceRate/targRate)*frames))'; %resample to imaging frame rate. This is the zero lag regressor.    
    
    % create full design matrix
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    dMat{iTrials} = false(frames, length(motorIdx));
    dMat{iTrials}(cIdx(:)) = true;
    dMat{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    dMat{iTrials}(end,2:end) = dMat{iTrials}(end-1,1:end-1); %replace with shifted version of previous timepoint

    if gaussShift > 1
        dMat{iTrials} = dMat{iTrials}(:,1:gaussShift:end);
    end
end
end