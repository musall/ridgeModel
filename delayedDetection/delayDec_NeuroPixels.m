function delayDec_NeuroPixels(Animal,reload)

if ~exist('reload','var')
    reload = false; %flag to reload video data and re-compute movement variables
end

%% some variables
nrBlocks = 16; %number of blocks for svd
nrDims = 200;
ridgeFolds = 10;
sRate = 30;
preStim = 0.5; %baseline before stim onset
stimDur = 2.5; %total duration for stimulus regressor/PSTH
fPath = [];
filterOrder = 50;
lowCut = 0.03; %lower cutoff for low-pass filter (for ACh)
highCut = 10; %higher cutoff for band-pass filter (for NE)
mPreTime = ceil(0.5 * sRate);  % precede motor events to capture preparatory activity in frames
mPostTime = ceil(3 * sRate);   % follow motor events for mPostStim in frames
motorIdx = [-(mPreTime: -1 : 1) 0 (1:mPostTime)]; %index for design matrix to cover pre- and post motor action

eyeThresh(2) = 0.075;
eyeThresh(1) = 0;
cPath = [pwd filesep 'Neuropixels' filesep]; %path to neuropixels data
fPath = [cPath filesep Animal filesep];

if strcmpi(Animal, 'NP9')
    videoData = 'NP9_006_Video_1.mp4';
    spikeData = 'NP9_006_spikes.mat';
    trialData = 'NP9_006_sync.mat';
    vidFrames = 76216;
elseif strcmpi(Animal, 'N14')
    videoData = 'N14_008_Video_1.mp4';
    spikeData = 'NP14_008_spikes.mat';
    trialData = 'NP14_008_sync.mat';
    vidFrames = 41334;
end

% regressor groups
regGroups = {'stim' 'pupil' 'nose' 'snout' 'motion' 'video' 'movement'}; %group names, second row are individual regressor labels
regGroups{2,1} = {'stim'};
regGroups{2,2} = {'fastPupil' 'slowPupil'};
regGroups{2,3} = {'nose'};
regGroups{2,4} = {'snout'};
regGroups{2,5} = {'ME'};
regGroups{2,6} = {'video'};
regGroups{2,7} = {'fastPupil' 'slowPupil' 'nose' 'snout' 'ME' 'video'};

%% load video data
if ~exist([fPath 'faceVars.mat'], 'file') || ~exist([fPath 'vidR.mat'], 'file') || ~exist([fPath 'absVidR.mat'], 'file') || reload
    Cnt = 0;
    v = VideoReader([fPath videoData]);
    cVideo = zeros(300,480,vidFrames, 'uint8');
    while hasFrame(v)
        Cnt = Cnt + 1;
        temp = readFrame(v);
        cVideo(:,:,Cnt) = arrayResize(temp(:,:,1),2);
        
        if rem(Cnt, 1000) == 0
            disp(Cnt)
        end
    end
    vidSize = size(cVideo);
    disp(Cnt);
end

%% check for facial features and create if missing
if ~exist([fPath 'faceVars.mat'], 'file') || reload
    [eyeVars, snoutMotion, noseMotion] = Behavior_getFaceMovements(cVideo, fPath, [], true);
    save([fPath 'faceVars.mat'],'eyeVars', 'snoutMotion', 'noseMotion');
else
    load([fPath 'faceVars.mat'],'eyeVars', 'snoutMotion', 'noseMotion');
end

%% create pupil variables
% remove outliers and filter pupil
pupil = mean(eyeVars.axes,2);
temp = (pupil - nanmean(pupil)) ./ nanstd(pupil);
idx = [0;diff(temp) > 0.25] | [0;diff(temp) < -0.25] | temp < -4 | temp > 4; %index for blinks / false reads
pupil(idx) = NaN;
pupilFill = fillgaps(pupil, 300); %fill NaN gaps between trials
pupilFill = [pupilFill;pupilFill(end-round(length(pupil) * 0.1):end)]; %add some padding to the end to avoid filter artefacts
pupilFill = zscore(pupilFill);

% high-pass butter filter for pupil data. this is to isolate adrenergic component
fastFilter  = fdesign.bandpass('N,F3dB1,F3dB2', filterOrder, lowCut, highCut, sRate);
fastFilter = design(fastFilter, 'butter');
fastPupil = fastFilter.filter(pupilFill); fastPupil = flipud(fastFilter.filter(flipud(fastPupil))); %filter twice to correct for phaseshift

% low-pass butter filter for pupil data. this is to isolate cholinergic component
slowFilter  = fdesign.lowpass('N,F3dB', filterOrder, lowCut, sRate);
slowFilter = design(slowFilter, 'butter');
slowPupil = slowFilter.filter(pupilFill); slowPupil = flipud(slowFilter.filter(flipud(slowPupil)));

%remove padding
pupilFill(end-round(length(pupil) * 0.1):end) = [];
fastPupil(end-round(length(pupil) * 0.1):end) = [];
slowPupil(end-round(length(pupil) * 0.1):end) = [];

%% check for video SVD and create if missing
if ~exist([fPath 'vidR.mat'], 'file') || reload
    cVideo = reshape(cVideo, [], vidSize(3)); %merge all pixels
    blockSize = ceil(vidSize(1:2)/sqrt(nrBlocks)); %size of each block
    ind = im2col(reshape(1:prod(vidSize(1:2)),vidSize(1:2)),blockSize,'distinct'); %build index for diffent image blocks
    
    tic
    vidV = NaN(nrDims, nrBlocks, size(cVideo,2));
    for x = 1 : size(ind,2)
        mov =  single(cVideo(ind(:,x),:));  %get current segment
        [~, s, Vr] = fsvd(mov, nrDims, 1, 0);
        vidV(:, x, :) = s * Vr'; %multiply S into V, so only U and V from here on
    end
    
    vidV = reshape(vidV, [], size(vidV,3)); %merge video dimensions
    [~, s, Vr] = fsvd(vidV, nrDims, 1, 0);
    vidR = s * Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    disp('First SVD complete'); toc
    save([fPath 'vidR.mat'],'vidR');
    cVideo = reshape(cVideo, vidSize); %back to original form
else
    load([fPath 'vidR.mat'],'vidR');
end

%% compute motion energy and repeat svd
if ~exist([fPath 'absVidR.mat'], 'file') || reload
    cVideo = reshape(cVideo, vidSize);
    cVideo = cat(3,cVideo(:,:,1),cVideo); %duplicate first frame
    cVideo = abs(diff(cVideo,[],3)); %compute absolute of temporal derivative (motion energy)
    cVideo = reshape(cVideo, [], vidSize(3)); %merge all pixels
    blockSize = ceil(vidSize(1:2)/sqrt(nrBlocks)); %size of each block
    ind = im2col(reshape(1:prod(vidSize(1:2)),vidSize(1:2)),blockSize,'distinct'); %build index for diffent image blocks
    
    vidV = NaN(nrDims, nrBlocks, size(cVideo,2));
    for x = 1 : size(ind,2)
        mov =  single(cVideo(ind(:,x),:));  %get current segment
        tic; [~, s, Vr] = fsvd(mov, nrDims, 1, 0); toc;
        vidV(:, x, :) = s * Vr'; %multiply S into V, so only U and V from here on
    end
    vidV = reshape(vidV, [], size(vidV,3)); %merge video dimensions
    [~, s, Vr] = fsvd(vidV, nrDims, 1, 0);
    absVidR = s * Vr';
    absVidR = bsxfun(@minus, absVidR, mean(absVidR, 2));
    absVidR = bsxfun(@rdivide, absVidR, std(absVidR, [], 2))'; %this is the motion energy regressor in the model
    disp('Second SVD complete'); toc
    save([fPath 'absVidR.mat'],'absVidR');
else
    load([fPath 'absVidR.mat'],'absVidR');
end
clear cVideo

%% get ephys data
load([fPath spikeData],'sp');
load([fPath trialData],'sync_data');
if strcmp(trialData, 'NP9_006_sync.mat')
    sync_data.photodiode(57) = [];
end

dataCut = min([round((sp.st(end)-sync_data.photodiode(1)) * sRate) size(vidR,1)]); %last datapoint to be used
spikeCluster = unique(sp.clu);
spikeTrace = zeros(dataCut, size(sp.cgs,2), 'single');
for iUnits = 1 : size(sp.cgs,2)
    cData = sp.st(sp.clu == spikeCluster(iUnits)) - sync_data.photodiode(1);
    cData(cData <= 0) = [];
    cData(cData > dataCut / sRate) = [];
    spikeTrace(:, iUnits) = histcounts(cData, 0 : 1/sRate : dataCut / sRate );
end
spikeTrace = smoothCol(spikeTrace, 5, 'gauss'); %do some smoothing
spikeTrace = bsxfun(@minus, spikeTrace, mean(spikeTrace)); %make zero-mean
save([fPath 'spikeTrace.mat'], 'spikeTrace');

%% cut behavioral regressors to size and make design matrices
vidR = vidR(1:dataCut, :);
absVidR = absVidR(1:dataCut, :);

% keep unfiltered pupil estimate in the model
temp = zeros(dataCut,2,'single');
temp(:,1) = zscore(pupilFill(1:dataCut));
temp(2:end,2) = diff(zscore(pupilFill(1:dataCut)));

slowPupilR = [temp makeDesignMat(slowPupil)];
fastPupilR = makeDesignMat(fastPupil);
noseMotionR = makeDesignMat(smooth(noseMotion, 10));
snoutMotionR = makeDesignMat(smooth(snoutMotion, 10));

%% build stim regressors
stimIDs = unique(sync_data.stimIDs);
stimCnt = length(unique(sync_data.stimIDs));
stimOn = sync_data.photodiode - sync_data.photodiode(1);
stimR = cell(1, stimCnt);

sMat = single(diag(ones(1,(stimDur*sRate))));
for iStim = 1 : stimCnt
    cIdx = sync_data.stimIDs == stimIDs(iStim);
    cIdx = round(stimOn(cIdx) * sRate);
    cIdx(cIdx + ((stimDur-1)*sRate) > size(spikeTrace,1)) = [];
    cIdx(cIdx - (sRate*preStim) < 0) = [];
    
    % create stim regressors
    stimR{iStim} = zeros(size(sMat,1), dataCut, 'single');
    for iRuns = 1 : length(cIdx)
        stimR{iStim}(:, cIdx(iRuns)- (sRate*preStim) : cIdx(iRuns) + ((stimDur-preStim)*sRate)-1) = sMat;
    end
end
stimR = cat(1, stimR{:});
stimR = bsxfun(@minus, stimR, mean(stimR, 2));
stimR = bsxfun(@rdivide, stimR, std(stimR, [], 2))';

%% create full design matrix
fullR = [stimR fastPupilR slowPupilR noseMotionR snoutMotionR absVidR vidR];
fullR = bsxfun(@minus, fullR, mean(fullR, 1));

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = {'stim' 'fastPupil' 'slowPupil' 'nose' 'snout' 'ME' 'video'};

%index to reconstruct different response kernels
regIdx = [
    ones(1,size(stimR,2))*find(ismember(regLabels,'stim')) ...
    ones(1,size(fastPupilR,2))*find(ismember(regLabels,'fastPupil')) ...
    ones(1,size(slowPupilR,2))*find(ismember(regLabels,'slowPupil')) ...
    ones(1,size(noseMotionR,2))*find(ismember(regLabels,'nose')) ...
    ones(1,size(snoutMotionR,2))*find(ismember(regLabels,'snout')) ...
    ones(1,size(absVidR,2))*find(ismember(regLabels,'ME')) ...
    ones(1,size(vidR,2))*find(ismember(regLabels,'video'))];

%% run full model fit for PSTHs
[~, dimBeta] = ridgeMML(spikeTrace, fullR, true); %make model fit
fullFit = fullR * dimBeta; %fit data
stimIdx = regIdx == find(ismember(regLabels,'stim'));
vidFit = fullR(:, ~stimIdx) * dimBeta(~stimIdx, :); %fit data
stimFit = fullR(:, stimIdx) * dimBeta(stimIdx, :); %fit data
save([fPath 'fullFit.mat'], 'fullFit', 'vidFit', 'stimFit');
save([fPath 'regData.mat'], 'fullR', 'dimBeta', 'regLabels', 'regIdx');
disp('Full model fit completed');

%% cross-validated models - create for different regressors
fullV =  crossValModel(regLabels); %cross-validated full model.
allV = NaN(size(fullV,1), size(fullV,2), 2, size(regGroups,2), 'single');
for iRegs = 1 : size(regGroups,2)
    allV(:, :, 1, iRegs) = crossValModel(regGroups{2,iRegs}(:)'); %regressor-alone model
    allV(:, :, 2, iRegs) = crossValModel(regLabels(~ismember(regLabels,regGroups{2,iRegs}(:)'))); %regressor-exclude model
    fprintf('CrossVal for %s complete (%g/%g)\n', regGroups{1,iRegs}, iRegs, size(regGroups,2));
end
save([fPath 'fullV.mat'], 'fullV', 'regLabels');
save([fPath 'allV.mat'], 'allV', 'regGroups');

%% nested functions
    function dataR = makeDesignMat(data)
        data = data(1:dataCut);
        temp = (data - prctile(data,1))./ nanstd(data); %minimum values are at 0, signal in standard deviation units
        [dMat, traceOut] = delayDec_analogToDesign(temp, prctile(temp,60), 1, sRate, sRate, motorIdx, 0);
        dataR = [zscore(single(traceOut)) cat(1,dMat{:})]; %rebuild continuous format
%         temp = [zscore(single(traceOut)) cat(1,dMat{:})]; %rebuild continuous format
%         [dMat, ~] = delayDec_analogToDesign(traceOut, prctile(traceOut,75), 1, sRate, sRate, motorIdx, 0);
%         dataR = [temp cat(1,dMat{:})]; %add high amplitude movements separately
    end

    function [Vm, cBeta, cLabels] =  crossValModel(cLabels)
        rng(1) % for reproducibility
        Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
        randIdx = randperm(size(spikeTrace,1)); %generate randum number index
        foldCnt = floor(size(spikeTrace,1) / ridgeFolds);
        cIdx = ismember(regIdx, find(ismember(regLabels,cLabels))); %get index for task regressors
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            if iFolds == 1
                [cRidge, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true); %get beta weights and ridge penalty for task only model
            else
                [~, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
            end
            Vm(~dataIdx, :) = (fullR(~dataIdx,cIdx) * cBeta); %predict remaining data
            
            if rem(iFolds,ridgeFolds/5) == 0
                fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            end
        end
    end
end

