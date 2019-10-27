function delayDec_RegressModel(cPath,Animal,Rec,dType)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('dType', 'var') || isempty(dType)
    dType = 'Widefield'; %default is widefield data
end

if strcmpi(dType,'twoP')
    sRate = 31;        % Sampling rate of imaging in Hz
    piezoLine = 5;     % channel in the analog data that contains data from piezo sensor
    stimLine = 4;      % channel in the analog data that contains stimulus trigger.
    
elseif strcmpi(dType,'Widefield')
    sRate = 30;        % Sampling rate of imaging in Hz
    piezoLine = 2;     % channel in the analog data that contains data from piezo sensor
    stimLine = 6;      % channel in the analog data that contains stimulus trigger.
end

%% general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
% sPath = ['/sonas-hs/churchland/hpc/home/space_managed_data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
preStimDur = ceil(1.8*sRate) / sRate; % Duration of trial before lever grab in seconds
postStimDur = ceil(4.5*sRate) / sRate; % Duration of trial after lever grab onset in seconds
frames = round((preStimDur + postStimDur) * sRate); %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%other variables
mPreTime = ceil(0.5*sRate) / sRate;  % precede motor events to capture preparatory activity in seconds
mPostTime = ceil(2*sRate) / sRate;   % follow motor events for mPostStim in seconds
motorIdx = [-((mPreTime * sRate): -1 : 1) 0 (1:(mPostTime * sRate))]; %index for design matrix to cover pre- and post motor action
tapDur = 0.1;      % minimum time of lever contact, required to count as a proper grab.
leverMoveDur = 0.25; %duration of lever movement. this is used to orthogonalize video against lever movement.
leverMoveDur = ceil(leverMoveDur * sRate); %convert to frames
ridgeFolds = 10;    %number of cross-validations for motor/task models
opMotorLabels = {'lLick' 'rLick' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'}; %operant motor regressors

bhvDimCnt = 200;    % number of dimensions from behavioral videos that are used as regressors.
% dims = 200;         % number of dimensions from V that will be used in the model - will be overwritten if analyzing twoP data instead.
gaussShift = 1;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.
trialSegments = {1:54 55:81 89:132 140:162 170:188}; %trial segments to determine motion.

[~, motorLabels] = delayDecRecordings; %get motor labels for motor-only model

%% load data
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds

if strcmpi(dType,'Widefield')
    if exist([cPath 'Vc.mat'],'file') ~= 2 %check if data file exists and get from server otherwise
        copyfile([sPath 'Vc.mat'],[cPath 'Vc.mat']);
        copyfile([sPath 'mask.mat'],[cPath 'mask.mat']);
        bhvFile = dir([sPath filesep Animal '_' Paradigm '*.mat']);
        copyfile([sPath bhvFile.name],[cPath bhvFile.name]);
    end
    
    load([cPath 'mask.mat'], 'mask')
    load([cPath 'Vc.mat'], 'Vc', 'U', 'trials')
    dims = size(Vc,1);
    
    Vc = Vc(1:dims,:,:);
    U = U(:,:,1:dims);
    
    % ensure there are not too many trials in Vc
    ind = trials > SessionData.nTrials;
    trials(ind) = [];
    Vc(:,:,ind) = [];
    if ~exist('bTrials','var')
        bTrials = trials;
    end

elseif strcmpi(dType,'twoP')
    
    load([cPath 'data'], 'data'); %load 2p data
    % ensure there are not too many trials in the dataset
    bTrials = data.trialNumbers;
    trials = bTrials;
    bTrials(~ismember(data.trialNumbers,data.bhvTrials)) = []; %don't use trials that have problems with trial onset times
    bTrials(SessionData.DidNotChoose(bTrials) | SessionData.DidNotLever(bTrials) | ~SessionData.Assisted(bTrials)) = []; %don't use unperformed/assisted trials
    
    data.dFOF(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.DS(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    data.analog(:,:,~ismember(data.trialNumbers,bTrials)) = [];
    
    Vc = data.dFOF; %Vc is now neurons x frames x trials
    dims = size(data.dFOF,1); %dims is now # of neurons instead
end

bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
trialCnt = length(bTrials);

%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'file') ~= 2 || ... %check if svd behavior exists on hdd and pull from server otherwise
   exist([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'file') ~= 2
    
    if ~exist([cPath 'BehaviorVideo' filesep], 'dir')
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],[cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'segInd1.mat'],[cPath 'BehaviorVideo' filesep 'segInd1.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'segInd2.mat'],[cPath 'BehaviorVideo' filesep 'segInd2.mat']);
    
    movFiles = dir([sPath 'BehaviorVideo' filesep '*Video_*1.mj2']);
    copyfile([sPath 'BehaviorVideo' filesep movFiles(1).name],[cPath 'BehaviorVideo' filesep movFiles(1).name]);
    movFiles = dir([sPath 'BehaviorVideo' filesep '*Video_*2.mj2']);
    copyfile([sPath 'BehaviorVideo' filesep movFiles(1).name],[cPath 'BehaviorVideo' filesep movFiles(1).name]);
    
    svdFiles = dir([sPath 'BehaviorVideo' filesep '*SVD*-Seg*']);
    for iFiles = 1:length(svdFiles)
        copyfile([sPath 'BehaviorVideo' filesep svdFiles(iFiles).name],[cPath 'BehaviorVideo' filesep svdFiles(iFiles).name]);
    end
end

load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'vidV'); %load behavior video data
V1 = vidV(:,1:bhvDimCnt); %behavioral video regressors
load([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'vidV'); %load abs motion video data
V2 = vidV(:,1:bhvDimCnt); % motion regressors

load([cPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'pTime', 'fPupil', 'sPupil', 'whisker', 'faceM', 'bodyM', 'nose', 'bTime'); %load pupil data
%check if timestamps from pupil data are shifted against bhv data
timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
timeCheck2 = (SessionData.TrialStartTime(1)) - (bTime{1}(1)); %time difference between first acquired frame and onset of first trial
if (timeCheck1 > 3590 && timeCheck1 < 3610) && (timeCheck2 > 3590 && timeCheck2 < 3610) %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
        bTime{iTrials} = bTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck1 > 30 || timeCheck1 < -30 || timeCheck2 > 30 || timeCheck2 < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end

if any(bTrials > length(pTime))
    warning(['There are insufficient trials in the pupil data. Rejected the last ' num2str(sum(bTrials > length(pTime))) ' trial(s)']);
    bTrials(bTrials > length(pTime)) = [];
    trialCnt = length(bTrials);
end

%% find events in BPod time - All timestamps are relative to stimulus onset event to synchronize to imaging data later
% pre-allocate vectors
lickL = cell(1,trialCnt);
lickR = cell(1,trialCnt);
leverIn = NaN(1,trialCnt);
levGrabL = cell(1,trialCnt);
levGrabR = cell(1,trialCnt);
levReleaseL = cell(1,trialCnt);
levReleaseR = cell(1,trialCnt);
water = NaN(1,trialCnt);

for iTrials = 1:trialCnt
        
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    
    try
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab(iTrials); %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
    
    %check for spout motion
    if isfield(bhv.RawEvents.Trial{iTrials}.States,'MoveSpout') 
        spoutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab(iTrials);
        
        %also get time when the other spout was moved out at 
        if bhv.Rewarded(iTrials)
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
        else
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.HardPunish(1) - stimGrab(iTrials);
        end
    else
        spoutTime(iTrials) = NaN;
        spoutOutTime(iTrials) = NaN;
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimGrab(iTrials);
    end
    
    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimGrab(iTrials); %first reset state causes lever to move in
        
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimGrab(iTrials);
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimGrab(iTrials);
    end
    
    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
    end
end

maxStimRegs = length(min(round((preStimDur + stimTime) * sRate)) : (preStimDur + postStimDur) * sRate); %maximal number of required stimulus regressors
maxSpoutRegs = length(min(round((preStimDur + spoutTime) * sRate)) : (preStimDur + postStimDur) * sRate); %maximal number of required stimulus regressors

%% build regressors - create design matrix based on event times
%basic time regressors
timeR = logical(diag(ones(1,frames)));

lGrabR = cell(1,trialCnt);
lGrabRelR = cell(1,trialCnt);
rGrabR = cell(1,trialCnt);
rGrabRelR = cell(1,trialCnt);
lLickR = cell(1,trialCnt);
rLickR = cell(1,trialCnt);
leverInR = cell(1,trialCnt);

lVisStimR = cell(1,trialCnt);
rVisStimR = cell(1,trialCnt);
lAudStimR = cell(1,trialCnt);
rAudStimR = cell(1,trialCnt);
spoutR = cell(1,trialCnt);
spoutOutR = cell(1,trialCnt);

rewardR = cell(1,trialCnt);
prevRewardR = cell(1,trialCnt);

ChoiceR = cell(1,trialCnt);

prevChoiceR = cell(1,trialCnt);
prevModR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
fastPupilR = cell(1,trialCnt);
slowPupilR = cell(1,trialCnt);

whiskR = cell(1,trialCnt);
noseR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);
piezoMoveR = cell(1,trialCnt);
faceR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);

%%
tic
for iTrials = 1:trialCnt
    %% vis/aud stim - regressors cover the remaining trial after stimulus onset
    stimIdx = round((preStimDur + stimTime(iTrials)) * sRate) : round((preStimDur + postStimDur) * sRate); %index for which part of the trial should be covered by stim regressors
       
    % vision
    lVisStimR{iTrials} = false(frames, maxStimRegs);
    rVisStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    % audio
    lAudStimR{iTrials} = false(frames, maxStimRegs);
    rAudStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        lVisStimR{iTrials} = lVisStimR{iTrials}(:,1:gaussShift:end);
        rVisStimR{iTrials} = rVisStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
    end
    
    %% spout regressors
    spoutIdx = round((preStimDur + spoutTime(iTrials)) * sRate) : round((preStimDur + postStimDur) * sRate); %index for which part of the trial should be covered by spout regressors
    spoutR{iTrials} = false(frames, maxSpoutRegs);
    spoutR{iTrials}(:, 1:length(spoutIdx)) = timeR(:, spoutIdx);
    
    spoutOutR{iTrials} = false(frames, 3);
    spoutOut = round((preStimDur + spoutOutTime(iTrials)) * sRate); %time when opposing spout moved out again
    if ~isnan(spoutOut) && spoutOut < (frames + 1)
        cInd = spoutOut : spoutOut + 2; cInd(cInd > frames) = [];
        temp = diag(ones(1,3));
        spoutOutR{iTrials}(cInd, :) = temp(1:length(cInd),:);
    end
    
    %% lick regressors
    lLickR{iTrials} = false(frames, length(motorIdx));
    rLickR{iTrials} = false(frames, length(motorIdx));
    
    for iRegs = 0 : length(motorIdx)-1
        licks = lickL{iTrials} - (mPreTime - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - (mPreTime - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    if gaussShift > 1
        % subsample regressors
        lLickR{iTrials} = lLickR{iTrials}(:,1:gaussShift:end);
        rLickR{iTrials} = rLickR{iTrials}(:,1:gaussShift:end);
    end      
    
    %% lever in
    leverInR{iTrials} = false(frames, leverMoveDur);
    leverShift = round((preStimDur + leverIn(iTrials))* sRate); %timepoint in frames when lever moved in, relative to lever grab
    
    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(leverShift : leverShift + leverMoveDur -1, :) = diag(ones(1, leverMoveDur));
        elseif (leverShift + leverMoveDur) > 0  %lever was moving before data was recorded but still moving at trial onset
            leverInR{iTrials}(1 : leverMoveDur + leverShift, :) = [zeros(leverMoveDur + leverShift, abs(leverShift)) diag(ones(1, leverMoveDur + leverShift))];
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end

    %% choice and reward
    stimShift = round((stimTime(iTrials))* sRate) - sRate; %timepoint in frames when the stimulus was presented. This is the shift relative to the expectation that the stimulus comes up 1s after grabing the lever.
    
    stimTemp = false(frames,frames);
    if stimShift > 0 %stim came later than 1s from lever grab
        stimTemp(:, 1: frames - stimShift) = timeR(:, stimShift+1:end);
    else %stim came earlier than 1s from lever grab
        stimTemp(:, abs(stimShift) + 1 : end) = timeR(:, 1: frames + stimShift);
    end
    stimTemp(:,end-4:end) = []; %don't use the last timepoint to avoid rank defficient design matrix

    rewardR{iTrials} = false(size(stimTemp));
    if bhv.Rewarded(iTrials) %rewarded
        rewardR{iTrials} = stimTemp; %trial was rewarded
    end
    
    % get L/R choices as binary design matrix
    ChoiceR{iTrials} = false(size(stimTemp));
    if bhv.ResponseSide(iTrials) == 1
        ChoiceR{iTrials} = stimTemp;
    end
          
    % previous trial regressors
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevModR{iTrials} = NaN(size(timeR(:,1:end-4)));
        
    else %for all subsequent trials
        % same as for regular choice regressors but for prevoious trial
        prevChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)-1) == 1
            prevChoiceR{iTrials} = timeR(:,1:end-4);
        end
             
        prevModR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.StimType(bTrials(iTrials)-1) == 1 || SessionData.StimType(bTrials(iTrials)-1) == 3
            prevModR{iTrials} = timeR(:,1:end-4); % if previous trial was vision
        end
        
        prevRewardR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.Rewarded(bTrials(iTrials)-1) %last trial was rewarded
            prevRewardR{iTrials} = timeR(:,1:end-4);
        end
    end
       
    if gaussShift > 1
        % subsample regressors
        rewardR{iTrials} = rewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);

        ChoiceR{iTrials} = ChoiceR{iTrials}(:,1:gaussShift:end);
        prevChoiceR{iTrials} = prevChoiceR{iTrials}(:,1:gaussShift:end);

        prevModR{iTrials} = prevModR{iTrials}(:,1:Shift:end);
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(frames, sRate);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((preStimDur + water(iTrials)) * sRate); %timepoint in frames when reward was given
        waterR{iTrials}(:, 1: size(timeR,2) - waterOn + 1) = timeR(:, waterOn:end);
    end
    
    if gaussShift > 1
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end); % subsample regressor
    end
    
    %% lever grabs
    cGrabs = levGrabL{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    lGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace
    
    cGrabs = levGrabR{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    rGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace
        
    %% pupil / whisk / nose / face / body regressors
    bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)}))); %framerate of face camera
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab(iTrials) - preStimDur);
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    idx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(idx);
    
    if trialTime(1) > 0 %check if there is missing time at the beginning of a trial
        warning(['Trial ' int2str(bTrials(iTrials)) ': Missing behavioral video frames at trial onset. Trial removed from analysis']);
        fastPupilR{iTrials} = NaN(frames, 1);
        slowPupilR{iTrials} = NaN(frames, 1);
        whiskR{iTrials} = NaN(frames, 1);
        noseR{iTrials} = NaN(frames, 1);
        faceR{iTrials} = NaN(frames, 1);
        bodyR{iTrials} = NaN(frames, 1);
        
    else
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end
        
        fastPupilR{iTrials} = Behavior_vidResamp(fPupil{bTrials(iTrials)}(idx), trialTime, sRate);
        fastPupilR{iTrials} = smooth(fastPupilR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        slowPupilR{iTrials} = Behavior_vidResamp(sPupil{bTrials(iTrials)}(idx), trialTime, sRate);
        slowPupilR{iTrials} =  smooth(slowPupilR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        whiskR{iTrials} = Behavior_vidResamp(whisker{bTrials(iTrials)}(idx), trialTime, sRate);
        whiskR{iTrials} = smooth(whiskR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        noseR{iTrials} = Behavior_vidResamp(nose{bTrials(iTrials)}(idx), trialTime, sRate);
        noseR{iTrials} = smooth(noseR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        faceR{iTrials} = Behavior_vidResamp(faceM{bTrials(iTrials)}(idx), trialTime, sRate);
        faceR{iTrials} = smooth(faceR{iTrials}(end - frames + 1 : end), 'rlowess');
        
        
        % body regressors
        bhvFrameRate = round(1/mean(diff(bTime{bTrials(iTrials)}))); %framerate of body camera
        trialTime = bTime{bTrials(iTrials)} - trialOn;
        idx = trialTime < trialDur; %don't use late frames
        trialTime = trialTime(idx);
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
        
        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end
        
        bodyR{iTrials} = Behavior_vidResamp(bodyM{bTrials(iTrials)}(idx), trialTime, sRate);
        bodyR{iTrials} = smooth(bodyR{iTrials}(end - frames + 1 : end), 'rlowess');
    end
    
    %% piezo sensor information
    if strcmpi(dType,'Widefield')
        if exist([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'file') ~= 2  %check if files exists on hdd and pull from server otherwise
            cFile = dir([sPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
            copyfile([sPath 'Analog_'  num2str(trials(iTrials)) '.dat'],[cPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
        end
        [~,Analog] = Widefield_LoadData([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'Analog'); %load analog data
        stimOn = find(diff(double(Analog(stimLine,:)) > 1500) == 1); %find stimulus onset in current trial
    elseif strcmpi(dType,'twoP')
        Analog = squeeze(data.analog(:,:,iTrials));
        stimOn = find(diff(double(Analog(stimLine,:)) > 1) == 1); %find stimulus onset in current trial
    end
    
    Analog(1,round(stimOn + ((postStimDur-stimTime(iTrials)) * 1000) - 1)) = 0; %make sure there are enough datapoints in analog signal
    temp = Analog(piezoLine,round(stimOn - ((preStimDur + stimTime(iTrials)) * 1000)) : round(stimOn + ((postStimDur - stimTime(iTrials))* 1000) - 1)); % data from piezo sensor. Should encode animals hindlimb motion.
    temp = smooth(double(temp), sRate*5, 'lowess')'; %do some smoothing
    temp = [repmat(temp(1),1,1000) temp repmat(temp(end),1,1000)]; %add some padding on both sides to avoid edge effects when resampling
    temp = resample(double(temp), sRate, 1000); %resample to imaging rate
    piezoR{iTrials} = temp(sRate + 1 : end - sRate)'; %remove padds again
    piezoR{iTrials} = piezoR{iTrials}(end - frames + 1:end); %make sure, the length is correct

    temp = abs(hilbert(diff(piezoR{iTrials})));
    piezoMoveR{iTrials} = [temp(1); temp]; %keep differential motion signal
    clear temp

    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% get proper design matrices for handle grab
lGrabR = cat(1,lGrabR{:});
lGrabR = delayDec_analogToDesign(lGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix

rGrabR = cat(1,rGrabR{:});
rGrabR = delayDec_analogToDesign(rGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix

%% rebuild analog motor regressors to get proper design matrices
temp = double(cat(1,fastPupilR{:}));
temp = (temp - prctile(temp,1))./ nanstd(temp); %minimum values are at 0, signal in standard deviation units
[dMat, traceOut] = delayDec_analogToDesign(temp, median(temp), trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(traceOut, prctile(traceOut,75), trialCnt, sRate, sRate, motorIdx, gaussShift);
fastPupilR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = delayDec_analogToDesign(double(cat(1,whiskR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(double(cat(1,whiskR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
whiskR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = delayDec_analogToDesign(double(cat(1,noseR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(double(cat(1,noseR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
noseR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

% [dMat, traceOut] = delayDec_analogToDesign(double(cat(1,piezoR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
% piezoR1 = [traceOut cat(1,dMat{:})]; %rebuild continuous format
[dMat, traceOut] = delayDec_analogToDesign(double(cat(1,piezoMoveR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [traceOut cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(double(cat(1,piezoMoveR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
piezoR = [temp cat(1,dMat{:})]; %add high amplitude movements separately
% piezoR = [piezoR1 temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = delayDec_analogToDesign(double(cat(1,faceR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(double(cat(1,faceR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
faceR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = delayDec_analogToDesign(double(cat(1,bodyR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = delayDec_analogToDesign(double(cat(1,bodyR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
bodyR = [temp cat(1,dMat{:})]; %add high amplitude movements separately
clear piezoR1 piezoR2 dMat traceOut temp

%% re-align behavioral video data and Vc to lever grab instead of stimulus onset
if strcmpi(dType,'Widefield')
    iiSpikeFrames = findInterictalSpikes(U, Vc, 2, false); %find interictal spikes
    Vc = interpOverInterictal(Vc, iiSpikeFrames); %interpolate over interictal spikes
end

V1 = reshape(V1,205,[],bhvDimCnt); %get to trial format
V2 = reshape(V2,205,[],bhvDimCnt); %get to trial format
vidR = V1(:,bTrials,:); clear V1 %get correct trials from behavioral video data.
moveR = V2(:,bTrials,:); clear V2 %get correct trials from behavioral video data.

% re-align video data
temp1 = NaN(dims,frames,trialCnt);
temp2 = NaN(frames,trialCnt,bhvDimCnt);
temp3 = NaN(frames,trialCnt,bhvDimCnt);
temp4 = NaN(2,frames,trialCnt);
shVal = sRate * 3 + 1; %expected position of stimulus relative to stim onset in s. Usually stimulus onset is after 3s.
for x = 1 : size(vidR,2)
    try
        temp1(:,:,x) = Vc(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        temp2(:,x,:) = vidR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        temp3(:,x,:) = moveR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        if strcmpi(dType,'twoP')
            temp4(:,:,x) = data.DS(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        end
    catch
        fprintf(1,'Could not align trial %d. Relative stim time: %fs\n', x, stimTime(x));
    end
end
Vc = reshape(temp1,dims,[]); clear temp1
vidR = reshape(temp2,[],bhvDimCnt); clear temp2
moveR = reshape(temp3,[],bhvDimCnt); clear temp3

if strcmpi(dType,'twoP')
    DS = reshape(temp4,2,[]); %keep image motion trace for 2p imaging
end
clear temp4

%% reshape regressors, make design matrix and indices for regressors that are used for the model
timeR = repmat(logical(diag(ones(1,frames))),trialCnt,1); %time regressor
timeR = timeR(:,1:end-4);

lGrabR = cat(1,lGrabR{:});
lGrabRelR = cat(1,lGrabRelR{:});
rGrabR = cat(1,rGrabR{:});
rGrabRelR = cat(1,rGrabRelR{:});

lLickR = cat(1,lLickR{:});
rLickR = cat(1,rLickR{:});
leverInR = cat(1,leverInR{:});
leverInR(:,sum(leverInR) == 0) = []; 

lVisStimR = cat(1,lVisStimR{:});
rVisStimR = cat(1,rVisStimR{:});
lAudStimR = cat(1,lAudStimR{:});
rAudStimR = cat(1,rAudStimR{:});
spoutR = cat(1,spoutR{:});
spoutOutR = cat(1,spoutOutR{:});
spoutR(:,sum(spoutR) == 0) = []; 
spoutOutR(:,sum(spoutOutR) == 0) = [];

rewardR = cat(1,rewardR{:});
prevRewardR = cat(1,prevRewardR{:});

ChoiceR = cat(1,ChoiceR{:});

prevChoiceR = cat(1,prevChoiceR{:});
prevModR = cat(1,prevModR{:});

waterR = cat(1,waterR{:});

slowPupilR = cat(1,slowPupilR{:});
slowPupilR(~isnan(slowPupilR(:,1)),:) = zscore(slowPupilR(~isnan(slowPupilR(:,1)),:));

%% compute average motion change in each segment for each trial
vidMove = nanmean(reshape([bodyR(:,1) faceR(:,1)],frames,[],2),3); %reshape to trials and average over face and body camera motion
for x = 1: length(trialSegments)-1 %this should be 4 segments in total (Baseline, Handle, Stim and Waiting period)
   trialMove(x,:) = mean(vidMove(trialSegments{x},:));  %compute average motion change in each segment for each trial
end

for iTrials = 1:size(trialMove,2)
%     BaselineMoveR{iTrials} = diag(repmat(trialMove(1,iTrials),frames,1));
%     HandleMoveR{iTrials} = diag(repmat(trialMove(2,iTrials),frames-trialSegments{2}(1),1),-trialSegments{2}(1));
%     StimulusMoveR{iTrials} = diag(repmat(trialMove(3,iTrials),frames-trialSegments{3}(1),1),-trialSegments{3}(1));
%     WaitMoveR{iTrials} = diag(repmat(trialMove(4,iTrials),frames-trialSegments{4}(1),1),-trialSegments{4}(1));
    
    BaselineMoveR{iTrials} = zeros(frames, 1, 'single');
    HandleMoveR{iTrials} = zeros(frames, 1, 'single');
    StimulusMoveR{iTrials} = zeros(frames, 1, 'single');
    WaitMoveR{iTrials} = zeros(frames, 1, 'single');
end

BaselineMoveR = cat(1,BaselineMoveR{:});
HandleMoveR = cat(1,HandleMoveR{:});
StimulusMoveR = cat(1,StimulusMoveR{:});
WaitMoveR = cat(1,WaitMoveR{:});

%% create full design matrix
fullR = [timeR ChoiceR rewardR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR lVisStimR rVisStimR lAudStimR rAudStimR ...
    prevRewardR prevChoiceR prevModR waterR piezoR whiskR noseR fastPupilR slowPupilR faceR bodyR BaselineMoveR HandleMoveR StimulusMoveR WaitMoveR moveR vidR];

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
recLabels = {
    'time' 'Choice' 'reward' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'lVisStim' 'rVisStim' ...
    'lAudStim' 'rAudStim' 'prevReward' 'prevChoice' 'prevMod' 'water' 'piezo' 'whisk' 'nose' 'fastPupil' 'slowPupil' 'face' 'body' 'BaselineMove' 'HandleMove' 'StimulusMove' 'WaitMove' 'Move' 'bhvVideo'};

%index to reconstruct different response kernels
recIdx = [
    ones(1,size(timeR,2))*find(ismember(recLabels,'time')) ...
    ones(1,size(ChoiceR,2))*find(ismember(recLabels,'Choice')) ...
    ones(1,size(rewardR,2))*find(ismember(recLabels,'reward')) ...
    ones(1,size(lGrabR,2))*find(ismember(recLabels,'lGrab')) ...
    ones(1,size(lGrabRelR,2))*find(ismember(recLabels,'lGrabRel')) ...
    ones(1,size(rGrabR,2))*find(ismember(recLabels,'rGrab')) ...
    ones(1,size(rGrabRelR,2))*find(ismember(recLabels,'rGrabRel')) ...
    ones(1,size(lLickR,2))*find(ismember(recLabels,'lLick')) ...
    ones(1,size(rLickR,2))*find(ismember(recLabels,'rLick')) ...
    ones(1,size(lVisStimR,2))*find(ismember(recLabels,'lVisStim')) ...
    ones(1,size(rVisStimR,2))*find(ismember(recLabels,'rVisStim')) ...
    ones(1,size(lAudStimR,2))*find(ismember(recLabels,'lAudStim')) ...
    ones(1,size(rAudStimR,2))*find(ismember(recLabels,'rAudStim')) ...
    ones(1,size(prevRewardR,2))*find(ismember(recLabels,'prevReward')) ...
    ones(1,size(prevChoiceR,2))*find(ismember(recLabels,'prevChoice')) ...
    ones(1,size(prevModR,2))*find(ismember(recLabels,'prevMod')) ...
    ones(1,size(waterR,2))*find(ismember(recLabels,'water')) ...
    ones(1,size(piezoR,2))*find(ismember(recLabels,'piezo')) ...
    ones(1,size(whiskR,2))*find(ismember(recLabels,'whisk')) ...
    ones(1,size(noseR,2))*find(ismember(recLabels,'nose')) ...
    ones(1,size(fastPupilR,2))*find(ismember(recLabels,'fastPupil')) ...
    ones(1,size(slowPupilR,2))*find(ismember(recLabels,'slowPupil')) ...
    ones(1,size(faceR,2))*find(ismember(recLabels,'face')) ...
    ones(1,size(bodyR,2))*find(ismember(recLabels,'body')) ...
    ones(1,size(BaselineMoveR,2))*find(ismember(recLabels,'BaselineMove')) ...
    ones(1,size(HandleMoveR,2))*find(ismember(recLabels,'HandleMove')) ...
    ones(1,size(StimulusMoveR,2))*find(ismember(recLabels,'StimulusMove')) ...
    ones(1,size(WaitMoveR,2))*find(ismember(recLabels,'WaitMove')) ...
    ones(1,size(moveR,2))*find(ismember(recLabels,'Move')) ...
    ones(1,size(vidR,2))*find(ismember(recLabels,'bhvVideo'))];

% check if enough stimuli of each modality were presented
%vision
cMod = ismember(recIdx,find(ismember(recLabels,{'lVisStim' 'rVisStim'})));
stimCheck = sum(fullR(:,cMod),1);
stimCheck = stimCheck < 10; %if 10 trials or less were presented, dont use them in the analysis
if ~any(~stimCheck)
    fullR = reshape(fullR, frames, [], size(fullR,2));
    fullR(:, bhv.StimType == 1, cMod) = NaN; %mark for deletion
    fullR = reshape(fullR, [], size(fullR,3));
    lVisStimR = []; rVisStimR = [];
end
%audio
cMod = ismember(recIdx,find(ismember(recLabels,{'lAudStim' 'rAudStim'})));
stimCheck = sum(fullR(:,cMod),1);
stimCheck = stimCheck < 10; %if 10 trials or less were presented, dont use them in the analysis
if ~any(~stimCheck)
    fullR = reshape(fullR, frames, [], size(fullR,2));
    fullR(:, bhv.StimType == 2, cMod) = NaN; %mark for deletion
    fullR = reshape(fullR, [], size(fullR,3));
end

% orthogonalize video against spout/handle movement and visual stimuli
if ~isempty(lVisStimR)
    cInd = find(ismember(recIdx, find(ismember(recLabels,'lVisStim'))));
    lVisStimR = fullR(:,cInd([1:ceil(0.6*sRate) ceil(1.1*sRate):ceil(1.7*sRate)])); %regressors for visual stimulus. Only use ON times of stimulus for QR(1.2s with 0.5s gap in between).
    cInd = find(ismember(recIdx, find(ismember(recLabels,'rVisStim'))));
    rVisStimR = fullR(:,cInd([1:ceil(0.6*sRate) ceil(1.1*sRate):ceil(1.7*sRate)])); %regressors for visual stimulus. Only use ON times of stimulus for QR(1.2s with 0.5s gap in between).
end
vidIdx = find(ismember(recIdx, find(ismember(recLabels,{'Move' 'bhvVideo'})))); %index for video regressors
trialIdx = ~isnan(mean(fullR(:,vidIdx),2)); %don't use trials that failed to contain behavioral video data
smallR = [leverInR spoutR spoutOutR lVisStimR rVisStimR];

for iRegs = 1 : length(vidIdx)
    Q = qr([smallR(trialIdx,:) fullR(trialIdx,vidIdx(iRegs))],0); %orthogonalize video against other regressors
    fullR(trialIdx,vidIdx(iRegs)) = Q(:,end); % transfer orthogonolized video regressors back to design matrix
end

% reject trials with broken regressors that contain NaNs
trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fprintf(1, 'Rejected %d/%d trials for NaN entries in regressors\n', sum(trialIdx)/frames,trialCnt);
fullR(trialIdx,:) = []; %clear bad trials

%reject regressors that are too sparse
idx = nansum(abs(fullR)) < 10; 
fullR(:,idx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(idx),length(idx));

%% save modified Vc
Vc(:,trialIdx) = []; %clear bad trials

if strcmpi(dType,'Widefield')
    save([cPath 'interpVc.mat'], 'Vc', 'frames');
elseif strcmpi(dType,'twoP')
    DS(:,trialIdx) = []; %clear bad trials
    save([cPath 'interpVc.mat'], 'Vc', 'DS', 'frames');
end

%% apply gaussian filter to design matrix if using sub-sampling
if gaussShift > 1
    [a,b] = size(fullR);
    
    % find non-continous regressors (contain values different from -1, 0 or 1)
    temp = false(size(fullR));
    temp(fullR(:) ~= 0 & fullR(:) ~= 1 & fullR(:) ~= -1 & ~isnan(fullR(:))) = true;
    regIdx = nanmean(temp) == 0; %index for non-continous regressors
    
    % do gaussian convolution. perform trialwise to avoid overlap across trials.
    trialCnt = a/frames;
    fullR = reshape(fullR,frames,trialCnt,b);
    for iTrials = 1:trialCnt
        fullR(:,iTrials,regIdx) = smoothCol(squeeze(fullR(:,iTrials,regIdx)),gaussShift*2,'gauss');
    end
    fullR = reshape(fullR,a,b);
end

%% clear individual regressors
clear stimR lGrabR lGrabRelR rGrabR rGrabRelR waterR lLickR rLickR ...
    lVisStimR rVisStimR lAudStimR rAudStimR rewardR prevRewardR visChoiceR audChoiceR ...
    prevChoiceR prevModR fastPupilR moveR piezoR whiskR noseR faceR bodyR

%% run QR and check for rank-defficiency
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize design matrix
% figure; plot(abs(diag(fullQRR))); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    error('Design matrix is rank-defficient')
end

%% run ridge regression in low-D
%run model. Zero-mean without intercept. only video qr.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'orgdimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath filesep 'orgregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','fullQRR','-v7.3');
Behavior_betaRebuild(cPath, 'org'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

mInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels(~ismember(motorLabels,opMotorLabels)))));
spontMotorR = fullR(:, mInd);
[spontMotorRidge, spontMotorBeta] = ridgeMML(Vc', spontMotorR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, spont-motor only model: %f\n', mean(spontMotorRidge));
Vm = (spontMotorR * spontMotorBeta)';
save([cPath 'orgVspontMotor.mat'], 'Vm', 'frames'); %save predicted data based on motor model

%% orthogonalize some regressors for clarity
% orthogonalize spontaneous from operant movement regressors
lInd = ismember(recIdx(~idx), find(ismember(recLabels,{'lLick', 'rLick'})));
hInd = ismember(recIdx(~idx), find(ismember(recLabels,{'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'})));
pInd = ismember(recIdx(~idx), find(ismember(recLabels,{'fastPupil', 'slowPupil'})));
wInd = ismember(recIdx(~idx), find(ismember(recLabels,'whisk')));
nInd = ismember(recIdx(~idx), find(ismember(recLabels,'nose')));
piInd = ismember(recIdx(~idx), find(ismember(recLabels,'piezo')));
fInd = ismember(recIdx(~idx), find(ismember(recLabels,'face')));
mInd = ismember(recIdx(~idx), find(ismember(recLabels,'Move')));
vInd = ismember(recIdx(~idx), find(ismember(recLabels,'bhvVideo')));

smallR = [fullR(:,lInd) fullR(:,hInd) fullR(:,pInd) fullR(:,wInd) fullR(:,nInd) fullR(:,piInd) fullR(:,fInd)  fullR(:,mInd) fullR(:,vInd)];
[Q, redQRR] = qr(smallR,0); clear smallR %orthogonalize spont. from operant movement

% replace original with orthogonalized regressors (only for spont. movements)
fullR(:,pInd) = Q(:,sum(lInd | hInd) + 1 : sum(lInd | hInd | pInd)); %pupil
fullR(:,wInd) = Q(:,sum(lInd | hInd | pInd) + 1 : sum(lInd | hInd | pInd | wInd)); %whisk
fullR(:,nInd) = Q(:,sum(lInd | hInd | pInd | wInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd)); %nose
fullR(:,piInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd)); %piezo
fullR(:,fInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd | piInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd)); %face
fullR(:,mInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd | mInd)); %motion energy
fullR(:,vInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd | mInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd | mInd | vInd)); %raw video

%% run model with orthogonalized spontaneous movement regressors. Zero-mean without intercept.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath 'regData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');
Behavior_betaRebuild(cPath); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

% compute cross-validated prediction for full model
[Vm, cBeta] = crossValModel(recLabels);
save([cPath 'interpVfull.mat'], 'Vm', 'frames', 'cBeta'); %save predicted data based on full model

%% run same model without choice.
tIdx = idx; tRecIdx = recIdx; tFullR = fullR;
cIdx = ~ismember(recIdx, find(ismember(recLabels,'Choice'))); %index to exclude choice
idx = idx(cIdx); 
recIdx = recIdx(cIdx);
fullR = fullR(:, cIdx(~tIdx));

[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'noChoicedimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath 'noChoiceregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');
% Behavior_betaRebuild(cPath); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

% run same model without choice and orthogonal to spont movements.
cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels(~ismember(motorLabels,opMotorLabels))))); %find spontaneous movements
smallR = [fullR(:,cInd) fullR(:,~cInd)];
[Q, redQRR] = qr(smallR,0); clear smallR %orthogonalize task vs spont movements
fullR = Q(:,sum(cInd)+1:end); %only use orthogonalized task/op movement regressors
cIdx = ismember(recIdx, find(~ismember(recLabels,motorLabels(~ismember(motorLabels,opMotorLabels))))); %find task/op movements
idx = idx(cIdx); 
recIdx = recIdx(cIdx);

[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'TaskNoChoiceNoSpontdimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath 'TaskNoChoiceNoSpontregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');
idx = tIdx; recIdx = tRecIdx; fullR = tFullR; %recreate full design matrx/indices

%% run motor/task/opMotor and spontMotor only models. Zero-mean without intercept.
cIdx = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels))); %get index for motor regressors
motorLabels = recLabels(sort(find(ismember(recLabels,motorLabels)))); %make sure motorLabels is in the right order

[Vmotor, motorBeta, motorR, motorIdx, motorRidge, motorLabels] = crossValModel(motorLabels);
fprintf('Mean ridge penalty for motor-only, zero-mean model: %f\n', mean(motorRidge));
save([cPath 'interpVmotor.mat'], 'Vmotor', 'frames'); %save predicted data based on motor model
save([cPath 'motorBeta.mat'], 'motorBeta', 'motorRidge');
save([cPath filesep 'motorregData.mat'], 'motorR','trialIdx', 'motorIdx', 'motorLabels','gaussShift','-v7.3');

[Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels] = crossValModel(recLabels(~ismember(recLabels,motorLabels)));
fprintf('Mean ridge penalty for task-only, zero-mean model: %f\n', mean(taskRidge));
save([cPath 'interpVtask.mat'], 'Vtask', 'frames'); %save predicted data based on motor model
save([cPath 'taskBeta.mat'], 'taskBeta', 'taskRidge');
save([cPath filesep 'taskregData.mat'], 'taskR','trialIdx', 'taskIdx', 'taskLabels','gaussShift','-v7.3');

[VopMotor, opMotorBeta, opMotorR, opMotorIdx, opMotorRidge, opMotorLabels] = crossValModel(opMotorLabels);
fprintf('Mean ridge penalty for opMotor-only, zero-mean model: %f\n', mean(opMotorRidge));
save([cPath 'interpVopMotor.mat'], 'VopMotor', 'frames'); %save predicted data based on motor model
save([cPath 'opMotorBeta.mat'], 'opMotorBeta', 'opMotorRidge');
save([cPath filesep 'opMotorregData.mat'], 'opMotorR','trialIdx', 'opMotorIdx', 'opMotorLabels','gaussShift','-v7.3');

[VspontMotor, spontMotorBeta, spontMotorR, spontMotorIdx, spontMotorRidge, spontMotorLabels] = crossValModel(motorLabels(~ismember(motorLabels,opMotorLabels)));
fprintf('Mean ridge penalty for spontMotor-only, zero-mean model: %f\n', mean(spontMotorRidge));
save([cPath 'interpVspontMotor.mat'], 'VspontMotor', 'frames'); %save predicted data based on motor model
save([cPath 'spontMotorBeta.mat'], 'spontMotorBeta', 'spontMotorRidge');
save([cPath filesep 'spontMotorregData.mat'], 'spontMotorR','trialIdx', 'spontMotorIdx', 'spontMotorLabels','gaussShift','-v7.3');

%% run video only model. Zero-mean without intercept.
idx = false(1,size(vidR,2));
recLabels = {'bhvVideo'};
recIdx = ones(1,size(vidR,2));
[ridgeVals, dimBeta] = ridgeMML(Vc', vidR(~trialIdx,:), true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for video-only zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'vidOnlydimBeta.mat'], 'dimBeta', 'ridgeVals');
save([cPath filesep 'vidOnlyregData.mat'], 'vidR', 'idx', 'trialIdx', 'recIdx', 'recLabels','gaussShift','-v7.3');
Behavior_betaRebuild(cPath, 'vidOnly'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

%% re-run model with orthogonalized video. Baseline-corrected with intercept.
Vc = reshape(Vc, dims, frames, []);
temp = squeeze(mean(mean(Vc(:,floor(1:sRate/2),:),3),2));
Vc = reshape(Vc, dims, []);
Vc = bsxfun(@minus, Vc, temp);

[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, false); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for baseline-corrected model: %f\n', mean(ridgeVals));
save([cPath 'offsetdimBeta.mat'], 'dimBeta', 'ridgeVals')

fullR = [ones(size(fullR,1), 1) fullR]; %add intercept to design matrix
recIdx = [1 recIdx+1]; %add zero for intercept
idx = [false idx]; %add zero for intercept
recLabels = ['offset' recLabels];
save([cPath filesep 'offsetregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','-v7.3');


%% nested functions
function [Vm, cBeta, cR, subIdx, cRidge, cLabels] =  crossValModel(cLabels)

cIdx = ismember(recIdx(~idx), find(ismember(recLabels,cLabels))); %get index for task regressors
cLabels = recLabels(sort(find(ismember(recLabels,cLabels)))); %make sure motorLabels is in the right order

%create new regressor index that matches motor labels
subIdx = recIdx(~idx);
subIdx = subIdx(cIdx);
temp = unique(subIdx);
for x = 1 : length(temp)
    subIdx(subIdx == temp(x)) = x;
end
cR = fullR(:,cIdx);

Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
randIdx = randperm(size(Vc,2)); %generate randum number index
foldCnt = floor(size(Vc,2) / ridgeFolds);
cBeta = cell(1,ridgeFolds);

for iFolds = 1:ridgeFolds
    dataIdx = true(1,size(Vc,2));
    
    if ridgeFolds > 1
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        if iFolds == 1
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
        else
            [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
        end
        Vm(:,~dataIdx) = (cR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
        end
    else
        [cRidge, cBeta{iFolds}] = ridgeMML(Vc', cR, true); %get beta weights for task-only model. 
        Vm = (cR * cBeta{iFolds})'; %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset instead');
    end
end
end
end

