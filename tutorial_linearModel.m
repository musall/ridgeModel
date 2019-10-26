% Demo code for how to use the linear encoding model used in the study 
% ‘Single-trial neural dynamics are dominated by richly varied movements’ 
% by Musall, Kaufman et al., 2019. 
% 
% This code shows how to build a design matrix based on task and movement
% events, runs the linear model, shows how to analyze the fitted beta weight
% and quantify cross-validated explained variance
%
% To run, download the widefield dataset from the paper from the CSHL repository, 
% under http://repository.cshl.edu/38599/. Place the data folder 'Widefield' in
% the same location as the folder 'MatlabCode', navigate to this location
% in Matlab and run the following lines to add the code to your path:

codePath = [pwd filesep 'MatlabCode' filesep]; %path to demo codes
addpath(codePath)

%% get some data from example recording
animal = 'mSM43'; %example animal
rec = '23-Nov-2017'; %example recording
fPath = [pwd filesep 'Widefield' filesep animal filesep rec filesep]; %path to demo recording

load([fPath 'opts2.mat'], 'opts'); % load some options
load('allenDorsalMapSM.mat', 'dorsalMaps'); % load allen atlas info
load([fPath 'Vc.mat'], 'U'); %load spatial components
load([fPath 'interpVc.mat'], 'Vc', 'frames'); %load adjusted temporal components
mask = squeeze(isnan(U(:,:,1)));
allenMask = dorsalMaps.allenMask;

%load behavioral data
bhvFile = dir([fPath '*SpatialDisc*.mat']);
load([fPath bhvFile.name]);

%% asign some basic options for the model
% There are three event types when building the design matrix.
% Event type 1 will span the rest of the current trial. Type 2 will span
% frames according to sPostTime. Type 3 will span frames before the event
% according to mPreTime and frames after the event according to mPostTime.

opts.sPostTime = ceil(6 * opts.frameRate);   % follow stim events for sPostStim in frames (used for eventType 2)
opts.mPreTime = ceil(0.5 * opts.frameRate);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
opts.mPostTime = ceil(2 * opts.frameRate);   % follow motor events for mPostStim in frames (used for eventType 3)
opts.framesPerTrial = frames; % nr. of frames per trial
opts.folds = 10; %nr of folds for cross-validation

%% get some events
load([fPath 'orgRegData.mat'], 'fullR', 'recLabels', 'recIdx', 'idx'); %load design matrix to isolate example events and video data (refer to the code 'delayDec_RegressModel' to see how this was generated in the paper)
recIdx(idx) = []; %reject non-used regressors
vidR = fullR(:,end-399:end); %last 400 PCs are video components

% task events
taskLabels = {'time' 'lVisStim' 'rVisStim' 'Choice' 'prevReward'}; %some task variables
taskEventType = [1 2 2 1 1]; %different type of events.
taskEvents(:,1) = fullR(:,find(recIdx == find(ismember(recLabels, taskLabels(1))),1)); %find time regressor. This happens every first frame in every trial.
taskEvents(:,2) = fullR(:,find(recIdx == find(ismember(recLabels,taskLabels(2))),1)); %find event regressor for left visual stimulus
taskEvents(:,3) = fullR(:,find(recIdx == find(ismember(recLabels,taskLabels(3))),1)); %find event regressor for right visual stimulus
taskEvents(:,4) = fullR(:,find(recIdx == find(ismember(recLabels,taskLabels(4))),1)); %find choice reressor. This is true when the animal responded on the left.
taskEvents(:,5) = fullR(:,find(recIdx == find(ismember(recLabels,taskLabels(5))),1)); %find previous reward regressor. This is true when previous trial was rewarded.

% movement events
moveLabels = {'lGrab' 'rGrab' 'lLick' 'rLick' 'nose' 'whisk'}; %some movement variables
moveEventType = [3 3 3 3 3 3]; %different type of events. these are all peri-event variables.
for x = 1 : length(moveLabels)
    moveEvents(:,x) = fullR(:,find(recIdx == find(ismember(recLabels, moveLabels(x))),1)+15); %find movement regressor.
end
clear fullR %clear old design matrix

% make design matrix
[taskR, taskIdx] = makeDesignMatrix(taskEvents, taskEventType, opts); %make design matrix for task variables
[moveR, moveIdx] = makeDesignMatrix(moveEvents, moveEventType, opts); %make design matrix for movement variables

fullR = [taskR, moveR, vidR]; %make new, single design matrix
moveLabels = [moveLabels, {'video'}];
regIdx = [taskIdx; moveIdx + max(taskIdx); repmat(max(moveIdx)+max(taskIdx)+1, size(vidR,2), 1)]; %regressor index
regLabels = [taskLabels, moveLabels];

%% run QR and check for rank-defficiency
rejIdx = false(1,size(fullR,2));
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end
% save([fPath filesep 'regData.mat'], 'fullR', 'regIdx', 'regLabels','fullQRR','-v7.3'); %save some model variables

%% fit model to imaging data
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
% save([fPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals'); %save beta kernels

%reconstruct imaging data and compute R^2
Vm = (fullR * dimBeta)';
corrMat = modelCorr(Vc,Vm,U) .^2; %compute explained variance
corrMat = arrayShrink(corrMat,mask,'split'); %recreate full frame
corrMat = alignAllenTransIm(corrMat,opts.transParams); %align to allen atlas
corrMat = corrMat(:, 1:size(allenMask,2));

%% check beta kernels
% select variable of interest. Must be included in 'regLabels'.
cVar = 'rVisStim';
% cVar = 'rGrab';
% cVar = 'whisk';

% find beta weights for current variable
cIdx = regIdx == find(ismember(regLabels,cVar));
U = reshape(U, [], size(Vc,1)); 
cBeta = U * dimBeta(cIdx, :)';
cBeta = reshape(cBeta, size(mask,1), size(mask,2), []); 
U = reshape(U, size(mask,1), size(mask,2), size(Vc,1)); 
compareMovie(cBeta)

%% run cross-validation
%full model - this will take a moment
[Vfull, fullBeta, ~, fullIdx, fullRidge, fullLabels] = crossValModel(fullR, Vc, regLabels, regIdx, regLabels, opts.folds);
save([fPath 'cvFull.mat'], 'Vfull', 'fullBeta', 'fullR', 'fullIdx', 'fullRidge', 'fullLabels'); %save some results

fullMat = modelCorr(Vc,Vfull,U) .^2; %compute explained variance
fullMat = arrayShrink(fullMat,mask,'split'); %recreate full frame
fullMat = alignAllenTransIm(fullMat,opts.transParams); %align to allen atlas
fullMat = fullMat(:, 1:size(allenMask,2));

%task model alone - this will take a moment
[Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels] = crossValModel(fullR, Vc, taskLabels, regIdx, regLabels, opts.folds);
save([fPath 'cvTask.mat'], 'Vtask', 'taskBeta', 'taskR', 'taskIdx', 'taskRidge', 'taskLabels'); %save some results

taskMat = modelCorr(Vc,Vtask,U) .^2; %compute explained variance
taskMat = arrayShrink(taskMat,mask,'split'); %recreate task frame
taskMat = alignAllenTransIm(taskMat,opts.transParams); %align to allen atlas
taskMat = taskMat(:, 1:size(allenMask,2));

%movement model alone - this will take a moment
[Vmove, moveBeta, moveR, moveIdx, moveRidge, moveLabels] = crossValModel(fullR, Vc, moveLabels, regIdx, regLabels, opts.folds);
save([fPath 'cvMove.mat'], 'Vmove', 'moveBeta', 'moveR', 'moveIdx', 'moveRidge', 'moveLabels'); %save some results

moveMat = modelCorr(Vc,Vmove,U) .^2; %compute explained variance
moveMat = arrayShrink(moveMat,mask,'split'); %recreate move frame
moveMat = alignAllenTransIm(moveMat,opts.transParams); %align to allen atlas
moveMat = moveMat(:, 1:size(allenMask,2));

%% show R^2 results
%cross-validated R^2
figure;
subplot(1,3,1);
mapImg = imshow(fullMat,[0 0.75]);
colormap(mapImg.Parent,'inferno'); axis image; title('cVR^2 - Full model');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        
subplot(1,3,2);
mapImg = imshow(taskMat,[0 0.75]);
colormap(mapImg.Parent,'inferno'); axis image; title('cVR^2 - Task model');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
  
subplot(1,3,3);
mapImg = imshow(moveMat,[0 0.75]);
colormap(mapImg.Parent,'inferno'); axis image; title('cVR^2 - Movement model');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
  
%unique R^2
figure;
subplot(1,2,1);
mapImg = imshow(fullMat - moveMat,[0 0.4]);
colormap(mapImg.Parent,'inferno'); axis image; title('deltaR^2 - Task model');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.

subplot(1,2,2);
mapImg = imshow(fullMat - taskMat,[0 0.4]);
colormap(mapImg.Parent,'inferno'); axis image; title('deltaR^2 - Movement model');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.

