function delayDec_test2pRegs(cPath,Animal,Rec,nonOrth,noMove)
% Code to compute predictive power in different regressor or regressor
% groups. This code will create different version of the design matrix,
% where all regressors that correspond to a given variable are shuffled in
% time to assess the importance of that variable for the model.
% There are different shuffling procedures here: 'shCurrent' shuffles all
% regressors of the current variable, 'shOther' shuffles all other
% regressors (resulting in a single variable model), 'shOtherMotor'
% shuffles all regressors that correspond to other movement variables.
% 'shOtherSpontMotor' shuffles all other uninstructed movement variables
% and 'shTaskOtherSpontMotor' shuffles all other uninstructed movement and
% all task variables. The last two modes weren't used for any analysis in
% the paper so subsequent analysis should be ok if they are skipped here.

if ~exist('nonOrth','var') || isempty(nonOrth) || ~nonOrth
    fileExt = ''; %this is to be able to use other kinds of models
else
    fileExt = 'org'; %this is for non-orthogonalized model
end

if ~exist('noMove','var') || isempty(noMove)
    noMove = false; %flag to only use data without XY movement
end
    
if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end
cPath = [cPath Animal filesep 'SpatialDisc' filesep Rec filesep]; %Widefield data path
tPath = [cPath 'predVariance' filesep]; %Path to save the results. Make subfolder to get a bit organized.
load([cPath 'interpVc.mat'],'Vc','DS','frames') %Vc that was used for the model
[dims, times] = size(Vc);

if ~isempty(fileExt)
    betaFile = dir([cPath fileExt '*imBeta.mat']);
    regFile = dir([cPath fileExt '*egData.mat']);
else
    betaFile = dir([cPath fileExt 'dimBeta.mat']);
    regFile = dir([cPath fileExt 'regData.mat']);
end
load([cPath betaFile.name],'ridgeVals'); %load model ridge penalty
load([cPath regFile.name],'fullR','recIdx','idx','recLabels') %load design matrix
 
%index for stimulus onset in all trials. This is needed later
stimRegs = sum(fullR(:,ismember(recIdx,find(ismember(recLabels,{'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'})))),2);
stimRegs = find([0;diff(stimRegs)] > 0.5) - 1;

if noMove
    fileExt = ['noMove' fileExt];
    [b, a] = butter(2, 0.01/31, 'high');
    DS = filtfilt(b, a, DS');
    DS = bsxfun(@minus, DS, median(DS));
    moveIdx = any(DS' > 2); %don't use data with motion above 2 pixels in any direction
else
    moveIdx = false(1,size(Vc,2));
end
Vc(:,moveIdx) = [];
fullR(moveIdx,:) = [];

%remove intercept if present
recIdx = recIdx(~idx);
if sum(fullR(:,1)) == size(fullR,1)
    fullR(:,1) = [];
    recLabels = recLabels(2:end);
    recIdx = recIdx(2:end) - 1;
    idx = idx(2:end);
end

%make sure Vc and fullR are zero-mean
Vc = bsxfun(@minus,Vc',nanmean(Vc, 2)')';
fullR = bsxfun(@minus,fullR,nanmean(fullR));

%indices for trial segments. Baseline and handle are based on time regressors, remaining segments are based on stimulus regressors.
segIdx = {1:19 57:75 1:19 35:53 59:77 84:102}; %set index so different segments have same #frames
[~, motorLabels, sensorLabels, cogLabels] = delayDecRecordings;

% assign extra regressor groups that are removed from model together
extraGroups = {'pupils' 'handles' 'licks' 'allMove' 'opMotor' 'spontMotor' 'motor' 'sensory' 'cognitive'};
extraGroups{2,1} = {'fastPupil' 'slowPupil'};
extraGroups{2,2} = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'};
extraGroups{2,3} = {'lLick' 'rLick'};
extraGroups{2,4} = {'BaselineMove' 'HandleMove' 'StimulusMove' 'WaitMove' 'Move'}; %all motor regressors without video dimensions
extraGroups{2,5} = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'}; %all operant motor regressors
extraGroups{2,6} = motorLabels(~ismember(motorLabels,extraGroups{2,5})); %all spontaneous motor regressors
extraGroups{2,7} = motorLabels; %all motor regressors
extraGroups{2,8} = sensorLabels;
extraGroups{2,9} = cogLabels;

ridgeFolds = 10;    %folds for cross-validation when assessing predicted variance
rng default %reset randum number generator
randIdx = randperm(size(Vc,2)); %generate randum number index if required
modLabels = {'shCurrent','shOther','shTaskCurrent','shOtherMotor'}; %shCurrent: Shuffle current regressor, shOther: Shuffle all other regressors, shTaskCurrent: Shuffle all task and current regressors, shOtherMotor: Shuffle all other motor regressors
taskLabels = [sensorLabels cogLabels]; %these are all task regressors

%motor regressors, used for shOtherMotor analysis
oMotorLabels = ['pupils' motorLabels(~ismember(motorLabels, extraGroups{2, ismember(extraGroups(1,:),'pupils')}))];
oMotorLabels = ['handles' oMotorLabels(~ismember(oMotorLabels, extraGroups{2, ismember(extraGroups(1,:),'handles')}))];
oMotorLabels = ['licks' oMotorLabels(~ismember(oMotorLabels, extraGroups{2, ismember(extraGroups(1,:),'licks')}))];
oMotorLabels = ['allMove' oMotorLabels(~ismember(oMotorLabels, extraGroups{2, ismember(extraGroups(1,:),'allMove')}))];
oMotorLabels = ['opMotor' oMotorLabels];
oMotorLabels = ['spontMotor' oMotorLabels];

%% test predictive power of individual regressors through cross-validation
Cnt = 0;
for iRegs = 0 : length(recLabels) + size(extraGroups,2) %four additional runs to produce motor (with and without video),sensory and cognitive sets
    
    Cnt = Cnt+1;
    fprintf('Current regressor is %d of %d\n', Cnt,length(recLabels) + size(extraGroups,2));
   
    for modRuns = 1:length(modLabels)
        
        %index for current regressor or group
        if iRegs <= length(recLabels)
            cIdx = recIdx == iRegs; % index for reduced model.
        else
            cIdx = ismember(recIdx, find(ismember(recLabels,extraGroups{2,iRegs - length(recLabels)}))); % index for extra regressor group
        end
        
        %control for shOtherMotor condition: Only use for selected subset of motor regressors.
        checker = true;
        if strcmpi(modLabels{modRuns}, 'shOtherMotor')
            if iRegs <= length(recLabels) && iRegs > 0
               checker = any(ismember(recLabels{iRegs},oMotorLabels));
            elseif iRegs > length(recLabels)
               checker = any(ismember(extraGroups{1,iRegs - length(recLabels)},oMotorLabels));
            end
        end
        
        if iRegs == 0
            cLabel = 'full';
        elseif iRegs <= length(recLabels)
            cLabel = recLabels{iRegs};
        else
            cLabel = extraGroups{1,iRegs - length(recLabels)};
        end
        
        if (strcmpi(modLabels{modRuns}, 'shCurrent') || sum(cIdx) > 0) && checker
            
            fakeR = fullR; %copy  design matrix to shuffle up some regressor set
            if strcmpi(modLabels{modRuns}, 'shCurrent') %this is to the shuffle-current regressor
                shIdx = cIdx;
            elseif strcmpi(modLabels{modRuns}, 'shOther') %this is to the shuffle remaining regressors
                shIdx = ~cIdx;
            elseif strcmpi(modLabels{modRuns}, 'shTaskCurrent') %this is to the shuffle current and all task regressors.
                shIdx = ismember(recIdx, find(ismember(recLabels,taskLabels)));
                shIdx = shIdx | cIdx;
            elseif strcmpi(modLabels{modRuns}, 'shOtherMotor') %this is to shuffle remaining motor regressors
                shIdx = ismember(recIdx, find(ismember(recLabels,taskLabels)));
                shIdx = ~(shIdx | cIdx);
            end
            
            %shuffle selected regressors
            for iCol = find(shIdx)
                fakeR(:,iCol) = fullR(randperm(size(fullR,1)),iCol);
            end
            
            %% run cross-validation
            Vm = NaN(size(Vc),'single');
            foldCnt = floor(size(Vc,2) / ridgeFolds);

            tic
            for iFolds = 1:ridgeFolds
                tic
                dataIdx = true(1,size(Vc,2)); %make sure to not use motion frames if those got rejected.
                
                if ridgeFolds > 1
                    if iFolds == ridgeFolds
                        dataIdx(randIdx(((iFolds - 1)*foldCnt) + 1 : end)) = false; %index for training data on the last fold.
                    else
                        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                    end
                    emptyRegs = mean(fakeR(dataIdx,:)) == 0; %check for empty regressors and make sure they are not used
                    
                    [~, betas{iFolds}] = ridgeMML(Vc(:,dataIdx)', fakeR(dataIdx,~emptyRegs), true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm(:,~dataIdx) = (fakeR(~dataIdx,~emptyRegs) * betas{iFolds})'; %predict remaining data
                    
                    if rem(iFolds,ridgeFolds/5) == 0
                        fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
                        toc
                    end
                else
                    [~, betas{iFolds}] = ridgeMML(Vc', fakeR, true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm = (fakeR * betas{iFolds})'; %predict remaining data
                    disp('Ridgefold is <= 1, fit to complete dataset');
                end
            end
            
            %% compute correlation between data and prediction over for all frames / neurons
            cMap = NaN(size(Vc,1), 1);
            for iCells = 1:size(Vc,1)
                cMap(iCells) = corr2(Vc(iCells,~isnan(Vm(iCells,:))), Vm(iCells,~isnan(Vm(iCells,:))));
            end
            
            %% compute correlation between data and prediction over specific trial segments
            segMovie = NaN(size(Vc,1),length(segIdx));
            for iSegs = 1:length(segIdx)
                
                if iSegs <= 2 %baseline and handle grab, use trial onset times
                    cIdx = repmat((0:frames:size(Vc,2)-1)',1,length(segIdx{iSegs}));
                else %stimulus and later segments, use stimulus onset times
                    cIdx = repmat(stimRegs,1,length(segIdx{iSegs}));
                end
                cIdx = bsxfun(@plus,cIdx,segIdx{iSegs});
                cIdx = cIdx(:);
                cIdx(cIdx > times) = [];

                logIdx = false(1, times); logIdx(cIdx) = true; % convert to logical index
                logIdx(moveIdx) = []; %remove movement times if present
                
                for iCells = 1:size(Vc,1)
                    segMovie(iCells, iSegs) = corr2(Vc(iCells,logIdx), Vm(iCells,logIdx));
                end
            end
            
            %% compute correlation between data and prediction for each frame in all trials
            cMovie = zeros(size(Vc,1),frames, 'single');
            for iFrames = 1:frames
                frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
                logIdx = false(1, times); logIdx(frameIdx) = true; % convert to logical index
                logIdx(moveIdx) = []; %remove movement times if present
                
                if sum(logIdx) > 10 %at least 5 trials should be left
                    for iCells = 1:size(Vc,1)
                        cMovie(iCells, iFrames) = corr2(Vc(iCells,logIdx), Vm(iCells,logIdx));
                    end
                else
                    cMovie(:, iFrames) = NaN(size(Vc,1), 1);
                end
            end
            fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));
            
            %% save results           
            if ~exist([tPath modLabels{modRuns}],'dir')
                mkdir([tPath modLabels{modRuns}]);
            end
            
            save([tPath modLabels{modRuns} filesep fileExt cLabel 'corr.mat'], 'cMap', 'cMovie', 'segMovie', 'iRegs', 'recLabels', '-v7.3');
%             save([tPath modLabels{modRuns} filesep fileExt cLabel 'Vm.mat'], 'Vm', 'fakeR', 'betas','-v7.3');
        end
    end
    fprintf('Finished. Current reg: %d\n', Cnt);
end
save([tPath fileExt 'extraGroups.mat'], 'extraGroups'); %save labels for extra groups
save([tPath fileExt 'oMotorLabels.mat'], 'oMotorLabels'); %save other motor labels where some regressors sets are combined