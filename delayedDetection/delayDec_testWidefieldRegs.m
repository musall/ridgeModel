function delayDec_testWidefieldRegs(cPath,Animal,Rec,nonOrth)
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

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end
cPath = [cPath Animal filesep 'SpatialDisc' filesep Rec filesep]; %Widefield data path
tPath = [cPath 'predVariance' filesep]; %Path to save the results. Make subfolder to get a bit organized.
[~, motorLabels, sensorLabels, cogLabels] = delayDecRecordings;

%% load some data and get into right format
load([cPath 'mask.mat'], 'mask')
load([cPath 'Vc.mat'],'U')
load([cPath 'interpVc.mat'],'Vc','frames') %Vc that was used for the model
dims = size(Vc,1);
U = U(:,:,1:dims);
U = arrayShrink(U,mask);

betaFile = dir([cPath fileExt 'dimBeta.mat']);
regFile = dir([cPath fileExt 'regData.mat']);
load([cPath betaFile.name],'ridgeVals'); %load model ridge penalty
load([cPath regFile.name],'fullR','recIdx','idx','recLabels') %load design matrix
recIdx = recIdx(~idx);

%indices for trial segments. Baseline and handle are based on time regressors, remaining segments are based on stimulus regressors.
segIdx = {1:18 55:72 1:18 34:51 57:74 81:98}; %set index so different segments have same #frames

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
% shCurrent: Shuffle current regressor, shOther: Shuffle all other regressors, 
% shOtherMotor: Shuffle all other motor regressors, shOtherSpontMotor: Shuffle all other spontaneous motor regressors, 
% shTaskOtherSpontMotor: Shuffle task and all other spontaneous movements.
modLabels = {'shCurrent','shOther','shOtherMotor','shOtherSpontMotor','shTaskOtherSpontMotor'}; 
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
        
        %control for shOtherMotor/shOtherSpontMotor condition: Only use for motor regressors.
        checker = true;
        if strcmpi(modLabels{modRuns}, 'shOtherMotor') || strcmpi(modLabels{modRuns}, 'shOtherSpontMotor') || strcmpi(modLabels{modRuns}, 'shTaskOtherSpontMotor')
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
            elseif strcmpi(modLabels{modRuns}, 'shOtherMotor') %this is to shuffle remaining motor regressors
                shIdx = ismember(recIdx, find(ismember(recLabels,taskLabels)));
                shIdx = ~(shIdx | cIdx);
            elseif strcmpi(modLabels{modRuns}, 'shOtherSpontMotor') %this is to the shuffle all spont. motor regressors
                shIdx = ismember(recIdx, find(ismember(recLabels,[taskLabels extraGroups{2,ismember(extraGroups(1,:),{'opMotor'})}]))); %keep task + operant movements
                shIdx = ~(shIdx | cIdx);
            elseif strcmpi(modLabels{modRuns}, 'shTaskOtherSpontMotor') %this is to the shuffle all spont. motor regressors
                shIdx = ismember(recIdx, find(ismember(recLabels,extraGroups{2,ismember(extraGroups(1,:),{'opMotor'})}))); %keep operant movements
                shIdx = ~(shIdx | cIdx);
            end
            
            %shuffle selected regressors
            for iCol = find(shIdx)
                fakeR(:,iCol) = fullR(randperm(size(fullR,1)),iCol);
            end
            
            %% run cross-validation
            Vm = zeros(size(Vc),'single');
            foldCnt = floor(size(Vc,2) / ridgeFolds);
            
            tic
            for iFolds = 1:ridgeFolds
                tic
                dataIdx = true(1,size(Vc,2));
                
                if ridgeFolds > 1
                    
                    dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                    
                    [~, betas] = ridgeMML(Vc(:,dataIdx)', fakeR(dataIdx,:), true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm(:,~dataIdx) = (fakeR(~dataIdx,:) * betas)'; %predict remaining data
                    
                    if rem(iFolds,ridgeFolds/5) == 0
                        fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
                        toc
                    end
                else
                    
                    [~, betas] = ridgeMML(Vc', fakeR, true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                    Vm = (fakeR * betas)'; %predict remaining data
                    disp('Ridgefold is <= 1, fit to complete dataset');
                end
            end
            
            
            %% compute correlation between data and prediction over for all frames
            if ispc
                Vc = gpuArray(Vc);
                Vm = gpuArray(Vm);
                U = gpuArray(U);
            end
            
            Vc = reshape(Vc,size(Vc,1),[]);
            Vm = reshape(Vm,size(Vm,1),[]);
            covVc = cov(Vc');  % S x S
            covVm = cov(Vm');  % S x S
            cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
            covP = sum((U * cCovV) .* U, 2)';  % 1 x P
            varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
            varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
            stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
            cMap = gather((covP ./ stdPxPy)');
            
            
            %% compute correlation between data and prediction over specific trial segments
            stimRegs = sum(fullR(:,ismember(regIdx,find(ismember(regLabels, sensorLabels)))),2); %index for stimulus onset in all trials
            stimRegs = find([0;diff(stimRegs)] > 0.5) - 1;
            
            segMovie = zeros(size(U,1),length(segIdx), 'single');
            for iSegs = 1:length(segIdx)
                if iSegs <= 2 %baseline and handle grab, use trial onset times
                    cIdx = repmat((0:frames:size(Vc,2)-1)',1,length(segIdx{iSegs}));
                else %stimulus and later segments, use stimulus onset times
                    cIdx = repmat(stimRegs,1,length(segIdx{iSegs}));
                end
                cIdx = bsxfun(@plus,cIdx,segIdx{iSegs});
                cIdx = cIdx(:);
                cIdx(cIdx > size(Vc,2)) = [];
                
                covVc = cov(Vc(:,cIdx)');  % S x S
                covVm = cov(Vm(:,cIdx)');  % S x S
                cCovV = bsxfun(@minus, Vm(:,cIdx), mean(Vm(:,cIdx),2)) * Vc(:,cIdx)' / (length(cIdx) - 1);  % S x S
                covP = sum((U * cCovV) .* U, 2)';  % 1 x P
                varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
                varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
                stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
                segMovie(:,iSegs) = gather((covP ./ stdPxPy))';
            end
            
            
            %% compute correlation between data and prediction for each frame in all trials
            cMovie = zeros(size(U,1),frames, 'single');
            for iFrames = 1:frames
                
                frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
                tic
                cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
                cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
                covVc = cov(cData');  % S x S
                covVm = cov(cModel');  % S x S
                cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
                covP = sum((U * cCovV) .* U, 2)';  % 1 x P
                varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
                varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
                stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
                cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
                clear cData cModel
                
                if rem(iFrames,round(frames/4)) == 0
                    fprintf(1, 'Current frame is %d out of %d\n', iFrames,frames);
                    toc
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