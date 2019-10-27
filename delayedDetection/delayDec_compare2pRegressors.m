reload = false;
fileExt = [];
cMod = 'all';
cPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\NNpaper_Dataset\2pData\Animals\';
sPath = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\NNpaper_Dataset\2pData\Animals\'; %local data path to save down results

%% select data sets
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = twoP_delayDecRecordings;
taskLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts
modTypes = {'shCurrent' 'shOther' 'shTaskCurrent' 'shOtherMotor'};
orgVars = {'fastPupil', 'slowPupil', 'whisk', 'nose', 'face', 'lLick', 'rLick', ...
    'pupils', 'licks', 'bhvVideo', 'Move', 'allMove', 'motorNoVideo', 'motor'}; %make sure to load original version during shOther condition to avoid false results from orthogonalization.

%% load data
check = true;
if ~reload %check for a saved dataset
    if ~exist(sPath,'dir')
        check = false;
    else
        try
            load([sPath fileExt 'corrMaps2p.mat'],'corrMaps');
            load([sPath fileExt 'meanRsq2p.mat'],'aRsq','aTimeRsq','nRecLabels');
            load([sPath fileExt 'segMovies2p.mat'],'segMovies');
            load([sPath fileExt 'fullCorrMaps2p.mat'],'fullCorrMaps', 'recSite', 'recDepth', 'recAnimal', 'recDate');
            load([sPath fileExt 'fullSegMovies2p.mat'],'fullSegMovies');
            load([sPath fileExt 'otherMsegMovies2p.mat'],'otherMsegMovies');
            load([sPath fileExt 'otherMcorrMaps2p.mat'],'otherMcorrMaps');
            load([sPath fileExt 'twoP_rebuild2p.mat'],'recV', 'semV', 'allBeta', 'recLabels', 'allRecIdx', 'dataPath', 'modIdx', 'sideIdx', 'alignIdx', 'baseLength', 'frames', 'stimTimes', 'cellCnt', 'otherMotorLabels');
            
        catch
            check = false;
        end
    end
    if ~check
        reload = true;
        disp('Could not find saved results. Load predVariance data instead. Maybe go get a coffe in the meantime.')
    end
end

%% go through predVariance folders if reloading that data
if reload
    animals = dataOverview(:,1);
    Cnt = 0;
    
    % check animal count
    if strcmpi(cMod,'all')
        animalCnt = size(dataOverview,1);
    else
        animalCnt = sum(ismember(dataOverview(:,2)',cMod));
    end
    
    cellCnt = 0;
    for iAnimals = 1:length(animals)
        regCnt = 0;
        
        if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'all')
            Cnt = Cnt + 1;
            fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep dataOverview{iAnimals,3} filesep];
            disp(fPath); 
            tic;
            try
            %load design matrix and check for regressors that were used in the model
            load([fPath 'regData.mat'],'recIdx','idx','recLabels')
            usedR{Cnt} = recLabels(unique(recIdx(~idx))); %regressors that were used in the model
            
            if Cnt == 1
                % get labels for all regressors
                load([fPath 'predVariance' filesep modTypes{1} filesep fileExt 'fullcorr.mat'],'recLabels'); %load current labels
                load([fPath 'predVariance' filesep 'extraGroups.mat'],'extraGroups'); %load extra group labels
                load([fPath 'predVariance' filesep 'oMotorLabels.mat'],'oMotorLabels'); %load other motor labels
                nRecLabels = [recLabels{:} extraGroups(1,:)]; %labels of all regressors
                nRecLabels = nRecLabels(~strcmpi(nRecLabels,'leverIn')); %don't use leverIn regressor
                
                % make recLabels where modality is swapped for expertise
                altRecLabels = nRecLabels;
                altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
                altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};
            end
            
            for iRuns = 1:length(modTypes)
                
                if iRuns == 1 %load full model data only once - running through different versions makes no difference here
                    load([fPath 'predVariance' filesep modTypes{iRuns} filesep fileExt 'fullcorr.mat'],'cMap','segMovie','recLabels'); %load current data
                    fprintf('%s_%d_%d\n',fPath,iAnimals,length(recLabels));
                    cMap = cMap.^2;
                    segMovie = segMovie.^2;
                    
                    fullCorrMaps(cellCnt + 1:cellCnt + size(cMap,1)) = cMap;
                    aRsq(1,Cnt,1) = nanmean(cMap(:));
                    
                    recSite(cellCnt + 1:cellCnt + size(cMap,1)) = repmat({dataOverview{iAnimals,7}}, 1, size(cMap,1)); %keep track of recording site
                    recDepth(cellCnt + 1:cellCnt + size(cMap,1)) = dataOverview{iAnimals,8}; %keep track of recording depth
                    recAnimal(cellCnt + 1:cellCnt + size(cMap,1)) = repmat({dataOverview{iAnimals,1}}, 1, size(cMap,1)); %keep track of recording depth
                    recDate(cellCnt + 1:cellCnt + size(cMap,1)) = repmat({dataOverview{iAnimals,3}}, 1, size(cMap,1)); %keep track of recording depth
                    if strcmp(dataOverview{iAnimals,2},'Visual')
                        recExp(cellCnt + 1:cellCnt + size(cMap,1)) = repmat('v', 1, size(cMap,1)); %keep track of visual expertise
                    else
                        recExp(cellCnt + 1:cellCnt + size(cMap,1)) = repmat('a', 1, size(cMap,1)); %keep track of auditory expertise
                    end
                    
                    % align segMovie to allen coordinates
                    fullSegMovies(cellCnt + 1:cellCnt + size(cMap,1), :) = segMovie;
                    aTimeRsq(1,:,Cnt,1) = nanmean(segMovie);
                    
                    load([fPath 'data.mat'],'data'); %load current data
                    cellRej(cellCnt + 1:cellCnt + size(cMap,1)) = data.rejIdx;
                    trialCnt(cellCnt + 1:cellCnt + size(cMap,1)) = length(data.bhvTrials);
                end
                
                for iRegs = 1:length(nRecLabels)
                    try
                        % get data from individual regressors
                        if ~ismember(nRecLabels{iRegs}, [usedR{Cnt} extraGroups(1,:)]) && strcmpi(modTypes{iRuns},'shOther') % if regressor was not used and other regressors were shuffled, set all results to zero
                            load([fPath 'predVariance' filesep 'shCurrent' filesep fileExt 'fullcorr.mat'],'cMap','segMovie'); %load full model instead of current regressor data
                            cMap = cMap - cMap; %set to 0
                            segMovie = segMovie - segMovie; %set to 0
                        elseif ismember(nRecLabels{iRegs}, orgVars) && strcmpi(modTypes{iRuns},'shOther') % use original for standalone video regressors
                            load([fPath 'predVariance' filesep modTypes{iRuns} filesep fileExt 'org' nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                        else
                            load([fPath 'predVariance' filesep modTypes{iRuns} filesep fileExt nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                        end
                        cMap = cMap.^2;
                        segMovie = segMovie.^2;
                    
                        cInd = iRegs;
                        if strcmpi(nRecLabels{iRegs},'visReward') || strcmpi(nRecLabels{iRegs},'audReward')
                            if (visExp(iAnimals) && strcmpi(nRecLabels{iRegs},'visReward')) || (audExp(iAnimals) && strcmpi(nRecLabels{iRegs},'audReward'))
                                cInd = find(ismember(altRecLabels, 'expReward')); %change index to point to expert instead of modality
                            else
                                cInd = find(ismember(altRecLabels, 'novReward')); %change index to point to novice instead of modality
                            end
                        end
                        
                        aTimeRsq(cInd + 1, :, iRuns) = nanmean(segMovie);
                        aRsq(cInd + 1, iRuns) = nanmean(cMap(:));
                        
                        if ~strcmpi(modTypes{iRuns},'shOtherMotor')
                            segMovies(cellCnt + 1:cellCnt + size(cMap,1), cInd, :, iRuns) = segMovie;
                            corrMaps(cellCnt + 1:cellCnt + size(cMap,1), cInd, iRuns) = cMap;
                        else
                            regCnt = regCnt + 1;
                            otherMsegMovies(cellCnt + 1:cellCnt + size(cMap,1), regCnt, :) = segMovie;
                            otherMcorrMaps(cellCnt + 1:cellCnt + size(cMap,1), regCnt) = cMap;
                            otherMotorLabels{regCnt} = nRecLabels{iRegs};
                        end
                    end
                end
            end
            cellCnt = cellCnt + size(cMap,1);
            cellPerAnimal(Cnt) = size(cMap,1);
            catch
                Cnt = Cnt - 1;
                disp('Cant load recording');
            end
            toc;
        end
        
        % check if shOtherMotor was loaded correctly
        if length(oMotorLabels) ~= sum(ismember(oMotorLabels,otherMotorLabels))
            error('oMotorLabels and otherMotorLabels are not the same. Problem with loading shOtherMotor variables?')
        end
    end
    clear segMovie cMap
        
    %% load reconstructed PSTH data
    [recV, semV, allBeta, recLabels, dataPath, modIdx, sideIdx, alignIdx, allRecIdx, baseLength, frames, stimTimes, cellCnt] = twoP_motorReconstruct(cMod); %get reconstructed V, used full model
    recV = cat(4,recV{:});
    semV = cat(4,semV{:});
    stimTimes = cat(2,stimTimes{:});
    stimTimes(stimTimes > 1.5) = [];
    recLabels = [recLabels {'task' 'motor' 'all'}];
    
    %% save down results
    if ~exist(sPath,'dir')
        mkdir(sPath);
    end
    save([sPath fileExt 'fullCorrMaps2p.mat'],'fullCorrMaps', 'recSite', 'recDepth', 'recDate', 'recAnimal');
    save([sPath fileExt 'fullSegMovies2p.mat'],'fullSegMovies');
    save([sPath fileExt 'corrMaps2p.mat'],'corrMaps');
    save([sPath fileExt 'segMovies2p.mat'],'segMovies', '-v7.3');
    save([sPath fileExt 'otherMsegMovies2p.mat'],'otherMsegMovies', '-v7.3');
    save([sPath fileExt 'otherMcorrMaps2p.mat'],'otherMcorrMaps');
    save([sPath fileExt 'meanRsq2p.mat'],'aRsq','aTimeRsq','nRecLabels');
    save([sPath fileExt 'twoP_rebuild2p.mat'],'recV', 'semV', 'allBeta', 'recLabels', 'allRecIdx', 'dataPath', 'modIdx', 'sideIdx', 'alignIdx', 'baseLength', 'frames', 'stimTimes', 'cellCnt', 'otherMotorLabels', '-v7.3');
end

%% make overview figure for overall task- and motor explained variance
rejIdx = (sum(abs(recV(:, :, 6, end)),2) > std(sum(abs(recV(:, :, 6, end)),2))*1)'; %exclude neurons with weird PSTH

% taskFull = corrMaps(:,strcmpi(nRecLabels,{'motor'}),strcmpi(modTypes,'shCurrent'))';
% spMotorFull = corrMaps(:,strcmpi(nRecLabels,{'spontMotor'}),strcmpi(modTypes,'shOther'))';
% opMotorFull = corrMaps(:,strcmpi(nRecLabels,{'opMotor'}),strcmpi(modTypes,'shOther'))';

taskRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'motor'}),strcmpi(modTypes,'shOther'))';
spMotorRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'spontMotor'}),strcmpi(modTypes,'shCurrent'))';
opMotorRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'opMotor'}),strcmpi(modTypes,'shCurrent'))';

% find cells that may suffer from overfitting and remove from analysis.
% Also remove cells with weird PSTH reconstruction (only 1 so far).
cIdx = spMotorRed < 0  | opMotorRed < 0 | taskRed < 0 | sum(abs(recV(:, :, 7, 1)), 2)' > 100;
fullCorrMaps(cIdx) = [];
corrMaps(cIdx, :, :) = [];
recSite = recSite(~cIdx);
recAnimal = recAnimal(~cIdx);
recDate = recDate(~cIdx);
recDepth(cIdx) = [];
recV(cIdx, :, :, :) = [];
fprintf('Removed %d/%d (%f) cells for having negative task/motor information\n', sum(cIdx), length(cIdx), 100*sum(cIdx)/length(cIdx))

taskFull = corrMaps(:,strcmpi(nRecLabels,{'motor'}),strcmpi(modTypes,'shCurrent'))';
spMotorFull = corrMaps(:,strcmpi(nRecLabels,{'spontMotor'}),strcmpi(modTypes,'shOther'))';
opMotorFull = corrMaps(:,strcmpi(nRecLabels,{'opMotor'}),strcmpi(modTypes,'shOther'))';

taskRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'motor'}),strcmpi(modTypes,'shOther'))';
spMotorRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'spontMotor'}),strcmpi(modTypes,'shCurrent'))';
opMotorRed = fullCorrMaps - corrMaps(:,strcmpi(nRecLabels,{'opMotor'}),strcmpi(modTypes,'shCurrent'))';

[~, cellIdx] = sort(fullCorrMaps, 'descend');
% [~, cellIdx] = sort(taskFull, 'descend');
% cellIdx = 1:cellCnt;

cellIdx(isnan(fullCorrMaps(cellIdx))) = [];

figure; hold on;
smoothFact = 3;
plot(smooth(fullCorrMaps(cellIdx),smoothFact),'color', [0.5 0.5 0.5], 'linewidth',4); axis square;
plot(smooth(spMotorFull(cellIdx),smoothFact),'k'); axis square;
plot(smooth(opMotorFull(cellIdx),smoothFact),'b'); axis square;
plot(smooth(taskFull(cellIdx),smoothFact),'g'); axis square;

plot(smooth(-spMotorRed(cellIdx),smoothFact),'k'); axis square;
plot(smooth(-opMotorRed(cellIdx),smoothFact),'b'); 
plot(smooth(-taskRed(cellIdx),smoothFact),'g'); axis square;

axis square;hline(0);
ylabel('crossVal R^2');
xlabel('Neurons')
xlim([0 length(cellIdx)]);


%% relative to full model plot
figure; hold on;
smoothFact = 5;
fullData = smooth(fullCorrMaps(cellIdx),smoothFact);
plot(smooth(opMotorFull(cellIdx),smoothFact) ./ fullData,'b'); axis square;
plot(smooth(taskFull(cellIdx),smoothFact) ./ fullData,'g'); axis square;
plot(smooth(spMotorFull(cellIdx),smoothFact) ./ fullData,'k'); axis square;
axis square;

axis square;hline(0);
ylabel('Relative cvR^2');
xlabel('Neurons')
% xlim([0 length(cellIdx)]);
xlim([0 10000]);

%% Task-dependent R2 reduction for all regressors
clear cData
opMotorLabels = {'handles' 'licks'}; %operant motor regressors
spMotorLabels = {'piezo' 'whisk' 'nose' 'body' 'bhvVideo' 'allMove' 'pupils'}; %spontaneous motor regressors
% cIdx = ismember(recSite,'RS');
cIdx = true(1, length(recSite));

clear cData cRegMod cRegRedM
for x = 1:length(nRecLabels)
    
    cReg = corrMaps(cIdx,x,ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
    cRegM = fullCorrMaps(cIdx)' - corrMaps(cIdx,x,ismember(modTypes,'shCurrent')); %this is the current regressor. This has non-redundant information.
    
    %normalize for animals
    animals = unique(recAnimal);
    for iAnimals = 1 : length(animals)
        aIdx = ismember(recAnimal,animals{iAnimals});
        cRegMod(iAnimals) = mean(cReg(aIdx));
        cRegRedM(iAnimals) = mean(cRegM(aIdx));
    end
    
    cData(1,x,:) = cRegMod;
    cData(2,x,:) = cRegRedM;
end

idx = zeros(1,length(nRecLabels));
idx(ismember(nRecLabels,opMotorLabels)) = 3;
idx(ismember(nRecLabels,cogLabels)) = 2;
idx(ismember(nRecLabels,sensorLabels)) = 1;
idx(ismember(nRecLabels,spMotorLabels)) = 4;

figure('renderer','painters')
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,idx>0,:))',nRecLabels(idx>0),5,subplot(2,2,1:2),[0 1 0],idx(idx>0),0.6,[-0.2 0.4]);
ax = regressorBoxPlot(squeeze(-cData(2,idx>0,:))',nRecLabels(idx>0),5,ax,[25 111 61]/255,idxGroup,0.6,[-0.2 0.4]);
legend(ax.Children(1:2),{'All information' 'Non-redundant information'},'location','northwest')
title('Single varibles: All (top) / Unique (bottom)');
ylabel('unique delta cv R^2');
xlim([0 25]); %make this a fixed value to ensure bars have always the same width

%% R2 reduction for regressor groups
clear cData
motorM = corrMaps(:,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shOther')); %this is the current regressor. This ha sall information.

taskM = corrMaps(:,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shCurrent')); %this is the current regressor. This ha sall information.
taskRedM = fullCorrMaps' - motorM; %this is the non-redundant task only model.

spMotorM = corrMaps(:,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shOther')); %this is the current regressor. This ha sall information.
spMotorRedM = fullCorrMaps' - corrMaps(:,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shCurrent')); %this is the non-redundant task only model.

opMotorM = corrMaps(:,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shOther')); %this is the current regressor. This ha sall information.
opMotorRedM = fullCorrMaps' - corrMaps(:,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shCurrent')); %this is the non-redundant task only model.

clear cData cError
cData(1,:,:) = [taskM spMotorM opMotorM fullCorrMaps']; %all reg information
cData(2,:,:) = [taskRedM spMotorRedM opMotorRedM fullCorrMaps']; %all reg information

figure
[ax, idxGroup] = regressorPlot(squeeze(cData(1,:,:)), {'Task' 'spontMove' 'opMove' 'full'}, 5, subplot(2,2,1:4), [0 1 0], {[1 3 2 4]},0.6,[-0.5 0.5]);
ax = regressorPlot(squeeze(-cData(2,:,:)), {'Task' 'spontMove' 'opMove' 'full'}, 5, ax, [25 111 61]/255, idxGroup,0.6,[-0.5 0.5]);
title('Regg groups: All (top) / Unique (bottom)'); xlim([0 30]); %make this a fixed value to ensure bars have always the same width
ylabel('total cv R^2');
hline(mean(fullCorrMaps),'k--')


%% R2 reduction for regressor groups: area specific
clear areaData areaError
% areas = unique(recSite); %different recording sites
areas = {'ALM' 'MM' 'V1' 'RS' 'S1'}; %different recording sites

for iAreas = 1 : length(areas) + 1
    
    if iAreas > 1
        cIdx = ismember(recSite,areas{iAreas - 1});
    else
        cIdx = true(1, length(recSite));
    end
    
    %cvR^2 - normalize for animals
    animals = unique(recAnimal); clear temp
    for iAnimals = 1 : length(animals)
        aIdx = cIdx & ismember(recAnimal,animals{iAnimals});
        temp(iAnimals,1) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shCurrent'))); %this is the task model
        temp(iAnimals,2) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shOther'))); %this is the operant motor model
        temp(iAnimals,3) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shOther'))); %this is the spont. motor model
        temp(iAnimals,4) = nanmean(fullCorrMaps(aIdx)); %this is the full model
    end
    areaData(iAreas, :, :, 1) = (temp);
    
    %dR^2 - normalize for animals
    animals = unique(recAnimal); clear temp
    for iAnimals = 1 : length(animals)
        aIdx = cIdx & ismember(recAnimal,animals{iAnimals});
        temp(iAnimals,1) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shOther'))); %this is the task model
        temp(iAnimals,2) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shCurrent'))); %this is the operant motor model
        temp(iAnimals,3) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shCurrent'))); %this is the spont. motor model
        temp(iAnimals,4) = nanmean(fullCorrMaps(aIdx)); %this is the full model
    end
    areaData(iAreas, :, :, 2) = (temp);
end

figure('renderer','painters'); 
areaLabels = ['Full' areas];
for x = 1 : size(areaData,1)
    subplot(2,max(size(areaData,1)/2),x);hold on;
    [ax, idxGroup] = regressorBoxPlot(areaData(x,:,:,1), {'Task' 'spontMove' 'opMove' 'full'}, 5, gca, [0 1 0], {[1 2 3 4]},0.6,[-0.5 0.5]);
    ax = regressorBoxPlot(squeeze(-areaData(x,:,:,2)), {'Task' 'spontMove' 'opMove' 'full'}, 5, ax, [25 111 61]/255, idxGroup,0.6,[-0.5 0.5]);
    axis square
    title(areaLabels{x});
end

ngroups = size(areaData, 1);
nbars = size(areaData, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, areaData(:,i,1), areaError(:,i,1), 'k', 'linestyle', 'none','LineWidth',2);
    errorbar(x, -areaData(:,i,2), areaError(:,i,2), 'k', 'linestyle', 'none','LineWidth',2);
end
ax = bar(areaData(:,:,1), 'FaceColor', 'g', 'barwidth', 0.6,'LineWidth',2);
bar(-areaData(:,:,2), 'FaceColor', [25 111 61]/255, 'barwidth', 0.6,'LineWidth',2)
set(ax(1).Parent,'XTick',1:ngroups)
set(ax(1).Parent,'xTickLabel',['Full' areas])
set(ax(1).Parent,'XTickLabelRotation',45)
axis square; grid on; title('Pred.Var / areas');

%% R2 reduction for regressor groups: depth specific
clear areaData areaError
depthLabels = {'Superficial' 'Infragranular'};
depths = [0 350 1000]; %define range of depths to group neurons
for iDepths = 1 : length(depths) - 1
    
    cIdx = recDepth >= depths(iDepths) & recDepth < depths(iDepths + 1);

    %cvR^2 - normalize for animals
    animals = unique(recAnimal); clear temp
    for iAnimals = 1 : length(animals)
        aIdx = cIdx & ismember(recAnimal,animals{iAnimals});
        temp(iAnimals,1) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shCurrent'))); %this is the task model
        temp(iAnimals,2) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shOther'))); %this is the operant motor model
        temp(iAnimals,3) = nanmean(corrMaps(aIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shOther'))); %this is the spont. motor model
        temp(iAnimals,4) = nanmean(fullCorrMaps(aIdx)); %this is the full model
    end
    areaData(iDepths, :, :, 1) = (temp);
    
    %dR^2 - normalize for animals
    animals = unique(recAnimal); clear temp
    for iAnimals = 1 : length(animals)
        aIdx = cIdx & ismember(recAnimal,animals{iAnimals});
        temp(iAnimals,1) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shOther'))); %this is the task model
        temp(iAnimals,2) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shCurrent'))); %this is the operant motor model
        temp(iAnimals,3) = nanmean(fullCorrMaps(aIdx)' - corrMaps(aIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shCurrent'))); %this is the spont. motor model
        temp(iAnimals,4) = nanmean(fullCorrMaps(aIdx)); %this is the full model
    end
    areaData(iDepths, :, :, 2) = (temp);
end

figure('renderer','painters'); 
for x = 1 : size(areaData,1)
    subplot(2,max(size(areaData,1)/2),x);hold on;
    [ax, idxGroup] = regressorBoxPlot(areaData(x,:,:,1), {'Task' 'spontMove' 'opMove' 'full'}, 5, gca, [0 1 0], {[1 2 3 4]},0.6,[-0.5 0.5]);
    ax = regressorBoxPlot(squeeze(-areaData(x,:,:,2)), {'Task' 'spontMove' 'opMove' 'full'}, 5, ax, [25 111 61]/255, idxGroup,0.6,[-0.5 0.5]);
    axis square
    grid on; title(depthLabels{x});
end


%% plot PSTH over all neurons for each area
clear areaData areaError
areas = {'ALM' 'MM' 'V1' 'RS' 'S1'}; %different recording sites

temp = recV(:, 1 : baseLength, 1, end) - recV(:, 1 : baseLength, 2, end);
baseRejIdx = ((abs(mean(recV(:, baseLength:end, 1, end),2)) < 0.01) | rejIdx')';
% baseRejIdx = ((abs(mean(temp,2)) > 0.05) | rejIdx)';
iMod = 1;

figure('units','normalized','outerposition',[0 0 1 1])
for iAreas = 1 : length(areas)
    cIdx = ismember(recSite,areas{iAreas}) ;
    subplot(ceil(length(areas)/2), 2, iAreas);   
    cStim = baseLength+2;
    
    stdshade(recV(cIdx & ~baseRejIdx, 1 : baseLength, iMod, end), 'k', (1 : baseLength) / 31, 0.5); hold on
    stdshade(recV(cIdx & ~baseRejIdx, cStim : end, iMod, end), 'k', (cStim : size(recV,2)) / 31, 0.5);
    
    stdshade(recV(cIdx & ~baseRejIdx, 1 : baseLength, iMod, end-3), 'r', (1 : baseLength) / 31, 0.5); hold on
    stdshade(recV(cIdx & ~baseRejIdx, cStim : end, iMod, end-3), 'r', (cStim : size(recV,2)) / 31, 0.5);
    
    stdshade(recV(cIdx & ~baseRejIdx, 1 : baseLength, iMod, end-2), 'g', (1 : baseLength) / 31, 0.5); hold on
    stdshade(recV(cIdx & ~baseRejIdx, cStim : end, iMod, end-2), 'g', (cStim : size(recV,2)) / 31, 0.5);
    
    stdshade(recV(cIdx & ~baseRejIdx, 1 : baseLength, iMod, end-1), 'b', (1 : baseLength) / 31, 0.5); hold on
    stdshade(recV(cIdx & ~baseRejIdx, cStim : end, iMod, end-1), 'b', (cStim : size(recV,2)) / 31, 0.5);
    
   
    title([areas{iAreas} ' - ' num2str(sum(cIdx & ~baseRejIdx))]); 
    axis square; xlim([0 6]);
    cStim = cStim / 31; vline([1.8 cStim cStim+0.6 cStim+1.1 cStim+1.7 cStim+2.7])
    ylim([-0.02 0.03])
    
end

%%
iMod = 6;
allVar = sum(abs(recV(:, :, iMod, end)),2);
spMotorVar = sum(abs(recV(:, :, iMod, end-1)),2);
opMotorVar = sum(abs(recV(:, :, iMod, end-2)),2);
taskVar = sum(abs(recV(:, :, iMod, end-3)),2);

cIdx = NaN(1, length(allVar));
cIdx(taskVar > opMotorVar & taskVar > spMotorVar) = 3;
cIdx(opMotorVar > taskVar & opMotorVar > spMotorVar) = 2;
cIdx(spMotorVar > taskVar & spMotorVar > opMotorVar) = 1;
cIdx(sum(abs(recV(:, :, iMod, end)),2) > std(sum(abs(recV(:, :, iMod, end)),2))*3) = NaN;

figure;
for x = 1:3
    
    [~, b] = sort(sum(abs(recV(cIdx == x, :, iMod, end-x)),2), 'descend');
    
    subplot(1,3,x);
    temp = sum(abs(recV(cIdx == x, :, iMod, end-3)),2); %task
    plot(temp(b),'r','linewidth',2); hold on
    temp = sum(abs(recV(cIdx == x, :, iMod, end-2)),2); %op motor
    plot(temp(b),'g','linewidth',2); hold on
    temp = sum(abs(recV(cIdx == x, :, iMod, end-1)),2); %sp motor
    plot(temp(b),'k','linewidth',2); hold on
    
end
fprintf('Modulation index: SpMajor: %d, OpMajor: %d, TaskMajor: %d\n', sum(cIdx == 1), sum(cIdx == 2), sum(cIdx == 3));

%% individual regressors / neurons
cReg = 'spontMotor';
clear cData;

cRegMod = corrMaps(:,ismember(nRecLabels,cReg),ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
cRegRedM = fullCorrMaps' - corrMaps(:,ismember(nRecLabels,cReg),ismember(modTypes,'shCurrent')); %this is the current regressor. This has non-redundant information.

cData(1,:) = cRegMod;
cData(2,:) = cRegRedM;

figure;
subplot(1,2,1);
plot(cData');
axis square
title(['Single regressor - ' fileExt]);
ylabel('crossVal R^2');
xlabel('Neurons');
ylim([0 1]);

subplot(1,2,2);
a = histogram(cData(1,:)); hold on
histogram(cData(2,:),a.BinEdges);
axis square
title(['Single regressor - ' fileExt]);
xlabel('crossVal R^2');
ylabel('Neurons')

% %% single cell betas
% lickBetas = cell(1, length(allBeta));
% handleBetas = cell(1, length(allBeta));
% visBetas = cell(1, length(allBeta));
% audBetas = cell(1, length(allBeta));
% for iSessions = 1 : length(allBeta)
%     cBeta = allBeta{iSessions};
%     cLick = mean(abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'lLick')), :)));
%     cLick(2,:) = mean(abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'rLick')), :)));
%     lickBetas{iSessions} = cLick;
%     
%     cHandle = mean(abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'lGrab')), :)));
%     cHandle(2,:) = mean(abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'rGrab')), :)));
%     handleBetas{iSessions} = cHandle;
%     
%     temp = abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'lVisStim')), :));
%     cVisStim = mean(temp(1:50, :));
%     temp = abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'rVisStim')), :));
%     cVisStim(2,:) = mean(temp(1:50, :));
%     visBetas{iSessions} = cVisStim;
%     
%     temp = abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'lAudStim')), :));
%     cAudStim = mean(temp(1:50, :));
%     temp = abs(cBeta(allRecIdx{iSessions} == find(ismember(recLabels, 'rAudStim')), :));
%     cAudStim(2,:) = mean(temp(1:50, :));
%     audBetas{iSessions} = cAudStim;
% end
% lickBetas = cat(2, lickBetas{:});    
% handleBetas = cat(2, handleBetas{:});    
% visBetas = cat(2, visBetas{:});    
% audBetas = cat(2, audBetas{:});   
% 
% %
% figure
% clear areaData depthData
% areas = unique(recSite); %different recording sites
% depths = [0 350 1000]; %define range of depths to group neurons
% plotCnt = 1 : (length(depths) - 1) * length(areas);
% plotCnt = reshape(plotCnt, length(areas), [])';
% Cnt = 0;
% for iAreas = 1 : length(areas)
%     for iDepths = 1 : length(depths) - 1
%         
%         cIdx = ismember(recSite,areas{iAreas}) & recDepth > depths(iDepths) & recDepth < depths(iDepths + 1);
% 
%         %make figure
%         try
%             [ax, idxGroup] = regressorPlot(squeeze(areaData{iAreas, iDepths}(1,:,:)), {'Task' 'Instructed' 'Uninstructed' 'Full'}, 5, subplot(length(depths) - 1, length(areas), plotCnt(Cnt)), [0 1 0], [1 1 1 1], 0.6);
%             regressorPlot(squeeze(-areaData{iAreas, iDepths}(2,:,:)), {'Task' 'Instructed' 'Uninstructed' 'Full'}, 5, ax, [25 111 61]/255, idxGroup, 0.6, [-0.5 0.5]);
%         catch
%             subplot(length(depths) - 1, length(areas), plotCnt(Cnt)); 
%             axis square;
%         end
%         title(['Regressor groups: ' areas{iAreas} ' - ' num2str(depths(iDepths+1)) ' (' num2str(sum(cIdx)) ' cells) - ' fileExt]);
%         ylabel('crossVal R^2'); axis square
%         xlim([0 5]); axis square
%         
%     end
% end


%%
% meanBase = mean((fullRec(:,1:17,7)),2);
% stdBase = mean((fullErr(:,[1:17 169:end],7)),2);
% meanLick = mean((fullRec(:,169:end,7)),2);
% glassLick = (meanLick - meanBase) ./ stdBase;

% figure
% histogram(glassLick,30); axis square
% title(['d" - ' fileExt]);
% ylabel('d"');
% 
% % make figure
% [a, b] = (sort(modIndex,'ascend')); %sort indices
% ind = round(a*200);
% modMat = false(size(modIndex,1),200);
% 
% for x = 1 : length(ind)
%     modMat(x, 1:ind(x)) = true;
% end
% 
% figure
% imagesc((modMat'));
% caxis([-0.5 1.5]);
% hline(100,'--w');
% vline(find(ind>100,1),'--w');
% colormap(colormap_blueblackred);
% axis square
% xlabel('Sorted neurons');
% ylabel('Modulation index'); hold on
% plot((ind),'w')
% 
% figure;
% cData = motorVar + taskVar;
% [a, b] = (sort(cData,'descend')); %sort indices
% plot(a,'k','linewidth',2); hold on
% plot(taskVar(b),'g','linewidth',2); axis square
% ylabel('Absolute trial-average modulation')
% xlabel('Sorted neurons');

%% trace examples
%% R2 reduction for regressor groups: area specific
figure
clear areaData depthData
% areas = unique(recSite); %different recording sites
areas = {'ALM' 'MM'}; %different recording sites
depths = [0 350 1000]; %define range of depths to group neurons
plotCnt = 1 : (length(depths) - 1) * length(areas);
plotCnt = reshape(plotCnt, length(areas), [])';
Cnt = 0;

for iAreas = 1 : length(areas)
    for iDepths = 1 : length(depths) - 1
        
        Cnt = Cnt + 1;
        cIdx = ismember(recSite,areas{iAreas}) & recDepth > depths(iDepths) & recDepth < depths(iDepths + 1);
        
        spMotorM = corrMaps(cIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shOther')); %this is the current regressor. This ha sall information.
        spMotorRedM = fullCorrMaps(cIdx)' - corrMaps(cIdx,ismember(nRecLabels,{'spontMotor'}),ismember(modTypes,'shCurrent')); %this is the non-redundant task only model.

        opMotorM = corrMaps(cIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shOther')); %this is the current regressor. This ha sall information.
        opMotorRedM = fullCorrMaps(cIdx)' - corrMaps(cIdx,ismember(nRecLabels,{'opMotor'}),ismember(modTypes,'shCurrent')); %this is the non-redundant task only model.

        taskM = corrMaps(cIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shCurrent')); %this is the current regressor. This ha sall information.
        taskRedM = fullCorrMaps(cIdx)' - corrMaps(cIdx,ismember(nRecLabels,{'motor'}),ismember(modTypes,'shOther')); %this is the non-redundant task only model.
        
        areaData{iAreas, iDepths}(1,:,:) = [taskM spMotorM opMotorM fullCorrMaps(cIdx)']; %all reg information
        areaData{iAreas, iDepths}(2,:,:) = [taskRedM spMotorRedM opMotorRedM fullCorrMaps(cIdx)']; %all reg information

        %make figure
        try
            [ax, idxGroup] = regressorPlot(squeeze(areaData{iAreas, iDepths}(1,:,:)), {'Task' 'Uninstructed' 'Instructed' 'Full'}, 5, subplot(length(depths) - 1, length(areas), plotCnt(Cnt)), [0 1 0], {[1 3 2 4]}, 0.6);
            regressorPlot(squeeze(-areaData{iAreas, iDepths}(2,:,:)), {'Task' 'Uninstructed' 'Instructed' 'Full'}, 5, ax, [25 111 61]/255, idxGroup, 0.6, [-0.5 0.5]);
        catch
            subplot(length(depths) - 1, length(areas), plotCnt(Cnt)); 
            axis square;
        end
        title(['Regressor groups: ' areas{iAreas} ' - ' num2str(depths(iDepths+1)) ' (' num2str(sum(cIdx)) ' cells) - ' fileExt]);
        ylabel('crossVal R^2'); axis square
        xlim([0 5]); axis square
        
    end
end

%% compute modulation index
avgType = 6; %all trials

spMotorRec = recV(:,:,:,end - 1);
opMotorRec = recV(:,:,:,end - 2);
taskRec = recV(:,:,:,end - 3);
motorRec = spMotorRec + opMotorRec;
taskOpRec = taskRec + opMotorRec;
taskSpRec = taskRec + spMotorRec;

taskVar = sum(abs(taskRec(:,:,avgType)),2);
motorVar = sum(abs(motorRec(:,:,avgType)),2);
taskOpVar = sum(abs(taskOpRec(:,:,avgType)),2);
taskSpVar = sum(abs(taskSpRec(:,:,avgType)),2);
spMotorVar = sum(abs(spMotorRec(:,:,avgType)),2);
opMotorVar = sum(abs(opMotorRec(:,:,avgType)),2);

taskIndex = (((taskVar-motorVar) ./ (taskVar + motorVar)) + 1) / 2;
spMotorIndex = (((spMotorVar-taskOpVar) ./ (spMotorVar + taskOpVar)) + 1) / 2;
opMotorIndex = (((opMotorVar-taskSpVar) ./ (opMotorVar + taskSpVar)) + 1) / 2;

taskModIndex = taskVar - recV(:,:,:,end);


h1 = figure;
subplot(3,1,1);
histogram(taskIndex,0:0.05:1, 'FaceColor', 'g'); ylim([0 3000]);
% vline(taskIndex([5849 5944 6005 5941 5986]),'r', 'task cell');
axis square; title(['Median modIdx: ' num2str(median(taskIndex)) '  - TaskIndex'])
xlabel('task index'); ylabel('cell count');

subplot(3,1,2);
histogram(spMotorIndex,0:0.05:1, 'FaceColor', 'k'); ylim([0 3000]);
% vline(spMotorIndex([5849 5944 6005 5941 5986]),'r', 'uninstructed cell');
axis square; title(['Median modIdx: ' num2str(median(spMotorIndex)) '  -  UninstructedIndex'])
xlabel('task index'); ylabel('cell count');

subplot(3,1,3);
histogram(opMotorIndex,0:0.05:1, 'FaceColor', 'b'); ylim([0 3000]);
% vline(opMotorIndex([5849 5944 6005 5941 5986]),'r', 'instructed cell');
axis square; title(['Median modIdx: ' num2str(median(opMotorIndex)) '  - InstructedIndex'])
xlabel('task index'); ylabel('cell count');

%% compute distribution of single variable contributions
cLabels = recLabels(ismember(recLabels,[taskLabels motorLabels])); %all task/motor variables
totalVar = sum(sum(abs(recV(:,:,avgType,ismember(recLabels,cLabels))),2),4); %total absolute modulation by all variables
variableVar = squeeze(sum(abs(recV(:,:,avgType,ismember(recLabels,cLabels))),2)); %absolute modulation for each variable

variableVar = bsxfun(@rdivide,variableVar,totalVar);
variableVar = cumsum(sort(variableVar,2,'descend'),2);


varPowerDist = NaN(length(cLabels)-6,1);
for x = 1 : length(cLabels)-6
    varPowerDist(x) = sum(variableVar(:,x) > 0.3 );
end
figure
% bar(varPowerDist); hold on;
bar(diff([0;varPowerDist]))
xlim([0 10]);

%% single cell examples
% cData = motorVar + taskVar;
% [~, cellIdx] = (sort(cData,'descend')); %sort indices
% [~, cellIdx] = sort(opMotorIndex, 'descend');
% [~, cellIdx] = sort(taskIndex, 'descend');
% [~, cellIdx] = sort(motorVar + taskVar, 'descend');
% [~, cellIdx] = sort(glassLick, 'descend');
[~, cellIdx] = sort(fullCorrMaps, 'descend');

figure
avgType = 6; %all trials
fullRec = recV(:,:,:,end);
spMotorRec = recV(:,:,:,end - 1);
opMotorRec = recV(:,:,:,end - 2);
taskRec = recV(:,:,:,end - 3);
fullErr = semV(:,:,:,end);
spMotorErr = semV(:,:,:,end - 1);
opMotorErr = semV(:,:,:,end - 2);
taskErr = semV(:,:,:,end - 3);

motorRec = recV(:,:,:,end - 1) + recV(:,:,:,end - 2);

% taskVar = sum(abs(taskRec(:,:,avgType)),2);
% spMotorVar = sum(abs(spMotorRec(:,:,avgType)),2);
% opMotorVar = sum(abs(opMotorRec(:,:,avgType)),2);
% motorVar = sum(abs(motorRec(:,:,avgType)),2);

% modIndex = ((taskVar-(opMotorVar)) ./ (taskVar + (opMotorVar)));
% modIndex = ((spMotorVar-(taskVar)) ./ (spMotorVar + (taskVar)));
% modIndex = ((taskVar-motorVar) ./ (taskVar + motorVar));

% [~, cellIdx] = sort((taskModIndex), 'descend');
% cellIdx(ismember(cellIdx, find(isnan(corrMaps(:,1,1))))) = [];

% for x = 6093
% for x = [5849 5944 6005 5941 5986]
% for x = [5973 3022 7361 9951 4625 2973 6444 6003 6442 5944 5947 6093 1711 3017 5957 5971 4877 11384 5986 5794 5817 5849 5853 6005 5851 5855 6008 5860 5872 5914 5941]
% for x = [3531 8195 8311 8000 12967 12979 1755 7910 8269 12978 7972 7771 8414 8609 2681 8263 8311 2897 8113 8289] % potential motor cell
% for x = [2946 2703 5801 2784 3652 3869 3511 6003 2861 2941 6596 2748 2695 6820 8117 3350] % potential task cells
for x = cellIdx'
    % cCell = cellIdx(8);
%     cCell = 6093;
%     cCell = 3460; %task cell1
%     cCell = 12979; %uninstructed cell1
%     cCell = 8263; %uninstructed cell1
%     cCell = 11972; %instructed cell1
    
    cCell = (x);
    
    
    clear cData
    for y = 1:length(nRecLabels)
        
        cRegMod = corrMaps(cCell,y,ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
        cRegRedM = fullCorrMaps(cCell)' - corrMaps(cCell,y,ismember(modTypes,'shCurrent')); %this is the current regressor. This has non-redundant information.
        
        cData(1,y) = cRegMod;
        cData(2,y) = cRegRedM;
        
    end
    
    if sum(isnan(cData(:))) ~= numel(cData)
    idx = zeros(1,length(nRecLabels));
    idx(ismember(nRecLabels,cogLabels)) = 1;
    idx(ismember(nRecLabels,sensorLabels)) = 2;
    idx(ismember(nRecLabels,otherMotorLabels)) = 3;
    % idx(cData(2,:) < 0.01) = 0;
    
    delete(subplot(2,2,1:2));
    ax = regressorPlot(squeeze(cData(1,idx>0,:)),nRecLabels(idx>0),5,subplot(2,2,1:2),[0 1 0],idx(idx>0),0.6, [0 max(squeeze(cData(1,idx>0,:)))]);
    ax = regressorPlot(squeeze(cData(2,idx>0,:)),nRecLabels(idx>0),5,ax,[25 111 61]/255,idx(idx>0),0.3, [0 max(squeeze(cData(1,idx>0,:)))]);
    legend(ax.Children([3 1]),{'All information' 'Non-redundant information'},'location','northwest')
    title(['Regressor comparison - ' recSite{cCell} ' - ' num2str(recDepth(cCell)) ' - ' fileExt]);
    ylabel('crossVal R^2');
    xlim([0 sum(idx>0)+1]);
    
    % subplot(2,8,9:11)
    ax1 = subplot(4,8,[17 18 25 26]);
    ax2 = subplot(4,8,[19 20 27 28]);
    ax3 = subplot(4,8,[21 22 29 30]);
    cla(ax1);cla(ax2);cla(ax3);
    
    amean{1} = fullRec(cCell,1:81,avgType); amean{2} = fullRec(cCell,83:end,avgType);
    astd{1} = fullErr(cCell,1:81,avgType); astd{2} = fullErr(cCell,83:end,avgType);
    F{1} = (1:length(amean{1}))/31;
    F{2} = (F{1}(end) + 0.1) + (1:length(amean{2}))/31;
    splitErrorshade(ax1, amean, [], F, [0.5 0.5 0.5]);
    splitErrorshade(ax2, amean, [], F, [0.5 0.5 0.5]);
    splitErrorshade(ax3, amean, [], F, [0.5 0.5 0.5]);
    
    amean{1} = opMotorRec(cCell,1:81,avgType); amean{2} = opMotorRec(cCell,83:end,avgType);
    astd{1} = opMotorErr(cCell,1:81,avgType); astd{2} = opMotorErr(cCell,83:end,avgType);
    F{1} = (1:length(amean{1}))/31;
    F{2} = (F{1}(end) + 0.1) + (1:length(amean{2}))/31;
    splitErrorshade(ax1, amean, astd, F, 'b');
    
    amean{1} = spMotorRec(cCell,1:81,avgType); amean{2} = spMotorRec(cCell,83:end,avgType);
    astd{1} = spMotorErr(cCell,1:81,avgType); astd{2} = spMotorErr(cCell,83:end,avgType);
    F{1} = (1:length(amean{1}))/31;
    F{2} = (F{1}(end) + 0.1) + (1:length(amean{2}))/31;
    splitErrorshade(ax2, amean, astd, F, 'k');
    
    amean{1} = taskRec(cCell,1:81,avgType); amean{2} = taskRec(cCell,83:end,avgType);
    astd{1} = taskErr(cCell,1:81,avgType); astd{2} = taskErr(cCell,83:end,avgType);
    F{1} = (1:length(amean{1}))/31;
    F{2} = (F{1}(end) + 0.1) + (1:length(amean{2}))/31;
    splitErrorshade(ax3, amean, astd, F, 'g');
    
    lines = findall(ax1.Children,'Type','line');
    ax1.XLim = [0 max(cat(2,lines(:).XData))];
    ax2.XLim = [0 max(cat(2,lines(:).XData))];
    ax3.XLim = [0 max(cat(2,lines(:).XData))];
    ax1.XTick = [];
    
    title(ax1,['CellNr: ' num2str(cCell) ' (Counter: ' num2str(x) ') ; R^2: ' num2str(fullCorrMaps(cCell))])
    ylabel(ax1,'dF / F');
    ylabel(ax2,'dF / F');
    xlabel(ax2, 'time (s)');
    
    
    ax = subplot(4,8,[23 24 31 32]); cla;
    cOpMotor = sum(abs(opMotorRec(cCell,:,avgType)),2);
    cSpMotor = sum(abs(spMotorRec(cCell,:,avgType)),2);
    cTask = sum(abs(taskRec(cCell,:,avgType)),2);
    bar(ax,[sum(abs(fullRec(cCell,:,avgType)),2) cOpMotor cSpMotor cTask])
    set(ax,'xTick',1:4)
    set(ax,'xTickLabel',{'Full' 'Instructed' 'Uninstructed' 'Task'})
    set(ax,'XTickLabelRotation',45)
    ax.TickLength = [0 0];%
    
    %total amount of PETH change from all regressors
    totalVar = sum(sum(abs(recV(cCell,:,avgType,ismember(recLabels,taskLabels) | ismember(recLabels,motorLabels))),2)); 

    %     title(['Mod. index ' num2str(modIndex(cCell))]);
    delete(subplot(2,8,1:8));
    ax1 = subplot(2,8,1:8); cla;
    cLabels = recLabels(ismember(recLabels,[taskLabels motorLabels]));
    [a, b] = sort(squeeze(sum(abs(recV(cCell,:,avgType,ismember(recLabels,cLabels))),2))', 'descend');
    bar(ax1,double(a)./totalVar);
    set(ax1,'xTick',1:length(a))
    set(ax1,'xTickLabel',cLabels(b))
    set(ax1,'XTickLabelRotation',45)
    ax1.TickLength = [0 0];%
    title('task');

    drawnow;
    
    pause
    end
end