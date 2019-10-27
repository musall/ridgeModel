%%
if ~exist('reload','var') || isempty(reload)
    reload = false;
end

%% select data sets
cMod = 'all';
[dataOverview, motorLabels, sensorLabels, cogLabels, ~, segLabels] = delayDecRecordings;
% modTypes = {'shCurrent' 'shOther' 'shTaskCurrent' 'shOtherMotor'};
modTypes = {'shCurrent' 'shOther' 'shOtherMotor' 'shOtherSpontMotor' 'shTaskOtherSpontMotor'};
orgVars = {'fastPupil', 'slowPupil', 'whisk', 'nose', 'face', 'lLick', 'rLick', 'piezo', ...
           'pupils', 'licks', 'bhvVideo', 'Move', 'allMove', 'motorNoVideo', 'motor', 'opMotor', 'spontMotor'}; %make sure to load original version during shOther condition to avoid false results from orthogonalization.
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts
dangerMice = ismember(dataOverview(:,1),{'mSM34' 'mSM36' 'mSM46' 'mSM43' 'mSM44'}); %index for danger mice. 46 has the most, 34 the least inter-ictal events.

%% general variables
cPath = [pwd filesep 'Widefield' filesep]; %path to widefield data
sPath = [pwd filesep 'Widefield' filesep]; %path to widefield data

opMotorLabels = {'handles' 'licks'}; %operant motor regressors
spMotorLabels = {'piezo' 'whisk' 'nose' 'body' 'bhvVideo' 'allMove' 'pupils'}; %spontaneous motor regressors

%% load data
check = true;
if ~reload %check for a saved dataset
    if ~exist(sPath,'dir')
        check = false;
    else
        try
            load([sPath cMod '_corrMaps.mat'],'corrMaps','allenMask','dorsalMaps');
            load([sPath cMod '_meanRsq.mat'],'aRsq','aTimeRsq','nRecLabels');
            load([sPath cMod '_segMovies.mat'],'segMovies');
            load([sPath cMod '_fullMovies.mat'],'fullMovies','taskMovies','opMotorMovies','spMotorMovies');
            load([sPath cMod '_fullCorrMaps.mat'],'fullCorrMaps', 'animals');
            load([sPath cMod '_fullSegMovies.mat'],'fullSegMovies');
            load([sPath cMod '_otherMsegMovies.mat'],'otherMsegMovies');
            load([sPath cMod '_otherMcorrMaps.mat'],'otherMcorrMaps', 'oMotorLabels', 'otherMotorLabels');

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
% check animal count
if strcmpi(cMod,'all')
    animalCnt = size(dataOverview,1);
else
    animalCnt = sum(ismember(dataOverview(:,2)',cMod));
end
animals = dataOverview(:,1);

if reload
    animals = dataOverview(:,1);
    Cnt = 0;
    
    for iAnimals = 1 : length(animals)
        regCnt = 0;
        
        if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'all')
            Cnt = Cnt + 1;
            fPath = [cPath animals{iAnimals} filesep dataOverview{iAnimals,3} filesep];
            disp(fPath); tic;
            
            %load design matrix and check for regressors that were used in the model
            load([fPath 'regData.mat'],'recIdx','idx','recLabels') 
            usedR{Cnt} = recLabels(unique(recIdx(~idx))); %regressors that were used in the model
            
            load([fPath 'Snapshot_1.mat']); %load snapshot
            load([fPath 'mask.mat']); %load mask
            allOpts(Cnt) = load([fPath 'opts2.mat']); %load allen alignment file
            
            if Cnt == 1
                % create mask from allen coordinates
                load('allenDorsalMap.mat')
                addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
                allenMask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
                [x1, y1] = size(allenMask);
                [x2, y2] = size(mask);
                allenMask = allenMask(1:min([x1 x2]), 1:min([y1 y2])); %cut allen mask to size
                [x1, y1] = size(allenMask);

                % get labels for all regressors
                load([fPath 'predVariance' filesep modTypes{1} filesep 'fullcorr.mat'],'recLabels', 'segMovie','cMovie'); %load current labels
                load([fPath 'predVariance' filesep 'extraGroups.mat'],'extraGroups'); %load extra group labels
                load([fPath 'predVariance' filesep 'oMotorLabels.mat'],'oMotorLabels'); %load other motor labels
                nRecLabels = [recLabels{:} extraGroups(1,:)]; %labels of all regressors
                nRecLabels = nRecLabels(~strcmpi(nRecLabels,'leverIn')); %don't use leverIn regressor

                % make recLabels where modality is swapped for expertise
                altRecLabels = nRecLabels;
                altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
                altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};

                % pre-allocate data matrices
                fullCorrMaps = NaN(sum(~allenMask(:)), animalCnt, 'single');
                fullSegMovies = NaN(sum(~allenMask(:)), animalCnt, size(segMovie,2), 'single');
                fullMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
                taskMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
                opMotorMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
                spMotorMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
                
                corrMaps = NaN(sum(~allenMask(:)), animalCnt, length(nRecLabels), 2, 'single');
                segMovies = NaN(sum(~allenMask(:)), animalCnt, length(nRecLabels), size(segMovie,2), 2, 'single');
                                
                otherMcorrMaps = NaN(sum(~allenMask(:)), animalCnt, length(oMotorLabels), 'single');
                otherMsegMovies = NaN(sum(~allenMask(:)), animalCnt, length(oMotorLabels), size(segMovie,2), 'single');
                                
                aRsq = NaN(length(nRecLabels) + 1, animalCnt, length(modTypes), 'single');
                aTimeRsq = NaN(length(nRecLabels) + 1, size(segMovie,2), animalCnt, length(modTypes), 'single');
            end
            
            
            for iRuns = 1:length(modTypes)
                
                if iRuns == 1 %load full model data only once - running through different versions makes no difference here
                    
                    load([fPath 'predVariance' filesep modTypes{iRuns} filesep 'fullcorr.mat'],'cMap','segMovie','cMovie','recLabels'); %load current data
                    
                    % align cMap to allen coordinates
                    cMap = arrayShrink(cMap.^2,mask,'split');
                    cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
                    cMap = arrayShrink(cMap(1:x1,1:y1),allenMask);
                    fullCorrMaps(:, Cnt) = cMap;
                    aRsq(1,Cnt,1) = nanmean(cMap(:));
                    
                    % align segMovie to allen coordinates
                    segMovie = arrayShrink(segMovie.^2,mask,'split');
                    segMovie = alignAllenTransIm(segMovie,allOpts(Cnt).opts.transParams);
                    segMovie = arrayShrink(segMovie(1:x1,1:y1,:),allenMask);
                    fullSegMovies(:, Cnt, :) = segMovie;
                    aTimeRsq(1,:,Cnt,1) = nanmean(segMovie);
                    
                    % align cMovie to allen coordinates
                    cMovie = arrayShrink(cMovie.^2,mask,'split');
                    cMovie = alignAllenTransIm(cMovie,allOpts(Cnt).opts.transParams);
                    cMovie = arrayShrink(cMovie(1:x1,1:y1,:),allenMask);
                    fullMovies(:, Cnt, :) = cMovie;

                end
                
                for iRegs = 1:length(nRecLabels)
                    try
                        %check if current regressor is task or motor model and load cMovie for that case
                        if (strcmpi(modTypes{iRuns},'shOther') && strcmpi(nRecLabels{iRegs},'motor')) || ... %task model. load movie
                                 (strcmpi(modTypes{iRuns},'shCurrent') && strcmpi(nRecLabels{iRegs},'opMotor')) || ... %operant motor model. load movie.
                                 (strcmpi(modTypes{iRuns},'shCurrent') && strcmpi(nRecLabels{iRegs},'spontMotor')) %spontaneous motor model. load movie
                                                         
                             load([fPath 'predVariance' filesep modTypes{iRuns} filesep nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie','cMovie'); %load current data
                             cMovie = arrayShrink(cMovie.^2,mask,'split');
                             cMovie = alignAllenTransIm(cMovie,allOpts(Cnt).opts.transParams);
                             cMovie = arrayShrink(cMovie(1:x1,1:y1,:),allenMask);
                             
                             if strcmpi(modTypes{iRuns},'shOther') && strcmpi(nRecLabels{iRegs},'motor')
                                 taskMovies(:, Cnt, :) = cMovie;
                             elseif strcmpi(modTypes{iRuns},'shCurrent') && strcmpi(nRecLabels{iRegs},'opMotor')
                                 opMotorMovies(:, Cnt, :) = cMovie;
                             elseif strcmpi(modTypes{iRuns},'shCurrent') && strcmpi(nRecLabels{iRegs},'spontMotor')
                                 spMotorMovies(:, Cnt, :) = cMovie;
                             end
                        else
                            if ~ismember(nRecLabels{iRegs}, [usedR{Cnt} extraGroups(1,:)]) && strcmpi(modTypes{iRuns},'shOther') % if regressor was not used and other regressors were shuffled, set all results to zero
                                load([fPath 'predVariance' filesep 'shCurrent' filesep 'fullcorr.mat'],'cMap','segMovie'); %load full model instead of current regressor data
                                cMap = cMap - cMap; %set to 0
                                segMovie = segMovie - segMovie; %set to 0
                            elseif (ismember(nRecLabels{iRegs}, orgVars) && strcmpi(modTypes{iRuns},'shOther')) || strcmpi(modTypes{iRuns},'shOtherMotor') || strcmpi(modTypes{iRuns},'shOtherSpontMotor') || strcmpi(modTypes{iRuns},'shTaskOtherSpontMotor') % use original for standalone regressors
                                load([fPath 'predVariance' filesep modTypes{iRuns} filesep 'org' nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                            else
                                load([fPath 'predVariance' filesep modTypes{iRuns} filesep nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                            end
                        end
                        
                        segMovie = arrayShrink(segMovie.^2,mask,'split');
                        segMovie = alignAllenTransIm(segMovie,allOpts(Cnt).opts.transParams);
                        segMovie = arrayShrink(segMovie(1:x1,1:y1,:),allenMask);
                                                
                        cMap = arrayShrink(cMap.^2,mask,'split');
                        cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
                        cMap = arrayShrink(cMap(1:x1,1:y1),allenMask);
                        
                        % change current regressor index if modality should point to animal experitise instead.
                        cInd = iRegs;
                        if strcmpi(nRecLabels{iRegs},'visReward') || strcmpi(nRecLabels{iRegs},'audReward')
                            if (visExp(iAnimals) && strcmpi(nRecLabels{iRegs},'visReward')) || (audExp(iAnimals) && strcmpi(nRecLabels{iRegs},'audReward'))
                                cInd = find(ismember(altRecLabels, 'expReward')); %change index to point to expert instead of modality
                            else
                                cInd = find(ismember(altRecLabels, 'novReward')); %change index to point to novice instead of modality
                            end
                        end
                        aTimeRsq(cInd + 1, :, Cnt, iRuns) = nanmean(segMovie);
                        aRsq(cInd + 1, Cnt, iRuns) = nanmean(cMap(:));
                            
                        if strcmpi(modTypes{iRuns},'shCurrent') || strcmpi(modTypes{iRuns},'shOther')
                            segMovies(:, Cnt, cInd, :, iRuns) = segMovie;
                            corrMaps(:, Cnt, cInd, iRuns) = cMap;
                        elseif strcmpi(modTypes{iRuns},'shOtherMotor') && any(ismember(oMotorLabels,nRecLabels{iRegs}))
                            regCnt = regCnt + 1;
                            otherMsegMovies(:, Cnt, regCnt, :) = segMovie;
                            otherMcorrMaps(:, Cnt, regCnt) = cMap;
                            otherMotorLabels{regCnt} = nRecLabels{iRegs}; %this is a control. should equal oMotorLabels in the end.
                        end
                    end
                end
            end
            toc;
        end
    end
    clear segMovie cMap
    
    %% save down results
    if ~exist(sPath,'dir')
        mkdir(sPath);
    end
    save([sPath cMod '_fullCorrMaps.mat'],'fullCorrMaps', 'animals');
    save([sPath cMod '_fullSegMovies.mat'],'fullSegMovies');
    save([sPath cMod '_fullMovies.mat'],'fullMovies','taskMovies','opMotorMovies','spMotorMovies', '-v7.3');
    save([sPath cMod '_corrMaps.mat'],'corrMaps','allenMask','dorsalMaps');
    save([sPath cMod '_segMovies.mat'],'segMovies', '-v7.3');
    save([sPath cMod '_otherMsegMovies.mat'],'otherMsegMovies', '-v7.3');
    save([sPath cMod '_otherMcorrMaps.mat'],'otherMcorrMaps', 'oMotorLabels', 'otherMotorLabels');
    save([sPath cMod '_meanRsq.mat'],'aRsq','aTimeRsq','nRecLabels');
    
end

%% check if shOtherMotor was loaded correctly
if length(oMotorLabels) ~= sum(ismember(oMotorLabels,otherMotorLabels))
    error('oMotorLabels and otherMotorLabels are not the same. Problem with loading shOtherMotor variables?')
end

%%
% fullRecLabels = [{'full'} nRecLabels]; %add full moodel label
segLabels = [{'all'} segLabels];
altRecLabels = nRecLabels;
altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};
fullRecLabels = [{'full'} altRecLabels]; %add full moodel label

%% crossvalidated R2 maps - full model
cMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');
cMovie = arrayShrink(nanmean(fullSegMovies,2),allenMask,'split');
cMovie = cat(3,cMap,squeeze(cMovie)); %combine full correlation map with different segments
% cRange = max(abs(cMap(:))).*0.75;
cRange = 0.5;

figure;
Cnt = 0;
for iSegs = [2 4 6 7]
    % for iSegs = 1:size(cMovie,3)
    Cnt = Cnt + 1;
    subplot(1,4,Cnt)
    mapImg = imshow(cMovie(:,:,iSegs),[0 cRange]);axis image;
    set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
    colormap(mapImg.Parent,inferno(256));
    title(['Full model, R^2 - ' segLabels{iSegs}])
end

%% make overview figure for overall task- and motor explained variance
fullMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');
taskMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'motor'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
spmotorMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'spontMotor'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
opmotorMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'opMotor'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');

% cRange = prctile(fullMap(:),99);
cRange = 0.2;

figure
subplot(1,3,1)
mapImg = imshow(taskMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(taskMap)); %make NaNs transparent.
title('Task model- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,3,2)
mapImg = imshow(opmotorMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(opmotorMap)); %make NaNs transparent.
title('Instructed motor model- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,3,3)
mapImg = imshow(spmotorMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(spmotorMap)); %make NaNs transparent.
title('Spontaneous motor model- R^2')
colormap(mapImg.Parent,inferno(256));

%% same figure but comparing visual and auditory experts
for iMod = 1:2
    if iMod == 1
        cIdx = visExp;
    else
        cIdx = ~visExp;
    end
    
    fullMap = arrayShrink(nanmean(fullCorrMaps(:, cIdx),2),allenMask,'split');
    taskMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'motor'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
    spmotorMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'spontMotor'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    opmotorMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'opMotor'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    
    % cRange = prctile(fullMap(:),99);
    cRange = 0.2;
    figure
    subplot(1,3,1)
    mapImg = imshow(taskMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(taskMap)); %make NaNs transparent.
    title('Task model- R^2')
    colormap(mapImg.Parent,inferno(256));
    
    subplot(1,3,2)
    mapImg = imshow(opmotorMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(opmotorMap)); %make NaNs transparent.
    title('Instructed motor model- R^2')
    colormap(mapImg.Parent,inferno(256));
    
    subplot(1,3,3)
    mapImg = imshow(spmotorMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(spmotorMap)); %make NaNs transparent.
    title('Spontaneous motor model- R^2')
    colormap(mapImg.Parent,inferno(256));
end

%% make movie for task- and motor explained variance
% cMovie = arrayShrink(squeeze(nanmean(fullMovies - taskMovies,2)), allenMask, 'split');
% cMovie = arrayShrink(squeeze(nanmean(fullMovies - opMotorMovies,2)), allenMask, 'split');
% cMovie = arrayShrink(squeeze(nanmean(fullMovies - spMotorMovies,2)), allenMask, 'split');
% cMovie = arrayShrink(squeeze(nanmean(fullMovies,2)), allenMask, 'split');
% compareMovie(cMovie,'outline',dorsalMaps.edgeOutlineSplit);

%% show predicted R2 over time
figure;
plot(squeeze(nanmean(nanmean(fullMovies))), 'k', 'linewidth', 2)
hold
plot(squeeze(nanmean(nanmean(fullMovies - taskMovies))), 'r', 'linewidth', 2)
plot(squeeze(nanmean(nanmean(fullMovies - opMotorMovies))), 'b', 'linewidth', 2)
plot(squeeze(nanmean(nanmean(fullMovies - spMotorMovies))), 'g', 'linewidth', 2)
axis square; legend({'full' 'task' 'opMotor' 'spMotor'});
xlabel('time(s)'); ylabel('predicted R^2');

%% single variable models 
handleMap = arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rGrab'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
stimMap = arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rVisStim'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
noseMap = arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'nose'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
pawMap = arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'piezo'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');

figure
subplot(2,2,1)
cRange = 0.1;
mapImg = imshow(handleMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(handleMap)); %make NaNs transparent.
title('Right handle- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,2)
cRange = 0.2;
mapImg = imshow(noseMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(noseMap)); %make NaNs transparent.
title('Nose - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,3)
cRange = 0.02;
mapImg = imshow(stimMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(stimMap)); %make NaNs transparent.
title('Right stim - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,4)
cRange = prctile(pawMap(:),99);
mapImg = imshow(pawMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(pawMap)); %make NaNs transparent.
title('Hindpaw - R^2')
colormap(mapImg.Parent,inferno(256));


%% single regressor effects - unique maps
fullMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');
handleMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rGrab'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
stimMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rVisStim'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
noseMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'nose'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
pawMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'piezo'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');

figure
subplot(2,2,1)
cRange = 0.035;
mapImg = imshow(handleMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(handleMap)); %make NaNs transparent.
title('Right handle- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,2)
% cRange = prctile(handleMap(:),99);
mapImg = imshow(noseMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(noseMap)); %make NaNs transparent.
title('Nose - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,3)
cRange = 0.02;
mapImg = imshow(stimMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(stimMap)); %make NaNs transparent.
title('Right stim - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(2,2,4)
% cRange = prctile(pawMap(:),99);
mapImg = imshow(pawMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(pawMap)); %make NaNs transparent.
title('Hindpaw - R^2')
colormap(mapImg.Parent,inferno(256));

%% same figure but comparing visual and auditory experts
for iMod = 1:2
    if iMod == 1
        cIdx = visExp;
    else
        cIdx = ~visExp;
    end
    
    fullMap = arrayShrink(nanmean(fullCorrMaps(:,cIdx),2),allenMask,'split');
    handleMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'rGrab'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    stimMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'rVisStim'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    noseMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'nose'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    pawMap = fullMap - arrayShrink(nanmean(corrMaps(:,cIdx,strcmpi(altRecLabels,{'piezo'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
    
    figure
    subplot(2,2,1)
    cRange = 0.035;
    mapImg = imshow(handleMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(handleMap)); %make NaNs transparent.
    title('Right handle- R^2')
    colormap(mapImg.Parent,inferno(256));
    
    subplot(2,2,2)
    % cRange = prctile(handleMap(:),99);
    mapImg = imshow(noseMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(noseMap)); %make NaNs transparent.
    title('Nose - R^2')
    colormap(mapImg.Parent,inferno(256));
    
    subplot(2,2,3)
    cRange = 0.015;
    mapImg = imshow(stimMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(stimMap)); %make NaNs transparent.
    title('Right stim - R^2')
    colormap(mapImg.Parent,inferno(256));
    
    
    subplot(2,2,4)
    % cRange = prctile(pawMap(:),99);
    mapImg = imshow(pawMap,[0 cRange]);axis image;
    set(mapImg,'AlphaData',~isnan(pawMap)); %make NaNs transparent.
    title('Hindpaw - R^2')
    colormap(mapImg.Parent,inferno(256));
end

%% R2 reduction for all regressors
clear cRegMod cRegRedM cData
cogLabels = [cogLabels {'expReward'} {'novReward'} {'lGrabRel'}];
fullM = aRsq(ismember(fullRecLabels,{'full'}), :, 1); %full model
cInd = dangerMice;
% cInd = true(1,length(animals));

cData  = NaN(2,length(fullRecLabels),sum(cInd));
for x = 1:length(fullRecLabels)
    
    cRegMod = aRsq(x, cInd, ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
    cRegRedM = fullM(cInd) - aRsq(x, cInd, ismember(modTypes,'shCurrent')); %this is the current regressor. This has non-redundant information.
        
    cData(1,x,:) = cRegMod;
    cData(2,x,:) = cRegRedM;
       
end

idx = zeros(1,length(fullRecLabels));
idx(ismember(fullRecLabels,opMotorLabels)) = 3;
idx(ismember(fullRecLabels,cogLabels)) = 2;
idx(ismember(fullRecLabels,sensorLabels)) = 1;
idx(ismember(fullRecLabels,spMotorLabels)) = 4;
% idx(isnan(cData(1,:,1))) = 0;

figure
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,idx>0,:))',fullRecLabels(idx>0),5,subplot(1,1,1),[0 1 0],idx(idx>0),0.6);
regressorBoxPlot(squeeze(-cData(2,idx>0,:))',fullRecLabels(idx>0),5,ax,[25 111 61]/255,idxGroup,0.6, [-0.1 0.35]);

title('Single varibles: All (top) / Unique (bottom)');
ylabel('cv R^2');
ylabel('unique delta cv R^2');
xlim([0 25]); %make this a fixed value to ensure bars have always the same width

%% overview figure for all unique R2 maps
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
idx = zeros(1,length(altRecLabels));
idx(ismember(altRecLabels,cogLabels)) = 1;
idx(ismember(altRecLabels,sensorLabels)) = 1;
idx(ismember(altRecLabels,otherMotorLabels(~ismember(otherMotorLabels,'face')))) = 2;
fullMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');

for iMods = 1 : 2
    figure('name',['iMods = ' int2str(iMods)])
    cIdx = find(idx == iMods);
    Cnt = 0;
    
    for x = cIdx
        
        Cnt = Cnt +1;
        cMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,x,strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
        allScale{iMods}(Cnt) = max(abs(prctile(cMap(:),[1 99])));

        subplot(ceil(length(cIdx) / 4), 4 ,Cnt);

        cRange = max(abs(prctile(cMap(:),[1 99])));
        mapImg = imshow(cMap,[0 cRange]);
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        colormap(mapImg.Parent,inferno(256));hold on;
        title(altRecLabels{x}); colorbar
        drawnow;
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
                
    end
end

%% R2 reduction for regressor groups
% cInd = ~dangerMice;
cInd = true(1,animalCnt);

fullM = aRsq(ismember(fullRecLabels,{'full'}), cInd, 1); %full model

spMotorM = aRsq(ismember(fullRecLabels,{'spontMotor'}), cInd, ismember(modTypes,'shOther')); %this is the spont. motor only model. Shuffle everything but motor regressors
spMotorRedM = fullM - aRsq(ismember(fullRecLabels,{'spontMotor'}), cInd, ismember(modTypes,'shCurrent')); %this is the non-redundant spont. motor only model. Shuffle all motor regressors and substract from full.

opMotorM = aRsq(ismember(fullRecLabels,{'opMotor'}), cInd, ismember(modTypes,'shOther')); %this is the operant motor only model. Shuffle everything but motor regressors
opMotorRedM = fullM - aRsq(ismember(fullRecLabels,{'opMotor'}), cInd, ismember(modTypes,'shCurrent')); %this is the non-redundant operant motor only model. Shuffle all motor regressors and substract from full.

taskM = aRsq(ismember(fullRecLabels,{'motor'}), cInd, ismember(modTypes,'shCurrent')); %this is the task only model. Shuffle only motor regressors to leave task regressors
taskRedM = fullM - aRsq(ismember(fullRecLabels,{'motor'}), cInd, ismember(modTypes,'shOther')); %this is the non-redundant task only model.

clear cData cError
cData(1,:,:) = [taskM' spMotorM' opMotorM' fullM']; %all reg information
cData(2,:,:) = [taskRedM' spMotorRedM' opMotorRedM' fullM']; %all reg information

figure
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,:,:)), {'Task' 'spontMove' 'opMove' 'full'}, 5, subplot(2,2,1:4), [0 1 0], {[1 3 2 4]},0.6,[-0.5 0.5]);
ax = regressorBoxPlot(squeeze(-cData(2,:,:)), {'Task' 'spontMove' 'opMove' 'full'}, 5, ax, [25 111 61]/255, idxGroup,0.6,[-0.5 0.5]);
title('Regg groups: All (top) / Unique (bottom)'); xlim([0 30]); %make this a fixed value to ensure bars have always the same width
ylabel('total cv R^2');
hline(mean(fullM),'k--')

%% Task-dependent R2 reduction for motor regressors
taskM = aRsq(ismember(fullRecLabels,{'motor'}), :, ismember(modTypes,'shCurrent')); %this is the task only model. Shuffle only motor regressors to leave task regressors

cData  = NaN(2,length(otherMotorLabels),animalCnt);
Cnt = 0;
for x = 1:length(otherMotorLabels)
    
        Cnt = Cnt +1;
        cRegMod = aRsq(ismember(fullRecLabels,otherMotorLabels{x}), :, ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
        cRegRedM = aRsq(ismember(fullRecLabels,otherMotorLabels{x}), :, ismember(modTypes,'shOtherMotor')) - taskM; %this is the current regressor - task model. This has non-task specific regressor information.
        
        cData(1,Cnt,:) = cRegMod - cRegRedM; %subtract unique information from full model. Remaining is task-redundant information.
        cData(2,Cnt,:) = cRegRedM;
        
end

idx = zeros(1,length(otherMotorLabels));
idx(ismember(otherMotorLabels,opMotorLabels)) = 1;
idx(ismember(otherMotorLabels,spMotorLabels)) = 2;
idx(ismember(otherMotorLabels,{'opMotor' 'spontMotor'})) = 3;

figure
[ax, allInd] = regressorBoxPlot(squeeze(cData(1,idx>0,:))',otherMotorLabels(idx>0),5,subplot(2,2,1:4),[0,191,255]/255,idx(idx>0),0.6);
ax = regressorBoxPlot(squeeze(-cData(2,idx>0,:))',otherMotorLabels(idx>0),5,ax,[0,0,139]/255,allInd,0.6,[-0.4 0.4]);
xlim([0 30]);
title('Regressor comparison - Task-dependancy');
ylabel('crossVal R^2');
xlim([0 30]);
