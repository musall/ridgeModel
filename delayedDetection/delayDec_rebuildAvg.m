%% basic variables
fPath = [pwd filesep 'Widefield' filesep]; %path to widefield data

trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2

% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, ~, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here

visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

 %get allen maps
load('allenDorsalMapSM.mat')
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
[x1, y1] = size(mask);
load([cPath dataOverview{1,1} filesep dataOverview{1,3} filesep 'snapshot_1.mat'])
load([cPath dataOverview{1,1} filesep dataOverview{1,3} filesep 'opts2.mat'])
snap = alignAllenTransIm(single(snap),opts.transParams);
[x2, y2] = size(snap);
mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
allenMask = dorsalMaps.allenMask;

segFrames = cumsum([54 30 105]);
xRange = [(-segFrames(2) : 1 : 0)./30, (1 : segFrames(end)-segFrames(2)-1) ./30]; %time vector for x-axis. Stimulus onset is at 0s.

%% load raw data
allU = NaN(sum(~allenMask(:)), 200, size(dataOverview,1), 'single'); %pre-allocate larger data array for all Us
allData = NaN(sum(~allenMask(:)), segFrames(end), size(dataOverview,1), 'single'); %pre-allocate larger data array for all Vs
dMaps =  NaN(sum(~allenMask(:)), 2, 2, size(dataOverview,1));

for iAnimals = 1:size(dataOverview,1)

    cFile = dir([cPath dataOverview{iAnimals,1} filesep dataOverview{iAnimals,3} filesep dataOverview{iAnimals,1} '_SpatialDisc*.mat']);
    load([cPath dataOverview{iAnimals,1} filesep dataOverview{iAnimals,3} filesep strtrim(cFile.name)]); %load behavior data
    load([cPath dataOverview{iAnimals,1} filesep dataOverview{iAnimals,3} filesep 'Vc.mat'],'U','Vc','bTrials');
    load([cPath dataOverview{iAnimals,1} filesep dataOverview{iAnimals,3} filesep 'snapshot_1.mat'])
    load([cPath dataOverview{iAnimals,1} filesep dataOverview{iAnimals,3} filesep 'opts2.mat'])
    
    bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dat
    Vc = Widefield_getRealignment(Vc, bhv, segFrames, opts);
    
    snap = alignAllenTransIm(single(snap),opts.transParams);
    U = alignAllenTransIm(U,opts.transParams);
    
    allSnap(:,iAnimals) = arrayShrink(snap(1:size(allenMask,1),1:size(allenMask,2)),allenMask);
    allU(:,:,iAnimals) = arrayShrink(U(1:size(allenMask,1),1:size(allenMask,2),:),allenMask);
        
    visIdx = bhv.StimType == 1;
    allData(:,:,iAnimals) = allU(:,:,iAnimals) * mean(Vc(:,:,visIdx),3); %average over all visual trials
    
    % get average responses for d' during stim1 and delay
    visMap(:,:,1) = squeeze(nanmean(Vc(:,segFrames(2):segFrames(2)+18, visIdx),2)); %visual average during stim1
    visMap(:,:,2) = squeeze(nanmean(Vc(:,segFrames(2) + 51 : segFrames(2) + 51, visIdx),2)); %visual average during delay

    audIdx = bhv.StimType == 2;
    audMap(:,:,1) = squeeze(nanmean(Vc(:,segFrames(2):segFrames(2)+18, audIdx),2)); %auditory average during stim1
    audMap(:,:,2) = squeeze(nanmean(Vc(:,segFrames(2) + 51 : segFrames(2) + 51, audIdx),2)); %auditory average during delay

    vLeftIdx = visIdx & bhv.CorrectSide == 1;
    visLeft(:,:,1) = squeeze(nanmean(Vc(:,segFrames(2):segFrames(2)+18, vLeftIdx),2)); %visual average during stim1
    visLeft(:,:,2) = squeeze(nanmean(Vc(:,segFrames(2) + 51 : segFrames(2) + 51, vLeftIdx),2)); %visual average during delay

    vRightIdx = visIdx & bhv.CorrectSide == 2;
    visRight(:,:,1) = squeeze(nanmean(Vc(:,segFrames(2):segFrames(2)+18, vRightIdx),2)); %visual average during stim1
    visRight(:,:,2) = squeeze(nanmean(Vc(:,segFrames(2) + 51 : segFrames(2) + 51, vRightIdx),2)); %visual average during delay
    
    % compute d'
    for x = 1 : 2
        %modality (vis - audio)
        varP(:,1) = (nansum((allU(:,:,iAnimals) * cov(visMap(:,:,x)')) .* allU(:,:,iAnimals), 2)); %variance map
        varP(:,2) = (nansum((allU(:,:,iAnimals) * cov(audMap(:,:,x)')) .* allU(:,:,iAnimals), 2)); %variance map
                
        dMaps(:,x,1,iAnimals) = allU(:,:,iAnimals) * (mean(visMap(:,:,x),2) - mean(audMap(:,:,x),2)); %response difference
        dMaps(:,x,1,iAnimals) = dMaps(:,x,1,iAnimals) ./ sqrt(sum(varP,2)/2); %normalize by standard deviation
        
        %sides (left - right)
        varP(:,1) = (nansum((allU(:,:,iAnimals) * cov(visLeft(:,:,x)')) .* allU(:,:,iAnimals), 2)); %variance map
        varP(:,2) = (nansum((allU(:,:,iAnimals) * cov(visRight(:,:,x)')) .* allU(:,:,iAnimals), 2)); %variance map
                
        dMaps(:,x,2,iAnimals) = allU(:,:,iAnimals) * (mean(visLeft(:,:,x),2) - mean(visRight(:,:,x),2)); %response difference
        dMaps(:,x,2,iAnimals) = dMaps(:,x,2,iAnimals) ./ sqrt(sum(varP,2)/2); %normalize by standard deviation
    end
    clear visMap audMap visLeft visRight
end


%% show averaged data from visual trials
cMovie = arrayShrink(nanmean(allData,3),mask,'split');
cMovie = cMovie - mean(cMovie(:,:,baseRange),3);
cRange = [-nanmean(abs(cMovie(:))*2) nanmean(abs(cMovie(:)))*2];
figure;
for iSegs = 1:5
    subplot(1,5,iSegs)
    mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['All mice. cVision: ' segLabels{iSegs+1}])
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% figure1 : show d' maps: modality and choice
cLabel = {'d" Vision/Audio - Stim' 'd" Choice - Stim' 'd" Vision/Audio - Stim' 'd" Choice - Delay'};
figure
Cnt = 0;
for iSegs = 1:2
    for iMod = 1:2
        Cnt = Cnt+1;
        subplot(2,2,Cnt)
        
        cMap = arrayShrink(nanmean(dMaps(:,iSegs,iMod,:),4),mask,'split');
        
        mapImg = imshow(cMap,[-0.5 0.5]);
        colormap(mapImg.Parent,viridis); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(cLabel{Cnt})
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(dorsalMaps.edgeOutlineSplit)
            plot(dorsalMaps.edgeOutlineSplit{x}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
    end
end