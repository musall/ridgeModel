function delayDec_NeuropixelsResults(Animal)

cPath = [pwd filesep 'Neuropixels' filesep]; %path to neuropixels data

% clusterLabels = {'V1', 'Cingulate', 'SC', 'Reticular'};
clusterLabels = {'VC', 'SC', 'RF'};
if strcmpi(Animal, 'all')
    fPath{1} = [cPath filesep 'NP9' filesep];
    fPath{2} = [cPath filesep 'N14' filesep];
    trialData{1} = 'NP9_006_sync.mat';
    trialData{2} = 'NP14_008_sync.mat';
    spikeData{1} = 'NP9_006_spikes.mat';
    spikeData{2} = 'NP14_008_spikes.mat';
    depthClusters{1} = {[2541 3500], [1371 2540], [0 1370]};
    depthClusters{2} = {[2011 3700], [871 2010], [0 870]};
elseif strcmpi(Animal, 'NP9')
    fPath{1} = [cPath filesep Animal filesep];
    trialData{1} = 'NP9_006_sync.mat';
    spikeData{1} = 'NP9_006_spikes.mat';
    depthClusters{1} = {[2541 3500], [1371 2540], [0 1370]};
elseif strcmpi(Animal, 'N14')
    fPath{1} = [cPath filesep Animal filesep];
    trialData{1} = 'NP14_008_sync.mat';
    spikeData{1} = 'NP14_008_spikes.mat';
    depthClusters{2} = {[2011 3700], [871 2010], [0 870]};
end
preStim = 0.5;
stimDur = 2.5;
sRate = 30;

%% load data for PSTHs
allSpikes = cell(1, length(fPath));
fullModel = cell(1, length(fPath));
vidModel = cell(1, length(fPath));
stimModel = cell(1, length(fPath));
for x = 1 : length(fPath)
    load([fPath{x} 'spikeTrace.mat'], 'spikeTrace');
    load([fPath{x} 'fullFit.mat'], 'fullFit', 'vidFit', 'stimFit');
    allSpikes{x} = spikeTrace;
    
    load([fPath{x} trialData{x}],'sync_data');
    stimIdx = round((sync_data.photodiode - sync_data.photodiode(1)) * sRate); %stimulus onset times
    stimIdx(stimIdx + ((stimDur-preStim)*sRate) > size(spikeTrace,1)) = [];
    stimIdx(stimIdx - sRate < 0) = [];
    stimIdx = bsxfun(@plus, repmat(stimIdx, 1, (stimDur*sRate)), -(preStim*sRate):((stimDur-preStim)*sRate)-1)';
    stimIdx = stimIdx(:);
    
    % create average over all stimuli
    cData = spikeTrace(stimIdx, :);
    cData = reshape(cData, (stimDur*sRate), [], size(cData,2));
    rawMean{x} = squeeze(mean(cData,2));
    
    cData = fullFit(stimIdx, :);
    cData = reshape(cData, (stimDur*sRate), [], size(cData,2));
    fullMean{x} = squeeze(mean(cData,2));
    
    cData = stimFit(stimIdx, :);
    cData = reshape(cData, (stimDur*sRate), [], size(cData,2));
    stimMean{x} = squeeze(mean(cData,2));    
        
    cData = vidFit(stimIdx, :);
    cData = reshape(cData, (stimDur*sRate), [], size(cData,2));
    vidMean{x} = squeeze(mean(cData,2));  
    
    % load spike data to get depth for each contact
    load([fPath{x} spikeData{x}],'sp');
    if strcmpi(fPath{x}, '\\grid-hs\churchland_nlsas_data\data\ashley_looming\PupilVideos\N14\')
        allDepths{x} = sp.clusterDepths' - 1020; %this recording needs to be corrected for some reason
    else
        allDepths{x} = sp.clusterDepths';
    end
end
rawMean = cat(2,rawMean{:});
fullMean = cat(2,fullMean{:});
stimMean = cat(2,stimMean{:});
vidMean = cat(2,vidMean{:});

%% make PSTH figure - mean over all neurons
figure;
subplot(2, 2, 1:2);
plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, mean(rawMean,2), 'linewidth', 4, 'color', 'k'); hold on
plot(1/sRate : 1/sRate : size(fullMean,1) / sRate, mean(fullMean,2), 'linewidth', 4, 'color', 'r'); axis square; hold off
xlim([0 stimDur])

subplot(2,2,3);
plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, mean(rawMean,2), 'linewidth', 4, 'color', [0.5 0.5 0.5]); axis square; hold on
plot(1/sRate : 1/sRate : size(stimMean,1) / sRate, mean(stimMean,2), 'linewidth', 4, 'color', 'g'); hold off
xlim([0 stimDur])

subplot(2,2,4);
plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, mean(rawMean,2), 'linewidth', 4, 'color', [0.5 0.5 0.5]); axis square; hold on
plot(1/sRate : 1/sRate : size(vidMean,1) / sRate, mean(vidMean,2), 'linewidth', 4, 'color', 'b'); hold off
xlim([0 stimDur])

%% show task/motor index histogram
vidVar = sum(abs(vidMean),1);
stimVar = sum(abs(stimMean),1);
taskIndex = (((stimVar-vidVar) ./ (stimVar + vidVar)) + 1) / 2;
motorIndex = (((vidVar-stimVar) ./ (stimVar + vidVar)) + 1) / 2;

h1 = figure;
subplot(2,1,1);
histogram(taskIndex,0:0.05:1, 'FaceColor', 'g'); 
ylim([0 75]);
axis square; title(['Median modIdx: ' num2str(median(taskIndex)) '  - TaskIndex'])
xlabel('task index'); ylabel('cell count');
% vline(taskIndex([5849 5944 6005 5941 5986]),'r', 'task cell');

subplot(2,1,2);
histogram(motorIndex,0:0.05:1, 'FaceColor', 'r'); 
ylim([0 75]);
axis square; title(['Median modIdx: ' num2str(median(motorIndex)) '  - MotorIndex'])
xlabel('Move index'); ylabel('cell count');
% vline(taskIndex([5849 5944 6005 5941 5986]),'r', 'task cell');

%% single cell examples
% [~, c] = sort(motorIndex,'descend');
% figure;
% for x = 1 : length(c)
%     
% subplot(2, 2, 1:2);
% plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, rawMean(:,c(x)), 'linewidth', 4, 'color', 'k'); hold on
% plot(1/sRate : 1/sRate : size(fullMean,1) / sRate, fullMean(:,c(x)), 'linewidth', 4, 'color', 'r'); axis square; hold off
% 
% subplot(2,2,3);
% plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, rawMean(:,c(x)), 'linewidth', 4, 'color', [0.5 0.5 0.5]); axis square; hold on
% plot(1/sRate : 1/sRate : size(stimMean,1) / sRate, stimMean(:,c(x)), 'linewidth', 4, 'color', 'g'); hold off
% 
% subplot(2,2,4);
% plot(1/sRate : 1/sRate : size(rawMean,1) / sRate, rawMean(:,c(x)), 'linewidth', 4, 'color', [0.5 0.5 0.5]); axis square; hold on
% plot(1/sRate : 1/sRate : size(vidMean,1) / sRate, vidMean(:,c(x)), 'linewidth', 4, 'color', 'b'); hold off
% 
% pause
% end

%% load data for predicted R2 using different models
cFull = cell(1,length(fPath));
cRegs = cell(1,length(fPath));
animalID = cell(1,length(fPath));
for x = 1 : length(fPath)
    load([fPath{x} 'fullV.mat'], 'fullV', 'regLabels');
    load([fPath{x} 'allV.mat'], 'allV', 'regGroups');
    Vm{x} = fullV; %
    animalID{x} = ones(size(Vm{x},2),1)*x;
    for y = 1:size(allSpikes{x},2)
        cFull{x}(y,1) = corr2(allSpikes{x}(:,y), fullV(:,y))^2;
        for z = 1:size(allV,4)
            cRegs{x}(y,1,z) = corr2(allSpikes{x}(:,y), allV(:,y,1,z))^2;
            cRegs{x}(y,2,z) = corr2(allSpikes{x}(:,y), allV(:,y,2,z))^2;
        end
    end
end
cFull = cat(1,cFull{:});
cRegs = cat(1,cRegs{:});
animalID = cat(1,animalID{:});

%% make R2 trace figure
stimReg = ismember(regGroups(1,:), 'stim');
moveReg = ismember(regGroups(1,:), 'movement');

figure
subplot(1,2,1)
[~, c] = sort(cFull,'descend');
plot(cFull(c), 'linewidth', 4, 'color', 'k'); hold on;
plot(cRegs(c, 1, moveReg), 'linewidth', 4, 'color', 'r');
plot(cRegs(c, 1, stimReg), 'linewidth', 4, 'color', 'g');
xlabel('Neurons');
ylabel('cross-val. R^2');
grid on; axis square
title('Predicted R^2 - video only');

subplot(1,2,2)
fullMean = [mean(cRegs(:, 1, stimReg)) mean(cRegs(:, 1, moveReg))];
fullError = [sem(cRegs(:, 1, stimReg)) sem(cRegs(:, 1, moveReg))];

errorbar(fullMean,fullError,'k-','linestyle','none','lineWidth',3); hold on
bar(fullMean,'FaceColor','g','EdgeColor','k','BarWidth',0.5,'LineWidth',2);

uniqueMean = [mean(cFull-cRegs(:, 2, stimReg)) mean(cFull-cRegs(:, 2, moveReg))];
uniqueError = [sem(cFull-cRegs(:, 2, stimReg)) sem(cFull-cRegs(:, 2, moveReg))];

errorbar(-uniqueMean,uniqueError,'k-','linestyle','none','lineWidth',3); hold on
bar(-uniqueMean,'FaceColor',[25 111 61]/255,'EdgeColor','k','BarWidth',0.5,'LineWidth',2);

ax = gca;
set(ax,'xTick',1:size(fullMean,2))
set(ax,'xTickLabel',{'Stimulus' 'Movement'})
set(ax,'XTickLabelRotation',45)
ax.TickLength = [0 0];
ylabel('cross-val. R^2 R^2'); ylim([-0.2 0.2]);
axis square

%% R^2 for different depths
figure
for iAreas = 1 : length(clusterLabels)
    % run analysis for selected neurons
    stimReg = ismember(regGroups(1,:), 'stim');
    moveReg = ismember(regGroups(1,:), 'movement');
    
    % find contacts at the right depth for each area
    cIdx = {}; clear cData uniqueData
    for x = 1 : length(fPath)
        if ~isempty(depthClusters{x}{iAreas})
            cIdx{x} = find(allDepths{x} > depthClusters{x}{iAreas}(1) & allDepths{x} < depthClusters{x}{iAreas}(2))';
        end
        if x == 2
            cIdx{x} = cIdx{x} + size(allSpikes{1},2);
        end
        cData(x,1) = mean(cRegs(cIdx{x}, 1, stimReg));
        cData(x,2) = mean(cRegs(cIdx{x}, 1, moveReg));
        
        uniqueData(x,1) = mean(cFull(cIdx{x}, 1) - cRegs(cIdx{x}, 1, moveReg));
        uniqueData(x,2) = mean(cFull(cIdx{x}, 1) - cRegs(cIdx{x}, 1, stimReg));
        disp(size(cIdx{x},2));
    end
                
    subplot(1,length(clusterLabels),iAreas)
    fullMean = [nanmean(cData(:,1),1) nanmean(cData(:,2),1)];
%     fullError = [sem(cData(:,1),1) sem(cData(:,2),1)];
%     errorbar(fullMean,fullError,'k-','linestyle','none','lineWidth',3); hold on
    bar(fullMean,'FaceColor','g','EdgeColor','k','BarWidth',0.5,'LineWidth',2); hold on;
    plot([1 1; 2 2]',cData,'og','MarkerFaceColor','w')
    
    uniqueMean = [nanmean(uniqueData(:,1),1) nanmean(uniqueData(:,2),1)];
%     uniqueError = [sem(uniqueData(:,1),1) sem(uniqueData(:,2),1)];
%     errorbar(uniqueMean,uniqueError,'k-','linestyle','none','lineWidth',3); hold on
    
    bar(-uniqueMean,'FaceColor',[25 111 61]/255,'EdgeColor','k','BarWidth',0.5,'LineWidth',2);
    plot([1 1; 2 2]',-uniqueData,'ob','Color',[25 111 61]/255,'MarkerFaceColor','w')

    ax = gca;
    set(ax,'xTick',1:size(fullMean,2))
    set(ax,'xTickLabel',{'Stimulus' 'Movement'})
    set(ax,'XTickLabelRotation',45)
    ax.TickLength = [0 0];
    ylabel('cross-val. R^2'); 
    ylim([-0.35 0.35]);
    axis square
    title(clusterLabels{iAreas});
end

%% single cell example
[~, c] = sort(cFull,'descend');
figure
% for x = 1 : length(c)
    x = 12;
    if c(x) <= size(allSpikes{1},2)
        cTime = (1/sRate : 1/sRate : size(allSpikes{1},1)/sRate)/60; %time in minutes
        plot(cTime, smooth(allSpikes{1}(:,c(x)),10),'linewidth',2,'color',[0.5 0.5 0.5]); hold on
        plot(cTime, smooth(Vm{1}(:,c(x)),10),'linewidth',2,'color','k'); hold off
    else
        cTime = (1/sRate : 1/sRate : size(allSpikes{2},1)/sRate)/60; %time in minutes
        plot(cTime, smooth(allSpikes{2}(:,c(x)-size(allSpikes{1},2)),10),'linewidth',2,'color',[0.5 0.5 0.5]); hold on
        plot(cTime, smooth(Vm{2}(:,c(x)-size(allSpikes{1},2)),10),'linewidth',2,'color','k'); hold off
    end
%     pause
% end
xlim([cTime(1) cTime(end)]); xlabel('time (min)');legend({'Raw data' 'Model'});
ylabel('Firing rate (Hz)');

%% correlation between neurons
[~, c] = sort(cFull,'descend');
% [~, c] = sort(cat(1,allDepths{:}),'ascend');
R1 = corrcoef(allSpikes{1}(:,c(c <= size(allSpikes{1},2))));
R2 = corrcoef(allSpikes{2}(:,c(c > size(allSpikes{1},2))-size(allSpikes{1},2)));

figure
subplot(1,2,1); imagesc(abs(R1));axis square; caxis([0 0.5]); title('NP9');xlabel('Neurons (sorted)');ylabel('Neurons (sorted)');
subplot(1,2,2); imagesc(abs(R2));axis square; caxis([0 0.5]); title('N14');xlabel('Neurons (sorted)');ylabel('Neurons (sorted)');

