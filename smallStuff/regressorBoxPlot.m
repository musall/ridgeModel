function [ax, allInd] = regressorBoxPlot(cData,recLabels,tickDist,ax,aColor,idx,bWidth,yWidth)

if ~exist('tickDist', 'var') || isempty(tickDist)
    tickDist = 10;
end

if ~exist('ax', 'var') || isempty(ax)
    ax = [];
    figure;
    ax = subplot(2,2,1:2);
end

aLine = 3;

if ~exist('aColor', 'var') || isempty(aColor)
    aColor = 'k';
    aLine = 3;
else
end

if ~exist('idx', 'var')
    idx = [];
end

if ~exist('bWidth', 'var')
    bWidth = 1;
end

if ~exist('yWidth', 'var')
    yWidth = [0 0.5];
end

showRec = false;
cData = squeeze(cData);
dMean =  nanmedian(cData,1);
dSem =  nansem(cData,1);
if isempty(idx)
    [~,idx] = sort(dMean,'ascend');
    allInd = [];
elseif iscell(idx)
    allInd = idx;
    idx = cat(2,allInd{:});
    showRec = true;
elseif ~islogical(idx)
    mods = unique(idx(idx > 0));
    for x = 1:length(unique(idx(idx > 0)))
        [~, mIdx{x}] = sort(dMean(idx == mods(x)),'ascend');
        allInd{x} = find(idx == mods(x));
        allInd{x} = allInd{x}(mIdx{x});
    end
    idx = cat(2,allInd{:});
    showRec = true;
end
axisLabel = recLabels(idx);
    
hold on
% % dPts = errorbar(dMean(idx),dSem(idx),'color','k','Marker','.','linestyle','none','lineWidth',aLine,'MarkerSize',10);
% dPts = errorbar(dMean(idx),dSem(idx),'k-','linestyle','none','lineWidth',aLine);
% bar(dMean(idx),'FaceColor',aColor,'EdgeColor','k','BarWidth',bWidth,'LineWidth',2);

boxplot(cData(:,idx), 'Color',aColor)

% if length(ax.Children) <= 2
    set(ax,'xTick',1:size(dMean,2))
    set(ax,'xTickLabel',axisLabel)
    set(ax,'XTickLabelRotation',45)
    ax.TickLength = [0 0];%
%     [~,b] = histcounts(dMean(idx),tickDist);
    %     set(ax,'yTick',(unique([fliplr(dPts.YData(1):-mean(diff(b)):ax.YLim(1)) dPts.YData(1):mean(diff(b)):ax.YLim(2)])));
    set(ax,'XLim',[0 length(axisLabel)+1])
    grid on;
% end
ax.YLim(1) = yWidth(1);
ax.YLim(2) = yWidth(2);

if showRec
    rColor = colormap('lines');
    a = get(gca,'YLim');
    rCnt = 0;
    for x = 1 : length(allInd)
        y(x) = rectangle('Position',[rCnt+0.5, a(1), length(allInd{x}), diff(a)],'EdgeColor',rColor(x,:),'LineWidth',3);
        rCnt = rCnt+length(allInd{x})+0.05;
    end
    uistack(y,'bottom');
end
