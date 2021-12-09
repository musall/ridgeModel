function [lineOut, fillOut] = stdshade(amatrix,alpha,acolor,F,smth)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  

if ne(size(F,1),1)
    F=F';
end

amean = nanmean(amatrix,1); %get man over first dimension
if smth > 1
    amean = boxFilter(nanmean(amatrix,1),smth); %use boxfilter to smooth data
end
% astd = nanstd(amatrix,[],1); % to get std shading
astd = nanstd(amatrix,[],1)/sqrt(size(amatrix,1)); % to get sem shading

if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

% check if current axis is on hold
check = false;
if ~ishold
    check=true;
end
hold on;

if any(isnan(amean)) %make multiple patches if there are nans
    nanIdx = find(isnan(amean));
    nanIdx = [0 nanIdx];
    
    for x = 1 : length(nanIdx)
        if x ~= length(nanIdx)
            cIdx = nanIdx(x)+1:nanIdx(x+1)-1; %data for current patch
        else
            cIdx = nanIdx(x)+1:length(amean); %data for current patch
        end
        
        if ~isempty(cIdx)
            cF = F(cIdx);
            cM = amean(cIdx);
            cE = astd(cIdx);
            fillOut(x) = fill([cF fliplr(cF)],[cM+cE fliplr(cM-cE)],acolor, 'FaceAlpha', alpha, 'linestyle','none');
        end
    end
else
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha, 'linestyle','none');
end
if alpha == 1; acolor='k'; end
lineOut = plot(F,amean, 'color', acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end
end


function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing

fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];

end

