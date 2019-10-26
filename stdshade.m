% function cLine = stdshade(amatrix,acolor,F,alpha,smth,varargin)
function stdshade(amatrix,alpha,acolor,F,smth,varargin)
% usage: stdshade(amatrix,acolor,F,alpha,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - amatrix is a data matrix (observations,datalength)
% - acolor defines the used color (default is red)
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
%
% Edit: Optionally added possibility to define amean and astd as additional
% input: stdshade(amatrix,acolor,F,alpha,smth,amean,astd)
%

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r';
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1;
end

if ne(size(F,1),1)
    F=F';
end

if length(varargin)==2
    amean=smooth(varargin{1},smth)';
    astd=smooth(varargin{2},smth)';
else
    if size(amatrix,1) == 1 || size(amatrix,2) == 1
        amean(1,:) = amatrix;
        astd = zeros(1,length(amean));
        F = 1 : length(amean);
    else
        amean=smooth(nanmean(amatrix,1),smth)';
        astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading
        astd(isnan(astd)) = 0;
    end
    %     astd=nanstd(amatrix); % to get std shading
end

lineMean = amean;
amean = fillmissing(amean, 'previous');

if exist('alpha','var')==0 || isempty(alpha)
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');
end

if ishold==0
    check=true; else check=false;
end

hold on;
cLine = plot(F,lineMean,'color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



