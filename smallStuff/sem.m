function avec=sem(avec,varargin)
% smusall 2010/1/18
% usage: sem(avec,dim)
% sem of a vector: std divided by sqrt of length, 
% for matrices, sem works on columns by default or can be adjusted by dim 

if size(varargin,1)==1
    dim=varargin{1};
else
    dim=1;
end

samples=size(avec,dim);
if samples==1 && size(varargin,1)==0
    dim=2;
    samples=size(avec,dim);
end

avec=nanstd(avec,0,dim)/sqrt(samples);
if isempty(avec)
    avec=NaN;
end
