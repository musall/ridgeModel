function dataOut = smoothCol(dataIn, dim, fWidth, fType, fLength)
% function to smooth columns of a data matrix 'dataIn'.
% dataOut = smoothCol(dataIn, dim, fWidth, fType, fLength)
% fType defines the filter type which is either a box or gaussian filter.
% When using a box filter 'fWidth' is the size of the moving average, when
% using a gaussian filter 'fWidth' defines its full width half maximum.
% dim is dimension that is smooth over (default is 1st dimension).
% Default filter is a 5pt box filter when only 'dataIn' is provided.
% ----------------------------------------
% To make the code run faster when using the gaussian filter, you can also
% use a shorter filter length. Default is 100*sigma.
%
% Simon Musall
% Cold Spring Harbor, 8/15/2019

if ~exist('fWidth','var') || isempty(fWidth)
    fWidth = 5; %default filter length is 5 dpoints
end

if ~exist('fType','var') || isempty(fType)
    fType = 'box'; %default is box filter
end

if ~exist('dim','var') || isempty(dim)
    dim = 1; %default is smoothing first dimension
end

if ~ismember(fType,{'box' 'gauss'})
    error('Unknown filter type. Use "box" or "gauss" for fType')
end

dSize = size(dataIn);
if dim > length(dSize)
    error('Dimension to be smoothed is above the dimensionality of the input data.');
end

% permute dataset to allow smoothing over other dimensions
dimOrder = 1: length(dSize);
if dim ~= 1
    dimOrder(dim) = [];
    dimOrder = [dim dimOrder];
    dataIn = permute(dataIn, dimOrder);
end

%add buffer to first dimension to avoid edge artefacts
cbegin = repmat(dataIn(1,:,:),[fWidth ones(1,length(dSize)-1)]); %start buffer
cend = repmat(dataIn(end,:,:),[fWidth ones(1,length(dSize)-1)]); %end buffer
dataIn = cat(1,cbegin,dataIn,cend);

dataIn = reshape(dataIn,size(dataIn,1), []); %merge other dimensions for smoothing

if nansum(abs(dataIn(:))) > 0 %check if there is any data available
    if ~strcmpi(class(dataIn),'double') %make sure data is double
        dataIn = double(dataIn);
    end
    
    if strcmpi(fType,'box') %box filter
        n = size(dataIn,1);
        cbegin = cumsum(dataIn(1:fWidth-2,:),1);
        cbegin = bsxfun(@rdivide, cbegin(1:2:end,:), (1:2:(fWidth-2))');
        cend = cumsum(dataIn(n:-1:n-fWidth+3,:),1);
        cend = bsxfun(@rdivide, cend(end:-2:1,:), (fWidth-2:-2:1)');
        dataOut = conv2(dataIn,ones(fWidth,1)/fWidth,'full'); %smooth trace with moving average
        dataOut = [cbegin;dataOut(fWidth:end-fWidth+1,:);cend];
        
    elseif strcmpi(fType,'gauss') %gaussian filter
        fSig = fWidth./(2*sqrt(2*log(2))); %in case of gaussian smooth, convert fwhm to sigma.
        if ~exist('fLength','var') || isempty(fLength)
            fLength = round(fSig * 100); %length of total gaussian filter
        end
        fLength = fLength-1+mod(fLength,2); % ensure kernel length is odd
        kernel = exp(-linspace(-fLength / 2, fLength / 2, fLength) .^ 2 / (2 * fSig ^ 2));
        dataOut = conv2(dataIn,kernel','same'); %smooth trace with gaussian
    end
else
    dataOut = dataIn;
end

dataOut = reshape(dataOut,[size(dataOut,1) dSize(dimOrder(2:end))]); %split dimensions again
dataOut = dataOut(fWidth+1:end-fWidth, : ,:); %remove buffers
if dim ~= 1
    dataOut = ipermute(dataOut, dimOrder);
end
