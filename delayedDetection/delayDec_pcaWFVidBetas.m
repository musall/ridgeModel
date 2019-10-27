function [proj, vidC] = delayDec_pcaWFVidBetas(cam, stdCam, meanCam, U, mask, transParams, dorsalMaps, varFracThresh, allPlots)
% Trying to find how video affects WF, sparsifying WF
% 
% PCA the camera-to-WF betas, to get common ways the behavioral video
% pixels influence the WF. This produces a handful of WF maps. Use varimax
% to sparsify the WF maps. To see how each of these maps corresponds to
% behavioral video, take the dot product of each WF map with each
% behavioral video pixel. To convert these dot products into measures of
% influence, multiply by the std devs of the pixels.

%% Parameters

% For scaling the video influence frames
maxFrac = 99;
minFrac = 0.1;

%% Reshape betas

Cr = reshape(cam, [], size(cam, 3));
Cr = Cr - mean(Cr, 2);


%% PCA

% PCA to reduce pixels and get strongest WF patterns
[~, scores, latent] = pca(Cr');

% Find dimensionality
nDim = find(cumsum(latent) / sum(latent) > varFracThresh, 1);
varRem = 1 - cumsum(latent) / sum(latent);


%% Compute sparsified WF maps -- project up and then use varimax

Ur = arrayShrink(U, mask, 'merge');
if size(Ur, 2) > size(scores, 1)
  Ur = Ur(:, 1:size(scores, 1));
end
proj = Ur * scores(:, 1:nDim);
[proj, R] = rotatefactors(proj, 'Normalize', 'off');    % sparsify basis


%% Show sparsified WF components with the Allen areas aligned
figure;
mapPlots = min([nDim 6]);
if allPlots
  for d = 1:mapPlots
    im = arrayShrink(proj(:, d), mask, 'split');
    [extreme, extremeI] = max(abs(im(:)));
    im = im * (sign(im(extremeI)) / extreme / 2);
    
    im = alignAllenTransImMasked(im, transParams, 0);
    im = arrayShrink(im(:, 1:size(dorsalMaps.allenMask,2)), dorsalMaps.allenMask,'merge');
    im = arrayShrink(im, dorsalMaps.allenMask,'split');
    
    subplot(2,3,d);
    mapImg = imshow(im,[-0.5 0.5]);
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
    colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
    hold(mapImg.Parent, 'on');

    for p = 1:length(dorsalMaps.edgeOutlineSplit)
      if strcmpi(dorsalMaps.sidesSplit{p}, 'R')
        plot(dorsalMaps.edgeOutlineSplit{p}(:, 2), dorsalMaps.edgeOutlineSplit{p}(:, 1), 'w-', 'LineWidth', 2);
      end
    end
    axis equal off; title(num2str(d));
  end
end


%% show weight projection on facecam
vidC = Cr * scores(:, 1:nDim) * R;
vidC = abs(reshape(vidC, [size(cam(:, :, 1)) size(vidC, 2)]));

backIm = rot90(meanCam(1:size(meanCam,1)/2,:), -1);
foreIm = rot90(sqrt(sum(cam(1:size(meanCam,1)/2,:,:) .^ 2, 3)) .* stdCam(1:size(meanCam,1)/2,:), -1);
theMax = prctile(foreIm(:), maxFrac);
lims = [theMax * minFrac, theMax];
hf = figure;
cBack = (mat2gray(backIm) .* lims(2)) +  lims(1);
backImg = imshow(cBack,[lims(1)*0.5 lims(2)]); hold on; colormap('gray');
freezeColors;
mapImg = imshow(foreIm,lims);
set(mapImg,'AlphaData',mapImg.CData > lims(1)); %make 0s transparent.
set(mapImg,'AlphaData',mat2gray(mapImg.CData)*10); %make 0s transparent.
colormap(mapImg.Parent,hot(256));hold on;
hold(mapImg.Parent, 'on');

%% Scree plot
if allPlots
  figure;
  hold on;
  plot([0 length(latent)], (1 - varFracThresh) * [1 1], 'k--', 'LineWidth', 1);
  plot(0:length(latent), [1; varRem], 'LineWidth', 2);
  plot(nDim, varRem(nDim), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
  xlim([0 20]);
  ylim([0 1]);
  set(gca, 'TickDir', 'out', 'Box', 'off');
  xlabel('Dimensions');
  ylabel('Fraction variance remaining');
end