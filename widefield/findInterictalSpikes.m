function iiSpikeFrames = findInterictalSpikes(U, V, padFrames, makePlot, iiMinProm, iiMinFrames, iiMaxFrames)
% iiSpikeFrames = findInterictalSpikes(U, V [, padFrames] [, makePlot] [, iiMinProm] [, iiMinFrames] [, iiMaxFrames])
% 
% Find interictal events. We do so by taking the mean over the whole brain,
% and looking for peaks of large prominence (height relative to flanking
% points) and specific width. By default, prominence must be >0.03 dF/F,
% and the width must be 2-7 frames.
% 
% As inputs, requires U and V. V may be dims x times x trials, or dims x
% times*trials.
% 
% The output is a vector of logicals of size 1 x times*trials. It is 1 for
% the duration of each peak that corresponds to an interictal event,
% including the minimum points on both sides and padFrames extra points on
% either side. Default padFrames = 2.
% 
% If makePlot = 1 (default), a plot is produced of peak prominence vs.
% width. A rectangle shows the interictal event identification window, with
% identified interictal events shown in red.


%% Optional arguments

if ~exist('padFrames', 'var')
  padFrames = 2;
end

if ~exist('makePlot', 'var')
  makePlot = 1;
end

if ~exist('iiMinProm', 'var')
  iiMinProm = 0.03;
end

if ~exist('iiMinFrames', 'var')
  iiMinFrames = 2;
end

if ~exist('iiMaxFrames', 'var')
  iiMaxFrames = 7;
end


%% Error checking

if ismatrix(V) && size(V, 1) > 2000 && size(V, 2) <= 2000
  warning('V is probably not the right shape! Expecting dims x times x trials, or dims x times*trials');
end


%% Compute the trace

if ndims(V) == 3
  Vr = reshape(V, size(V, 1), []);
else
  Vr = V;
end

mask = isnan(U(:, :, 1));
Ur = arrayShrink(U, mask, 'merge');
meanTrace = mean(Ur, 1) * Vr;


%% Find the interictal events

[~, peakLocs, widths, prominences] = findpeaks(meanTrace);

iiEvents = (prominences >= iiMinProm & widths >= iiMinFrames & widths <= iiMaxFrames);


%% If requested, make the figure

if makePlot
  figure;
  hold on;
  
  plot(widths(~iiEvents), prominences(~iiEvents), 'b.');
  plot(widths(iiEvents), prominences(iiEvents), 'r.');
  
  iiMinProm = min([iiMinProm max(prominences)]);
  rectangle('Position', [iiMinFrames, iiMinProm, iiMaxFrames-iiMinFrames, max(prominences)-iiMinProm]);
  xlabel('Peak width (frames)');
  ylabel('Peak prominence (dF/F)');
  title('Whole-brain mean');
end


%% Produce the output
% a logical trace of the interictal events, padded

iiSpikeFrames = false(1, size(Vr, 2));

if any(iiEvents) % if any event was detected
    iiPeakLocs = peakLocs(iiEvents);
    
    % The edges of peaks are where the trace goes from having a negative
    % derivative to a positive one.
    d = sign(diff(meanTrace));
    if any(d == 0)
        warning('Algorithm problem: some flat spots. Fixable, but not implemented');
    end
    d2 = diff(d);
    edges = find(d2 > 0) + 1;
    edges(end + 1) = length(meanTrace);
    
    % Loop through each peak edge, and look at the peak following that edge. If
    % it's one of the interictal spikes, mark it.
    p = 1;
    for e = 1:length(edges) - 1
        % Find the interictal event following this edge
        while iiPeakLocs(p) < edges(e)
            p = p + 1;
            % When we run out of interictal events, bail out
            if p > length(iiPeakLocs)
                break;
            end
        end
        if p > length(iiPeakLocs)
            break;
        end
        
        % If this peak is an interictal event, mark it and the padFrames frames
        % on either side
        if iiPeakLocs(p) >= edges(e) && iiPeakLocs(p) <= edges(e + 1)
            iiSpikeFrames(max(edges(e) - padFrames, 1) : min(edges(e+1) + padFrames, length(meanTrace))) = 1;
        end
    end
end
