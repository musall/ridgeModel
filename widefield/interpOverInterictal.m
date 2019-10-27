function V = interpOverInterictal(V, iiSpikeFrames)
% V = interpOverInterictal(V, iiSpikeFrames)
% 
% For use with widefield data that contains interictal spikes. This
% function excises the times around interictal events, and patches over
% them using the Matlab built-in fillgaps.
% 
% Interictal events should be found using findInterictalSpikes and the
% result supplied as iiSpikeFrames.
% 
% A bit slow, ~75s on Matt's laptop.


% How many frames (around each interictal event) to use in order to patch
% over the hole in the data. Supplied as an argument to fillgaps.
fillWin = 20;


siz = size(V);
V = reshape(V, siz(1), []);

V(:, iiSpikeFrames) = NaN;
V = double(reshape(V, siz));

for d = 1:size(V, 1)
  for tr = 1:size(V, 3)
    V(d, :, tr) = fillgaps(V(d, :, tr), fillWin);
  end
end

V = single(V);
