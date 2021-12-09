function data = Behavior_vidResamp(data, trialTime, targRate)
%quick function to adjust single trial video data if it was taken at
%irregular framerate or is too short. Resamples from current sampling rate
%to target rate and ensures to return the requested number of frames. If
%trialtimes is a single number it is treated to indicate the current
%sampling rate of the video data.

if isvector(data)
    data = reshape(data,length(data),1); %make sure input is column vector
end

if length(trialTime) == 1
    trialTime = 1/trialTime : 1/trialTime : size(data,1)/trialTime;
end

lengthCheck = length(trialTime) - size(data,1);
if lengthCheck > 0
    data = [data; NaN(lengthCheck,size(data,2))]; %fill up missing frames and interpolate
    data = fillgaps(data, length(data)); %use interpolation to fill missing frames
elseif lengthCheck < 0
    error('Less timestamps as frames provided.');
end

trialTime = trialTime - min(trialTime);
timePad = 1/targRate : 1/targRate : 1/targRate*21;
trialTime = trialTime' + timePad(end) + mean(diff(trialTime));
trialTime = [timePad trialTime' timePad+trialTime(end)];
data = [repmat(data(1,:),21,1); data; repmat(data(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
data = resample(double(data), trialTime, targRate); %resample to set target framerate

if isvector(data)
    data = data(22:end-21); %remove pads.
else
    data = data(22:end-21,:); %remove pads.
end