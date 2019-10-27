[dataOverview, motorLabels, sensorLabels, cogLabels, ~, segLabels] = delayDecRecordings;
cPath = [pwd filesep 'Widefield' filesep]; %path to widefield data
load allenDorsalMapSM

varFracThresh = 0.9;
allPlots = 1;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
iAnimals = 14; % for figure, used mSM43 23-Nov-2017 (iAnimals = 14), with varFracThresh = 0.9

fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
load([fPath 'betaBhvVid'], 'cam1', 'cam2', 'stdCam1', 'stdCam2', 'meanCam1', 'meanCam2'); %load video data
load([fPath 'interpVc'], 'Vc'); %load model raw data
load([fPath 'Vc'], 'U');
load([fPath 'mask'], 'mask');
load([fPath 'opts2'], 'opts');

cam = cat(1, cam1, cam2);
stdCam = cat(1, stdCam1, stdCam2);
meanCam = cat(1, meanCam1, meanCam2);
[maps, vid] = delayDec_pcaWFVidBetas(cam, stdCam, meanCam, U, mask, opts.transParams, dorsalMaps, varFracThresh, allPlots);

