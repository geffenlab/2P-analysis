function [raster_out,x,fr] = getToRasters(cmSessions,mouse,dataLoc,cm,ss,preEv,postEv)


load([dataLoc mouse '\' cmSessions(cm).name])
date1 = cell2mat(cellMatching.session{ss}(2)); 
exptNo = cell2mat(cellMatching.session{ss}(3));
dataLoc2 = ['E:\dataAnalysed\' mouse '\' date1 mouse '_tifStacks\'];
load([dataLoc2 exptNo '\F_' mouse '_' date1 mouse '_tifStacks_plane1_proc.mat'])


% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)


% load events
stimLoc=['C:\data\' mouse '\' date1 mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);

% load stim info
stimLoc=['C:\data\' mouse '\' date1 mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
load([stimLoc files(1).name])

stimDur = stimInfo.tDur/1000;
ITI = stimInfo.ITI/1000;


% Make raster
Wn = 14/fr/2;
[b,a] = butter(5,Wn,'low');
preEv = floor(preEv*fr); % 
postEv = floor(postEv*fr); % 
raster = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,1);
rawRast = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,0);
preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1)-1);

% Make spike 'traces' and raster
spikeTraces = zeros(size(npilSubTraces));
for ii=1:size(spikeTraces,1)
    spikeTraces(ii,st{ii})=sa{ii};
end
spikeRaster = makeCaRaster(spikeTraces,eventsOn,preEv,postEv,0);

% order it by frequency and attenuation
p = stimInfo.index(stimInfo.order,:);
order = repmat(p',1,floor(length(eventsOn)/length(stimInfo.order)))';
[x,ind]=sortrows(order);
raster = raster(ind,:,:);

index = cellMatching.index(~isnan(cellMatching.index(:,ss)),ss);
index2 = find(~isnan(cellMatching.index(:,ss)));
raster_out = NaN(size(raster,1),size(raster,2),length(cellMatching.index));
for ii = 1:length(index2)
    raster_out(:,:,index2(ii)) = raster(:,:,index(ii));
end
% rawRast = rawRast(ind,:,:);
% spikeRaster = spikeRaster(ind,:,:);