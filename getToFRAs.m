function [spFRA_out,FRA_out,uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,cm,ss,method)


load([dataLoc mouse '\fearCellMatching\' cmSessions(cm).name])
date1 = cell2mat(cellMatching.session{ss}(2)); % start with the pre-fear conditioning recording
exptNo = cell2mat(cellMatching.session{ss}(3));
dataLoc2 = ['E:\dataAnalysed\' mouse '\' date1 mouse '_tifStacks\'];
proc = dir([dataLoc2 exptNo '\*proc.mat']);
load([dataLoc2 exptNo '\' proc.name])


% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces

% load events
stimLoc=['E:\dataAnalysed\' mouse '\' date1 mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);

% process traces
[bf,af] = butter(2,10/fr/2,'low'); % 10 is good for low pass
for ii=1:size(npilSubTraces,1)
    npilSubTraces(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
%     npst_z(ii,:) = zscore(npilSubTraces(ii,:));
end
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)




if strcmp(mouse,'K056') && strcmp(date1,'20170330')
    eventsOn = eventsOn(1:end-4);
end

% load stim info
stimLoc=['E:\dataAnalysed\' mouse '\' date1 mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
if ~isempty(files)
    load([stimLoc files(1).name]);
else
    files  = dir([stimLoc '*exptInfo.mat']);
    load([stimLoc files(1).name]);
    stimInfo = exptInfo.stimInfo{1};
end

stimDur = stimInfo.tDur/1000;
ITI = stimInfo.ITI/1000;


% Make raster
preEv = floor(1*fr); % 1 second before stim
postEv = ceil((stimDur+ITI)*fr); % 3 seconds post stim
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
rawRast = rawRast(ind,:,:);
spikeRaster = spikeRaster(ind,:,:);

% Work out FRAs
window = [preEv+1,preEv+round(fr*1)];
FRA = makeCaFRA(raster,window,x(:,1),x(:,2),method);
spike_FRA = makeSpikeFRA(spikeRaster,window,x(:,1),x(:,2));

FRA_out = NaN(size(FRA,1),size(FRA,2),length(cellMatching.index));
index = cellMatching.index(~isnan(cellMatching.index(:,ss)),ss);
index2 = find(~isnan(cellMatching.index(:,ss)));
for ii = 1:length(index2)
    FRA_out(:,:,index2(ii)) = FRA(:,:,index(ii));
end

spFRA_out = NaN(size(spike_FRA,1),size(spike_FRA,2),length(cellMatching.index));
index = cellMatching.index(~isnan(cellMatching.index(:,ss)),ss);
index2 = find(~isnan(cellMatching.index(:,ss)));
for ii = 1:length(index2)
    spFRA_out(:,:,index2(ii)) = spike_FRA(:,:,index(ii));
end

uT = unique(x(:,1));
uA = unique(x(:,2));