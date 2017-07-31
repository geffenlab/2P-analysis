% make data files

% load suite2P processed file
clear
mouse = 'K048';
date = '20170315';
exptNo = '2';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
proc = dir([dataLoc exptNo '\*proc.mat']);
load([dataLoc exptNo '\' proc.name])

data.mouse = mouse;
data.exptDate = date;
data.exptNo = exptNo;
data.suite2P_procFile = [dataLoc exptNo '\' proc.name];

% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)

data.rawTraces = traces;
data.npilTraces = npilTraces;
data.npilCoeffs = npilCoeffs;
data.npilSubTraces = npilSubTraces;
data.spikeTimes = st;
data.spikeAmps = sa;

% load events
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

data.stimInfoLoc = [stimLoc files(1).name];

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);
data.frameRate = fr;
data.stimOn = eventsOn;

% load stim info
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
if ~isempty(files)
    load([stimLoc files(1).name])
else
    stimInfo = events.exptInfo.stimInfo{1};
end

stimDur = stimInfo.tDur/1000;
ITI = stimInfo.ITI/1000;
data.stimDur = stimDur;
data.ISI = ITI;


% Make raster
preEv = floor(1*fr); % 1 second before stim
postEv = ceil((stimDur+ITI)*fr); % x seconds post stim
data.nFramesPreStim = preEv;
data.nFramesPostStim = postEv;
raster = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,1);
data.raster_DFF0 = raster;


% get stim presentation info
p = stimInfo.index(stimInfo.order,:);
% p = index(toneOrder);
order = repmat(p',1,floor(length(eventsOn)/length(p)))';
data.stimIndex = order;


% Is response significant?
wind = round(2*fr); % 2 seconds
pr = zeros(1,length(n)); h = pr;
for ii=1:length(n)  
    [pr(ii),h(ii)]=signrank(mean(raster(:,preEv:preEv+wind,ii),2),mean(raster(:,1:preEv-1,ii),2));       
end
data.sigResp = h;
exptNo = sprintf('%02d',str2double(exptNo));

save(['E:\dataAnalysed\Kath\bathellier\' date '_' mouse '_' exptNo '.mat'],'data')
