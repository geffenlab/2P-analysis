% FRA tone analysis

% load processed file
clear
mouse = 'K056';
date = '20170424';
exptNo = '1';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat']) % load the processed file
saveSig = 1;

% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces => CALCIUM TRACES
% st = {dat.stat(n).st}; % spike times
% sa = {dat.stat(n).c}; % spike amplitudes (relative)

% load events
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\']; % CHANGE TO WHERE YOUR FILES ARE
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr); % WHEN STIMULI COME ON IN FRAMES

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


% order it by frequency and attenuation
p = stimInfo.index(stimInfo.order,:);
order = repmat(p',1,floor(length(eventsOn)/length(p)))'; % ORDER OF STIM PRESENTATION






