% FRA tone analysis

% load processed file
clear
mouse = 'K056';
date = '20170417';
exptNo = '1';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
proc = dir([dataLoc exptNo '\*proc.mat']);
load([dataLoc exptNo '\' proc.name])
saveSig = 1;

% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)
k = {dat.stat(n).kernel};

% load events
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);
% eventsOn = eventsOn(1:end-3);
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

% Make raster
preEv = floor(1*fr); % 1 second before stim
postEv = ceil((stimDur+ITI)*fr); % x seconds post stim

% filter the traces?
% [bf,af] = butter(2,[0.5/fr/2 10/fr/2]); % 10 is good for low pass
[bf,af] = butter(2,10/fr/2,'low'); % 10 is good for low pass
% [bf,af] = butter(2,0.5/fr/2,'high');
for ii=1:size(npilSubTraces,1)
    npst_f(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
    npst_g(ii,:) = SmoothGaus(double(npilSubTraces(ii,:)),1);
    npst_z(ii,:) = zscore(npilSubTraces(ii,:));
end

F0 = min(npst_g,[],2);
df_F0 = (npst_g-F0)./abs(F0);
raster = makeCaRaster(npst_g,eventsOn,preEv,postEv,1);

o = stimInfo.index(stimInfo.order,:);
order = repmat(o',1,floor(length(eventsOn)/length(o)))';
uStim = unique(order,'rows');
uFreq = unique(uStim(:,1));
uAtten = unique(uStim(:,2));
raster  = raster(1:length(order),:,:);

%% significant?
wind = round(2*fr); % 2 seconds

for ii=1:length(n)
    
    [p(ii),h(ii)]=signrank(mean(raster(:,preEv:preEv+wind,ii),2),mean(raster(:,1:preEv-1,ii),2));
        
end

%% selective?
clc
uT = unique(order(:,1));
if size(order,2)>1
    uA = unique(order(:,2));
else
    uA=1;
end

ha = zeros(length(uT),length(uA),length(n));
for ii=1:length(n)
    if h(ii)==1
        disp(ii)
            resp = mean(raster(:,preEv:preEv+wind,ii),2);
           [p,tbl,stats] =  kruskalwallis(resp,order(:,1),'on');
           if p<0.05
               multcompare(stats)
           end
           pause()
           close all
    end
end




%%
for ii=1:size(raster,3)
    if h(ii)==1
    plot(squeeze(mean(raster(:,:,ii),1)))
    title(num2str(ii))
    pause()
    end
end

%%
cells = 1:30;
figure('position',[ 1000         175         637        1163])
time = 1/fr/fr:1/fr:length(npilSubTraces)/fr;
for ii=1:length(cells)
    norm = (df_F0(cells(ii),:)-min(df_F0(cells(ii),:)))/max((df_F0(cells(ii),:)-min(df_F0(cells(ii),:))));
    plot(time,norm+(1.5*(ii-1)),'b','LineWidth',0.5)
    hold on
    xlabel('Time(s)')
    hold on
    
end
% for jj=1:length(eventsOn)
%          plot([time(eventsOn(jj)) time(eventsOn(jj))], [0 (sum(ha)+1)*1.5],'r-')
% end
axis tight
 xlim([0 360])
% ylim([0   68.9102])
set(gca,'YTick',0.75:1.5:sum(ha)*1.5,'YTickLabels',1:sum(ha))
ylabel('cell number')

