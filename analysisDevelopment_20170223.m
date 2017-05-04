% load processed file
clear
mouse = 'K056';
date = '20170411';
exptNo = '2';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

%%
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces



% c = {dat.stat(n).c}; % calcium traces
% st = {dat.stat(n).st}; % deconvolved spikes


% npilCoef =[dat.cl.dcell{n-length(dat.stat0)}]; npilCoef = [npilCoef.B];
% baseline = npilCoef(2,:); npilCoef = npilCoef(3,:);
% cat = bsxfun(@minus,traces,baseline')-bsxfun(@times,npilTraces,npilCoef');
%
% F0 = mean(cat(:,1:300),2);
%
% dff = bsxfun(@rdivide,bsxfun(@minus,cat,F0),F0*100);
%
% spikeInfo = {dat.cl.dcell(n-length(dat.stat0))};

%% Load the stimulus info

% stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([dataLoc exptNo '\*exptInfo.mat']);
load([stimLoc files.name])

% Work out frame rate
fs=400000;
dt = diff(y(:,1));
ev = find(dt>0.5)+1;
evT = diff(ev);
evR = evT(evT>200);
evR = evR/fs;
plot(evR,'.')
fr = 1/mean(evR); % For now use mean but could interpolate traces to make frames equal length
data.averageFR = fr;
data.expt = files(1).name;
data.frameTimes = cumsum(evR);
if length(npilSubTraces)<length(data.frameTimes)
    data.frameTimes = data.frameTimes(1:length(npilSubTraces));
else
    npilSubTraces = npilSubTraces(:,1:length(data.frameTimes));
end
data.caTraces = npilSubTraces;

% fr = 29.534;
% Work out events
s=y(:,2); % voltage trace of sound
conversionFactor = (fs/fr);
ds = [s(1:conversionFactor:length(s)); zeros(size(traces,2)-length(1:conversionFactor:length(s)),1)]';
eventsOn = [find(diff(ds) > 0.3)+1];
eventsOff = find(diff(ds) < -0.3)+1;
data.stimOn = eventsOn';

%% Make raster
Wn = 20/fr;
[b,a] = butter(5,Wn,'low');
preEv = floor(1*fr); % 1 second before stim
postEv = ceil(5*fr); % 3 seconds post stim
% postEv = ceil(8*fr);
% stimDur = mode(eventsOff(2:end)-eventsOn(1:end));
% eOff = eventsOn+stimDur;
ind=1; raster = zeros(size(npilTraces,1),preEv+postEv+1,length(eventsOn));
for ii=1:length(eventsOn) % ignore 1st event because will be the start recording event - events from then on are relevant
    raster(:,:,ind) = npilSubTraces(:,eventsOn(ii)-preEv: eventsOn(ii)+postEv);
    raster(:,:,ind) = (raster(:,:,ind)-mean(squeeze(raster(:,1:preEv-1,ind)),2))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
%     raster(:,:,ind) = (raster(:,:,ind))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
    for jj=1:size(raster,1)
        raster(jj,:,ii) = filtfilt(b,a,double(raster(jj,:,ii)));
    end
    ind = ind+1;
end
preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1));
raster = permute(raster,[3,2,1]);

%% spike raster
spikeTraces = zeros(size(npilSubTraces));
for ii=1:size(spikeTraces,1)
    spikeTraces(ii,st{ii})=c{ii};
end

ind=1; spikeRaster = zeros(size(npilTraces,1),preEv+postEv+1,length(eventsOn));
for ii=1:length(eventsOn) % ignore 1st event because will be the start recording event - events from then on are relevant
%     spikeRaster(:,:,ind) = spikeTraces(:,eventsOn(ii)-preEv: eventsOn(ii)+postEv);
    raster(:,:,ind) = (raster(:,:,ind)-mean(squeeze(raster(:,1:preEv-1,ind)),2))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
%     spikeRaster(:,:,ind) = (spikeRaster(:,:,ind))./std(squeeze(spikeRaster(:,1:preEv-1,ind)),[],2);
%     for jj=1:size(spikeRaster,1)
%         spikeRaster(jj,:,ii) = filtfilt(b,a,double(spikeRaster(jj,:,ii)));
%     end
    ind = ind+1;
end

spikeRaster = permute(spikeRaster,[3,2,1]);


%% Load stimulus info
stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
load([stimLoc files.name])

% data.stimInfo = stimInfo;
if exist('attnOrder','var')
    index(:,2) =70-index(:,1);
    order = index(attnOrder,2);
    order = repmat(order,4,1);
    data.stimOrder = order;
elseif exist('toneOrder','var')
    order = index(toneOrder,1);
    order = repmat(order,4,1);
    data.stimOrder = order;
elseif exist('params','var')
    order = repmat([1,2]',length(eventsOn)/2,1);
elseif exist('stimInfo','var')
    p = stimInfo.index(stimInfo.order,:);
    order = repmat(p',1,length(eventsOn)/length(stimInfo.order))';
else
    order = repmat(stimInfo.index(stimInfo.order,:),2,1); 
end
[x,ind]=sortrows(order(:,1));
raster = raster(ind,:,:);
spikeRaster = spikeRaster(ind,:,:);

%% basic plot
Wn = 20/fr;
[b,a] = butter(5,Wn,'low');
doPlot = 1;
for ii=1:size(raster,3)
    
    t = double(mean(squeeze(raster(:,:,ii)),1));
    ts = double(std(squeeze(raster(:,:,ii)),[],1));
    h(ii) = ttest2(t(1:preEv)',t(preEv+1:preEv+preEv)','alpha',0.05/5);
    if doPlot
        subplot(1,2,1)
    plot(-preEv/fr:1/fr:(postEv)/fr,t,'LineWidth',2)
    time = -preEv/fr:1/fr:(postEv)/fr;
    axis tight
    hold on
    patch([time fliplr(time)],[(t+ts) fliplr(t-ts)],'b','FaceAlpha',0.3,'edgealpha',0)
    plot([0 0],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
    plot([0.2 0.2],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
%     plot([0.05 0.05],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
%     plot([2.5 2.5],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
%     plot([5 5],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
    xlabel('Time (s)')
    ylabel('\DeltaF/F_0')
    subplot(1,2,2)
    imagesc(raster(:,:,ii));
    % subplot(1,2,2)
    % sc = ones(1,10);
    % for jj=1:size(raster,1)
    %     plot((-preEv/fr:1/fr:(postEv)/fr),convn(raster(jj,:,ii)+((jj-1)*2),sc,'same'))
    %     hold on
    % end
    pause()
    clf
    end
end
ha = h;
BF = ones(1,length(h));
%% FOR RLF
cs = get(groot,'DefaultAxesColorOrder');
uAtt = unique(order);
leg = {'100 dB','90 dB','80 dB','70 dB'};
plotFig = 0;
psth=[]; sempsth=[];
for jj = 1:size(raster,3)
    for ii=1:length(uAtt)
        rows = find(x==uAtt(ii));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);     
    end
      h(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.noiseDur/1000 stimInfo.noiseDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        legend(leg)
        axis tight
        subplot(1,2,2)
        imagesc(raster(:,:,jj))
        pause()
        clf
    end
end

ha = any(h);


% best level
for ii=1:size(psth,3)
    peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
    BF(ii) = uAtt(peakResp==max(peakResp));
end

%% FOR quick tonotopy
uT = unique(order(:,1));
leg = {'5 kHz','8.891 kHz','15.811 kHz','28.117 kHz','50 kHz'};
cs = get(groot,'DefaultAxesColorOrder');
cs=[cs;cs;cs;cs];
plotFig = 1;
psth=[]; sempsth=[];
for jj = 1:size(raster,3)
    for ii=1:length(uT)
        rows = find(x==uT(ii));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);     
    end
       h(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
%        h_ttest(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.tDur/1000 stimInfo.tDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
%         legend(leg)
        axis tight
        subplot(1,2,2)
        imagesc(raster(:,:,jj))
        pause()
        clf
    end
end

ha = any(h);

% best frequency
for ii=1:size(psth,3)
    peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
    BF(ii) = uT(peakResp==max(peakResp));
end


%% For contrast gain

cs = get(groot,'DefaultAxesColorOrder');
uStim = unique(order);

leg = {'low','high'};
plotFig = 1;
psth=[]; sempsth=[];
for jj = 1:size(raster,3)
    for ii=1:length(uStim)
        rows = find(x==uStim(ii));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
%         mr(ii,:,jj) = mean(raster(rows,preEv+1:preEv*9,jj),2);
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);  
         spike_psth(ii,:,jj) = (mean(spikeRaster(rows,:,jj)));
        spike_sempsth(ii,:,jj) = std(spikeRaster(rows,:,jj))./sqrt(length(rows)-1);   
    end
   
       h(jj) = ttest2(preRast(jj,:),raster(1,preEv+1:preEv*9,jj))
       difT = psth(1,:,jj)-psth(2,:,jj);
       if h(jj)==1
           dir(jj) = (mean(difT(1:preEv))-mean(difT(preEv+1:preEv+preEv))>0)+1
       end
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(2,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
%         plot(time,psth(1,:,jj)-psth(2,:,jj),'-k','LineWidth',2)
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
         axis tight
          legend(leg)
         
        subplot(2,2,3)
        for ii=1:size(spike_psth,1)
        plot(time,smooth(spike_psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
%         plot(time,spike_psth(1,:,jj)-spike_psth(2,:,jj),'-k','LineWidth',2)
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(spike_psth(ii,:,jj)+spike_sempsth(ii,:,jj))' flipud(smooth(spike_psth(ii,:,jj)-spike_sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = spike_psth(:,:,jj); b = spike_sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        
%         plot([stimInfo.noiseDur/1000 stimInfo.noiseDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
       
        axis tight
        subplot(2,2,2)
        imagesc(raster(:,:,jj))
        
        subplot(2,2,4)
        imagesc(spikeRaster(:,:,jj))
        pause()
        clf
    end
end

ha = any(h);
ha=h;

% best level
% for ii=1:size(psth,3)
% %     peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
% %     minResp = min(psth(:,preEv+1:preEv*2,ii),[],2);
%     BF(ii) = uAtt(peakResp==max(peakResp));
% end

BF = dir;
%% for contrast gain spikes
cs = get(groot,'DefaultAxesColorOrder');
uStim = unique(order);

leg = {'low','high'};
plotFig = 1;
psth=[]; sempsth=[];
for jj = 1:size(spikeRaster,3)
    for ii=1:length(uStim)
        rows = find(x==uStim(ii));
        psth(ii,:,jj) = (mean(spikeRaster(rows,:,jj)));
%         mr(ii,:,jj) = mean(raster(rows,preEv+1:preEv*9,jj),2);
        sempsth(ii,:,jj) = std(spikeRaster(rows,:,jj))./sqrt(length(rows)-1);     
    end
%        h(jj) = ttest2(preRast(jj,:),spikeRaster(1,preEv+1:preEv*9,jj))
%        difT = psth(1,:,jj)-psth(2,:,jj);
%        if h(jj)==1
%            dir(jj) = (mean(difT(1:preEv))-mean(difT(preEv+1:preEv+preEv))>0)+1
%        end
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
%         plot(time,psth(1,:,jj)-psth(2,:,jj),'-k','LineWidth',2)
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
%         plot([stimInfo.noiseDur/1000 stimInfo.noiseDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        legend(leg)
        axis tight
        subplot(1,2,2)
        imagesc(spikeRaster(:,:,jj))
        pause()
        clf
    end
end

ha = any(h);
ha=h;

% best level
% for ii=1:size(psth,3)
% %     peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
% %     minResp = min(psth(:,preEv+1:preEv*2,ii),[],2);
%     BF(ii) = uAtt(peakResp==max(peakResp));
% end

BF = dir;


%% For FRA
uStim = unique(order(:,1));
leg = {'5 kHz','6.74 kHz','9.086 kHz','12.247 kHz','16.51 kHz','22.255 kHz','30 kHz'};
cs = get(groot,'DefaultAxesColorOrder');
plotFig = 1;
psth=[]; sempsth=[];
for jj = 1:size(raster,3)
    for ii=1:length(uStim)
        rows = find(x(:,1)==uStim(ii));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);     
    end
       h(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
%       h_ttest(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.tDur/1000 stimInfo.tDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        legend(leg)
        axis tight
        subplot(1,2,2)
        imagesc(raster(:,:,jj))
        pause()
        clf
    end
end

ha = any(h);

% best frequency
for ii=1:size(psth,3)
    peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
    BF(ii) = uT(peakResp==max(peakResp));
end






%% for FM tones
cs = get(groot,'DefaultAxesColorOrder');
cs = [cs;cs;cs];
uStim = unique(order);
leg = uStim;
plotFig = 1;
psth=[]; sempsth=[]; h=[];
for jj = 1:size(raster,3)
    for ii=1:length(uStim)
        rows = find(x==uStim(ii));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);     
    end
      h(:,jj) = ttest(psth(:,1:preEv,jj)',psth(:,preEv+1:preEv+preEv,jj)','alpha',0.05/5);
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
        plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        hold on
        end
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        a = psth(:,:,jj); b = sempsth(:,:,jj);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.stimDur/1000 stimInfo.stimDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
%         legend(leg)
        axis tight
        subplot(1,2,2)
        imagesc(raster(:,:,jj))
        pause(0.5)
        clf
    end
end

ha = any(h);


% best level
for ii=1:size(psth,3)
    peakResp = max(psth(:,preEv+1:preEv*2,ii),[],2);
    BF(ii) = uStim(peakResp==max(peakResp));
end


%%
% figure
% for ii=1:4:length(n)
%     for jj=1:4
%         subplot(2,2,jj)
% %         plot(1:7776,traces(ii+(jj-1),1:7776)-mean(traces(ii+(jj-1),1:7776))+100,'.')
%         t = traces(ii+(jj-1),:)-mean(traces(ii+(jj-1),:))+100;
%         patch([1:length(traces) length(traces):-1:1],[t fliplr(t)],'b','FaceAlpha',0,'edgeAlpha',1)
% %         subplot(2,4,jj+4)
%         hold on
%         spikeTimes = spikeInfo{1}{ii+(jj-1)}.st;
%         spikeAmps = spikeInfo{1}{ii+(jj-1)}.c;
%         y = zeros(1,length(dff)); y(spikeTimes) = spikeAmps*100;
%         h = bar(1:length(traces),y(1:length(traces)));
%         h.FaceAlpha = 0.5; h.FaceColor = 'm'; h.EdgeAlpha = 0;
%         h.LineWidth = 1; h.BarWidth = 30;
% %         xlim([0 7776])
%     end
%     pause()
%     clf
% end

%%

f1 = figure;
co = colormap ('jet');
colormap gray
im = dat.mimg(:,:,2);
im=brighten(im,1);
imagesc(im);
uStim = unique(order(:,1));
% ci = floor(length(co)/length(uT));
cs = get(groot,'defaultAxesColorOrder');
cs=[cs;cs;cs;cs];
for ii = 1:length(n)
    img = zeros(1,size(im,1)*size(im,2));
    img(dat.stat(n(ii)).ipix)=1;
    img = reshape(img,size(im,1),size(im,2));
    c = regionprops(img,'Centroid');
    cent(ii,:) = round(c(end).Centroid);
    celln = num2str(ii);
    % t = text(c(1)-5,c(2)+7,celln); t.Color = [1 1 0];
    s = [bwboundaries(img)];
    boundaries(ii) = s(1);
    hold on
    if ha(ii)==1
        plot(boundaries{ii}(:,2),boundaries{ii}(:,1),'Color',cs(find(BF(ii)==uStim),:),'LineWidth',1.5)
    else
        plot(boundaries{ii}(:,2),boundaries{ii}(:,1),'Color','k','LineWidth',1.5)
    end
end
img = reshape(img,size(im,1),size(im,2));
axis off
axis equal
data.cellCentroids = cent;
data.cellOutlines = boundaries;
data.image = im;

%

%% plot imagesc with coloured blobs where cells are

img = zeros(1,size(im,1)*size(im,2));
for ii = 1:length(n)
    if ha(ii)==1
        img(dat.stat(n(ii)).ipix)=BF(ii);
    else
        img(dat.stat(n(ii)).ipix)=0;
    end
end
img = reshape(img,size(im,1),size(im,2));


imagesc(img)
h = colorbar;
axis off
axis equal




% %%
% figure
% fr = 30;
% ind = 1;
% t1 = 1;
% for ii = 1:size(cat,1)
%     if mod(ii,10)==0
%         savefig(['/Users/kath/Documents/ImagingData/K019AnlaysisForMariaGrant/traces_' sprintf('%d',t1)])
%         t1 = t1+1;
%         close all
%         figure
%         ind = 1;
%     end
%     plot(1/30:1/30:300,(SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10))/max((SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10)))+ind*2,'k')
%     hold on
%     ind = ind+1;
%     disp(ii)
% end
%
% figure
% for ii = 1:size(cat,1)
%     plot(1/30:1/30:300,(SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10))/max((SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10)))+ii*2,'k')
%     hold on
% end
%  savefig(['/Users/kath/Documents/ImagingData/K019AnlaysisForMariaGrant/alltraces.fig'])
%


