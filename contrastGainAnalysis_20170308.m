%% Load processed data and extract cell traces
clear
mouse = 'K048';
date = '20170307';
exptNo = '5';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:);
npilTraces = dat.FcellNeu{1}(n,:);
npilCoeffs = [dat.stat(n).neuropilCoefficient]';
npilSubTraces = traces-(npilTraces.*npilCoeffs);

data.caTraces = npilSubTraces;

c = {dat.stat(n).c};
st = {dat.stat(n).st};

%% Load the event traces

stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*_events.mat']);
load([stimLoc files(1).name])

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

%% Events and Conditions
% fr = 29.534;
% Work out events
s=y(:,2); % voltage trace of sound
s = [0; s];
onS = find(diff(s) > .3);
offS = find(diff(s) < -.3);
dt = offS - onS;
blockTS = onS(dt == 4001)./fs;
toneTS = onS(dt == 10000)./fs;
startTS = onS(dt < 4000)./fs;

% get nearest frames for each event time
blockFs = knnsearch(data.frameTimes,blockTS);
toneFs = knnsearch(data.frameTimes,toneTS);
startFs = knnsearch(data.frameTimes,startTS);

% block type index (0 = low, 1 = high)
blockInd = repmat([0 1],1,length(blockFs)/2);

% load stimulus info to get target times
file = dir([stimLoc '*_stim.mat']);
stim = load([stimLoc file(1).name]);

% get stimulus order for each block
filesUsed = 1:4;
toneIndex = [];
for i = filesUsed
    [~,ind] = sort(stim.index{i}(:,5));
    toneIndex = [toneIndex; stim.index{i}(ind,6)];
end
toneIndex(isnan(toneIndex)) = 0;



%% Make raster for transitions
Wn = 20/fr;
[b,a] = butter(5,Wn,'low');
preEv = floor(1*fr); % 1 second before stim
postEv = ceil(10*fr); % 3 seconds post stim
% postEv = ceil(8*fr);
% stimDur = mode(eventsOff(2:end)-eventsOn(1:end));
% eOff = eventsOn+stimDur;
ind=1; raster = zeros(size(npilTraces,1),preEv+postEv+1,length(blockFs));
for ii=1:length(blockFs) % ignore 1st event because will be the start recording event - events from then on are relevant
    raster(:,:,ind) = npilSubTraces(:,blockFs(ii)-preEv: blockFs(ii)+postEv);
    %     raster(:,:,ind) = (raster(:,:,ind)-mean(squeeze(raster(:,1:preEv-1,ind)),2))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
    raster(:,:,ind) = (raster(:,:,ind))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
    for jj=1:size(raster,1)
        raster(jj,:,ii) = filtfilt(b,a,double(raster(jj,:,ii)));
    end
    ind = ind+1;
end
preRast = npilSubTraces(:,blockFs(1)-preEv:blockFs(1));
raster = permute(raster,[3,2,1]);

spikeTraces = zeros(size(npilSubTraces));
for ii=1:size(spikeTraces,1)
    spikeTraces(ii,st{ii})=c{ii};
end

ind=1; spikeRaster = zeros(size(npilTraces,1),preEv+postEv+1,length(blockFs));
for ii=1:length(blockFs) % ignore 1st event because will be the start recording event - events from then on are relevant
    spikeRaster(:,:,ind) = spikeTraces(:,blockFs(ii)-preEv: blockFs(ii)+postEv);
    %     raster(:,:,ind) = (raster(:,:,ind)-mean(squeeze(raster(:,1:preEv-1,ind)),2))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
    %     spikeRaster(:,:,ind) = (spikeRaster(:,:,ind))./std(squeeze(spikeRaster(:,1:preEv-1,ind)),[],2);
    %     for jj=1:size(spikeRaster,1)
    %         spikeRaster(jj,:,ii) = filtfilt(b,a,double(spikeRaster(jj,:,ii)));
    %     end
    ind = ind+1;
end

% reshape raster
spikeRaster = permute(spikeRaster,[3,2,1]);

% sort by block condition
conds = [blockInd' toneIndex];
[x,ind]=sortrows(conds);
raster = raster(ind,:,:);
spikeRaster = spikeRaster(ind,:,:);



%% For contrast gain
close all
cs = get(groot,'DefaultAxesColorOrder');
uStim = unique(conds,'rows');

leg = {'low','high'};
plotFig = 1;h=[];
psth=[]; sempsth=[];
spike_psth = [];
spike_sempsth = [];
for jj = 1:size(raster,3)
    for ii=1:length(uStim)
        rows = find(ismember(x,uStim(ii,:),'rows'));
        psth(ii,:,jj) = (mean(raster(rows,:,jj)));
        %         mr(ii,:,jj) = mean(raster(rows,preEv+1:preEv*9,jj),2);
        sempsth(ii,:,jj) = std(raster(rows,:,jj))./sqrt(length(rows)-1);
        spike_psth(ii,:,jj) = (mean(spikeRaster(rows,:,jj)));
        spike_sempsth(ii,:,jj) = std(spikeRaster(rows,:,jj))./sqrt(length(rows)-1);
    end
    
    if plotFig==1
        time = -preEv/fr:1/fr:(postEv)/fr;
        %         subplot(2,2,1)
        %         for ii=1:size(psth,1)
        %         plot(time,smooth(psth(ii,:,jj))','-','LineWidth',2);
        %         hold on
        %         end
        % %         plot(time,psth(1,:,jj)-psth(2,:,jj),'-k','LineWidth',2)
        %         hold on
        %         for ii=1:size(psth,1)
        %             patch([time fliplr(time)],[smooth(psth(ii,:,jj)+sempsth(ii,:,jj))' flipud(smooth(psth(ii,:,jj)-sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        %         end
        %         %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        %         a = psth(:,:,jj); b = sempsth(:,:,jj);
        %         plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        %          axis tight
        %           legend(leg)
        %
        %         subplot(2,2,3)
        %         for ii=1:size(spike_psth,1)
        %         plot(time,smooth(spike_psth(ii,:,jj))','-','LineWidth',2);
        %         hold on
        %         end
        % %         plot(time,spike_psth(1,:,jj)-spike_psth(2,:,jj),'-k','LineWidth',2)
        %         hold on
        %         for ii=1:size(psth,1)
        %             patch([time fliplr(time)],[smooth(spike_psth(ii,:,jj)+spike_sempsth(ii,:,jj))' flipud(smooth(spike_psth(ii,:,jj)-spike_sempsth(ii,:,jj)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        %         end
        %         %      plot(time,psth(:,:,jj)','-','LineWidth',2)
        %         a = spike_psth(:,:,jj); b = spike_sempsth(:,:,jj);
        %         plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        %
        % %         plot([stimInfo.noiseDur/1000 stimInfo.noiseDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        
        offsets = unique(uStim(:,2));
        subplot(2,2,1)
        imagesc(raster(x(:,1)==0,:,jj))
        subplot(2,2,3)
        ind = find(uStim(:,1)==0);
        mx = max(max(psth(:,:,jj)));
        mn = min(min(psth(:,:,jj)));
        colors = [0 0 0];
        colors = [colors; [linspace(.2,1,8)' repmat(.2,8,2)]];
        hold on
        for i = 1:length(ind)
            plot(time,smooth(psth(ind(i),:,jj))+(i*.5*(mx-mn)),'Color',colors(i,:),'LineWidth',2);
            if i > 1
                plot([offsets(i) offsets(i)+.05],...
                    [max(psth(ind(i),:,jj)) max(psth(ind(i),:,jj))]+(i*.25*(mx-mn)),...
                    'k','LineWidth',2);
            end
        end
        %ylim([mn - .1*mn mx + .1*mx]);
        hold off
        axis tight
        
        legend
        
        
        
        subplot(2,2,2)
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


%%

f1 = figure;
co = colormap ('jet');
colormap gray
im = dat.mimg(:,:,2);
im=brighten(im,1);
imagesc(im);
uStim = unique(order);
% ci = floor(length(co)/length(uT));
BF=BF-1;
cs = get(groot,'defaultAxesColorOrder');
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


