% FRA tone analysis

% load processed file
clear
mouse = 'K048';
date = '20170315';
exptNo = '2';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

%% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)

%% load events
stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);

%% load stim info
stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
load([stimLoc files(1).name])

stimDur = stimInfo.tDur/1000;
ITI = stimInfo.ITI/1000;

%% Make raster
Wn = 14/fr/2;
[b,a] = butter(5,Wn,'low');
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


%% Work out which cells respond to sounds
uT = unique(order(:,1));
uA = unique(order(:,2));
resp=zeros(length(n),length(uT));
for jj=1:length(n)
    for ii=1:length(uT)
        stdPre=zeros(1,length(uA)); auc=stdPre;
        for kk=1:length(uA)
            rows = x(:,1)==uT(ii) & x(:,2)==uA(kk);
            stdPre(:,kk) = std(preRast(jj,:),[],2); % pre stimulus          
            auc(:,kk) = mean(mean((rawRast(rows,preEv:preEv*2,jj)')));
        end        
        r1 = sum((auc)>stdPre);
       resp(jj,ii) = r1;
    end
end

ha = sum(resp,2)>0;

%% Work out FRAs
window = [preEv+1,preEv+round(fr*1)];
FRA = makeCaFRA(raster,window,x(:,1),x(:,2));
spike_FRA = makeSpikeFRA(spikeRaster,window,x(:,1),x(:,2));
filt = fspecial('gaussian',3,0.75);
for ii=1:length(n); subplot(1,2,1); imagesc(conv2(FRA(:,:,ii),filt,'same'));hold on; subplot(122); imagesc(conv2(spike_FRA(:,:,ii),filt,'same'));colormap gray; pause();clf;end
% % shuffle the FRAs - pretty sure this doesn't tell you anything useful
% nShuf = 2000;
% for ii=1:length(n)
%     for jj=1:nShuf
%         r = randperm(size(raster,1));
%         FRA_shuf(:,:,jj) = makeCaFRA(raster(r,:,ii),window,x(:,1),x(:,2));
%     end
%     for xx = 1:size(FRA,1)
%         for   yy=1:size(FRA,2)
%             FRA_p(xx,yy,ii) = sum(sort(squeeze(FRA_shuf(xx,yy,:)))>=FRA(xx,yy,ii))/nShuf;
%         end
%     end
% end

%% 1) Plot mean by frequency across attenuations
figure
uT = unique(order(:,1));
leg = num2str(uT);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
plotFig = 1;

index = find(ha==1);
for ff = 44:sum(ha)
    jj = index(ff);
    psth = zeros(length(uT),size(raster,2));
    sempsth = psth;
    for ii=1:length(uT)
        rows = find(x(:,1)==uT(ii));
        psth(ii,:) = (mean(raster(rows,:,jj)));
        sempsth(ii,:) = std(raster(rows,:,jj))./sqrt(length(rows));
    end
    
    if plotFig==1
        time = (-preEv-1)/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
            plot(time,smooth(psth(ii,:))','-','LineWidth',2,'Color',cs(ii,:));
            hold on
        end
        hold on
        for ii=1:size(psth,1)
            patch([time fliplr(time)],[smooth(psth(ii,:)+sempsth(ii,:))' flipud(smooth(psth(ii,:)-sempsth(ii,:)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        end
        a = psth(:,:); b = sempsth(:,:);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.tDur/1000 stimInfo.tDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        legend(leg)
        axis tight
        subplot(1,2,2)
        colormap('gray')
        imagesc(raster(:,:,jj))
        disp(jj)
        pause()
        clf
    end
    
end


%% 2) Plot mean response for different attenuations at each frequency
intWindow = floor(1*fr);
uT = unique(order(:,1));
uA = unique(order(:,2));
attns = [100,90,80,70,60,50];
leg = num2str(attns');
cs = colormap('summer');  close all;
cs=cs(1:floor(length(cs)/length(uA)):floor(length(cs)/length(uA))*length(uA),:);
plotFig = 1;
psth=[]; sempsth=[];
f(1) = figure; f(2) = figure;
h_attn=[];
for jj=69:length(n)
    for ii=1:length(uT)
        psth = zeros(length(uA),size(raster,2))';
        sempsth = psth;
        for kk=1:length(uA)
            rows = find(x(:,1)==uT(ii) & x(:,2)==uA(kk));
            psth(:,kk) = mean(raster(rows,:,jj));
            sempsth(:,kk) = std(raster(rows,:,jj))./sqrt(length(rows)-1);
        end
        
        if plotFig==1 && ha(jj)==1
            figure(f(1))
            time = (-preEv-1)/fr:1/fr:5.35; time = time(1:end-1);
            subplot(3,5,ii)
            for kk=1:size(psth,2)
                plot(time,smooth(psth(:,kk))','-','LineWidth',2,'Color',cs(kk,:));
                hold on
            end
            hold on
            for kk=1:size(psth,2)
                patch([time fliplr(time)],[smooth(psth(:,kk)+sempsth(:,kk))' flipud(smooth(psth(:,kk)-sempsth(:,kk)))'],cs(kk,:),'FaceAlpha',0.2,'edgealpha',0)
            end
            plot([0 0],[min(psth(:)-sempsth(:))-0.1 max(psth(:)+sempsth(:))+0.1],'r--','LineWidth',1);
            plot([stimInfo.tDur/1000 stimInfo.tDur/1000],[min(psth(:)-sempsth(:))-0.1 max(psth(:)+sempsth(:))+0.1],'r--','LineWidth',1);
            if ii==1; legend(leg); end
            axis tight
            ylim([-2 5])
            xlabel('time(s)')
            ylabel('mean \DeltaF/F_0')
            title([num2str(uT(ii)) ' Hz'])
           
        end
        
        auc(ii,:) = trapz(psth(preEv+1:preEv+intWindow,:));
        
    end
    if plotFig==1 && ha(jj)==1
%     SetFigure(45,30)
%      print(['E:\dataAnalysed\K048\figures\'  mouse '_cell' sprintf('%02d',cellMatching.index(zz,1)) '_respByFreq_03'],'-dpng')
    end
    if plotFig==1 && ha(jj)==1
        figure(f(2))
        colormap('gray')
%         filt = [1,1,1;1,1,1;1,1,1];
%         filt = fspecial('gaussian', [3 3], 1);
%         auc = conv2(auc,filt,'same');
        imagesc(auc')
        caxis([-30 60]);
        set(gca,'XTick',1:length(uT),'XTickLabel',uT,'YTick',1:length(attns),'YTickLabel',attns)
        xtickangle(45)
        xlabel('Frequency (Hz)')
        ylabel('Intensity (dB SPL)')
%         title([num2str(cellMatching.index(zz,1)) '=' num2str(jj)])
%         SetFigure(15,12)
%          print(['E:\dataAnalysed\K048\figures\'  mouse '_cell' sprintf('%02d',cellMatching.index(zz,1)) '_FRA_03'],'-dpng')
         pause()
         clf(f(1))
         clf(f(2))
    end
end


%% plot figs
f1 = figure;
im = dat.mimg(:,:,2);
im=brighten(im,1);
colormap gray
imagesc(im);
uStim = unique(order(:,1));
% ci = floor(length(co)/length(uT));
cs = get(groot,'defaultAxesColorOrder');
cs=[cs;cs;cs;cs];
for ii =  1:length(n)
    img = zeros(1,size(im,1)*size(im,2));
    img(dat.stat(n(ii)).ipix)=1;
    img = reshape(img,size(im,1),size(im,2));
    c = regionprops(img,'Centroid');
    cent(ii,:) = round(c(end).Centroid);
    celln = num2str(ii);
    t = text(c.Centroid(1)-5,c.Centroid(2)+7,celln); t.Color = [1 1 0];
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



