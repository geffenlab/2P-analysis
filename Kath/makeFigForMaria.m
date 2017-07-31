% FRA tone analysis

% load processed file
clear
mouse = 'K056';
date = '20170424';
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
    npst_z(ii,:) = zscore(npilSubTraces(ii,:));
end

raster = makeCaRaster(npst_f,eventsOn,preEv,postEv,1);
rawRast = makeCaRaster(npst_f,eventsOn,preEv,postEv,0);
preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1)-1);

% Make spike 'traces' and raster
spikeTraces = zeros(size(npilSubTraces,1),size(npilSubTraces,2)*2);
for ii=1:size(spikeTraces,1)
%     spikeTraces(ii,st{ii})=sa{ii};
    kernel = dat.stat(ii).kernel';
    for jj=1:length(st{ii})
        spikeTraces(ii,st{ii}(jj):st{ii}(jj)+length(kernel)-1) = spikeTraces(ii,st{ii}(jj):st{ii}(jj)+length(kernel)-1) + kernel*sa{ii}(jj);
    end
end
spikeTraces = spikeTraces(:,1:length(npilSubTraces));
spikeRaster = makeCaRaster(spikeTraces,eventsOn,preEv,postEv,0);

% order it by frequency and attenuation
p = stimInfo.index(stimInfo.order,:);
% p = index(toneOrder);
order = repmat(p',1,floor(length(eventsOn)/length(p)))';
[x,ind]=sortrows(order,[1 2]);
raster = raster(ind,:,:);
rawRast = rawRast(ind,:,:);
spikeRaster = spikeRaster(ind,:,:);


% Work out which cells respond to sounds
uT = unique(order(:,1));
if size(order,2)>1
    uA = unique(order(:,2));
else
    uA=1;
end

ha = zeros(length(uT),length(uA),length(n));
for jj=1:length(n)
    for tt = 1:length(uT)
        for aa=1:length(uA)
            rows = x(:,1)==uT(tt) & x(:,2)==uA(aa);
            mPre = mean(rawRast(rows,1:preEv-1,jj),1);
            stdPre = std(mPre);
%             mPost = max(mean(rawRast(rows,preEv:preEv*3,jj),1));
            mPost = mean(mean(rawRast(rows,preEv:preEv*3,jj),2));
            c = (mPost-mean(mPre))/stdPre;
            ha(tt,aa,jj)=c>4; % 3 for <10 filter
        end
    end
    
end

ha = squeeze(sum(sum(ha),2))>0;
disp([num2str(sum(ha)) ' significantly responding cells out of ' num2str(size(npilSubTraces,1))])

if saveSig==1
    sig = ha;
    save(['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\sigResp.mat'],'sig');
end

BF = zeros(1,size(raster,3));
for ff = 1:size(raster,3)
    jj=(ff);
    psth = zeros(length(uT),size(raster,2));
    sempsth = psth;
    for ii=1:length(uT)
        rows = find(x(:,1)==uT(ii));
        psth(ii,:) = (mean(raster(rows,:,jj)));
        sempsth(ii,:) = std(raster(rows,:,jj))./sqrt(length(rows));
    end
    
    BF(ff) = uT(find(max(max(psth,[],2))==max(psth,[],2)));
end

%% find most active cells

nSp = cellfun(@sum,st)';
[nSp_sorted,ind] = sortrows(nSp,-1);
n2plot = 30;
dur = 200; % seconds



%%
npst_sorted =npst_f(ind,:);% npilSubTraces;%%npst_f(ind,:);
ha_sorted = ha;%ha(ind);
cells = find(ha_sorted==1);
figure('position',[ 427   812   278   413])
time = 1/fr/fr:1/fr:length(npst_sorted)/fr;
for ii=12:12+n2plot
    norm = (npst_sorted(cells(ii),:)-min(npst_sorted(cells(ii),:)))/max((npst_sorted(cells(ii),:)-min(npst_sorted(cells(ii),:))));
    plot(time,norm+(1*(ii-1)),'k','LineWidth',1)
    hold on
    xlabel('Time(s)')
    hold on
    
    
end
% for jj=1:length(eventsOn)
%          plot([time(eventsOn(jj)) time(eventsOn(jj))], [0 (sum(ha)+1)*1.5],'r-')
% end
axis tight
z = 200;
 xlim([z z+dur])
% ylim([0   68.9102])
set(gca,'YTick',1:10:sum(ha))
set(gca,'YTickLabels',-10:10:sum(ha))
ylabel('cell number')


%% Plot spike rasters


cells = find(ha==1);
figure('position',[ 1000         175         637        1163])
time = 1/fr/fr:1/fr:length(npilSubTraces)/fr;
c={'r','k'};
for ii=12:12+n2plot
    t_spike = (st{cells(ii)}/fr)';
    n_cell = repmat(ii-1,length(t_spike),1)';
   line([t_spike; t_spike], [n_cell; n_cell+1],'color',c{rem(ii,2)+1});
   hold on 
end
% for jj=1:length(eventsOn)
%          plot([time(eventsOn(jj)) time(eventsOn(jj))], [0 (sum(ha)+1)*1.5],'r-')
% end
axis tight
z = 200;
 xlim([z z+dur])
% ylim([0   68.9102])
set(gca,'YTick',1:10:sum(ha))
set(gca,'YTickLabels',-10:10:sum(ha))
ylabel('cell number')
xlabel('time (s)')


%% 1) Plot mean by frequency across attenuations
figure
set(gcf, 'position',[385 448 1774 890])
uT = unique(order(:,1));
leg = num2str(uT);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
plotFig =1;


index = find(ha==1);
for ff = 1:length(index)%1:size(raster,3)%
    jj = index(ff);
    % jj=(ff);
    psth = zeros(length(uT),size(raster,2));
    sempsth = psth;
    for ii=1:length(uT)
        rows = find(x(:,1)==uT(ii));
        psth(ii,:) = (mean(raster(rows,:,jj)));
        sempsth(ii,:) = std(raster(rows,:,jj))./sqrt(length(rows));
    end
    
    %     BF(ff) = uT(find(max(max(psth,[],2))==max(psth,[],2)));
    
    if plotFig==1
        
        time = (-preEv-1)/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
            plot(time,psth(ii,:)','-','LineWidth',2,'Color',cs(ii,:));
            %             plot(time,mean(psth(ii,:))','-','LineWidth',2,'Color',cs(ii,:));
            hold on
        end
        set(gca,'fontSize',16)
        ylabel('\DeltaF/F_0')
        xlabel('Time (s)')
        hold on
        %         for ii=1:size(psth,1)
        %              patch([time fliplr(time)],[smooth(psth(ii,:)+sempsth(ii,:))' flipud(smooth(psth(ii,:)-sempsth(ii,:)))'],cs(ii,:),'FaceAlpha',0.2,'edgealpha',0)
        %         end
        a = psth(:,:); b = sempsth(:,:);
        plot([0 0],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        plot([stimInfo.tDur/1000 stimInfo.tDur/1000],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
        legend(leg)
        axis tight
        subplot(1,2,2)
        colormap('gray')
        surf(time, (1:size(raster, 1)), flipud(squeeze(raster(:,:,jj))))
        shading flat; view(2);
        axis tight
        xlabel('Time (s)')
        hold on
%         freezeColors
        plot([0 0], [1 size(raster,1)],'y-','LineWidth',1)
        nn = round(size(raster,1)/length(uT));
        csflip = flipud(cs);
        for ii=1:round(size(raster,1)/nn)
%             line([time(1) time(time==0)],[ii*n ii*n],'Color',csflip(ii,:),'LineWidth',1)
            patch([time(1) time(5) time(5) time(1)], [(ii-1)*nn+1 (ii-1)*nn+1 ii*nn ii*nn],csflip(ii,:),'edgealpha',0,'facealpha',1)
        end
        set(gca,'fontSize',16)
        set(gca,'YTick',[nn/2:nn:size(raster,1)],'YTickLabel',flipud(round(uT/1000,1)))
        ylabel('Frequency (kHz)')
        
        %  plot(mean(raster(:,:,jj)),'Linewidth',2)
        disp(jj)
        pause()
        clf
    end
    
end


    
%% plot figs

cs = colormap('hsv');close(gcf)
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
f1 = figure;
im = dat.mimg(:,:,2);
im=brighten(im,1);
colormap gray
imagesc(im);
uStim = unique(order(:,1));
hold on
for ii =  1:size(raster,3)
    img = zeros(1,size(im,1)*size(im,2));
    img(dat.stat(n(ii)).ipix)=1;
    img = reshape(img,size(im,1),size(im,2));
    c = regionprops(img,'Centroid');
    cent(ii,:) = round(c(end).Centroid);
    celln = num2str(ii);
%     t = text(c.Centroid(1)-3,c.Centroid(2),celln); t.Color = [1 1 0];
    s = [bwboundaries(img)];
    boundaries(ii) = s(1);
%     if ha(ii)==1
%         plot(boundaries{ii}(:,2),boundaries{ii}(:,1),'Color',cs(BF(ii)==uStim,:),'LineWidth',1.5)
%     else
        plot(boundaries{ii}(:,2),boundaries{ii}(:,1),'Color','k','LineWidth',1.5)
%     end
end
img = reshape(img,size(im,1),size(im,2));
axis off
axis equal
data.cellCentroids = cent;
data.cellOutlines = boundaries;
data.image = im;

