% FRA tone analysis

% load processed file
clear
mouse = 'K070';
date = '20170728';
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
% st = dat.sp{1}(n,:); % spike times
% sa = {dat.stat(n).c}; % spike amplitudes (relative)

% load events
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*events.mat']);
if isempty(files)
    files = dir([stimLoc '*exptInfo.mat']);
end
events = load([stimLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr);

% load stim info
stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*stimInfo.mat']);
if ~isempty(files)
    load([stimLoc files(1).name])
else
        stimInfo = events.exptInfo.stimInfo;
end

stimDur = stimInfo{1}.tDur/1000;
ITI = stimInfo{1}.ITI/1000;

% Make raster
preEv = floor(1*fr); % 1 second before stim
postEv = ceil((stimDur+ITI)*fr); % x seconds post stim
eventsOn = eventsOn(eventsOn+postEv<length(traces));
% filter the traces?
% [bf,af] = butter(2,[0.5/fr/2 10/fr/2]); % 10 is good for low pass
[bf,af] = butter(2,20/fr/2,'low'); % 10 is good for low pass
% [bf,af] = butter(2,0.5/fr/2,'high');
for ii=1:size(npilSubTraces,1)
    npst_f(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
    npst_z(ii,:) = zscore(npilSubTraces(ii,:));
end

raster = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,1);
rawRast = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,0);
preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1)-1);

noResp = squeeze(sum(sum(raster)));
raster(:,:,isnan(noResp))=[];
rawRast(:,:,isnan(noResp))=[];
preRast(isnan(noResp),:)=[];
npilSubTraces(isnan(noResp),:)=[];
n(isnan(noResp)) = [];


%% Make spike 'traces' and raster
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

%% order it by frequency and attenuation
p=[];
if length(stimInfo)==1
    p = stimInfo{1}.index(stimInfo{1}.order,:);
else
    for ii=1:length(stimInfo)
        p=[p ; stimInfo{ii}.index(stimInfo{ii}.order,:)];
    end
end
% p = index(toneOrder);
if length(p)>size(raster,1)
    p = p(1:size(raster,1),:);
end
order = repmat(p',1,floor(length(eventsOn)/length(p)))';
[x,ind]=sortrows(order,[1 2]);
raster = raster(ind,:,:);
rawRast = rawRast(ind,:,:);
% spikeRaster = spikeRaster(ind,:,:);


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
            mPost = mean(mean(rawRast(rows,preEv:preEv*2,jj),2));
            c = (mPost-mean(mPre))/stdPre;
            ha(tt,aa,jj)=c>4; % 3 for <10 filter
        end
    end
    
end

wind = round(1*fr); % 2 seconds

for ii=1:length(n)
    
    [p(ii),h(ii)]=signrank(mean(rawRast(:,preEv:preEv+wind,ii),2),mean(rawRast(:,1:preEv-1,ii),2));
        
end

% ha = squeeze(sum(sum(ha),2))>0;
ha = h;
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



%%

cells = find(ha==1);
figure('position',[ 1000         175         637        1163])
time = 1/fr/fr:1/fr:length(npilSubTraces)/fr;
for ii=1:length(cells)
    norm = (npilSubTraces(cells(ii),:)-min(npilSubTraces(cells(ii),:)))/max((npilSubTraces(cells(ii),:)-min(npilSubTraces(cells(ii),:))));
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


%% Plot spike rasters


cells = find(ha==1);
figure('position',[ 1000         175         637        1163])
time = 1/fr/fr:1/fr:length(npilSubTraces)/fr;
for ii=1:length(cells)
    t_spike = (st{cells(ii)}/fr)';
    n_cell = repmat(ii-1,length(t_spike),1)';
   line([t_spike; t_spike], [n_cell; n_cell+1],'color','k');
   hold on
    
end
% for jj=1:length(eventsOn)
%          plot([time(eventsOn(jj)) time(eventsOn(jj))], [0 (sum(ha)+1)*1.5],'r-')
% end
axis tight
 xlim([0 360])
% ylim([0   68.9102])
set(gca,'YTick',1:sum(ha),'YTickLabels',1:sum(ha))
ylabel('cell number')
   


%% Work out FRAs
window = [preEv+1,preEv+round(fr*1)];
FRA = makeCaFRA(raster,window,x(:,1),x(:,2),'trapz');
% spike_FRA = makeSpikeFRA(spikeRaster,window,x(:,1),x(:,2));
filt = fspecial('gaussian',3,0.75);
index = find(ha==1);
figure; set(gcf,'position',[ 510         509        1087         562]);
attnInd = 70-uA;
for ff = 1:length(index)
    ii = index(ff); 
    disp(index(ff)); 
%     subplot(1,2,1); 
%     imagesc(conv2(FRA(:,:,ii),filt,'same'));
    imagesc(FRA(:,:,ii));
    ylabel('Intensity (dB)')
    set(gca,'YTickLabel',attnInd,'FontSize',16)
    xlabel('Frequency (Hz)')
    set(gca,'XTick',1:length(uT),'XTickLabel',uT,'XTickLabelRotation',45)
    title('Peak Calcium response')
    hold on; 
%     subplot(122); 
%     imagesc(conv2(spike_FRA(:,:,ii),filt,'same'));
%     imagesc(spike_FRA(:,:,ii));
%      ylabel('Intensity (dB)')
%     set(gca,'YTickLabel',attnInd,'FontSize',16)
%     xlabel('Frequency (Hz)')
%     set(gca,'XTick',1:length(uT),'XTickLabel',uT,'XTickLabelRotation',45)
%     title('Deconvolved spikes')
%     colormap gray; 
    pause();
    clf;
end

%% 1) Plot mean by frequency across attenuations
figure
set(gcf, 'position',[385 448 1774 890])
uT = unique(order(:,1));
leg = num2str(uT);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
plotFig =1;


index = find(ha==1);
% index = 1:length(ha);
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
        plot([stimDur stimDur],[min(a(:)-b(:))-0.1 max(a(:)+b(:))+0.1],'r--','LineWidth',1);
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

%% Plot frequency tuning by attenuation
figure
set(gcf, 'position',[385 448 1774 890])
uT = unique(order(:,1));
uA = unique(order(:,2));
leg = num2str(uT);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
plotFig =1;


index = find(ha==1);
% index = 1:length(ha);
for ff = 1:length(index)%1:size(raster,3)%
    jj = index(ff);
    
    psth = zeros(length(uT),length(uA),size(raster,2));
    sempsth = psth;
    for ii=1:length(uT)
        for kk = 1:length(uA)
            rows = find(x(:,1)==uT(ii) & x(:,2)==uA(kk));
            psth(ii,kk,:) = (trapz(raster(rows,:,jj)));
            sempsth(ii,kk,:) = std(raster(rows,:,jj))./sqrt(length(rows));
        end
    end
    
    plot(uT,max(squeeze(psth(:,1,:))'))
    hold on
    plot(uT,max(squeeze(psth(:,2,:))'))
    plot(uT,max(squeeze(psth(:,3,:))'))
    
    pause()
    
    clf
    
end

%%
figure
ind = find(ha==1);
for ff=1:sum(ha)
    jj = ind(ff);
    plot(mean(raster(:,:,jj)))
    pause()
    clf
end

%% 2) Plot mean response for different attenuations at each frequency
intWindow = floor(1*fr);
uT = unique(order(:,1));
uA = unique(order(:,2));
attns = [90,80,70,60,50];
leg = num2str(attns');
cs = colormap('summer');  close all;
cs=cs(1:floor(length(cs)/length(uA)):floor(length(cs)/length(uA))*length(uA),:);
plotFig = 1;
psth=[]; sempsth=[];
f(1) = figure; f(2) = figure;
h_attn=[];
ind = find(ha==1);
for ff=1:sum(ha)
    jj = ind(ff);
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
            time = (-preEv-1)/fr:1/fr:(postEv+1)/fr; time = time(1:end-1);
            subplot(3,5,ii)
            for kk=1:size(psth,2)
                plot(time,smooth(psth(:,kk))','-','LineWidth',2,'Color',cs(kk,:));
                hold on
            end
            hold on
            %             for kk=1:size(psth,2)
            %                 patch([time fliplr(time)],[smooth(psth(:,kk)+sempsth(:,kk))' flipud(smooth(psth(:,kk)-sempsth(:,kk)))'],cs(kk,:),'FaceAlpha',0.2,'edgealpha',0)
            %             end
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


%% 1) Plot mean by frequency across attenuations
figure
uT = unique(order(:,1));
leg = num2str(uT);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uT)):floor(length(cs)/length(uT))*length(uT),:);
plotFig = 1;

index = find(ha==1);
for ff = 1:sum(ha)%size(raster,3)%
    jj = index(ff);
    % jj=ff;
    psth = zeros(length(uT),size(raster,2));
    sempsth = psth;
    for ii=1:length(uT) 
        rows = find(x(:,1)==uT(ii));
        psth(ii,:) = (mean(spikeRaster(rows,:,jj)));
        sempsth(ii,:) = std(spikeRaster(rows,:,jj))./sqrt(length(rows));
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
        imagesc(spikeRaster(:,:,jj))
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

%% 1) Plot mean by spike rasters
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
    psth = zeros(length(uT),size(spikeRaster,2));
    sempsth = psth;
    for ii=1:length(uT)
        rows = find(x(:,1)==uT(ii));
        psth(ii,:) = (mean(spikeRaster(rows,:,jj)));
        sempsth(ii,:) = std(spikeRaster(rows,:,jj))./sqrt(length(rows));
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
        surf(time, (1:size(spikeRaster, 1)), flipud(squeeze(spikeRaster(:,:,jj))))
        shading flat; view(2);
        axis tight
        xlabel('Time (s)')
        hold on
%         freezeColors
        plot([0 0], [1 size(spikeRaster,1)],'y-','LineWidth',1)
        nn = round(size(spikeRaster,1)/length(uT));
        csflip = flipud(cs);
        for ii=1:round(size(spikeRaster,1)/nn)
%             line([time(1) time(time==0)],[ii*n ii*n],'Color',csflip(ii,:),'LineWidth',1)
            patch([time(1) time(5) time(5) time(1)], [(ii-1)*nn+1 (ii-1)*nn+1 ii*nn ii*nn],csflip(ii,:),'edgealpha',0,'facealpha',1)
        end
        set(gca,'fontSize',16)
        set(gca,'YTick',[nn/2:nn:size(spikeRaster,1)],'YTickLabel',flipud(round(uT/1000,1)))
        ylabel('Frequency (kHz)')
        
        %  plot(mean(raster(:,:,jj)),'Linewidth',2)
        disp(jj)
        pause()
        clf
    end
    
end

