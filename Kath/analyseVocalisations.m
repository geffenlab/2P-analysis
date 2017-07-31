%% load processed file
clear
saveSig = 1;
mouse = 'K056';
date = '20170421';
exptNo = '1';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

% get Ca data
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:); % raw fluorescence
badFrames = find(dat.ops.badframes==1);
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times
sa = {dat.stat(n).c}; % spike amplitudes (relative)

% Load the stimulus info
files = dir([dataLoc exptNo '\*exptInfo.mat']);
load([dataLoc exptNo '\' files.name])

% Get stim info
stimInfo = exptInfo.stimInfo{1};
stimDur = stimInfo.IVO*stimInfo.vocrep; % stim duration
ISI = stimInfo.ISI; % inter-stimulus interval

% Get events
fr = 1/(mean(diff(frameOn/fs))); % mean frame rate
eventsOn = floor(eventOn/fs*fr);
eventsOn = eventsOn(1:stimInfo.vocrep:length(eventsOn));

% Make raster
Wn = 14/fr/2;
[b,a] = butter(5,Wn,'low');
preEv = floor(1*fr); % 1 second before stim
postEv = ceil((stimDur+ISI)*fr); % x seconds post stim

% filter the traces?
[bf,af] = butter(2,[0.5/fr/2 25/fr/2]);
for ii=1:size(npilSubTraces,1)
    npst_f(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
    npst_z(ii,:) = zscore(npilSubTraces(ii,:));
end

raster = makeCaRaster(npst_f,eventsOn,preEv,postEv,1);
rawRast = makeCaRaster(npilSubTraces,eventsOn,preEv,postEv,0);
preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1)-1);

% Make spike 'traces' and raster
spikeTraces = zeros(size(npilSubTraces));
for ii=1:size(spikeTraces,1)
    spikeTraces(ii,st{ii})=sa{ii};
end
spikeRaster = makeCaRaster(spikeTraces,eventsOn,preEv,postEv,0);

% order it by frequency and attenuation
p = stimInfo.index(stimInfo.order)';
types = [1,2,3,1,3,2,4,4,1,5,5];
p(:,2) = types(stimInfo.order);
order = repmat(p',1,floor(length(eventsOn)/length(p)))';
%5=complex,6=chevron,8=complex,...
% 14=FMup,25=FMup,34=chevron, 86 = FMdown,101= FMdown,107=complex,
% 71&120=flat

[x,ind]=sortrows(order,[2 1]);
raster = raster(ind,:,:);
rawRast = rawRast(ind,:,:);
spikeRaster = spikeRaster(ind,:,:);


% Work out which cells respond to sounds
uV = unique(order(:,1));

ha = zeros(length(uV),length(n));
for jj=1:length(n)
    for tt = 1:length(uV)
        rows = x(:,1)==uV(tt);
        mPre = mean(rawRast(rows,1:preEv-1,jj),1);
        stdPre = std(mPre);
        mPost = mean(mean(rawRast(rows,preEv:preEv+stimDur+fr,jj),2));
        c = (mPost-mean(mPre))/stdPre;
        ha(tt,jj)=c>2;
    end   
end

 ha = sum(ha)>0;

if saveSig==1
    sig = ha;
    save(['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\sigResp.mat'],'sig');
end

BV = zeros(1,size(raster,3));
for ff = 1:size(raster,3)
    jj=(ff);
    psth = zeros(length(uV),size(raster,2));
    sempsth = psth;
    for ii=1:length(uV)
        rows = find(x(:,1)==uV(ii));
        psth(ii,:) = (mean(raster(rows,:,jj)));
        sempsth(ii,:) = std(raster(rows,:,jj))./sqrt(length(rows));
    end
    
    BV(ff) = uV(find(max(max(psth,[],2))==max(psth,[],2)));
end

disp([num2str(sum(ha)) ' significantly responding cells out of ' num2str(size(npilSubTraces,1))])

%% plot the data
figure
set(gcf, 'position',[385 448 1774 890])
uV = unique(order(:,1));
leg = num2str(uV);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uV)):floor(length(cs)/length(uV))*length(uV),:);
plotFig =1;


% ha = ones(1,size(raster,3));
index = find(ha==1);
for ff = 1:length(index)%1:size(raster,3)%
    jj = index(ff);
    % jj=(ff);
    psth = zeros(length(uV),size(raster,2));
    sempsth = psth;
    for ii=1:length(uV)
        rows = find(x(:,1)==uV(ii));
        psth(ii,:) = (mean(raster(rows,:,jj)));
        sempsth(ii,:) = std(raster(rows,:,jj))./sqrt(length(rows));
    end
    
    %     BF(ff) = uV(find(max(max(psth,[],2))==max(psth,[],2)));
    
    if plotFig==1
        time = (-preEv-1)/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
            plot(time(1:round(preEv*4)),psth(ii,1:round(preEv*4))','-','LineWidth',2,'Color',cs(ii,:));
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
        nn = round(size(raster,1)/length(uV));
        csflip = flipud(cs);
        for ii=1:round(size(raster,1)/nn)
%             line([time(1) time(time==0)],[ii*n ii*n],'Color',csflip(ii,:),'LineWidth',1)
            patch([time(1) time(5) time(5) time(1)], [(ii-1)*nn+1 (ii-1)*nn+1 ii*nn ii*nn],csflip(ii,:),'edgealpha',0,'facealpha',1)
        end
        set(gca,'fontSize',16)
        set(gca,'YTick',[nn/2:nn:size(raster,1)],'YTickLabel',flipud(round(uV/1000,1)))
        ylabel('Frequency (kHz)')
        
        %  plot(mean(raster(:,:,jj)),'Linewidth',2)
        disp(jj)
        pause()
        clf
    end
    
end

%% plot the spikes
figure
set(gcf, 'position',[385 448 1774 890])
uV = unique(order(:,1));
leg = num2str(uV);
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/length(uV)):floor(length(cs)/length(uV))*length(uV),:);
plotFig =1;


ha = ones(1,size(spikeRaster,3));
index = find(ha==1);
for ff = 1:length(index)%1:size(spikeRaster,3)%
    jj = index(ff);
    % jj=(ff);
    psth = zeros(length(uV),size(spikeRaster,2));
    sempsth = psth;
    for ii=1:length(uV)
        rows = find(x(:,1)==uV(ii));
        psth(ii,:) = (mean(spikeRaster(rows,:,jj)));
        sempsth(ii,:) = std(spikeRaster(rows,:,jj))./sqrt(length(rows));
    end
    
    %     BF(ff) = uV(find(max(max(psth,[],2))==max(psth,[],2)));
    
    if plotFig==1
        time = (-preEv-1)/fr:1/fr:(postEv)/fr;
        subplot(1,2,1)
        for ii=1:size(psth,1)
            plot(time(1:round(preEv*4)),psth(ii,1:round(preEv*4))','-','LineWidth',2,'Color',cs(ii,:));
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
        surf(time, (1:size(spikeRaster, 1)), flipud(squeeze(spikeRaster(:,:,jj))))
        shading flat; view(2);
        axis tight
        xlabel('Time (s)')
        hold on
%         freezeColors
        plot([0 0], [1 size(spikeRaster,1)],'y-','LineWidth',1)
        nn = round(size(spikeRaster,1)/length(uV));
        csflip = flipud(cs);
        for ii=1:round(size(spikeRaster,1)/nn)
%             line([time(1) time(time==0)],[ii*n ii*n],'Color',csflip(ii,:),'LineWidth',1)
            patch([time(1) time(5) time(5) time(1)], [(ii-1)*nn+1 (ii-1)*nn+1 ii*nn ii*nn],csflip(ii,:),'edgealpha',0,'facealpha',1)
        end
        set(gca,'fontSize',16)
        set(gca,'YTick',[nn/2:nn:size(spikeRaster,1)],'YTickLabel',flipud(round(uV/1000,1)))
        ylabel('Frequency (kHz)')
        
        %  plot(mean(spikeRaster(:,:,jj)),'Linewidth',2)
        disp(jj)
        pause()
        clf
    end
    
end


