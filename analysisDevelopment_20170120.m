% load processed file
clear
mouse = 'K053';
date = '20170223';
exptNo = '1';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

%%
n = find([dat.stat.iscell]==1);
traces = dat.Fcell{1}(n,:);
npilTraces = dat.FcellNeu{1}(n,:);
npilCoeffs = [dat.stat(n).neuropilCoefficient]';
npilSubTraces = traces-(npilTraces.*npilCoeffs);

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

stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
files = dir([stimLoc '*.xml']);
voltageFile = files(~cellfun(@isempty,strfind({files.name},'Voltage'))); 
recFile = files(cellfun(@isempty,strfind({files.name},'Voltage')));
% Import xml file and get (relative) time stamp
% [s]=xml2struct([stimLoc recFile.name]);
% tn = {s.Children.Name};
% la = strcmp(tn,'Sequence');
% tnn = {s.Children(la).Children.Name};
% laa = find(strcmp(tnn,'Frame')==1);
% for jj=1:length(laa)
%     laaa = strcmp({s.Children(la).Children(laa(jj)).Attributes.Name},'relativeTime');
%     ta(jj) = str2num(s.Children(la).Children(laa(jj)).Attributes(laaa).Value);
% end

%%
% dt=diff(ta);
% fr=1/mean(dt); % frame rate
% fr=30.35;
% fr=30.4;
fr = 29.4238;

% Import the voltage trace
soundFile=dir([stimLoc '*.csv']);
sound=importdata([stimLoc soundFile.name],',');

vt.t=sound.data(:,1); % sound time 
vt.s=sound.data(:,2); % voltage trace of sound
vt.fs = 40000; % sample rate
conversionFactor = (vt.fs/fr);
vt.ds = [vt.s(1:conversionFactor:length(vt.s)); zeros(size(traces,2)-length(1:conversionFactor:length(vt.s)),1)]';
% vt.ds(1)=1;
% b=1;
% for ii=2:length(ta)
%     b(ii) = floor(ta(ii)*vt.fs);
%     vt.ds(ii) = vt.s(b(ii)+1);
% end
vt.eventsOn = [find(diff(vt.ds) > 0.3)+1];
vt.eventsOff = find(diff(vt.ds) < -0.3)+1;

%  vt.eventsOn = vt.eventsOn(1:end-1);

% Frames
% fr = sound.data(:,3);
% dt = diff(fr);
% ev = find(dt>0.1);
% evT = diff(ev);
% evR = evT(evT>100);


%% Plot
preEv = floor(1*fr); % 1 second before stim
postEv = ceil(3*fr); % 3 seconds post stim
stimDur = mode(vt.eventsOff(2:end)-vt.eventsOn(1:end)); vt.eventsOff = vt.eventsOn+stimDur;
ind=1; raster = zeros(size(npilTraces,1),preEv+postEv+1,length(vt.eventsOn));
for ii=1:length(vt.eventsOn) % ignore 1st event because will be the start recording event - events from then on are relevant
   raster(:,:,ind) = traces(:,vt.eventsOn(ii)-preEv: vt.eventsOn(ii)+postEv);
   raster(:,:,ind) = (raster(:,:,ind)-mean(squeeze(raster(:,1:preEv-1,ind)),2))./std(squeeze(raster(:,1:preEv-1,ind)),[],2);
   ind = ind+1;
end

raster = permute(raster,[3,2,1]);
[x,ind]=sortrows(order');
raster = raster(ind,:,:);

for ii=1:size(raster,3)
     subplot(1,2,1)
    t = smooth(mean(squeeze(raster(:,:,ii)),1),10);
    ts = smooth(std(squeeze(raster(:,:,ii)),[],1),10);
plot(-preEv/fr:1/fr:(postEv)/fr,t,'LineWidth',2)
time = -preEv/fr:1/fr:(postEv)/fr;
axis tight
hold on
patch([time fliplr(time)],[(t+ts)' flipud(t-ts)'],'b','FaceAlpha',0.3,'edgealpha',0)
plot([0 0],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
plot([1 1],[min(t-ts) max(t+ts)],'r--','LineWidth',1)
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
%% For Tones

stimLoc=['C:\data\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
stimFile = dir([stimLoc '*.mat']);

stimInfo = load([stimLoc stimFile.name]);
stimInfo.f = stimInfo.index(stimInfo.toneOrder,:);
stimInfo.f = repmat(stimInfo.f,4,1);

% stimInfo.f = stimInfo.f(2:96); %raster = raster(2:end,:,:); % for debugging
[z,ind]=sortrows(stimInfo.f,[1,2]);
rasSorted = raster(ind,:,:);

uF = unique(stimInfo.f(:,1));
for ii=1:size(raster,3)
    subplot(1,2,1)
    for jj=1:length(uF)
        psth(jj,:) = smooth(mean(squeeze(rasSorted(stimInfo.f==uF(jj),:,ii))));
    end
plot((0:length(psth)-1)*fr, psth')
hold on
plot([preEv preEv+stimDur]*fr,[0 0],'r-','LineWidth',3)
subplot(1,2,2)
imagesc(rasSorted(:,:,ii))

pause()
clf
end





%%
for jj=1:length(npilTraces)
for ii=1:length(vt.eventsOn)
plot([vt.eventsOn(ii) vt.eventsOff(ii)],[zeros(length(vt.eventsOn),1) zeros(length(vt.eventsOn),1)],'r-')
hold on
end
plot(smooth(npilSubTraces(jj,:),10))
pause()
clf

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
% f1 = figure;
% colormap gray
% % imagesc(scaleK019);
% for ii = 1:length(n)
% img = zeros(1,size(dat.img1.V,1)*size(dat.img1.V,2));
% img(dat.stat(n(ii)).ipix)=1;
% img = reshape(img,size(dat.img1.V,1),size(dat.img1.V,2));
% c = regionprops(img,'Centroid');
% cent(ii,:) = round(c.Centroid);
% celln = num2str(ii);
% % t = text(c(1)-5,c(2)+7,celln); t.Color = [1 1 0];
% stats = [bwboundaries(img)];
% hold on
% % plot(stats{1}(:,2),stats{1}(:,1),'y')
% end
% axis off
% axis equal
% 
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


