% load processed file
clear
load('E:\dataAnalysed\K019\20161027K019_tifStacks\3\F_K019_20161027K019_tifStacks_plane1_Nk25_proc.mat')

%%
n = find(dat.cl.iscell==1);
traces = dat.F.Fcell{1}(n,:);
npilTraces = dat.F.FcellNeu{1}(n,:);
npilCoef =[dat.cl.dcell{n-length(dat.stat0)}]; npilCoef = [npilCoef.B]; 
baseline = npilCoef(2,:); npilCoef = npilCoef(3,:);
cat = bsxfun(@minus,traces,baseline')-bsxfun(@times,npilTraces,npilCoef');

F0 = mean(cat(:,1:300),2);

dff = bsxfun(@rdivide,bsxfun(@minus,cat,F0),F0*100);

spikeInfo = {dat.cl.dcell(n-length(dat.stat0))};


% ii=5
% pix = zeros(1,500*490);
% ipix = dat.stat(n(ii)).ipix;
% pix(ipix) = ii;
% imagesc(reshape(pix,500,490))
%%
figure
for ii=1:4:length(n)
    for jj=1:4
        subplot(2,2,jj)
%         plot(1:7776,traces(ii+(jj-1),1:7776)-mean(traces(ii+(jj-1),1:7776))+100,'.')
        t = traces(ii+(jj-1),:)-mean(traces(ii+(jj-1),:))+100;
        patch([1:length(traces) length(traces):-1:1],[t fliplr(t)],'b','FaceAlpha',0,'edgeAlpha',1)
%         subplot(2,4,jj+4)
        hold on
        spikeTimes = spikeInfo{1}{ii+(jj-1)}.st;
        spikeAmps = spikeInfo{1}{ii+(jj-1)}.c;
        y = zeros(1,length(dff)); y(spikeTimes) = spikeAmps*100;
        h = bar(1:length(traces),y(1:length(traces)));
        h.FaceAlpha = 0.5; h.FaceColor = 'm'; h.EdgeAlpha = 0;
        h.LineWidth = 1; h.BarWidth = 30;
%         xlim([0 7776])
    end
    pause()
    clf
end

%%
f1 = figure;
colormap gray
imagesc(scaleK019);
for ii = 1:length(n)
img = zeros(1,509*506);
img(dat.stat(n(ii)).ipix)=1;
img = reshape(img,509,506);
c = regionprops(img,'Centroid');
c = round(c.Centroid);
celln = num2str(ii);
% t = text(c(1)-5,c(2)+7,celln); t.Color = [1 1 0];
stats = [bwboundaries(img)];
hold on
plot(stats{1}(:,2),stats{1}(:,1),'y')
end
axis off
axis equal

%%
figure
fr = 30; 
ind = 1;
t1 = 1;
for ii = 1:size(cat,1)
    if mod(ii,10)==0
        savefig(['/Users/kath/Documents/ImagingData/K019AnlaysisForMariaGrant/traces_' sprintf('%d',t1)])
        t1 = t1+1;
        close all
        figure
        ind = 1;
    end
    plot(1/30:1/30:300,(SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10))/max((SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10)))+ind*2,'k')
    hold on
    ind = ind+1;   
    disp(ii)
end

figure
for ii = 1:size(cat,1)
    plot(1/30:1/30:300,(SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10))/max((SmoothGaus(cat(ii,1:fr*300)-mean(cat(ii,1:fr*300)),10)))+ii*2,'k')
    hold on
end
 savefig(['/Users/kath/Documents/ImagingData/K019AnlaysisForMariaGrant/alltraces.fig'])



