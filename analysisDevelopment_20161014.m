% load processed file
clear
load('E:\dataAnalysed\K019\20161027K019_tifStacks\3\F_K019_20161027K019_tifStacks_plane1_Nk25_proc.mat')

%%
n = find(dat.cl.iscell==1);
traces = dat.F.Fcell{1}(n,:);
npilTraces = dat.F.FcellNeu{1}(n,:);
npilCoef =[dat.cl.dcell{n-length(dat.stat0)}]; npilCoef = [npilCoef.B]; 
baseline = npilCoef(2,:); npilCoef = npilCoef(3,:);
cat = traces-baseline'-(npilTraces.*npilCoef');

F0 = mean(cat(:,1:300),2);

dff = (cat-F0)./F0*100;

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
        patch([1:length(traces) length(traces):-1:1],[t fliplr(t)],'b','FaceAlpha',1,'edgeAlpha',1)
%         subplot(2,4,jj+4)
        hold on
        spikeTimes = spikeInfo{1}{ii+(jj-1)}.st;
        s = spikeTimes<7777;
        spikeAmps = spikeInfo{1}{ii+(jj-1)}.c;
        y = zeros(1,length(dff)); y(spikeTimes) = spikeAmps*100;
        h = bar(1:length(traces),y(1:length(traces)));
        h.FaceAlpha = 0.5; h.FaceColor = 'm';
        h.LineWidth = 1; h.BarWidth = 30;
%         xlim([0 7776])
    end
    pause()
    clf
end

