%clustering

clear
mouse = 'K056';
date = '20170426';
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
b = [dat.stat.noiseLevel];


fr = 29.8743; % mean frame rate

% filter the traces?
% [bf,af] = butter(2,[0.5/fr/2 10/fr/2]); % 10 is good for low pass
[bf,af] = butter(2,10/fr/2,'low'); % 10 is good for low pass
% [bf,af] = butter(2,0.5/fr/2,'high');
for ii=1:size(npilSubTraces,1)
    npst_f(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
    npst_g(ii,:) = SmoothGaus(double(npilSubTraces(ii,:)),1);
    npst_z(ii,:) = zscore(npst_g(ii,:));
end

F0 = min(npst_g,[],2);
df_F0 = (npst_g-F0)./abs(F0);


%% deconvolved Ca trace

for ii=1:length(n)
    a = zeros(1,size(npst_z,2));
    a(st{ii}) = sa{ii};
    dc(ii,:) = conv(a,k{ii});
    deca(ii,:) = dc(ii,1:size(npst_z,2))+(npilTraces(ii,:)*npilCoeffs(ii))+b(ii);
end

% Plot the traces and deconvolved calcium signals
% figure('position',[5        1173        2554         165])
% for ii=1:length(n)
%     
% plot(traces(ii,:)/max(traces(ii,:))-mean(traces(ii,:)/max(traces(ii,:))))
% hold on
% plot(dc(ii,:)/max(dc(ii,:)))
% pause()
% clf
% end

for ii=1:length(n)
    for jj=1:length(n)
        ed(ii,jj) = mean(sqrt((dc(ii,:)-dc(jj,:)).^2));
    end
end




