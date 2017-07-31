% Work out freezing vs changes in frequency tuning

clear
sortAll = false;
mouse = 'K056';

%% Work out changes of index for calcium data
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\fearCellMatching\*cellMatch*']);
% cmSessions=cmSessions([6,8,9,10]);
load([dataLoc mouse '\fearCellMatching\' cmSessions(1).name])
load([dataLoc mouse '\' cell2mat(cellMatching.session{1}(2)) mouse '_tifStacks\' cell2mat(cellMatching.session{1}(3)) '\sigResp.mat']);
sigCell = sig;
method = 'peak'; % do you want to look at the peak Calcium or the area under the curve ('trapz')
attns = 1; % which attenuations do you want to look at?
session=1;
[s_FRA{1},FRA{1},uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,1,session,method);
m_FRA(:,:,1) = squeeze(mean(FRA{1}(attns,:,:),1));
sm_FRA(:,:,1) = squeeze(mean(s_FRA{1}(attns,:,:),1));
BF = zeros(length(FRA{1}),length(cmSessions)+1);
SBF = zeros(length(FRA{1}),length(cmSessions)+1);
for ii=1:length(m_FRA)
    BF(ii,1) = find(max(m_FRA(:,ii))==m_FRA(:,ii));
    sm = find(max(sm_FRA(:,ii))==sm_FRA(:,ii));
    SBF(ii,1) = sm(randperm(length(sm),1));
end
session=2;
for ff = 1:length(cmSessions)
    [s_FRA{ff+1},FRA{ff+1}] = getToFRAs(cmSessions,mouse,dataLoc,ff,session,method);
    m_FRA(:,:,ff+1) = squeeze(mean(FRA{ff+1}(attns,:,:),1));
    sm_FRA(:,:,ff+1) = squeeze(mean(s_FRA{ff+1}(attns,:,:),1));
    for ii=1:length(m_FRA)
        fm = find(max(m_FRA(:,ii,ff+1))==m_FRA(:,ii,ff+1));
        sfm = find(max(sm_FRA(:,ii,ff+1))==sm_FRA(:,ii,ff+1));
        if ~isnan(fm)
            BF(ii,ff+1) = fm;
            SBF(ii,ff+1) = sfm;
        else
            BF(ii,ff+1) = NaN;
            SBF(ii,ff+1) = NaN;
        end
    end
end

nSessionsToTrack = 2; % minimum number of sessions needed pre fear conditioning
a = squeeze(sum(m_FRA,1));
nPreFC = 5; % sessions before fear conditioning
preFRAs = m_FRA(:,(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack),1:nPreFC);
postFRAs = m_FRA(:,(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack),nPreFC+1:end);
sig = sig(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack);
BFs = BF(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack,:);
% preFRAs = preFRAs(:,sig,:);
% postFRAs = postFRAs(:,sig,:);

csplus = 15000; csminus = 11400;
for ii=1:size(preFRAs,2)
    pre = nanmean(preFRAs(:,ii,:),3);
    pre = SmoothGaus(pre,2)';
    Rpre_plus(ii,1) = interp1(uT,pre,csplus);
    Rpre_minus(ii,1) = interp1(uT,pre,csminus);
    Rpre_higher(ii,1) = nanmean(pre(uT>csplus));
    Rpre_lower(ii,1) = nanmean(pre(uT<csminus));
    for jj=1:size(postFRAs,3)
        post = SmoothGaus(postFRAs(:,ii,jj),2)';
        Rpost_plus(ii,jj) = interp1(uT,post,csplus);
        Rpost_minus(ii,jj) = interp1(uT,post,csminus);
        Rpost_higher(ii,jj) = nanmean(post(uT>csplus));
        Rpost_lower(ii,jj) = nanmean(post(uT<csminus));
    end
end



for ii=1:size(Rpost_plus,2)
    ioc_plus(:,ii) = (Rpost_plus(:,ii)-Rpre_plus)./(Rpost_plus(:,ii)+Rpre_plus);
    ioc_minus(:,ii) = (Rpost_minus(:,ii)-Rpre_minus)./(Rpost_minus(:,ii)+Rpre_minus);
    ioc_higher(:,ii) = (Rpost_higher(:,ii)-Rpre_higher)./(Rpost_higher(:,ii)+Rpre_higher);
    ioc_lower(:,ii) = (Rpost_lower(:,ii)-Rpre_lower)./(Rpost_lower(:,ii)+Rpre_lower);
end


dioc = ioc_plus-ioc_minus;


%% Work out changes of index for behavioural data
freqs = [5 10 11.4 15 20]*10^3;
fdata = [20170427	20	4	20.3	9.4	64.4	20.1
20170428	20	0	30.4	3.3	33.9	29
20170429	20	0.9	0	0	19.4	14.7
20170501	20	5.1	7.1	19.1	29.9	1.5
20170502	20	10	58.8	21.9	39.5	4];


bioc = (fdata(:,freqs==csplus)-fdata(:,freqs==csminus))./(fdata(:,freqs==csplus)+fdata(:,freqs==csminus));
% bioc = (fdata(:,freqs==csplus)-fdata(:,freqs==csminus));




