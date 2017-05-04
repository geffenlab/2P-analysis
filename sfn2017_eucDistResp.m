%% Plot ordered best frequency of sessions
clear
% close all

sortAll = false;

mouse = 'K056';
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

%% keep cells without NaNs
% a = squeeze(sum(m_FRA,1));
% fc = 6; % sessions after fear conditioning
% preFRAs = m_FRA(:,(nansum(~isnan(a(:,1:fc-1)),2)>2),1:fc-1);
% postFRAs = m_FRA(:,(nansum(~isnan(a(:,1:fc-1)),2)>2),fc:end);
% sig = sig(nansum(~isnan(a(:,1:fc-1)),2)>3);
% BFs = BF(nansum(~isnan(a(:,1:fc-1)),2)>3,:);


%% Euclidean distance

csplus = 15000;
csminus = 11400;

csplus_resp = zeros(size(m_FRA,2),1); csminus_resp = csplus_resp;
for ii=1:size(m_FRA,3)
    csplus_resp(:,ii) = interp1(uT,m_FRA(:,:,ii),csplus);
    csminus_resp(:,ii) = interp1(uT,m_FRA(:,:,ii),csminus);
end

cspr_mean = nanmean(csplus_resp(:,1:5),2);
csmr_mean = nanmean(csminus_resp(:,1:5),2);
for ii=1:size(m_FRA,1)
    for jj=1:size(m_FRA,3)
        cspr = csplus_resp(~isnan(m_FRA(ii,:,jj)'),jj);
        csmr = csminus_resp(~isnan(m_FRA(ii,:,jj)'),jj);
        mf = m_FRA(ii,~isnan(m_FRA(ii,:,jj)'),jj)';
%         ed_csp(ii,jj) = sum(sqrt((cspr_mean(~isnan(m_FRA(ii,:,jj)'))-mf).^2));
%         ed_csm(ii,jj) = sum(sqrt((csmr_mean(~isnan(m_FRA(ii,:,jj)'))-mf).^2));
        ed_csp(ii,jj) = sum(sqrt((cspr-mf).^2));
        ed_csm(ii,jj) = sum(sqrt((csmr-mf).^2));
    end
end

% plot eucDist
figure
cs = colormap('cool'); close all;
cs=cs(1:floor(length(cs)/size(m_FRA,3)):floor(length(cs)/size(m_FRA,3))*size(m_FRA,3),:);
for ii=1:length(cmSessions)+1
    t{ii} = sprintf('session%02d',ii);
end
for ii=1:size(m_FRA,3)
    plot(uT,ed_csp(:,ii),'-','Color',cs(ii,:));
    hold on
end
set(gca,'xscale','log','XTick',uT)
axis tight
title('CS plus')
legend(t)

figure

for ii=1:size(m_FRA,3)
    plot(uT,ed_csm(:,ii),'-','Color',cs(ii,:));
    hold on
end
set(gca,'xscale','log','XTick',uT)
axis tight
title('CS minus')
legend(t)












