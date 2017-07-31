%% Plot ordered best frequency of sessions
clear
% close all

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
cellIndex = cellMatching.index;
[s_FRA{1},FRA{1},uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,1,session,method);
m_FRA(:,:,1) = squeeze(mean(FRA{1}(attns,:,:),1));
sm_FRA(:,:,1) = squeeze(mean(s_FRA{1}(attns,:,:),1));
BF = zeros(length(FRA{1}),length(cmSessions)+1);
SBF = zeros(length(FRA{1}),length(cmSessions)+1);
for ii=1:length(m_FRA)
    sm = SmoothGaus(m_FRA(:,ii),2)';
    BF(ii,1) = find(max(sm)==sm);
    sm = find(max(sm_FRA(:,ii))==sm_FRA(:,ii));
    SBF(ii,1) = sm(randperm(length(sm),1));
end
session=2;
for ff = 1:length(cmSessions)
    load([dataLoc mouse '\fearCellMatching\' cmSessions(ff).name]);
    load([dataLoc mouse '\' cell2mat(cellMatching.session{2}(2)) mouse '_tifStacks\' cell2mat(cellMatching.session{2}(3)) '\sigResp.mat']);
    cellIndex(:,ff+1) = cellMatching.index(:,2);
    a = find(~isnan(cellMatching.index(:,2)));
    b = cellMatching.index(~isnan(cellMatching.index(:,2)),2);
    sigCell(a,ff+1) = sig(b);
    [s_FRA{ff+1},FRA{ff+1}] = getToFRAs(cmSessions,mouse,dataLoc,ff,session,method);
    m_FRA(:,:,ff+1) = squeeze(mean(FRA{ff+1}(attns,:,:),1));
    sm_FRA(:,:,ff+1) = squeeze(mean(s_FRA{ff+1}(attns,:,:),1));
    for ii=1:length(m_FRA)
        sm = SmoothGaus(m_FRA(:,ii,ff+1),2)';
        fm = find(max(sm)==sm);
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

%% sort

nSessionsToTrack = 1; % minimum number of sessions needed pre fear conditioning
a = squeeze(sum(m_FRA,1));
nPreFC = 5; % sessions before fear conditioning
preFRAs = m_FRA(:,(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack),1:nPreFC);
postFRAs = m_FRA(:,(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack),nPreFC+1:end);
sigCell = sigCell(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack,:);
BFs = BF(nansum(~isnan(a(:,1:nPreFC)),2)>=nSessionsToTrack,:);
s = sum(sigCell,2)>0;
preFRAs = preFRAs(:,s,:);
postFRAs = postFRAs(:,s,:);

type = 'CELLS SIG IN AT LEAST ONE SESSION';

%% plot best frequencies
% for ii=1:length(uT)
%     BFs(BFs==ii)=uT(ii);
% end
% 
% pre = nanmean(BFs(:,1:5),2);
% post = nanmean(BFs(:,6:end),2);
% pre(isnan(post))=[];
% post(isnan(post))=[];
% figure
% plot(pre,post,'.')
% xlabel('Pre-fear conditioning best frequency (Hz)')
% ylabel('Post-fear conditioning best frequency (Hz)')
% hold on
% plot([0 max(uT)],[0 max(uT)],'r--')
% 
% figure
% plot(abs(pre-15000),abs(post-15000),'.')
%%
% 
% for ii=1:size(preFRAs,2)
%     
%     pre = nanmean(preFRAs(:,ii,:),3);
%   
%     pre = SmoothGaus(pre,2)';
%     prestd = nanstd(preFRAs(:,ii,:),[],3);
%     prestd = SmoothGaus(prestd,2)';
%     post = nanmean(postFRAs(:,ii,:),3);
%     post = SmoothGaus(post,2)';
%     poststd = nanstd(postFRAs(:,ii,:),[],3);
%     poststd = SmoothGaus(poststd,2)';
%     
%     plot(uT,pre,'-k','LineWidth',2)
%     hold on
%     plot(uT,post,'-r','LineWidth',2)
%     legend('pre-fear conditioning','post-fear conditioning')
%     patch([uT' flipud(uT)'],[(pre+prestd)' flipud(pre-prestd)'],'k','facealpha',0.5,'edgealpha',0)
%     patch([uT' flipud(uT)'],[(post+poststd)' flipud(post-poststd)'],'r','facealpha',0.5,'edgealpha',0)
%     set(gca,'xscale','log','XTick',uT)
%     axis tight
%     pause()
%     clf
% end

%% Calculate changes in response pre and post fear conditioning (MEAN ACROSS ALL SESSIONS)
csplus = 15000;
csminus = 11400;
f = normpdf(-2:2,0,0.75);
for ii=1:size(preFRAs,2)   
    pre = nanmean(preFRAs(:,ii,:),3);
    pre = SmoothGaus(pre,2)';
    rcsplus(ii,1) = interp1(uT,pre,csplus);
    rcsminus(ii,1) = interp1(uT,pre,csminus);
    lower(ii,1) = mean(pre(uT<csminus));
    higher(ii,1) = mean(pre(uT>csplus));
    prestd = nanstd(preFRAs(:,ii,:),[],3);
    prestd = SmoothGaus(prestd,2)';
    post = nanmean(postFRAs(:,ii,:),3);
    post = SmoothGaus(post,2)';
    rcsplus(ii,2) = interp1(uT,post,csplus);
    rcsminus(ii,2) = interp1(uT,post,csminus);
    lower(ii,2) = mean(post(uT<csminus));
    higher(ii,2) = mean(post(uT>csplus));
    poststd = nanstd(postFRAs(:,ii,:),[],3);
    poststd = SmoothGaus(poststd,2)';
end

%%  plot Changes in cs plus peak 
figure
plot(rcsplus(:,1),rcsplus(:,2),'.')
[h,p]=ttest(rcsplus(:,1),rcsplus(:,2));
hold on
maxf = max(rcsplus(:)); minf = min(rcsplus(:));
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': CS plus response'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('pre-fear conditioning')
ylabel('post-fear conditioning')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'CSplusChanges.fig'])

%%  plot Changes in cs minus peak 
figure
plot(rcsminus(:,1),rcsminus(:,2),'.')
[h,p]=ttest(rcsminus(:,1),rcsminus(:,2))
hold on
maxf = max(rcsminus(:)); minf = min(rcsminus(:));
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': CS minus response'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('pre-fear conditioning')
ylabel('post-fear conditioning')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'CSminusChanges.fig'])

%% Plot change in lower than cs minus
figure
plot(lower(:,1),lower(:,2),'.')
[h,p]=ttest(lower(:,1),lower(:,2))
hold on
maxf = max(lower(:)); minf = min(lower(:));
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': lower than CS- response'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('pre-fear conditioning')
ylabel('post-fear conditioning')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'lowerCSminusChanges.fig'])

%% Plot change in higher than cs minus
figure
plot(higher(:,1),higher(:,2),'.')
[h,p]=ttest(higher(:,1),higher(:,2))
hold on
maxf = max(higher(:)); minf = min(higher(:));
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': higher than CS+ response'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('pre-fear conditioning')
ylabel('post-fear conditioning')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'higherCSplusChanges.fig'])


%% Plot change in cs minus against change in cs plus
figure
plot(rcsminus(:,2)-rcsminus(:,1),rcsplus(:,2)-rcsplus(:,1),'.')
[h,p]=ttest(rcsminus(:,2)-rcsminus(:,1),rcsplus(:,2)-rcsplus(:,1));
hold on
maxf = max([rcsminus(:,2)-rcsminus(:,1);rcsplus(:,2)-rcsplus(:,1)]); 
minf = min([rcsminus(:,2)-rcsminus(:,1);rcsplus(:,2)-rcsplus(:,1)]);
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': change in CS- against change in CS+'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('CS minus')
ylabel('CS plus')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'changeCS-vsCS+.fig'])

%% Plot change in lower than CS- against change higher than CS+
figure
plot(lower(:,2)-lower(:,1),higher(:,2)-higher(:,1),'.')
[h,p]=ttest(lower(:,2)-lower(:,1),higher(:,2)-higher(:,1));
hold on
maxf = max([lower(:,2)-lower(:,1);higher(:,2)-higher(:,1)]); 
minf = min([lower(:,2)-lower(:,1);higher(:,2)-higher(:,1)]);
plot([minf maxf],[minf maxf],'r--')
axis equal
title([type ': change in lower than CS- against change higher than CS+'])
text(minf*0.9, maxf*0.9,sprintf('paired t-test: p = %0.4f',p))
xlabel('lower than CS-')
ylabel('higher than CS+')
savefig(['E:\dataAnalysed\' mouse '\sfnFigures\' type 'changeLowerVsHigher.fig'])

%%


figure
plot(rcsplus(:,1)-rcsminus(:,1),rcsplus(:,2)-rcsminus(:,2),'.')
[h,p]=ttest(rcsplus(:,1)-rcsminus(:,1),rcsplus(:,2)-rcsminus(:,2));
hold on
maxf = max([rcsplus(:,1)-rcsminus(:,1);rcsplus(:,2)-rcsminus(:,2)]); 
minf = min([rcsplus(:,1)-rcsminus(:,1);rcsplus(:,2)-rcsminus(:,2)]);
plot([minf maxf],[0 0],'r--')
axis equal
title(sprintf('Change in difference between CS plus and CS minus response (paired t-test: p = %2d)',p))
xlabel('pre-fear conditioning')
ylabel('post-fear conditioning')

figure
plot(rcsminus(:,2)-rcsminus(:,1),rcsplus(:,2)-rcsplus(:,1),'.')
[h,p]=ttest(rcsminus(:,2)-rcsminus(:,1),rcsplus(:,2)-rcsplus(:,1));
hold on
maxf = max([rcsminus(:,2)-rcsminus(:,1);rcsplus(:,2)-rcsplus(:,1)]); 
minf = min([rcsminus(:,2)-rcsminus(:,1);rcsplus(:,2)-rcsplus(:,1)]);
plot([minf maxf],[minf maxf],'r--')
axis equal
title(sprintf('Change in CS plus vs CS minus response (paired t-test: p = %d)',p))
xlabel('CS minus')
ylabel('CS plus')

%%

figure
plot(nanmean(preFRAs((uT==csplus),:,:),3),nanmean(postFRAs((uT==csplus),:,:),3),'x')
axis equal
maxf = max([nanmean(preFRAs((uT==csplus),:,:),3),nanmean(postFRAs((uT==csplus),:,:),3)]);
xlim([0 maxf])
ylim([0 maxf])
hold on
plot([0 maxf],[0 maxf],'r--')

figure
plot(nanmean(preFRAs((uT==csminus),:,:),3),nanmean(postFRAs((uT==csminus),:,:),3),'x')
axis equal
maxf = max([nanmean(preFRAs((uT==csminus),:,:),3),nanmean(postFRAs((uT==csminus),:,:),3)]);
xlim([0 maxf])
ylim([0 maxf])
hold on
plot([0 maxf],[0 maxf],'r--')

%%
figure
plot(nanmean(preFRAs((uT==csplus),:,:),3)-nanmean(preFRAs((uT==csminus),:,:),3)...
    ,nanmean(postFRAs((uT==csplus),:,:),3)-nanmean(postFRAs((uT==csminus),:,:),3),'x')

axis equal
maxf = max([nanmean(preFRAs((uT==csminus),:,:),3),nanmean(postFRAs((uT==csminus),:,:),3)]);
xlim([0 maxf])
ylim([0 maxf])
hold on
plot([0 maxf],[0 maxf],'r--')
%% plot frequency tuning of responsive neurons over days

cs = colormap('cool'); close all;
cs=cs(1:floor(length(cs)/size(m_FRA,3)):floor(length(cs)/size(m_FRA,3))*size(m_FRA,3),:);
for ii=1:size(m_FRA,2)
    
    for jj=1:size(m_FRA,3)
        
%         plot(m_FRA(:,ii,jj)/max(m_FRA(:,ii,jj)),'LineWidth',2,'Color',cs(jj,:))
          plot(m_FRA(:,ii,jj),'LineWidth',2,'Color',cs(jj,:))
        hold on
        
    end
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    ylabel('Mean peak \DeltaF/F_0')
    xlabel('Frequency (Hz)')
    axis tight
    legend(t)
    pause()
    clf
    
end

figure
for ii=1:size(m_FRA,3)
    plot(mean(m_FRA(:,:,ii),2),'LineWidth',2,'Color',cs(ii,:))
    hold on
end
set(gca,'XTick',1:length(uT),'XTickLabel',uT)
xtickangle(45)
ylabel('Mean peak \DeltaF/F_0')
xlabel('Frequency (Hz)')
axis tight
legend(t)



%% Sort spikes

sm_FRA = sm_FRA(:,a,:);
% sum_smFRA = sum(sum(sm_FRA,1)==0,3);
% sm_FRA=sm_FRA(:,sum_smFRA==0,:);
SBF = SBF(a,:);
% SBF(sum_smFRA==1)=[];

if sortAll
    for ii=1:size(sm_FRA,3)
        [~,ind]=sort(SBF(:,ii));
        sm_FRA(:,:,ii) = sm_FRA(:,ind,ii);
    end
else
    [~,ind]=sort(SBF(:,1));
    sm_FRA(:,:,:) = sm_FRA(:,ind,:);
end

% csplus = 15000;
% csminus = 11400;
% nnplus = knnsearch(uT,csplus);
% nnminus = knnsearch(uT,csminus);
% 
% for ii=1:size(sm_FRA,3)
%     [~,ind] = sort(sm_FRA(nnplus,:,ii));
%     sm_FRA(:,:,ii) = sm_FRA(:,ind,ii);
% end
%%
cs = colormap('hsv');
cs=cs(1:floor(length(cs)/(length(cmSessions)+1)):floor(length(cs)/(length(cmSessions)+1))*(length(cmSessions)+1),:);
for ii=1:length(cmSessions)+1
    t{ii} = sprintf('session%02d',ii);
end
f1 = figure; f2 = figure;
for ii=1:size(sm_FRA,3)
    figure(f1)
    n_FRA = sm_FRA(:,:,ii)';
    n_FRA = n_FRA./max(n_FRA,[],2);
%     for jj=1:size(n_FRA,1)
%         temp = smooth([n_FRA(jj,1),n_FRA(jj,1), n_FRA(jj,:), n_FRA(jj,end), n_FRA(jj,end)],3);
%         n_FRA(jj,:) = temp(3:end-2);
%     end
    subplot(1,size(sm_FRA,3),ii)
    imagesc(n_FRA)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
    title(t{ii})
    figure(f2)
    %     subplot(1,size(sm_FRA,3),ii)
    hold on
    plot(uT,nanmean(n_FRA,1),'Color',cs(ii,:))
    axis tight
    %     set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    legend(t)
    %     figure(f3)
    %     hold on
    %     subplot(1,size(sm_FRA,3),ii)
    ub = unique(SBF(:,ii));
    for jj=1:length(ub)
        n(ii,jj) = sum(SBF(:,ii)==ub(jj));
    end
    %     bar(n)
    
end
figure
bar(n')
set(gca,'XTick',1:length(uT),'XTickLabel',uT)
xtickangle(45)
ylabel('Number of cells')
legend(t)


%%
for ii=1:length(cmSessions)+1
    t{ii} = sprintf('session%02d',ii);
end
cs = colormap('cool'); close all;
cs=cs(1:floor(length(cs)/size(sm_FRA,3)):floor(length(cs)/size(sm_FRA,3))*size(sm_FRA,3),:);
for ii=1:size(sm_FRA,2)
    subplot(8,9,ii)
    for jj=1:size(sm_FRA,3)
        
        plot(sm_FRA(:,ii,jj)/max(sm_FRA(:,ii,jj)),'LineWidth',2,'Color',cs(jj,:))
        hold on
        
    end
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    %         xtickangle(45)
    %         ylabel('Mean peak \DeltaF/F_0')
    %         xlabel('Frequency (Hz)')
    set(gca,'XTickLabel',[],'YTickLabel',[])
    axis tight
    %         legend(t)
    %         pause()
    %         clf
    
end

%%
figure
for ii=1:size(sm_FRA,3)
    plot(mean(sm_FRA(:,:,ii),2),'LineWidth',2,'Color',cs(ii,:))
    hold on
end
set(gca,'XTick',1:length(uT),'XTickLabel',uT)
xtickangle(45)
ylabel('Mean peak \DeltaF/F_0')
xlabel('Frequency (Hz)')
axis tight
legend(t)


%%
t = {'pre-FC','24hr post-FC','3 days post-FC','24hr post 2nd FC','3 days post 2nd FC','5 days post 2nd FC'};
f1 = figure; f2 = figure; f3 = figure;
for ii=2:size(sm_FRA,3)
    figure(f1)
    n_FRA = sm_FRA(:,:,ii)'-sm_FRA(:,:,1)';
    %     n_FRA = n_FRA./max(n_FRA,[],2);
    for jj=1:size(n_FRA,1)
        temp = smooth([n_FRA(jj,1),n_FRA(jj,1), n_FRA(jj,:), n_FRA(jj,end), n_FRA(jj,end)],3);
        n_FRA(jj,:) = temp(3:end-2);
    end
    subplot(1,size(sm_FRA,3)-1,ii-1)
    imagesc(n_FRA)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
    title(t{ii})
    figure(f2)
    subplot(1,size(sm_FRA,3)-1,ii-1)
    plot(uT,mean(n_FRA,1))
    axis tight
    %     set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    figure(f3)
    subplot(1,size(sm_FRA,3)-1,ii-1)
    ub = unique(SBF(:,ii));
    for jj=1:length(ub)
        n(jj) = sum(SBF(:,ii)==ub(jj));
    end
    bar(n)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    ylabel('Number of cells')
end









