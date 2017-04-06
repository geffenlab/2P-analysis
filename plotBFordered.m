%% Plot ordered best frequency of sessions
clear
close all

sortAll = true;

mouse = 'K048';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
cmSessions(4:5)=[];
method = 'peak'; % do you want to look at the peak Calcium or the area under the curve ('trapz')
attns = 1:6; % which attenuations do you want to look at?
session=1;
[s_FRA{1},FRA{1},uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,1,session,method);
m_FRA(:,:,1) = squeeze(mean(FRA{1}(attns,:,:),1));
sm_FRA(:,:,1) = squeeze(mean(s_FRA{1}(attns,:,:),1));
BF = zeros(length(FRA{1}),length(cmSessions)+1);
SBF = zeros(length(FRA{1}),length(cmSessions)+1);
for ii=1:length(m_FRA)
    BF(ii,1) = find(max(m_FRA(:,ii))==m_FRA(:,ii));
    SBF(ii,1) = find(max(sm_FRA(:,ii))==sm_FRA(:,ii));
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

%% sort

a = squeeze(sum(m_FRA,1));
a = sum(isnan(a),2)==0;
% a = 1:size(m_FRA,2);
m_FRA = m_FRA(:,a,:);
BF = BF(a,:);

if sortAll
    for ii=1:size(m_FRA,3)
        [~,ind]=sort(BF(:,ii));
        m_FRA(:,:,ii) = m_FRA(:,ind,ii);
    end
else
    [~,ind]=sort(BF(:,1));
    m_FRA(:,:,:) = m_FRA(:,ind,:);
end

%%


t = {'pre-FC','24hr post-FC','3 days post-FC','24hr post 2nd FC','3 days post 2nd FC','5 days post 2nd FC'};
f1 = figure; f2 = figure; f3 = figure;
for ii=1:size(m_FRA,3)
    figure(f1)
    n_FRA = m_FRA(:,:,ii)';
    minFRA = min(n_FRA,[],2);
    n_FRA = n_FRA-minFRA;
    maxFRA = max(n_FRA,[],2);
    n_FRA = n_FRA./maxFRA;
    subplot(1,size(m_FRA,3),ii)
    imagesc(n_FRA)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
    title(t{ii})
    figure(f2)
    subplot(1,size(m_FRA,3),ii)
    plot(uT,mean(n_FRA,1))
    axis tight
%     set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    figure(f3)
    subplot(1,size(m_FRA,3),ii)
    ub = unique(BF(:,ii));
    for jj=1:length(ub)
       n(jj) = sum(BF(:,ii)==ub(jj));
    end
    bar(n)
      set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    ylabel('Number of cells')
end


% figure
% plot(mean(m_FRA(:,:,1),2))
% hold on
% plot(mean(m_FRA(:,:,2),2))
% plot(mean(m_FRA(:,:,2),2)-mean(m_FRA(:,:,1),2),'k')
% set(gca,'XTick',1:length(uT),'XTickLabel',uT)
% xtickangle(45)
% xlabel('Frequency (Hz)')
% legend('pre','post')
% % ylabel('Cell number')



%% Sort spikes

sm_FRA = sm_FRA(:,a,:);
SBF = SBF(a,:);

if sortAll
    for ii=1:size(m_FRA,3)
        [~,ind]=sort(SBF(:,ii));
        sm_FRA(:,:,ii) = sm_FRA(:,ind,ii);
    end
else
    [~,ind]=sort(SBF(:,1));
    sm_FRA(:,:,:) = sm_FRA(:,ind,:);
end
%%
t = {'pre-FC','24hr post-FC','3 days post-FC','24hr post 2nd FC','3 days post 2nd FC','5 days post 2nd FC'};
f1 = figure; f2 = figure; f3 = figure;
for ii=1:size(sm_FRA,3)
    figure(f1)
    n_FRA = sm_FRA(:,:,ii)';
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
    subplot(1,size(sm_FRA,3),ii)
    plot(uT,mean(n_FRA,1))
    axis tight
%     set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    figure(f3)
    subplot(1,size(sm_FRA,3),ii)
    ub = unique(SBF(:,ii));
    for jj=1:length(ub)
       n(jj) = sum(SBF(:,ii)==ub(jj));
    end
    bar(n)
      set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    ylabel('Number of cells')
end













