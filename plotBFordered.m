%% Plot ordered best frequency of sessions


clear
mouse = 'K048';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
session=1;
[~,FRA{1},uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,1,session);
m_FRA(:,:,1) = squeeze(mean(FRA{1}(:,:,:),1));
BF = zeros(1,length(FRA{1}));
for ii=1:length(m_FRA)
    BF(ii) = find(max(m_FRA(:,ii))==m_FRA(:,ii));
end
session=2;
for ff = 1:length(cmSessions)
    [~,FRA{ff+1}] = getToFRAs(cmSessions,mouse,dataLoc,ff,session);
    m_FRA(:,:,ff+1) = squeeze(mean(FRA{ff+1}(:,:,:),1));
end

[~,ind]=sort(BF);
m_FRA = m_FRA(:,ind,:);

a = squeeze(sum(m_FRA,1));
a = sum(isnan(a),2)==0;

t = {'pre-FC','24hr post-FC','3 days post-FC','24hr post 2nd FC','3 days post 2nd FC','5 days post 2nd FC'};
figure
for ii=1:size(m_FRA,3)
    n_FRA = m_FRA(:,a,ii)';
%     n_FRA = (n_FRA-min(n_FRA,[],2));
%     n_FRA = n_FRA./max(n_FRA,[],2);
    subplot(1,size(m_FRA,3),ii)
%         imagesc(m_FRA(:,a,ii)');
%     for jj=1:length(n_FRA)
%         ni(jj) = find(max(n_FRA(jj,:))==n_FRA(jj,:));
%     end
%     [~,ind]=sort(ni);
%     n_FRA = n_FRA(ind,:);
    
%      filt = fspecial('gaussian',3,0.75);
%      n_FRA = conv2(n_FRA,filt(1,:),'same');
    imagesc(n_FRA)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
    title(t{ii})
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



















