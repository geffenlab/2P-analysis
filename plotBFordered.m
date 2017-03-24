%% Plot ordered best frequency of sessions


clear
mouse = 'K048';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
session=1;
[FRA,uT,uA] = getToFRAs(cmSessions,mouse,dataLoc,1,session);
m_FRA(:,:,1) = squeeze(mean(FRA(1:3,:,:),1));
BF = zeros(1,length(FRA));
for ii=1:length(m_FRA)
    BF(ii) = find(max(m_FRA(:,ii))==m_FRA(:,ii));
end
session=2;
for ff = 1:length(cmSessions)   
    [FRA] = getToFRAs(cmSessions,mouse,dataLoc,ff,session);
    m_FRA(:,:,ff+1) = squeeze(mean(FRA(1:3,:,:),1));
end

[~,ind]=sort(BF);
m_FRA = m_FRA(:,ind,:);

a = squeeze(sum(m_FRA,1));
a = sum(isnan(a),2)==0;


figure
for ii=1:size(m_FRA,3)
    subplot(1,size(m_FRA,3),ii)
    imagesc(m_FRA(:,a,ii)');
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
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



















