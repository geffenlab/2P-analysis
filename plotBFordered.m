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

%% sort

a = squeeze(sum(m_FRA,1));
a = sum(isnan(a),2)==0 & sig==1;
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

% csplus = 15000;
% csminus = 11400;
% nnplus = knnsearch(uT,csplus);
% nnminus = knnsearch(uT,csminus);
% 
% for ii=1:size(m_FRA,3)
%     [~,ind] = sort(m_FRA(nnplus,:,ii));
%     m_FRA(:,:,ii) = m_FRA(:,ind,ii);
% end

%% plot calcium
for ii=1:length(cmSessions)+1
    t{ii} = sprintf('session%02d',ii);
end
cs = colormap('cool'); close all;
cs=cs(1:floor(length(cs)/size(m_FRA,3)):floor(length(cs)/size(m_FRA,3))*size(m_FRA,3),:);
n=[];
f1 = figure('position',[942, 373, 1176, 965]); f2 = figure; f3 = figure;
for ii=1:size(m_FRA,3)
    figure(f1)
    n_FRA = m_FRA(:,:,ii)';
    %     minFRA = min(n_FRA,[],2);
    %     n_FRA = n_FRA-minFRA;
    maxFRA = max(n_FRA,[],2);
    n_FRA = n_FRA./maxFRA;
    subplot(1,size(m_FRA,3),ii)
    %     imagesc(m_FRA(:,:,ii)')
%     for jj=1:size(n_FRA,1)
%         temp = smooth([n_FRA(jj,1),n_FRA(jj,1), n_FRA(jj,:), n_FRA(jj,end), n_FRA(jj,end)],3);
%         n_FRA(jj,:) = temp(3:end-2);
%     end
    imagesc(n_FRA)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    ylabel('Cell number')
    colormap gray
    hold on
    title(t{ii})
    figure(f2)
    %     subplot(1,size(m_FRA,3),ii)
    plot(uT,mean(n_FRA,1),'color',cs(ii,:))
    axis tight
    %     set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    xlabel('Frequency (Hz)')
    legend(t)
    hold on
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

%% plot frequency tuning of responsive neurons over days

% t = {'pre-FC','24hr post-FC','3 days post-FC','24hr post 2nd FC','3 days post 2nd FC','5 days post 2nd FC'};
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









