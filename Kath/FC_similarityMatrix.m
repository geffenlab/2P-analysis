%% Plot ordered best frequency of sessions
clear

sortAll = true;

mouse = 'K048';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
cmSessions(4:5)=[];
session=1;

[raster{1},index{1},fr] = getToRasters(cmSessions,mouse,dataLoc,1,session,1,3);

session=2;
for ff = 1:length(cmSessions)
    [raster{ff+1},index{ff+1}] = getToRasters(cmSessions,mouse,dataLoc,ff,session,1,3);
end


%% Work out similarity matrix
for ff = 1:size(raster{1},3)
    aup = zeros(length(uT),length(raster));
    for ii=1:length(raster)
        
        uT = unique(index{ii}(:,1));
        uA = unique(index{ii}(:,2));
        
        
        psth = zeros(length(uT),size(raster{ii}(:,floor(fr):floor(fr*4),ff),2));
        
        for jj=1:length(uT)
            rows = find(index{ii}(:,1)==uT(jj));
            psth(jj,:) = (mean(raster{ii}(rows,floor(fr):floor(fr*4),ff)));
            aup(jj,ii) = trapz(psth(jj,:));
        end
        
        for jj=1:15
            for kk = 1:15
                m(jj,kk,ii) = dot(psth(jj,:),psth(kk,:));
            end
        end
    end
    
    %     figure
    for jj=1:size(m,3)
        subplot(2,4,jj)
        imagesc(m(:,:,jj))
        colormap gray
        set(gca,'XTick',1:length(uT),'XTickLabel',uT)
        xtickangle(45)
        hold on
        if jj==1
            set(gca,'YTick',1:length(uT),'YTickLabel',uT)
        end
        subplot(2,4,jj+4)
        plot(aup(:,jj).^2)
        set(gca,'XTick',1:length(uT),'XTickLabel',uT)
        xtickangle(45)
        axis tight
    end
    
    pause()
    clf
    
end

%% across cells
a=[];
for ii=1:4
    a(:,ii) = sum(sum(raster{ii}));
end
a = sum(isnan(a),2)==0;

aup = zeros(length(uT),length(raster));
for ii=1:length(raster)
    raster{ii} = raster{ii}(:,:,a);
    uT = unique(index{ii}(:,1));
    uA = unique(index{ii}(:,2));
    
    
    psth = zeros(length(uT),size(raster{ii}(:,floor(fr):floor(fr*4),ff),2));
    
    for jj=1:length(uT)
        rows = find(index{ii}(:,1)==uT(jj));
        psth(jj,:) = nanmean(nanmean(raster{ii}(rows,floor(fr):floor(fr*4),:),1),3);
        aup(jj,ii) = trapz(psth(jj,:));
    end
    
    for jj=1:15
        for kk = 1:15
            m(jj,kk,ii) = dot(psth(jj,:),psth(kk,:));
        end
    end
end

%     figure
for jj=1:size(m,3)
    subplot(2,4,jj)
    imagesc(m(:,:,jj))
    colormap gray
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    hold on
    if jj==1
        set(gca,'YTick',1:length(uT),'YTickLabel',uT)
    end
    subplot(2,4,jj+4)
    plot(aup(:,jj).^2)
    set(gca,'XTick',1:length(uT),'XTickLabel',uT)
    xtickangle(45)
    axis tight
end



% axis equal




