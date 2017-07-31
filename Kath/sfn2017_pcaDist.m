%% Plot ordered best frequency of sessions
clear
% close all

sortAll = false;

mouse = 'K056';
dataLoc = ['E:\dataAnalysed\' mouse '\network\'];
load([dataLoc 'sfnData.mat'])


%% keep cells without NaNs
% a = squeeze(sum(m_FRA,1));
% fc = 6; % sessions after fear conditioning
% preFRAs = m_FRA(:,(nansum(~isnan(a(:,1:fc-1)),2)>2),1:fc-1);
% postFRAs = m_FRA(:,(nansum(~isnan(a(:,1:fc-1)),2)>2),fc:end);
% sig = sig(nansum(~isnan(a(:,1:fc-1)),2)>3);
% BFs = BF(nansum(~isnan(a(:,1:fc-1)),2)>3,:);
dd=1
uT = unique(data(dd).BF);
index = ~isnan(sum([data.index],2));
for ii=1:length(data)
    ind(:,ii) = data(ii).index(index);
end
csplus = 15000;
csminus = 11400;

for ii = 1:length(data)
    csplus_resp = interp1(uT,data(ii).FRA,csplus);
    csminus_resp = interp1(uT,data(ii).FRA,csminus);
    FRA{ii} = [data(ii).FRA(:,ind(:,ii))' csplus_resp(ind(:,ii))' csminus_resp(ind(:,ii))'];
end

freqs = [uT csplus csminus];

%% PCA analysis projected into first session PCs
dd=1;
fra = [FRA{dd}];
[eigenvec,tData,eigenval,tsquared,explained,mu] = pca(fra');
V = eigenvec(:,1:4);

% project data
p_FRA = cell(1,length(data));
for ii=1:length(data)
    for jj=1:size(FRA{ii},2)
        p_FRA{ii}(:,:,jj) = FRA{ii}(:,jj).*V;
    end
end

csp = find(freqs==csplus);
dist = zeros(length(data),length(freqs));
for ii=1:length(data)
    for jj=1:size(p_FRA{ii},3)
        dist(ii,jj) = sum( sqrt( ( p_FRA{ii}(:,csp)-p_FRA{ii}(:,jj) ).^2 ));
    end
end

% make plots
imagesc(dist)
colorbar
colormap gray
set(gca,'XTickLabels',freqs)
set(gca,'XTick',1:length(freqs),'XTickLabels',freqs)

figure
plot(mean(dist,2))



%% PCA analysis plot into pre-fear conditioning
dd=1:5; 
fra = [FRA{dd}];
[eigenvec,tData,eigenval,tsquared,explained,mu] = pca(fra');
V = eigenvec(:,1:4);

% project data
p_FRA = cell(1,length(data));
for ii=1:length(data)
    for jj=1:size(FRA{ii},2)
        p_FRA{ii}(:,:,jj) = FRA{ii}(:,jj).*V;
    end
end

csp = find(freqs==csplus);
dist = zeros(length(data),length(freqs));
for ii=1:length(data)
    for jj=1:size(p_FRA{ii},3)
        dist(ii,jj) = mean( sqrt( ( p_FRA{ii}(:,csp)-p_FRA{ii}(:,jj) ).^2 ));
    end
end
% make plots
figure
imagesc(dist)
colorbar
colormap gray
set(gca,'XTickLabels',freqs)
set(gca,'XTick',1:length(freqs),'XTickLabels',freqs)

figure
plot(mean(dist,2))


%% PCA analysis plot into ALL
dd=1:length(data); 
fra = [FRA{dd}];
[eigenvec,tData,eigenval,tsquared,explained,mu] = pca(fra');
V = eigenvec(:,1:10);

% project data
p_FRA = cell(1,length(data));
for ii=1:length(data)
    for jj=1:size(FRA{ii},2)
        p_FRA{ii}(:,:,jj) = FRA{ii}(:,jj).*V;
    end
end

csp = find(freqs==csplus);
dist = zeros(length(data),length(freqs));
for ii=1:length(data)
    for jj=1:size(p_FRA{ii},3)
        dist(ii,jj) = mean( sqrt( ( p_FRA{ii}(:,csp)-p_FRA{ii}(:,jj) ).^2 ));
    end
end
% make plots
figure
imagesc(dist)
colorbar
colormap gray
set(gca,'XTickLabels',freqs)
set(gca,'XTick',1:length(freqs),'XTickLabels',freqs)

figure
plot(mean(dist,2))


%% PCA analysis project into post-fear conditioning
dd=6:10; 
fra = [FRA{dd}];
[eigenvec,tData,eigenval,tsquared,explained,mu] = pca(fra');
V = eigenvec(:,1:4);

% project data
p_FRA = cell(1,length(data));
for ii=1:length(data)
    for jj=1:size(FRA{ii},2)
        p_FRA{ii}(:,:,jj) = FRA{ii}(:,jj).*V;
    end
end

csp = find(freqs==csplus);
dist = zeros(length(data),length(freqs));
for ii=1:length(data)
    for jj=1:size(p_FRA{ii},3)
        dist(ii,jj) = mean( sqrt( ( p_FRA{ii}(:,csp)-p_FRA{ii}(:,jj) ).^2 ));
    end
end
% make plots
figure
imagesc(dist)
colorbar
colormap gray
set(gca,'XTickLabels',freqs)
set(gca,'XTick',1:length(freqs),'XTickLabels',freqs)

figure
plot(dist)


%% PCA analysis project into itself
dist = zeros(length(data),length(freqs));
for dd=1:length(data) 
fra = [FRA{dd}];
[eigenvec,tData,eigenval,tsquared,explained,mu] = pca(fra');
V = eigenvec(:,1:4);

% project data
p_FRA = cell(1,length(data));
for ii=1:length(data)
    for jj=1:size(FRA{ii},2)
        p_FRA{ii}(:,:,jj) = FRA{ii}(:,jj).*V;
    end
end

csp = find(freqs==csplus);

for ii=dd
    for jj=1:size(p_FRA{ii},3)
        dist(ii,jj) = mean( sqrt( ( p_FRA{ii}(:,csp)-p_FRA{ii}(:,jj) ).^2 ));
    end
end
end
% make plots
figure
imagesc(dist)
colorbar
colormap gray
set(gca,'XTickLabels',freqs)
set(gca,'XTick',1:length(freqs),'XTickLabels',freqs)

figure
plot(dist)






