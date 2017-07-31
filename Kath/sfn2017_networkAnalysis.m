%% Look at Rick's data
clear
calciumLoc = 'E:\dataAnalysed\K056\';
dataLoc = 'E:\dataAnalysed\K056\network\';
load([dataLoc 'sfnData.mat'])
communities = dir([dataLoc '*.csv']);
cellList  = load([dataLoc 'cellList.mat']); cellList = cellList.cells;
calFold = dir(calciumLoc); calFold(~[calFold.isdir])=[];
dates = {'0425','0426','0427','0428'};
exptNo = {'2_3',[],'2','2'};
exptDates = {data.date};
for ii=1:length(dates)
    calciumFolders(ii) = calFold(~cellfun(@isempty,strfind({calFold.name},dates{ii})));
    if sum(~cellfun(@isempty,strfind(exptDates,dates{ii})))~=0
        edu(ii) = find(~cellfun(@isempty,strfind(exptDates,dates{ii}))==1);
    end
end
edu(edu==0)=[];
data = data(edu);
calciumFolders(cellfun(@isempty,exptNo))=[];
cellList(:,cellfun(@isempty,exptNo))=[];
a = exptNo;
exptNo(cellfun(@isempty,exptNo))=[];
%%
clc
for cc=5
    disp(' ')
    disp(['Number of communities: ' num2str(cc+1)])
    dd=1;
    c = csvread([dataLoc communities(cc).name]); % load community data
    c(:,logical([cellfun(@isempty,a),1]))=[];
    f = dir([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\*_proc.mat']);
    load([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\' f.name]);
    
    x=[];
    x(:,1) = cellList(:,dd);
    x(:,2) = data(dd).BF(x(:,1));
    x(:,3) = data(dd).sig(x(:,1));
    x(:,4) = c(:,dd);
    
    uC = unique(x(:,4));
    nC=[]; nSig=[];
    for ii=1:length(uC)
        nC(ii) = sum(x(:,4)==uC(ii));
        nSig(ii) = sum(x(x(:,4)==uC(ii),3));
    end
    
    nShufs = 2000;
    shuf_s = zeros(nShufs,length(nC));
    for ii=1:length(nC)
        for jj=1:nShufs
            shuf_s(jj,ii) = mean(x(randperm(length(x),nC(ii)),3));
        end
    end
    
    ss=[];
    for ii=1:length(nC)
        ss(:,ii) = [(nSig(ii)/nC(ii))<prctile(shuf_s(:,ii),2.5),(nSig(ii)/nC(ii))>prctile(shuf_s(:,ii),97.5)];
    end
    sig_resp{cc}=ss;
    
    for ii=1:length(nC)
        meanbf(ii) = mean( x( x(:,4)==uC(ii) & x(:,3)==1,2));
    end
    
    nShufs = 2000;
    shuf_bf = zeros(nShufs,length(nC));
    bfs = NaN(1,length(x));
    bfs(x(:,3)==1) = x(x(:,3)==1,2);
    for ii=1:length(nC)
        for jj=1:nShufs
            shuf_bf(jj,ii) = nanmean(bfs(randperm(length(x),nC(ii))));
        end
    end
    
    ss=[];
     for ii=1:length(nC)
        ss(:,ii) = [meanbf(ii)<prctile(shuf_bf(:,ii),2.5),meanbf(ii)>prctile(shuf_bf(:,ii),97.5)];
    end
    sig_bf{cc}=ss;
end

%% plot tuning curves

cs = colormap('hsv');
uC = unique(x(:,4)); % unique communities
cs=cs(1:floor(length(cs)/length(uC)):floor(length(cs)/length(uC))*length(uC),:);
uT = unique(x(:,2));
for ii=1:length(x)
    fr = data(dd).FRA(:,x(ii,1));
    if x(ii,3)==1
        plot(uT,fr,'Color',cs(uC==x(ii,4),:),'LineWidth',2);
    else
        patch([uT' flipud(uT)'],[(fr)' flipud(fr)'],[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0);
    end
    hold on
end
set(gca,'xscale','log','XTick',uT)
axis tight
figure
FRA = data(dd).FRA(:,x(:,1));
for ii=1:length(uC)
    plot(uT,mean(FRA(:,x(:,4)==uC(ii) & x(:,3)==1),2))
    hold on
end
plot(uT,mean(FRA(:,x(:,3)==0),2),'color',[0.5 0.5 0.5])

c5 = FRA(:,x(:,4)==5 & x(:,3)==1);
cr = corr(c5)
cn5 = FRA(:,x(:,4)~=5 & x(:,3)==1);
crn = corr(cn5)
figure
histogram(cr(~eye(length(cr))),20)
hold on
histogram(crn(~eye(length(crn))),20)

%% plot figures
cc=5;
c = csvread([dataLoc communities(cc).name]); % load community data
c(:,logical([cellfun(@isempty,exptNo),1]))=[];
f = dir([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\*_proc.mat']);
load([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\' f.name]);

x=[];
x(:,1) = cellList(:,dd);
x(:,2) = data(dd).BF(x(:,1));
x(:,3) = data(dd).sig(x(:,1));
x(:,4) = c(:,dd);
ipix = {dat.stat(logical([dat.stat.iscell])).ipix};
f1 = figure;
im = dat.mimg(:,:,2);
im=brighten(im,1);
colormap gray
imagesc(im);
freezeColors
% plot communities
cs = colormap('hsv');
uC = unique(c(:,1)); % unique communities
cs=cs(1:floor(length(cs)/length(uC)):floor(length(cs)/length(uC))*length(uC),:);
hold on
for jj =  1:length(x)
    cell = x(jj,1);
    img = zeros(1,size(im,1)*size(im,2));
    img(ipix{cell})=1;
    img = reshape(img,size(im,1),size(im,2));
    cent = regionprops(img,'Centroid');
    celln = num2str(jj);
    %         t = text(cent.Centroid(1)-3,cent.Centroid(2),celln); t.Color = [1 1 0];
    s = [bwboundaries(img)];
    boundaries(jj) = s(1);
    %         if x(jj,3)==1
    %             plot(boundaries{jj}(:,2),boundaries{jj}(:,1),'Color',cs(c(jj,ii)==uC,:),'LineWidth',1.5)
    patch(boundaries{jj}(:,2),boundaries{jj}(:,1),cs(x(jj,4)==uC,:),'edgealpha',0,'facealpha',0.65)
    %         else
    %             plot(boundaries{jj}(:,2),boundaries{jj}(:,1),'Color','k','LineWidth',1.5)
    %         end
end
axis off
axis equal

cs = colormap('jet');
uStim = unique(x(:,2));
cs=cs(1:floor(length(cs)/length(uStim)):floor(length(cs)/length(uStim))*length(uStim),:);
hold on
for jj =  1:length(x)
    cell = x(jj,1);
    if x(jj,3)==1
        plot(boundaries{jj}(:,2),boundaries{jj}(:,1),'Color',cs(find(x(jj,2)==uStim),:),'LineWidth',2)
    else
        plot(boundaries{jj}(:,2),boundaries{jj}(:,1),'Color','k','LineWidth',1.5)
    end
end
img = reshape(img,size(im,1),size(im,2));
axis off
axis equal



















