% Look at Rick's data
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


for cc=1:length(communities)
    dd=1;
    c = csvread([dataLoc communities(cc).name]); % load community data
    c(:,logical([cellfun(@isempty,exptNo),1]))=[];
    f = dir([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\*_proc.mat']);
    load([calciumLoc calciumFolders(dd).name '\' exptNo{dd} '\' f.name]);
       
    x=[];
    x(:,1) = cellList(:,dd);
    x(:,2) = data(dd).BF(x(:,1));
    x(:,3) = data(dd).sig(x(:,1));
    x(:,4) = c(:,dd);
   
    
    
    
    
    
    
    %% plot figures
    
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
        img(dat.stat(cell).ipix)=1;
        img = reshape(img,size(im,1),size(im,2));
        cent = regionprops(img,'Centroid');
        celln = num2str(jj);
        t = text(cent.Centroid(1)-3,cent.Centroid(2),celln); t.Color = [1 1 0];
        s = [bwboundaries(img)];
        boundaries(jj) = s(1);
%         if x(jj,3)==1
%             plot(boundaries{jj}(:,2),boundaries{jj}(:,1),'Color',cs(c(jj,ii)==uC,:),'LineWidth',1.5)
            patch(boundaries{jj}(:,2),boundaries{jj}(:,1),cs(x(jj,4)==uC,:),'edgealpha',0,'facealpha',0.75)
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
    
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    