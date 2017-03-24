% plot cells over multiple sessions:
clear
mouse = 'K048';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
ii=1
load([dataLoc mouse '\' cmSessions(ii).name])
date1 = cellMatching.session{1}(2); % start with the pre-fear conditioning recording
exptNo = cellMatching.session{1}(3);
dataLoc = ['E:\dataAnalysed\' mouse '\' date1{1} mouse '_tifStacks\'];
load([dataLoc exptNo{1} '\F_' mouse '_' date1{1} mouse '_tifStacks_plane1_proc.mat'])
n1 = find([dat.stat.iscell]==1);
im = dat.mimg(:,:,2);
im = uint16(im);
f1=figure;
colormap gray
imb = brighten(double(im),1);
imagesc(imb)
dat1 = dat;

% get centroids and boundaries
for ii = 1:length(n1)
    img = zeros(1,size(im,1)*size(im,2));
    img(dat1.stat(n1(ii)).ipix)=1;
    img = reshape(img,size(im,1),size(im,2));
    c = regionprops(img,'Centroid');
    if ~isempty(c)
        cent(ii,:) = round(c(end).Centroid);
    end
    s = [bwboundaries(img)];
    boundaries(ii) = s(1);
end

date2 = cellMatching.session{2}(2); % start with the pre-fear conditioning recording
exptNo = cellMatching.session{2}(3);
dataLoc = ['E:\dataAnalysed\' mouse '\' date2{1} mouse '_tifStacks\'];
load([dataLoc exptNo{1} '\F_' mouse '_' date2{1} mouse '_tifStacks_plane1_proc.mat'])
n2 = find([dat.stat.iscell]==1);
im2 = dat.mimg(:,:,2);
im2 = uint16(im2);
f2 = figure;
colormap gray
im2b = brighten(double(im2),1);
imagesc(im2b)
dat2 = dat;

tform = cellMatching.transfrom;

Roriginal = imref2d(size(im));
recovered = imwarp(im2,tform,'OutputView',Roriginal);

cent2 = zeros(length(n2),2); boundaries2 = cell(1,length(n2));
for ii = 1:length(n2)
    img = zeros(1,size(im2,1)*size(im2,2));
    img(dat2.stat(n2(ii)).ipix)=1;
    img = reshape(img,size(im2,1),size(im2,2));
    recovered2 = imwarp(img,tform,'OutputView',Roriginal);
    c = regionprops(recovered2,'Centroid');
    s = [bwboundaries(recovered2)];
    if ~isempty(c)
        cent2(ii,:) = round(c(end).Centroid);
        boundaries2(ii) = s(1);
    end
end


f3=figure;
blank = zeros(size(im));
imshowpair(blank,recovered)


%% work out the closest neighbours to the 2nd session
for ii=1:length(n1)
    dist = sqrt(sum((cent2-cent(ii,:)).^2,2));
    a = find(min(dist)==dist);
    if length(a)>1
        nn(ii) = a(randperm(length(a),1));
    else
        nn(ii)=a;
    end
    nnd(ii) = dist(nn(ii));
end

% plot boundaries of selected cells
thresh = 6;
index = find(nnd<thresh);
index2 = nn(nnd<thresh);
% f5 = figure;
% colormap gray
% imagesc(brighten(double(recovered),1))
hold on
for ii=1:length(index2)
    plot(boundaries{index(ii)}(:,2),boundaries{index(ii)}(:,1),'Color','g','LineWidth',1.5)
    plot(boundaries2{index2(ii)}(:,2),boundaries2{index2(ii)}(:,1),'Color','m','LineWidth',1.5)
%           pause()
end
plot(cent(index,1),cent(index,2),'gx')
plot(cent2(index2,1),cent2(index2,2),'mx')
