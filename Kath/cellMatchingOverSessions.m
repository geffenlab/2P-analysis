% Make cell mapping

clear
close all
mouse = 'K070';
date1 = '20170728'; % start with the pre-fear conditioning recording
exptNo = '1';
expt1 = {mouse,date1,exptNo};
dataLoc = ['E:\dataAnalysed\' mouse '\' date1 mouse '_tifStacks\'];
load([dataLoc exptNo '\F_' mouse '_' date1 mouse '_tifStacks_plane1_proc.mat'])
n1 = find([dat.stat.iscell]==1);
im = dat.mimg(:,:,2);
im_norm1 =dat.mimg(:,:,2);
im_norm1 = uint16(im_norm1);
f1=figure;
colormap gray
imb = brighten(double(im),1);
imb = brighten(double(imb),0.25);
imagesc(imb)
% imagesc(im)
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


%% Now open the file you want to map with it:

date2 = '20170730'; % start with the pre-fear conditioning recording
exptNo = '2';
expt2 = {mouse,date2,exptNo};
dataLoc = ['E:\dataAnalysed\' mouse '\' date2 mouse '_tifStacks\'];
proc = dir([dataLoc exptNo '\*proc.mat']);
load([dataLoc exptNo '\' proc.name])
n2 = find([dat.stat.iscell]==1);
im2 = dat.mimg(:,:,2);
im_norm2 = dat.mimg(:,:,2);
im_norm2 = uint16(im_norm2);
f2 = figure;
colormap gray
im2b = brighten(double(im2),1);
im2b = brighten(double(im2b),0.25);
imagesc(im2b)
% imagesc(im2)
dat2 = dat;




%
% correlationOutput = normxcorr2(im2(:,:,1), im(:,:,1));
%
% [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
% [ypeak, xpeak] = ind2sub(size(correlationOutput),maxIndex(1));
% corr_offset = [(xpeak-size(im2,2)) (ypeak-size(im2,1))];


%
% %
figure(f1)
fixedPoints = ginput();
hold on
plot(fixedPoints(:,1),fixedPoints(:,2),'x')
hold on
for ii=1:length(fixedPoints)
    t = text(fixedPoints(ii,1)-5,fixedPoints(ii,2)+7,num2str(ii)); t.Color = [1 1 0];
end
figure(f2)
movingPoints = ginput(length(fixedPoints));

%% transform
method = 'projective';
tform = fitgeotrans(movingPoints,fixedPoints,method);
tFormInv = invert(tform);
Tinv = tFormInv.T;
ss = Tinv(2,1);
sc = Tinv(1,1);
scale_recovered = sqrt(ss*ss + sc*sc);
theta_recovered = atan2(ss,sc)*180/pi;

Roriginal = imref2d(size(im));
recovered = imwarp(im2,tform,'OutputView',Roriginal);

f3=figure;
imshowpair(im,recovered)


%% Get transformed centroids and boundaries
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

% have a look at the centroid locations
f4 = figure;
colormap gray
imagesc(im)
hold on
plot(cent(:,1),cent(:,2),'gx')
plot(cent2(:,1),cent2(:,2),'mx')

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
thresh = 20;
index = 1:length(n1);
index2 = NaN(size(index));
index2(nnd<thresh) = nn(nnd<thresh);
while sum(~isnan(index2))> length(unique(index2(~isnan(index2))))
    thresh = thresh-0.1;
    index = 1:length(n1);
    index2 = NaN(size(index));
    index2(nnd<thresh) = nn(nnd<thresh);
end
f5 = figure;
colormap gray
imagesc(im_norm1)
hold on
for ii=1:length(index2)
    if ~isnan(index2(ii))
        plot(boundaries{index(ii)}(:,2),boundaries{index(ii)}(:,1),'Color','g','LineWidth',1.5)
        plot(boundaries2{index2(ii)}(:,2),boundaries2{index2(ii)}(:,1),'Color','m','LineWidth',1.5)
        text(cent(ii,1),cent(ii,2),num2str(ii),'Color','m');
        text(cent2(index2(ii),1),cent2(index2(ii),2),num2str(ii),'Color','y')
        pause()
    end
end
plot(cent(index,1),cent(index,2),'gx')
plot(cent2(:,1),cent2(:,2),'mx')

% Check there are no duplicates
if length(unique(index))<length(index) || length(unique(index2))<length(index)
    disp('DUPLICATES!!!!')
end

disp([num2str(sum(~isnan(index2))) ' neurons matched'])

%% save stuff

cellMatching.index = [index',index2'];
cellMatching.session = {expt1,expt2};
cellMatching.transfrom = tform;
cellMatching.refFixedPoints = fixedPoints;
cellMatching.refMovingPoints = movingPoints;
cellMatching.imRegMethod = method;
cellMatching.threshold = thresh;
cellMatching.ROIs = {boundaries,boundaries2};
cellMatching.centroids = {cent,cent2};


save(['E:\dataAnalysed\' mouse '\' mouse '_' date1 '_' date2 '_cellMatch_FRA.mat'],'cellMatching')



