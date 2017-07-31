clear

% SET VARIABLES
mouse = 'K056';
dataLoc = 'E:\dataAnalysed\';
cmSessions = dir([dataLoc mouse '\*cellMatch*']);
cm=8; % which session do you want to plot? CHANGE HERE


%%
session=1;
load([dataLoc mouse '\' cmSessions(cm).name])
date = cell2mat(cellMatching.session{session}(2)); % start with the pre-fear conditioning recording
exptNo = cell2mat(cellMatching.session{session}(3));
dataLoc2 = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc2 exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

% Get first session image
im =dat.mimg(:,:,2);

% Get 2nd session image
session=2;
date = cell2mat(cellMatching.session{session}(2)); % start with the pre-fear conditioning recording
exptNo = cell2mat(cellMatching.session{session}(3));
dataLoc2 = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
load([dataLoc2 exptNo '\F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat'])

im2 =dat.mimg(:,:,2);
% im_norm1 = uint16(im2);
f1=figure;
colormap gray
imb = brighten(double(im2),0.25);

% Transform the 2nd image
Roriginal = imref2d(size(im));
recovered = imwarp(imb,cellMatching.transfrom,'OutputView',Roriginal);

% Get the centroids
a = find(~isnan(cellMatching.index(:,2)));
b = cellMatching.index(~isnan(cellMatching.index(:,2)),2);

% plot the overlaid images
imshowpair(im,recovered)

% Plot the cell outlines
figure('position',[530          92        1287        1178])
imagesc(recovered)
colormap gray
hold on

for ii=1:length(a)
    plot(cellMatching.ROIs{1}{a(ii)}(:,2),cellMatching.ROIs{1}{a(ii)}(:,1),'y');
    plot(cellMatching.ROIs{2}{b(ii)}(:,2),cellMatching.ROIs{2}{b(ii)}(:,1),'m');
    text(cellMatching.centroids{1}(a(ii),1),cellMatching.centroids{1}(a(ii),2),num2str(cellMatching.index(a(ii))),'Color','y');
%     text(cellMatching.centroids{1}(ii,1),cellMatching.centroids{1}(ii,2),num2str(ii),'Color','y');
%     pause()
end

