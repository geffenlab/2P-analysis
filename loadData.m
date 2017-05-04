%% Load data set and get basic parameters
%% Info about K048 recordings
% 20170315 = pre-fear conditioning - ALL CELLS ARE ALIGNED TO THIS CELL INDEX
% 20170316 = post first fear conditioning
% 20170318 = 3 days post first fear conditioning/24 hours pre 2nd fear conditioning
% 20170320 = 24 hours post 2nd fear conditioning

%%
clear
mouse = 'K048'; % mouse number
date = '20170315'; % date of recording
exptNo = '2'; % folder that data is contained in
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\']; % location of data folder
load([dataLoc 'F_' mouse '_' date mouse '_tifStacks_plane1_proc.mat']) % load processed data
%% STRUCTURE OF RESULTS (dat) FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.stat.iscell
% st are the deconvolved spike times (in frames)
% c  are the deconvolved amplitudes
% kernel is the estimated kernel
% dat.ops.badframes indicates which frames were not aligned well (probably
% because of excessive z-motion artifact), might be good to remove those
% frames and interpolate the trace which I haven't done here

%% get Ca data
n = find([dat.stat.iscell]==1); % ROIs that have been defined as cells (initially automatic then checked by me manually)
traces = dat.Fcell{1}(n,:); % raw fluorescence
npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
st = {dat.stat(n).st}; % spike times in frames
sa = {dat.stat(n).c}; % spike amplitudes (relative)

%% load/define events and frames
files = dir([dataLoc '*exptInfo.mat']);
if isempty(files)
    files = dir([dataLoc '*events.mat']);
end
events = load([dataLoc files(1).name]);

fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
eventsOn = floor(events.eventOn/events.fs*fr); % in frames

%% load stim info if there is a separate file for it else just define from
files = dir([dataLoc '*stimInfo.mat']);
if ~isempty(files)
    load([dataLoc files(1).name])
else
    stimInfo = events.exptInfo.stimInfo{1};
end
if ~isfield(stimInfo,'tDur')
    toneDuration = stimInfo.stimDur/1000; % duration of each tone in seconds
    ITI = stimInfo.ISI/1000; % inter tone interval in seconds
else
    toneDuration = stimInfo.tDur/1000; % duration of each tone in seconds
    ITI = stimInfo.ITI/1000; % inter tone interval in seconds
end
stimIndex = stimInfo.index(stimInfo.order,:);
stimIndex = repmat(stimIndex',1,floor(length(eventsOn)/length(stimIndex)))';
eventsOn = eventsOn(1:length(stimIndex)); % for 20170315 the end of the recording cut off the 5th repeat of the stimulus so there are only 4 reps for this recording

%% Get cell index if not from recording 1
if ~strcmp(date,'20170315')
    files = dir(['E:\dataAnalysed\' mouse '\*cellMatch*.mat']);
    load([dataLoc files(1).name])
    cellIndex = cellMatching.index;
    % column 2 of cellIndex tells you which cell number of this recording
    % is equivalent to the cell numbers from expt 1    
end


%% Get cell outlines/boundaries of the ROIs and centroids
im = dat.mimg(:,:,2);

centroids = zeros(length(n),2); boundaries = cell(length(n),1);
for ii = 1:length(n)
    img = zeros(1,size(im,1)*size(im,2));
    img(dat.stat(n(ii)).ipix)=1;
    img = reshape(img,size(im,1),size(im,2));
    c = regionprops(img,'Centroid');
    if ~isempty(c)
        centroids(ii,:) = round(c(end).Centroid);
    end
    s = [bwboundaries(img)];
    boundaries(ii) = s(1);
end

% Plot outlines
imb = brighten(double(im),1);
imb = brighten(double(imb),0.75);
% imb = dat.mimg_proc(:,:,2); % the processed image from which cells are defined
figure
imagesc(imb)
colormap gray
hold on
for ii=1:length(n)
        plot(boundaries{ii}(:,2),boundaries{ii}(:,1),'Color','y','LineWidth',1.5)
        text(centroids(ii,1),centroids(ii,2),num2str(ii),'Color','m');
end


%% plot Ca traces
figure('position',[1000,175,637,1163])
time = 1/fr/fr:1/fr:length(npilSubTraces)/fr;
for ii=1:length(n)
    norm = (npilSubTraces(ii,:)-min(npilSubTraces(ii,:)))/max((npilSubTraces(ii,:)-min(npilSubTraces(ii,:))));
    plot(time,norm+(1.5*(ii-1)),'k','LineWidth',0.5)
    hold on
    
    hold on 
end

% for jj=1:length(eventsOn) % ******* UNCOMMENT if want to overlay stimulus onsets ******
%          plot([time(eventsOn(jj)) time(eventsOn(jj))], [0 (sum(ha)+1)*1.5],'r-')
% end

axis tight
set(gca,'YTick',0.75:1.5:length(n)*1.5,'YTickLabels',1:length(n))
ylabel('cell number')
xlabel('Time(s)')



