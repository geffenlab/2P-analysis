% FRA tone analysis
clear
cmdataLoc = 'E:\dataAnalysed\K056\fearCellMatching\';
cmfiles = dir([cmdataLoc '*cellMatch*.mat']);

for ff=1:length(cmfiles)
    load([cmdataLoc cmfiles(ff).name])
    mouse = cellMatching.session{2}{1};
    date = cellMatching.session{2}{2};
    exptNo = cellMatching.session{2}{3};
    dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
    proc = dir([dataLoc exptNo '\*proc.mat']);
    load([dataLoc exptNo '\' proc.name])
    saveSig = 1;
    
    % get Ca data
    n = find([dat.stat.iscell]==1);
    traces = dat.Fcell{1}(n,:); % raw fluorescence
    npilTraces = dat.FcellNeu{1}(n,:); % neuropil estimates
    npilCoeffs = [dat.stat(n).neuropilCoefficient]'; % neuropil coefficients
    npilSubTraces = traces-(npilTraces.*npilCoeffs); % subtract neuropil from traces
    st = {dat.stat(n).st}; % spike times
    sa = {dat.stat(n).c}; % spike amplitudes (relative)
    
    % load events
    stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
    files = dir([stimLoc '*events.mat']);
    if isempty(files)
        files = dir([stimLoc '*exptInfo.mat']);
    end
    events = load([stimLoc files(1).name]);
    
    fr = 1/(mean(diff(events.frameOn/events.fs))); % mean frame rate
    eventsOn = floor(events.eventOn/events.fs*fr);
    % load stim info
    stimLoc=['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
    files = dir([stimLoc '*stimInfo.mat']);
    if ~isempty(files)
        load([stimLoc files(1).name])
    else
        stimInfo = events.exptInfo.stimInfo{1};
    end
    
    stimDur = stimInfo.tDur/1000;
    ITI = stimInfo.ITI/1000;
    
    % Make raster
    preEv = floor(1*fr); % 1 second before stim
    postEv = ceil((stimDur+ITI)*fr); % x seconds post stim
    
    % filter the traces?
    % [bf,af] = butter(2,[0.5/fr/2 10/fr/2]); % 10 is good for low pass
    [bf,af] = butter(2,10/fr/2,'low'); % 10 is good for low pass
    % [bf,af] = butter(2,0.5/fr/2,'high');
    npst_f=[];
    for ii=1:size(npilSubTraces,1)
        npst_f(ii,:) = filtfilt(bf,af,double(npilSubTraces(ii,:)));
%         npst_z(ii,:) = zscore(npilSubTraces(ii,:));
    end
    
    raster = makeCaRaster(npst_f,eventsOn,preEv,postEv,1);
    rawRast = makeCaRaster(npst_f,eventsOn,preEv,postEv,0);
    preRast = npilSubTraces(:,eventsOn(1)-preEv:eventsOn(1)-1);
    
    % Make spike 'traces' and raster
    spikeTraces = zeros(size(npilSubTraces,1),size(npilSubTraces,2)*2);
    for ii=1:size(spikeTraces,1)
        %     spikeTraces(ii,st{ii})=sa{ii};
        kernel = dat.stat(ii).kernel';
        for jj=1:length(st{ii})
            spikeTraces(ii,st{ii}(jj):st{ii}(jj)+length(kernel)-1) = spikeTraces(ii,st{ii}(jj):st{ii}(jj)+length(kernel)-1) + kernel*sa{ii}(jj);
        end
    end
    spikeTraces = spikeTraces(:,1:length(npilSubTraces));
    spikeRaster = makeCaRaster(spikeTraces,eventsOn,preEv,postEv,0);
    
    % order it by frequency and attenuation
    p = stimInfo.index(stimInfo.order,:);
    order = repmat(p',1,floor(length(eventsOn)/length(p)))';
    [x,ind]=sortrows(order,[1 2]);
    raster = raster(ind,:,:);
    rawRast = rawRast(ind,:,:);
    spikeRaster = spikeRaster(ind,:,:);
    
    
    % Work out which cells respond to sounds
    uT = unique(order(:,1));
    if size(order,2)>1
        uA = unique(order(:,2));
    else
        uA=1;
    end
    
    ha = zeros(length(uT),length(uA),length(n));
    for jj=1:length(n)
        for tt = 1:length(uT)
            for aa=1:length(uA)
                rows = x(:,1)==uT(tt) & x(:,2)==uA(aa);
                mPre = mean(rawRast(rows,1:preEv-1,jj),1);
                stdPre = std(mPre);
                %             mPost = max(mean(rawRast(rows,preEv:preEv*3,jj),1));
                mPost = mean(mean(rawRast(rows,preEv:preEv*3,jj),2));
                c = (mPost-mean(mPre))/stdPre;
                ha(tt,aa,jj)=c>4; % 3 for <10 filter
            end
        end
        
    end
    
    ha = squeeze(sum(sum(ha),2))>0;
    disp([num2str(sum(ha)) ' significantly responding cells out of ' num2str(size(npilSubTraces,1))])
    
    if saveSig==1
        sig = ha;
        save(['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\sigResp.mat'],'sig');
    end
    
    BF = zeros(1,size(raster,3));
    FRA = zeros(length(uT),size(raster,3));
    for jj=1:length(n)
        for tt = 1:length(uT)
            for aa=1:length(uA)
                rows = x(:,1)==uT(tt) & x(:,2)==uA(aa);
                mPre = mean(rawRast(rows,1:preEv-1,jj),1);
                stdPre = std(mPre);
                mPost = mean(rawRast(rows,preEv:preEv*3,jj));
                resp = (mPost-mean(mPre))/stdPre;
                FRA(tt,jj) = max(resp);
            end
        end
        BF(jj) = uT(max(FRA(:,jj))==FRA(:,jj));
    end
    
    index = cellMatching.index(:,2);
    
    
    data(ff+1).mouse = mouse;
    data(ff+1).date = date;
    data(ff+1).exptNo = exptNo;
    data(ff+1).FRA = FRA;
    data(ff+1).BF = BF;
    data(ff+1).index = index;
    data(ff+1).sig = sig;
end
    
save(['E:\dataAnalysed\' mouse '\sfnData.mat'],'data')