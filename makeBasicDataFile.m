% function makeBasicDataFile(mouse)


% load processed file
mouse = 'MT002';
date = '20180315';
exptNo = '1';
dataLoc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\' exptNo '\'];
proc = dir([dataLoc '\*proc.mat']);

f2a{1} = [proc.folder filesep proc.name];

for ff = 1:length(f2a)
%     clearvars -except mouse f2a ff mouseFolder
    
    disp(f2a{ff})
    % load events/experimental data
    bs = strfind(f2a{ff},'\');
    stimLoc = f2a{ff}(1:bs(5));
    exptNo = f2a{ff}(bs(4)+1:bs(5)-1);
    %     mouse = f2a{ff}(bs(2)+1:bs(3)-1);
    date = f2a{ff}(bs(3)+1:bs(3)+8);
    files = dir([stimLoc '*events.mat']);
    if isempty(files)
        files = dir([stimLoc '*exptInfo.mat']) 
    end
    if ~isempty(files)
        ev = load([stimLoc files(1).name]);
        fr = 1/(mean(diff(ev.frameOn/ev.fs))); % mean frame rate
        events.eventsOn = floor(ev.eventOn/ev.fs*fr);
        events.eventsOff = floor(ev.eventOff/ev.fs*fr);
        
        if isfield(ev,'exptInfo')
            exptInfo = ev.exptInfo;
            if isfield(exptInfo,'stimInfo')
            else
                si = dir([stimLoc '*stimInfo.mat']);
                si = load([stimLoc si.name]);
                exptInfo.stimInfo = si.stimInfo;
            end
        else
            si = dir([stimLoc '*stimInfo.mat']);
            si = load([stimLoc si.name]);
            exptInfo.stimInfo = si.stimInfo;
        end
        
        exptInfo.recDate = date;
        exptInfo.folderNumber = exptNo;
        exptInfo.dataLoc = f2a{ff};
        exptInfo.fr = fr;
        exptInfo.mouse = mouse;
        
        stimInfo = exptInfo.stimInfo;
        
        if iscell(stimInfo)
            stimInfo = stimInfo{1};
        end
        if isfield(stimInfo,'stimInfo')
            stimInfo = stimInfo.stimInfo;
        end
            
        
        %% SAVE the basic data with a useful title
        
        us = strfind(files.name,mouse);
        bl = strfind(files.name,'block');
        if ~isempty(bl)
            exptName = files.name(us+5:bl-2);
        else
            fra = strfind(files.name,'FRA');
            if ~isempty(fra)
                exptName = '2P_FRA';
            else
                exptName = files.name(us+13:end-11);
            end
        end
        
        % mouse_dateOfRecording_exptName_number
        fn = [dataLoc exptInfo.mouse '_' exptInfo.recDate...
            '_' exptName '_' sprintf('%02d',str2double(exptNo)) '.mat'];
        if ~exist(fn,'file')
            load(f2a{ff})
            % get Ca data
            calcium.n = find([dat.stat.iscell]==1);
            checkTraces = dat.Fcell{1}(calcium.n,:); % raw fluorescence
            rmCells = isnan(sum(checkTraces,2));
            clear checkTraces
            calcium.n(rmCells) = [];
            calcium.rawTraces = dat.Fcell{1}(calcium.n,:); % raw fluorescence
            calcium.npilTraces = dat.FcellNeu{1}(calcium.n,:); % neuropil estimates
            calcium.npilCoeffs = [dat.stat(calcium.n).neuropilCoefficient]'; % neuropil coefficients
            calcium.npilSubTraces = calcium.rawTraces-(calcium.npilTraces.*calcium.npilCoeffs); % subtract neuropil from traces
            
            spatialInfo.im = dat.mimg(:,:,2);
            
            for ii = 1:length(calcium.n)
                img = zeros(1,size(spatialInfo.im,1)*size(spatialInfo.im,2));
                img(dat.stat(calcium.n(ii)).ipix)=1;
                img = reshape(img,size(spatialInfo.im,1),size(spatialInfo.im,2));
                spatialInfo.ipix{ii} = dat.stat(calcium.n(ii)).ipix;
                c = regionprops(img,'Centroid');
                spatialInfo.centroid(ii,:) = round(c(end).Centroid);
                b = bwboundaries(img);
                spatialInfo.ROIs{ii} = b{1};
            end
%             
%             if isfield(dat,'sp')
%                 if iscell(dat.sp)
                    spikes.raster = dat.sp{1}(calcium.n,:); % spike times
                    spikes.kernel = dat.kernels{1};
                    spikes.deconvTraces = dat.ca{1}(calcium.n,:);
%                     spikes.kernel = dat.kernels;
%                     spikes.raster = dat.sp(calcium.n,:);
%                     spikes.deconvTraces = dat.ca(calcium.n,:);
%                 else
%                     spikes.raster = dat.sp(calcium.n,:); % spike times
%                 end
%             else
%                 spikes.raster = zeros(size(calcium.rawTraces));
%                 for ii = 1:size(spikes.raster,1)
%                     spikes.raster(ii,dat.stat(ii).st) = dat.stat(ii).c;
%                 end
%                 spikes.kernel = dat.stat(1).kernel;
%             end
            
            
            
            
            save(fn,'calcium','spikes','exptInfo','stimInfo','events','spatialInfo')
        end
        
    end
end
