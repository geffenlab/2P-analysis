% separate blocks sorted together:

mouse = 'K056';
date = '20170502';
fold_loc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
folders = dir(fold_loc);
folders(strcmp({folders.name},'.'))=[]; folders(strcmp({folders.name},'..'))=[];

% define folders:
us = strfind(folders.name,'_');
f{1} = folders.name(1:us(1)-1);
for ii=2:length(us)
    f{ii} = folders.name(us(ii-1)+1:us(ii)-1);
end
f{end+1} = folders.name(us(end)+1:end);

%% Move the xml files ...
for ii=1:length(f)
    fLoc = ['C:\data\' mouse '\' '\' date mouse '_tifStacks\' f{ii} '\'];
    if ~isdir([fold_loc f{ii} '\']); mkdir([fold_loc f{ii} '\']); end
    files = dir(fLoc);
    files(~cellfun(@isempty,strfind({files.name},'.tif')))=[];
    files(~cellfun(@isempty,strfind({files.name},'.db')))=[];
    files([files.isdir])=[];
    for mm = 1:length(files)
        movefile([fLoc files(mm).name],[fold_loc f{ii} '\']);
    end
end

%% move exptInfo into correct folders
ei = dir([fold_loc folders.name '\*exptInfo.mat']);
if ~isempty(ei)
    t = zeros(1,length(ei));
    for ii=1:length(ei)
        us = strfind(ei(ii).name,'_');
        t(ii) = str2double(ei(ii).name(us(1)+1:us(1)+4)); % times of expt files
    end
else
    t = zeros(1,length(f));
    for ii=1:length(f)
        exf = dir([fold_loc f{ii} '\*exptInfo.mat']);
        us = strfind(exf.name,'_');
        t(ii) = str2double(exf.name(us(1)+1:us(1)+4)); % times of expt files
    end
end

for ii=1:length(f)
    xml = dir([fold_loc f{ii} '\*xml']);
    recInfo = readXmlFile([fold_loc f{ii} '\' xml.name]);
    save([fold_loc f{ii} '\recInfo.mat'],'recInfo')
    ri(ii) = recInfo;
end

nFrames = [ri.nFrames];
time = {ri.timeOfRec};
for ii=1:length(time)
    c = strfind(time{ii},':');
    t2(ii) = str2double(time{ii}([1:c(1)-1,c(1)+1:c(2)-1])); % times of recordings = order of folders
end

for ii=1:length(t)
    ind(ii) = knnsearch(t',t2(ii));
end

proc = dir([fold_loc folders.name '\*_proc.mat']);
load([fold_loc folders.name '\' proc.name]);
datAll = dat;
for ii=1:length(f)     
    if ~isempty(ei)
        ei_file = ei(ind(ii)).name;
        movefile([fold_loc folders(2).name '\' ei_file],[fold_loc f{ii} '\'])
    end
    dat = datAll; dat.Fcell = dat.Fcell(ii); dat.FcellNeu = dat.FcellNeu(ii);
    save([fold_loc f{ii} '\' proc.name],'dat');
end

%% Add recInfo to exptInfo

for ii=1:length(f)
    eit = dir([fold_loc f{ii} '\*_exptInfo.mat']);
    load([fold_loc f{ii} '\' eit.name])
    rci = dir([fold_loc f{ii} '\recInfo.mat']);
    load([fold_loc f{ii} '\' rci.name])
    exptInfo.recInfo = recInfo;
    delete([fold_loc f{ii} '\' rci.name])
    save([fold_loc f{ii} '\' eit.name],'exptInfo','-append')
end
    












