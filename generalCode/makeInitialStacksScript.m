clear

% Enter info here
mouse = 'K019';
dateOfRecording = '20160914';
dataFolder = [dateOfRecording mouse];

% open and/or make excel file for mouse
dataAnalysisComp = '\\DESKTOP-GK8OVIP\data\';
f = dir([dataAnalysisComp mouse '\']);
csvExist = any(~cellfun(@isempty,strfind({f.name},'.txt')));
csvName = [dataAnalysisComp mouse '\' mouse '_experiments.txt'];
if ~csvExist
    headers = {'date','recording name','folder number'};
    fid = fopen(csvName,'a');
    for ii=1:length(headers)
        fprintf(fid,'%s\t',headers{ii});
    end
else
    fid = fopen(csvName,'a');
end


%% make stacks
dataPath = ['E:\' mouse '\' dataFolder '\'];
folders = dir(dataPath);
folders = folders([folders.isdir]);
folders(strcmp({folders.name},'.'))=[]; % get rid of stupid windows directories
folders(strcmp({folders.name},'..'))=[]; % get rid of stupid windows directories
folders(~cellfun(@isempty,regexp({folders.name},'.tifStacks')))=[]; % get rid of one you are saving to
tic
ind = 1;
for ii =1:length(folders)
    disp(['Folder ' num2str(ii) '/' num2str(length(folders))])
    tiffs = dir([dataPath folders(ii).name '\*.tif']);
    if length(tiffs)>1
        saveFileName =  folders(ii).name;
        %         dash = strfind(folders(ii).name,'-');
        %         fne = folders(ii).name(dash(1)+1:end);
        saveFolder = ['G:\TEMP\' dataFolder '_tifStacks\' sprintf('%02d',ind) '\'];
        tiffStackMaker([dataPath folders(ii).name],'.ome',saveFileName,saveFolder);
        fprintf(fid,'\n %s\t%s\t%02d',dateOfRecording,saveFileName,ind);
        ind = ind+1;
        otherData = dir([dataPath folders(ii).name]);
        %       otherData = otherData(~[otherData.isdir]);
        otherData(strcmp({otherData.name},'.'))=[];
        otherData(strcmp({otherData.name},'..'))=[];
        otherData(~cellfun(@isempty,regexp({otherData.name},'.tif')))=[];
        otherData(~cellfun(@isempty,regexp({otherData.name},'CYCLE')))=[];
        otherData(~cellfun(@isempty,regexp({otherData.name},'Cycle')))=[];
        for jj=1:length(otherData)
            if otherData(jj).isdir
                copyfile([dataPath folders(ii).name '\' otherData(jj).name], [saveFolder otherData(jj).name])
            else
                copyfile([dataPath folders(ii).name '\' otherData(jj).name], saveFolder)
            end
        end
        
    end
    
    
    toc
end

for ii =1:length(folders)
    disp(['Folder ' num2str(ii) '/' num2str(length(folders))])
    tiffs = dir([dataPath folders(ii).name '\*.tif']);
    if length(tiffs)==1
        saveFileName =  folders(ii).name;
        saveFolder = ['G:\TEMP\' dataFolder '_tifStacks\singleImages\'];
        copyfile([dataPath folders(ii).name], saveFolder)
        
    end
end
fclose(fid);
%% Now transfer the files to the analysis computer...
disp('Moving files')
if ~isdir(['\\DESKTOP-GK8OVIP\data\' mouse '\' dataFolder '_tifStacks\'])
    mkdir(['\\DESKTOP-GK8OVIP\data\' mouse '\' dataFolder '_tifStacks\']);
end
movefile(['G:\TEMP\' dataFolder '_tifStacks'], ['\\DESKTOP-GK8OVIP\data\' mouse '\'])


disp('finished');







