%% Makes tif stacks 2000 frames long and then transfers them to the 2P analysis computer via the 10 gigabit ethernet cable.
% If the movefile stops working (can happen after windows update - I hate
% Windows), go to network ans sharing center. Temporarily disable the
% internet connection then click on the 'unidentified network' at the top
% of the window, make sure that sharing and network discovery is enabled.
% To reenable internet, go to device manager and enable the network from
% there. PROBLEM (HOPEFULLY) SOLVED.

clear

% Enter info here
mouse = {'K087'};
dateOfRecording = '20171122';

for mm = 1:length(mouse)
    
dataFolder = [dateOfRecording mouse{mm}];

% open and/or make excel file for mouse
dataAnalysisComp = '\\DESKTOP-GK8OVIP\data\';
f = dir([dataAnalysisComp mouse{mm} '\']);
csvExist = any(~cellfun(@isempty,strfind({f.name},'.csv')));
csvName = [dataAnalysisComp mouse{mm} '\' mouse{mm} '_experiments.csv'];
if ~csvExist
    headers = {'date','recording name','folder number'};
    if ~isdir([dataAnalysisComp mouse{mm} '\']); mkdir([dataAnalysisComp mouse{mm} '\']); end
    fid = fopen(csvName,'w');
    fprintf(fid, '%s,', headers{1,1:end-1}) ;
    fprintf(fid, '%s', headers{1,end}) ;
    fclose('all');
end

fid = fopen(csvName,'a');

%% make stacks
% dataPath = ['G:\tempDataStorage\' mouse{mm} '\' dataFolder '\'];
dataPath = ['E:\' mouse{mm} '\' dataFolder '\'];
folders = dir(dataPath);
folders = folders([folders.isdir]);
folders(strcmp({folders.name},'.'))=[]; % get rid of stupid windows directories
folders(strcmp({folders.name},'..'))=[]; % get rid of stupid windows directories
folders(~cellfun(@isempty,regexp({folders.name},'.tifStacks')))=[]; % get rid of one you are saving to
tic
ind = 1;
for ii = 1:length(folders)
    disp(['Folder ' num2str(ii) '/' num2str(length(folders))])
    tiffs = dir([dataPath folders(ii).name '\*.tif']);
    if length(tiffs)>1
        saveFileName =  folders(ii).name;
        %         dash = strfind(folders(ii).name,'-');
        %         fne = folders(ii).name(dash(1)+1:end);
        saveFolder = ['G:\TEMP\' dataFolder '_tifStacks\' sprintf('%d',ind) '\'];
        tiffStackMaker([dataPath folders(ii).name],'.ome',saveFileName,saveFolder);
%         fprintf(fid,'\n%s,%s,%02d',dateOfRecording,saveFileName,ind);
        
        otherData = dir([dataPath folders(ii).name]);
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
    ind = ind+1;
    
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
fclose('all');
%% Now transfer the files to the analysis computer...
disp('Moving files to 2P analysis computer')
if ~isdir(['\\DESKTOP-GK8OVIP\data\' mouse{mm} '\' dataFolder '_tifStacks\'])
    mkdir(['\\DESKTOP-GK8OVIP\data\' mouse{mm} '\' dataFolder '_tifStacks\']);
end
% [Status, Msg] = FileRename(['G:\TEMP\' dataFolder '_tifStacks'], ['\\DESKTOP-GK8OVIP\data\' mouse{mm} '\']);
movefile(['G:\TEMP\' dataFolder '_tifStacks'], ['\\DESKTOP-GK8OVIP\data\' mouse{mm} '\'])
disp('Finished moving data to 2P analysis computer')
% delete raw files and move tifs to temporary drive
disp('Moving raw movies to temp storage on G drive')
% [Status, Msg] = FileRename(dataPath(1:end-1), ['G:\tempDataStorage\' mouse{mm} '\']);
movefile(dataPath(1:end-1), ['G:\tempDataStorage\' mouse{mm} '\']);

disp(['finished mouse ' sprintf('%02d',mm)]);

end

disp('Finished all')



