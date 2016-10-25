clear
mouse = 'K024';
dataFolder = '20161025K024';
dataPath = ['E:\' mouse '\' dataFolder '\'];
folders = dir(dataPath);
folders = folders([folders.isdir]); 
folders(strcmp({folders.name},'.'))=[]; % get rid of stupid windows directories
folders(strcmp({folders.name},'..'))=[]; % get rid of stupid windows directories
folders(strcmp({folders.name},'tifStacks'))=[]; % get rid of one you are saving to
tic
ind = 1;
for ii = 1:length(folders)
    disp(['Folder ' num2str(ii) '/' num2str(length(folders))])
    tiffs = dir([dataPath folders(ii).name '\*.tif']);    
    if length(tiffs)>1
        saveFileName =  folders(ii).name;
        saveFolder = [dataPath dataFolder '_tifStacks\' num2str(ind) '\'];
        tiffStackMaker([dataPath folders(ii).name],'.ome',saveFileName,saveFolder)
        ind = ind+1;
    end
    toc
end
toc
    
%% Now transfer the files to the analysis computer
copyfile([dataPath dataFolder '_tifStacks\'], ['\\DESKTOP-GK8OVIP\data\' mouse '\'])
    
    
    
    
    
    
    
    