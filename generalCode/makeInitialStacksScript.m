clear

%% Enter info here
mouse = 'K048';
dataFolder = '20170205K048';

%% make stacks
dataPath = ['F:\' mouse '\' dataFolder '\'];
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
        saveFolder = [dataPath dataFolder '_tifStacks\' num2str(ind) '\'];
        tiffStackMaker([dataPath folders(ii).name],'.ome',saveFileName,saveFolder)
        ind = ind+1;
      otherData = dir([dataPath folders(ii).name]);
      otherData = otherData(~[otherData.isdir]);
      otherData(~cellfun(@isempty,regexp({otherData.name},'.tif')))=[];
      for jj=1:length(otherData)
        copyfile([dataPath folders(ii).name '\' otherData(jj).name], saveFolder)
      end

    end
    
  
    toc
end
toc
    
%% Now transfer the files to the analysis computer...
if ~isdir(['\\DESKTOP-GK8OVIP\data\' mouse '\' dataFolder '_tifStacks\'])
    mkdir(['\\DESKTOP-GK8OVIP\data\' mouse '\' dataFolder '_tifStacks\']);
end
movefile([dataPath dataFolder '_tifStacks'], ['\\DESKTOP-GK8OVIP\data\' mouse '\'])


disp('finished');
    
    
    
    
    
    
    
    