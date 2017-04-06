function moveAnalysedFiles
da = 'E:\dataAnalysed\';
folders = dir(da);
folders(strcmp({folders.name},'.'))=[];
folders(strcmp({folders.name},'..'))=[];
folders(~[folders.isdir])=[];

for ii=1:length(folders)
    subdir = dir([da folders(ii).name]);
    subdir(strcmp({subdir.name},'.'))=[];
    subdir(strcmp({subdir.name},'..'))=[];
    
    for jj=1:length(subdir)
        ssdir = dir([da folders(ii).name '\' subdir(jj).name]);
        ssdir(strcmp({ssdir.name},'.'))=[];
        ssdir(strcmp({ssdir.name},'..'))=[];
        
        for kk=1:length(ssdir)
            if isdir(['C:\data\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name])
                files = dir(['C:\data\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name]);
                files(~cellfun(@isempty,strfind({files.name},'.tif')))=[];
                files(~cellfun(@isempty,strfind({files.name},'.db')))=[];
                files([files.isdir])=[];
                for mm = 1:length(files)
                    movefile(['C:\data\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\' files(mm).name],...
                        [da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\']);
                end
            end
        end
    end
end



