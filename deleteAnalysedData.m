function deleteAnalysedData

%% Check all info has been moved to the analysed data file...

moveAnalysedFiles

%% Delete analysed data ONLY if it is backed up

da = 'C:\data\';
folders = dir(da);
folders(strcmp({folders.name},'.'))=[];
folders(strcmp({folders.name},'..'))=[];
folders(strcmp({folders.name},'retigaData'))=[];
folders(strcmp({folders.name},'_tempreg'))=[];
folders(~[folders.isdir])=[];

for ii=1:length(folders)
    subdir = dir([da folders(ii).name]);
    subdir(strcmp({subdir.name},'.'))=[];
    subdir(strcmp({subdir.name},'..'))=[];
    subdir(~[subdir.isdir])=[];
    
    for jj=1:length(subdir)
        ssdir = dir([da folders(ii).name '\' subdir(jj).name]);
        ssdir(strcmp({ssdir.name},'.'))=[];
        ssdir(strcmp({ssdir.name},'..'))=[];
        ssdir(~[ssdir.isdir])=[];
        
        for kk=1:length(ssdir)
            if isdir(['E:\dataAnalysed\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name])
                if isdir(['G:\dataAnalysed\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name]) && isdir(['I:\dataAnalysed\' folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name])
                    files = dir([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name]);
                    files(strcmp({files.name},'.'))=[];
                    files(strcmp({files.name},'..'))=[];
                    for mm = 1:length(files)
                        if files(mm).isdir
                            subfiles = dir([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\' files(mm).name]);
                            subfiles(strcmp({subfiles.name},'.'))=[];
                            subfiles(strcmp({subfiles.name},'..'))=[];
                            for nn=1:length(subfiles)
                                delete([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\' files(mm).name '\' subfiles(nn).name]);
                            end
                            rmdir([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\' files(mm).name]);
                        else
                            delete([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name '\' files(mm).name]);
                        end
                    end
                    rmdir([da folders(ii).name '\' subdir(jj).name '\' ssdir(kk).name])
                    check = dir([da folders(ii).name '\' subdir(jj).name]);
                    check(strcmp({check.name},'.'))=[];
                    check(strcmp({check.name},'..'))=[];
                    if isempty(check)
                        rmdir([da folders(ii).name '\' subdir(jj).name])
                    end
                end
            end
        end
    end
end