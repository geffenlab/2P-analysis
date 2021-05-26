% separate blocks sorted together:
clear
mouse = {'XD036'};
date = {'20210422'};

for zz = 1:length(mouse)
    for dd = 1:length(date)
        fold_loc = ['E:\dataAnalysed\' mouse{zz} '\' date{dd} mouse{zz} '_tifStacks\'];
        folders = dir(fold_loc);
        folders(contains({folders.name},'.'))=[];
        folders(contains({folders.name},'_'))=[]; 
        folders(contains({folders.name},'suite2p')) = [];
        
        for uu = 1:length(folders)
            
            
            % Move the xml files ...
            fLoc = ['C:\data\' mouse{zz} '\' date{dd} mouse{zz} '_tifStacks\' folders(uu).name '\'];
            files = dir(fLoc);
            files(~cellfun(@isempty,strfind({files.name},'.tif')))=[];
            files(~cellfun(@isempty,strfind({files.name},'.db')))=[];
            files([files.isdir])=[];
            for mm = 1:length(files)
                movefile([fLoc files(mm).name],[fold_loc folders(uu).name '\']);
            end
            
            % get xml data
            xml = dir([fold_loc folders(uu).name '\*xml']);
            fprintf('Reading xml file...\n')
            recInfo = readXmlFile_v2_20170730([fold_loc folders(uu).name '\' xml.name]);
            save([fold_loc folders(uu).name '\recInfo.mat'],'recInfo')
            ri = recInfo;
            
            
            colon = strfind(ri.date,':');
            space = strfind(ri.date,' ');
            %             tr(ii) = str2double([ri(ii).date(space+1:colon(1)-1),ri(ii).date(colon(1)+1:colon(2)-1)]);
            tr = datenum(ri.date,'yyyy-mm-dd HH:MM:SS');
            
            %% move exptInfo into correct folders
            
            %         exptInfoFold = [fold_loc folders(uu).name '\*exptInfo.mat'];
            %         ei = dir(exptInfoFold);
            %         if ~isempty(ei)
            %             t = zeros(1,length(ei));
            %             for ii=1:length(ei)
            %                 us = strfind(ei(ii).name,'_');
            %                 t(ii) = str2double(ei(ii).name(us(1)+1:us(1)+4)); % times of expt files
            %             end
            %         else
            
            exptInfoFold = ['C:\data\' mouse{zz} '\' date{dd} mouse{zz} '_tifStacks\*exptInfo.mat'];
            exf = dir(exptInfoFold);
            t = zeros(1,length(exf));
            for ii = 1:length(t)
                us = strfind(exf(ii).name,'_');
                %                 t(ii) = str2double(exf(ii).name(us(1)+1:us(1)+4)); % times of expt files
                t(ii) = datenum(exf(ii).name(1:us(2)-1),'yymmdd_HHMMSS');
            end
            
            x = etime(datevec(t),datevec(tr));
            tm = knnsearch(x,0);
            %                 tm(ii) = knnsearch(t',tr(ii));
            movefile(['C:\data\' mouse{zz} '\' date{dd} mouse{zz} '_tifStacks\' exf(tm).name], [fold_loc folders(uu).name '\'])
            
            
            
            
            
            %% Split the processed suite 2P file
            
            proc = dir([fold_loc folders(uu).name '\*_proc.mat']);
            load([fold_loc folders(uu).name '\' proc.name]);
            nframes = [0,dat.ops.Nframes];
            bf = dat.ops.badframes;
            badFrames = cell(1,length(nframes)-1);
            for ii = 2:length(nframes)
                badFrames{ii-1} = bf(sum(nframes(1:ii-1))+1:sum(nframes(1:ii)));
            end
            
            datAll = dat;
            dat = datAll;
            dat.badFrames = badFrames;
            save([fold_loc folders(uu).name '\' proc.name],'dat');
            
            
            %% Add recInfo to exptInfo
            
            eit = dir([fold_loc folders(uu).name '\*exptInfo.mat']);
            load([fold_loc folders(uu).name '\' eit.name])
            rci = dir([fold_loc folders(uu).name '\recInfo.mat']);
            load([fold_loc folders(uu).name '\' rci.name])
            exptInfo.recInfo = recInfo;
            delete([fold_loc folders(uu).name '\' rci.name])
            save([fold_loc folders(uu).name '\' eit.name],'exptInfo','-append')
            
        end
        
        clearvars -except mouse date
    end
end

disp('Finished separating blocks');






